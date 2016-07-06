#' Treatment effect estimation
#'
#' Estimate treatment effect and median of the composite endpoint from using
#' imputed data
#'
#' @inheritParams plotCompleters
#' @inheritParams plotImputed
#'
#' @param offset A constant value to be added to survival days for reporting the
#'     median value
#'
#' @param quantiles Quantiles of the composite endpoint to be reported
#'
#' @param ... Options for ranking subjects using the composite endpoint
#' \describe{
#' \item{cut.z}{Clinically meaningful difference in the functional outcome}
#' \item{cut.surv}{Clinically meaningful difference in survival time}}
#'
#' @return A class \code{IDEM.RST} list contains
#' \describe{
#' \item{theta}{A dataset with columns \code{Delta0}, \code{Delta1}, \eqn{\hat{\theta}}}
#' \item{quantiles}{ A dataset with columns \code{Delta}, \code{Trt}, \code{Quantiles}}}
#'
#' @examples
#' \dontrun{
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- fit.model(abc, lst.var);
#' rst.imp <- get.imp.all(abc, rst.fit, lst.var, deltas=c(-0.25,0,0.25),
#'                    normal=TRUE, iter=300, n.imp=2, thin=10, p.scale=10);
#' rst.est <- get.theta.quantiles(rst.imp, lst.var,
#'                                quantiles=c(0.25,0.5,0.75));}
#'
#' @export
#'
get.theta.quantiles <- function(imp.data,
                                lst.var,
                                deltas=NULL,
                                offset=10000,
                                quantiles=0.5,
                                ...) {

    stopifnot(any(class(imp.data) == get.const("IMP.CLASS")));

    vtrt     <- lst.var$trt;
    duration <- lst.var$duration;

    ##trt arms
    atrt <- sort(unique(imp.data[, vtrt]));

    ##not imputed subjects
    sub.noimp   <- imp.data[is.na(imp.data$IMP),];
    ready.noimp <- get.data.ready(sub.noimp, lst.var, atrt);

    ##imputed
    if (all(is.na(imp.data$IMP))) {
        n.imp <- 1;
    } else {
        n.imp <- max(imp.data$IMP, na.rm=TRUE);
    }

    if (is.null(deltas)) {
        deltas <- sort(unique(imp.data$DELTA[which(!is.na(imp.data$DELTA))]));
    }

    rst.median <- NULL;
    rst.rank   <- NULL;
    for (i in 1:n.imp) {
        tmp.lst <- rep(list(NULL), length(deltas));
        for (j in 1:length(deltas)) {
            cur.data     <- subset(imp.data, imp.data$DELTA == deltas[j] & imp.data$IMP == i);
            cur.ready    <- get.data.ready(cur.data, lst.var, atrt);
            tmp.lst[[j]] <- rep(list(NULL), length(atrt));

            for (k in 1:length(atrt)) {
                tmp.lst[[j]][[k]] <- rbind(ready.noimp[[k]], cur.ready[[k]]);
                c.med             <- get.median(tmp.lst[[j]][[k]], duration, quantiles=quantiles, ...);
                rst.median        <- rbind(rst.median,
                                           cbind(i, deltas[j], atrt[k], c.med));
            }
        }

        for (t1 in 1:length(deltas)) {
            for (t2 in 1:length(deltas)) {
                cur.rank <- c.rankall(tmp.lst[[t1]][[1]],
                                      tmp.lst[[t2]][[2]],
                                      duration, ...);
                rst.rank <- rbind(rst.rank,
                                  c(i, deltas[t1], deltas[t2], cur.rank));
            }
        }
    }

    colnames(rst.median) <- c("Imputation", "Delta", "TRT", "Q", "QuantSurv", "QuantY");
    colnames(rst.rank)   <- c("Imputation", "Delta0", "Delta1", "Theta");

    ##median of median
    dfmedian                <- data.frame(rst.median);
    inx                     <- which(dfmedian$MedianSurv <= lst.var$duration);
    dfmedian[inx, 'QuantY'] <- -offset + dfmedian[inx, 'QuantSurv'];
    median.median           <- sqldf('select Delta, Trt, Q,
                                      median(QuantY) as Quant
                                      from dfmedian
                                      group by Delta, TRT, Q');

    ##average rank
    dfrank  <- data.frame(rst.rank);
    avg.rank <- sqldf('select Delta0, Delta1, avg(Theta) as Theta
                       from dfrank group by Delta1, Delta0');

    rst <- list(quantiles=median.median,
                theta=avg.rank,
                raw.theta=rst.rank,
                raw.quantiles=rst.median);

    class(rst) <- get.const("RST.CLASS");
    rst
}


#' Boostrap analysis
#'
#' @inheritParams plotCompleters
#' @inheritParams get.imp.all
#' @inheritParams get.theta.quantiles
#'
#' @param n.boot Number of bootstrap samples
#' @param n.cores Number of cores for parallel computation
#' @param ... parameters for imputation. See \code{\link{get.imp.all}}
#'
#' @return A class \code{IDEM.BOOT} list
#'
#' @examples
#'
#' \dontrun{
#' lst.var  <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                  endp=c("Y2"), unitTime="days",
#'                  trt.label = c("UC+SBT", "SAT+SBT"),
#'                  cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.boot <- get.bs.all(n.boot = 10, n.cores = 5, data.all = abc, lst.var = lst.var,
#'                        deltas = c(-0.25, 0, 0.25), quantiles = c(0.25,0.5,0.75),
#'                        normal=TRUE, iter=300, n.imp=2, thin=10, p.scale=10);}
#'
#' @export
#'
get.bs.all <- function(n.boot,
                       data.all,
                       lst.var,
                       deltas = 0,
                       quantiles = 0.5,
                       n.cores = 1,
                       update.progress=NULL,
                       ...) {

    ##number of cores
    n.cores <- min(n.cores, parallel::detectCores()-1);

    if ("PROGRESS" %in% toupper(class(update.progress)))
        update.progress$set(value=1, detail=paste(""));

    rst <- parallel::mclapply(1:n.boot,
                              function(x) {
        get.boot.single(data.all,
                        lst.var,
                        deltas=deltas,
                        quantiles=quantiles,
                        ...);
    }, mc.cores=n.cores);

    class(rst) <- get.const("BOOT.CLASS");
    rst
}


#' Hypothesis testing
#'
#' Hypothesis testing using the estimation for the original dataset and
#' Summarize Boostrap analysis results
#'
#' @param rst.org A class \code{IDEM.RST} result list from
#'     \code{\link{get.theta.quantiles}} using the original data
#'
#' @param rst.boot A class \code{IDEM.BOOT} result list from
#'     \code{\link{get.theta.quantiles}} using the original data
#'
#' @param quantiles Quantiles for extracting bootstrap confidence intervals
#'
#' @inheritParams list.vars
#'
#' @return A class \code{IDEM.TEST} containing two datasets
#'
#' \describe{
#'
#' \item{theta}{ With columns
#'
#' \itemize{
#'
#' \item \code{Delta0}: Sensitivity parameter for control arm,
#' \item \code{Delta1}: Sensitivity parameter for intervention arm
#' \item \code{Theta}: Estimated \eqn{\theta}
#' \item \code{SD}: Standard deviation
#' \item \code{PValue}: p-value
#' }}
#'
#' \item{quantiles}{With columns
#'
#' \itemize{
#'
#' \item \code{Delta}:Sensitivity parameter
#' \item \code{TRT}:Treatment arm
#' \item \code{Q}: Quantiles of the composite endpoint to be estimated
#' \item \code{Quant}: Estimation
#' \item \code{LB}: Lower bound of the specified confidence interval
#' \item \code{UB}: Upper bound of the specified confidence interval
#' }}
#' }
#'
#' @examples
#' \dontrun{
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- fit.model(abc, lst.var);
#' rst.imp <- get.imp.all(abc, rst.fit, lst.var, deltas=c(-0.25,0,0.25),
#'                        normal=TRUE, iter=300, n.imp=2, thin=10, p.scale=10);
#' rst.est <- get.theta.quantiles(rst.imp, lst.var,
#'                                quantiles=c(0.25,0.5,0.75));
#' rst.boot <- get.bs.all(n.boot = 10, n.cores = 5, data.all = abc,
#'                        lst.var = lst.var, deltas = c(-0.25, 0, 0.25),
#'                        quantiles = c(0.25,0.5,0.75), normal=TRUE,
#'                        iter=300, n.imp=2, thin=10, p.scale=10);
#'
#' rst.final <- get.overall.rst(rst.est, rst.boot);}
#' @export
#'
get.overall.rst <- function(rst.org, rst.boot, quantiles=c(0.025, 0.975)) {

    stopifnot(class(rst.org)  == get.const("RST.CLASS") &
              class(rst.boot) == get.const("BOOT.CLASS"));

    ##bootstrap
    meta  <- rep(list(NULL),2);
    rst   <- rep(list(NULL),2);
    for (i in 1:length(rst.boot)) {
        cur.rank <- rst.boot[[i]]$theta;
        cur.med  <- rst.boot[[i]]$quantiles;

        if (1 == i) {
            meta[[1]] <- cur.rank[, -ncol(cur.rank)];
            meta[[2]] <- cur.med[, -ncol(cur.med)];
        }

        rst[[1]] <- cbind(rst[[1]], cur.rank[,ncol(cur.rank)]);
        rst[[2]] <- cbind(rst[[2]], cur.med[,ncol(cur.med)]);
    }

    rsd <- cbind(meta[[1]], SD=apply(rst[[1]], 1, sd));
    rqs <- cbind(meta[[2]], t(apply(rst[[1]], 1, quantile, quantiles)));
    colnames(rqs)[4:5] <- c("LB", "UB");

    ##original
    orank <- rst.org$theta;
    omed  <- rst.org$quantiles;

    ##rank
    rstrank <- sqldf("select a.*, b.sd
                      from orank a
                      left join rsd b on (a.Delta0 = b.Delta0 and a.Delta1 = b.Delta1)");

    rstrank$PValue<- apply(rstrank,
                           1,
                           function(x) { 2*min(pnorm(x[3],0,x[4]),
                                               1-pnorm(x[3],0,x[4]))});
    ##quantiles
    rstquan <- sqldf("select a.*, b.lb, b.ub
                      from omed a
                      left join rqs b on (a.Delta = b.Delta and
                                          a.TRT   = b.TRT   and
                                          a.Q     = b.Q)");


    rtn <- list(theta=rstrank, quantiles=rstquan);
    class(rtn) <- get.const("TEST.CLASS");

    rtn
}


