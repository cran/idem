#' Treatment effect estimation
#'
#' Estimate treatment effect and median of the composite endpoint from using
#' imputed data
#'
#' @inheritParams imPlotCompleters
#' @inheritParams imPlotImputed
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
#' \item{list.var}{List of parameters}
#' \item{theta}{A dataset with columns \code{Delta0}, \code{Delta1}, \eqn{\hat{\theta}}}
#' \item{quantiles}{ A dataset with columns \code{Delta}, \code{Trt}, \code{Quantiles}}
#' \item{survivor}{A dataset for survivors with columns \code{Delta0}, \code{Delta1}, \code{Mean0},
#' \code{Mean1}, \code{Diff}}}
#'
#' @examples
#' \dontrun{
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- imFit(abc, lst.var);
#' rst.imp <- imImpAll(abc, rst.fit, lst.var, deltas=c(-0.25,0,0.25),
#'                     normal=TRUE, iter=300, n.imp=2, thin=10, p.scale=10);
#' rst.est <- imEstimate(rst.imp, quantiles=c(0.25,0.5,0.75));}
#'
#' @export
#'
imEstimate <- function(imp.rst,
                       quantiles=0.5,
                       ...) {

    if (is.null(imp.rst))
        return(NULL);

    stopifnot(any(class(imp.rst) == get.const("IMP.CLASS")));

    lst.var  <- imp.rst$lst.var;
    imp.data <- imp.rst$complete;
    deltas   <- imp.rst$deltas;
    n.imp    <- imp.rst$n.imp;

    vtrt     <- lst.var$trt;
    duration <- lst.var$duration;
    vsurv    <- lst.var$surv;

    ##trt arms
    atrt <- sort(unique(imp.data[, vtrt]));

    ##not imputed subjects
    sub.noimp   <- imp.data[is.na(imp.data$IMP),,drop=FALSE];
    ready.noimp <- get.data.ready(sub.noimp, lst.var, atrt);

    rst.median   <- NULL;
    rst.rank     <- NULL;
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
    dfmedian       <- data.frame(rst.median);
    dfmedian$Quant <- get.comp(dfmedian$QuantY, dfmedian$QuantSurv, lst.var$duration);
    dfmedian       <- sqldf('select * from dfmedian order by Delta, Trt, Q, Quant');
    inx.median     <- ceiling(n.imp/2);
    median.median  <- dfmedian[seq(inx.median, nrow(dfmedian), by=n.imp),
                               c("Delta", "TRT", "Q", "QuantY", "QuantSurv")];

    inx.2 <- which(median.median$QuantSurv <= lst.var$duration);
    if (length(inx.2) > 0) {
        median.median[inx.2,  "QuantY"]    <- NA;
        median.median[-inx.2, "QuantSurv"] <- NA;
    }

    ##average rank
    dfrank  <- data.frame(rst.rank);
    avg.rank <- sqldf('select Delta0, Delta1, avg(Theta) as Theta
                       from dfrank group by Delta1, Delta0');

    ##survivor functional means
    dalive       <- imp.data[imp.data[,vsurv] > duration, , drop=FALSE];
    rst.survivor <- NULL;
    if (nrow(dalive) > 0) {
        ##there exist survivors
        endp    <- get.const("TXT.ENDP");
        txt.sql <- paste("select delta, trt, avg(endp) as Mean",
                         "from (select distinct delta, ", vtrt, " as trt, ID,",
                         "avg(", endp, ") as endp",
                         "from dalive group by delta, ID) group by delta, trt",
                         sep=" ");
        survtrt    <- sqldf(txt.sql);
        survtrt.t0 <- survtrt[seq(1, nrow(survtrt)-1, 2), c("delta", "Mean")];
        survtrt.t1 <- survtrt[seq(2, nrow(survtrt), 2), c("delta", "Mean")];
        ndelta     <- nrow(survtrt.t0);
        rst.survivor <- cbind(survtrt.t0[ rep(1:ndelta, each = ndelta),],
                              survtrt.t1[ rep(1:ndelta, ndelta),]);

        colnames(rst.survivor) <- c("Delta0", "Mean0", "Delta1", "Mean1");
        rst.survivor <- data.frame(rst.survivor);
        rst.survivor$Diff <- rst.survivor$Mean1 - rst.survivor$Mean0;
    }

    rst <- list(lst.var=lst.var,
                quantiles=median.median,
                theta=avg.rank,
                survivor=rst.survivor,
                raw.theta=rst.rank,
                raw.quantiles=rst.median);

    class(rst) <- get.const("RST.CLASS");
    rst
}


#' Boostrap analysis
#'
#' @inheritParams imEstimate
#'
#' @param n.boot Number of bootstrap samples
#' @param n.cores Number of cores for parallel computation
#' @param update.progress Parameter reserved for run \code{idem} in GUI mode
#'
#' @return A class \code{IDEM.BOOT} list with length \code{n.boot+1}. Each item
#'     in the list is a class \code{{IDEM.RST}} list (see
#'     \code{\link{imEstimate}}). The first item correspons to the
#'     estimation result on the original dataset.
#'
#' @examples
#'
#' \dontrun{
#' lst.var  <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                  endp=c("Y2"), unitTime="days",
#'                  trt.label = c("UC+SBT", "SAT+SBT"),
#'                  cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- imFitModel(abc, lst.var);
#' rst.imp <- imImpAll(abc, rst.fit, deltas=c(-0.25,0,0.25),
#'                     normal=TRUE, chains = 4, iter = 2000, warmup = 1000);
#' rst.boot <- imBs(rst.imp, n.boot = 10, n.cores = 5, quantiles = c(0.25,0.5,0.75));}
#'
#' @export
#'
imBs <- function(imp.rst,
                 n.boot = 100,
                 n.cores = 1,
                 update.progress=NULL,
                 quantiles=0.5
                 ) {

    if (is.null(imp.rst))
        return(NULL);

    stopifnot(any(class(imp.rst) == get.const("IMP.CLASS")));
    stopifnot(!is.null(imp.rst$org.data));

    ##original result
    rst.org <- imEstimate(imp.rst, quantiles=quantiles);

    data.all  <- imp.rst$org.data;
    lst.var   <- imp.rst$lst.var;
    deltas    <- imp.rst$deltas;
    n.imp     <- imp.rst$n.imp;
    stan.par  <- imp.rst$stan.par;
    normal    <- imp.rst$normal;

    n.cores   <- min(n.cores, parallel::detectCores()-1);

    if ("PROGRESS" %in% toupper(class(update.progress)))
        update.progress$set(value=1, detail=paste(""));

    rst <- parallel::mclapply(1:n.boot,
                              function(x) {
                         if ("PROGRESS" %in% toupper(class(update.progress)))
                             update.progress$set(value=x/n.boot,
                                                 detail=paste("Bootstrap", x, sep=" "));

                         get.boot.single(data.all,
                                         lst.var,
                                         deltas=deltas,
                                         n.imp=n.imp,
                                         normal=normal,
                                         stan.par=stan.par,
                                         quantiles=quantiles);
                     }, mc.cores=n.cores);

    rst        <- c(list(rst.org), rst);
    class(rst) <- get.const("BOOT.CLASS");
    rst
}


#' Hypothesis testing
#'
#' Hypothesis testing using the estimation for the original dataset and
#' Summarize Boostrap analysis results
#'
#' @param bs.rst A class \code{IDEM.BOOT} result list from
#'     \code{\link{imBs}} for bootstrap analysis
#'
#' @param quantiles Quantiles for extracting bootstrap confidence intervals
#'
#' @return A class \code{IDEM.TEST} containing two datasets
#'
#' \describe{
#'
#' \item{list.var}{List of parameters}
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
#'                   endp=c("Y2"), unitTime="days",
#'                   trt.label = c("UC+SBT", "SAT+SBT"),
#'                   cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit   <- imFitModel(abc, lst.var);
#' rst.imp   <- imImpAll(abc, rst.fit, deltas=c(-0.25,0,0.25),
#'                       normal=TRUE, chains = 4, iter = 2000, warmup = 1000);
#' rst.boot  <- imBs(rst.imp, n.boot = 10, n.cores = 5);
#' rst.final <- imTest(rst.boot);}
#'
#' @export
#'
imTest <- function(bs.rst, quantiles=c(0.025, 0.975)) {

    stopifnot(class(bs.rst) == get.const("BOOT.CLASS"));

    rst.org  <- bs.rst[[1]];
    rst.boot <- bs.rst[-1];

    ##bootstrap
    meta  <- rep(list(NULL), 2);
    rst   <- rep(list(NULL), 2);
    for (i in 1:length(rst.boot)) {
        cur.rank <- rst.boot[[i]]$theta;
        cur.surv <- rst.boot[[i]]$survivor;

        if (1 == i) {
            meta[[1]] <- cur.rank[, -ncol(cur.rank)];
            meta[[2]] <- cur.surv[, -ncol(cur.surv)];
        }

        rst[[1]] <- cbind(rst[[1]], cur.rank[,ncol(cur.rank)]);
        rst[[2]] <- cbind(rst[[2]], cur.surv[,ncol(cur.surv)]);
    }
    rsd   <- cbind(meta[[1]], SD=apply(rst[[1]], 1, sd));
    rsurv <- cbind(meta[[2]], SD=apply(rst[[2]], 1, sd));

    rst.qs  <- array(NA, dim=c(nrow(rst.boot[[1]]$quantiles),
                               length(rst.boot),
                               2));
    for (i in 1:length(rst.boot)) {
        cur.med     <- rst.boot[[i]]$quantiles;
        rst.qs[,i,] <- as.matrix(cur.med[, c("QuantY", "QuantSurv")]);
    }

    rqs  <- NULL;
    inxs <- ceiling(quantile(1:length(rst.boot), quantiles));
    for (i in 1:nrow(rst.boot[[1]]$quantiles)) {
        cur.qs  <- rst.qs[i,,];
        cur.qsc <- get.comp(cur.qs[,1], cur.qs[,2]);
        cur.ord <- order(cur.qsc);

        cur.q   <- NULL;
        cur.qi  <- NULL;
        for (j in 1:length(inxs)) {
            cur.i <- cur.ord[inxs[j]];
            if (is.na(cur.qs[cur.i,1])) {
                cur.q  <- c(cur.q,  cur.qs[cur.i, 2]);
                cur.qi <- c(cur.qi, 1);
            } else {
                cur.q  <- c(cur.q,  cur.qs[cur.i, 1]);
                cur.qi <- c(cur.qi, 0);
            }
        }
        names(cur.q)  <- paste("BSQ", quantiles*100, sep = "");
        names(cur.qi) <- paste(names(cur.q), "_Surv", sep = "");
        cur.rst       <- c(cur.q, cur.qi);
        rqs           <- rbind(rqs, cur.rst);
    }
    row.names(rqs) <- NULL;
    rqs            <- data.frame(rqs);


    ##original
    orank <- rst.org$theta;
    omed  <- rst.org$quantiles;
    osurv <- rst.org$survivor;

    ##rank
    rstrank <- sqldf("select a.*, b.sd
                      from orank a
                      left join rsd b on (a.Delta0 = b.Delta0 and a.Delta1 = b.Delta1)");

    rstrank$PValue<- apply(rstrank,
                           1,
                           function(x) { 2*min(pnorm(x[3],0,x[4]),
                                               1-pnorm(x[3],0,x[4]))});
    ##quantiles
    rstquan        <- cbind(omed, rqs);

    ##survivors
    rstsurv <- sqldf("select a.delta0, a.delta1, a.diff, b.sd
                      from osurv a
                      left join rsurv b on (a.Delta0 = b.Delta0 and a.Delta1 = b.Delta1)");

    rstsurv$PValue<- apply(rstsurv,
                           1,
                           function(x) { 2*min(pnorm(x[3],0,x[4]),
                                               1-pnorm(x[3],0,x[4]))});

    rtn <- list(theta=rstrank,
                quantiles=rstquan,
                survivor=rstsurv,
                lst.var=rst.org$lst.var);
    class(rtn) <- get.const("TEST.CLASS");
    rtn
}


