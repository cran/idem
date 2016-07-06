
##------------------------------------------------------
##
##           MODEL FITTING
##
##------------------------------------------------------

#' Model fitting
#'
#' Fit linear imputation models to the observed data from complete survivors for
#' each treatment arm at each time point
#'
#' @inheritParams plotCompleters
#'
#' @return A class \code{idem.fit} list of modeling fitting results from each model with the following
#'     items
#'\describe{
#'   \item{lm}{results from function \code{lm}}
#'   \item{formula}{model formula}
#'   \item{coef}{model coefficients}
#'   \item{res}{residuals}
#'   \item{h}{bandwidth of residuals for kernel density estimation}
#' }
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", outcome=c("Y1","Y2"),
#'                 y0=NULL, cov=c("AGE"), bounds=c(0,100));
#'
#' rst.fit <- fit.model(abc, lst.var)
#'
#' @export
#'
#'
fit.model <- function(data.all=NULL, lst.var=NULL) {

    if (is.null(data.all)|is.null(lst.var))
        return(NULL);

    vtrt     <- lst.var$trt;
    voutcome <- lst.var$outcome;
    vy0      <- lst.var$y0;
    vcov     <- lst.var$cov;

    ##transfer data
    data.all <- get.transfer.all(data.all, lst.var);

    ##completers
    is.comp  <- apply(data.all[, voutcome],
                      1,
                      function(x){all(!is.na(x))});

    ##treatment arms
    a.trt <- get.trt(data.all[,vtrt]);
    rst   <- list(NULL);
    for (i in 1:length(a.trt)) {
        cur.data  <- data.all[which(a.trt[i] == data.all[,vtrt] & is.comp),];
        cur.rst   <- list(NULL);
        cur.names <- NULL;
        for (j in 1:length(voutcome)) {
            if (1 == j) {
                prev.y <- NULL;
            } else {
                prev.y <- voutcome[1:(j-1)];
            }
            cur.f    <- paste(voutcome[j], "~", paste(c(prev.y, vy0, vcov), collapse="+"));
            cur.lm   <- lm(as.formula(cur.f), data=cur.data);
            cur.coef <- get.coef(cur.lm);
            cur.res  <- residuals(cur.lm);
            cur.h    <- get.band.h(cur.res);
            cur.rst[[j]] <- list(lm=cur.lm,
                                 formula=cur.f,
                                 coef=cur.coef,
                                 res=cur.res,
                                 h=cur.h);
        }
        rst[[i]] <- cur.rst;
    }

    names(rst) <- a.trt;
    class(rst) <- get.const("FIT.CLASS");;
    rst;
}


#' Impute missing data
#'
#' Impute missing data for all the subjects or a small sample of the subjects
#'
#' @inheritParams plotCompleters
#'
#' @param fit.all A class \code{idem.fit} results of linear regression.
#'     See \code{\link{fit.model}}.
#'
#' @param endponly Logical variable indicating whether clinical outcomes not
#'     used in calculating the final clinical outcome will be imputed. The
#'     default is FALSE, indicating that all missing clinical outcomes will be
#'     imputed sequentially

#' @param trace.n Number of subjects to impute in order to access convergence of
#'     the MCMC chain. If \code{trace.n} is 0, all subjects will be imputed
#'
#' @param deltas Vector of imputation sensitivity parameters
#'
#' @param update.progress Parameter reserved for run \code{idem} in GUI mode
#'
#' @param imputeNone If \code{TRUE}, return subjects that do not need imputation
#'
#' @param ... Options for MCMC sampling including
#'     \describe{
#'     \item{normal}{Logical variable indicating whether normality assumption should be made for the residuals}
#'     \item{n.imp}{number of imputed missing values for each subject}
#'     \item{iter}{burn-in in the MCMC chain}
#'     \item{thin}{thinning after burn-in}
#'     \item{p.scale}{scale factor for the variances of the candidate
#'                    normal distribution in random-walk Metroplis-Hasting
#'     sampling} }
#'
#' @return
#'
#' If \code{imputeNone} is TRUE, return a dataset with the original data for the
#' subset of subjects who died at the end of the study or had no missing outcomes.
#'
#' If \code{trace.n >0}, return a list of MCMC sampling chains of the imputed
#' missing values for \code{trace.n} randomly selected subjects that need
#' imputation.
#'
#' Otherwise, return a class \code{IDEM.IMP} dataset with the original data for the subset of subjects
#' who died at the end of the study or had no missing outcomes and the
#' \code{n.imp} imputed missing outcomes for subjects who need missing value
#' imputation.
#'
#' @examples
#'
#' \dontrun{
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- fit.model(abc, lst.var);
#' rst.imp <- get.imp.all(abc, rst.fit, lst.var, deltas=c(-0.25,0,0.25),
#'                        normal=TRUE, iter=300, n.imp=2, thin=1, p.scale=10);}
#'
#' @export
#'
get.imp.all <- function(data.all,
                        fit.all,
                        lst.var,
                        endponly=TRUE,
                        trace.n=0,
                        deltas=0,
                        ...,
                        update.progress=NULL,
                        imputeNone=FALSE) {

    f.addcols <- function(dset) {
        cbind('ID'=1:nrow(dset),
              'DELTA'=0,
              'IMP'=NA,
              dset);
    }

    stopifnot(class(fit.all) == get.const("FIT.CLASS"));

    if (is.null(deltas)) {
        deltas <- 0;
    }

    voutcome <- NULL;
    eoutcome <- NULL;
    vsurv    <- NULL;
    duration <- NULL;
    endfml   <- NULL;
    tmp.endp <- NULL;

    ##get parameters in current enviroment
    get.para(lst.var, environment());

    ##save org voutcome
    data.all[, paste(get.const("ORG.PREFIX"), voutcome, sep="")] <- data.all[, voutcome];

    ##chk subjects need imputation
    if (endponly) {
        vendp <- eoutcome;
    } else {
        vendp <- voutcome;
    }

    need.imp <- apply(data.all,
                      1,
                      function(x) {
                         rst <- (x[vsurv] > duration) & any(is.na(x[vendp]));
    });
    need.imp <- which(need.imp);

    ##get complete data. return if no imputation needed
    if (0 == length(need.imp) | imputeNone) {
        if (imputeNone & length(need.imp) > 0) {
            rst <- data.all[-need.imp,];
        } else {
            rst <- data.all;
        }

        eval(parse(text=paste("tmp.endp <- with(rst, {", endfml,"})")));
        rst[[get.const("TXT.ENDP")]] <- tmp.endp;
        rst <- f.addcols(rst);
        return(rst);
    } else if (nrow(data.all) == length(need.imp)) {
        data.comp <- NULL;
    } else {
        data.comp <- data.all[-need.imp,];
        data.comp <- f.addcols(data.comp);
        data.comp <- data.frame(data.comp);
    }

    ##sd as scale
    sd.y <- NULL;
    for (i in 1:length(fit.all)) {
        cur.sd <- NULL;
        for (j in 1:length(fit.all[[i]])) {
            cur.sd <- c(cur.sd, fit.all[[i]][[j]]$coef[1]);
            ##cur.sd <- c(cur.sd, 1);
        }
        sd.y <- rbind(sd.y, cur.sd);
    }
    rownames(sd.y) <- names(fit.all);

    ##trace only several subjects
    if (trace.n > 0) {
        need.imp <- sample(need.imp, min(trace.n, length(need.imp)));
        rst      <- list(imp.sub=data.all[need.imp,]);
        rst$mcmc <- rep(list(NULL), length(need.imp));
    } else {
        rst <- data.comp;
    }

    ##actual imputation
    n.imp <- length(need.imp);
    for (i in 1:n.imp) {
        cur.rst <- imp.single(as.matrix(data.all[need.imp[i],]),
                              fit.all, lst.var, sd.y,
                              trace=trace.n>0,
                              deltas=deltas,
                              ...);

        if (trace.n > 0) {
            rst$mcmc[[i]] <- cur.rst;
        } else {
            cur.rst <- cbind('ID'=nrow(data.comp)+i, cur.rst);
            rst     <- rbind(rst, cur.rst);
        }

        if ("PROGRESS" %in% toupper(class(update.progress))) {
            update.progress$set(value=i/n.imp,
                                detail=paste(i, " out of ", n.imp, sep=""));
        } else {
            print(i);
        }
    }

    ##compute endp
    if (trace.n <=0) {
        eval(parse(text=paste("tmp.endp <- with(rst, {", endfml,"})")));
        rst[[get.const("TXT.ENDP")]] <- tmp.endp;
        class(rst) <- c(class(rst), get.const("IMP.CLASS"));
    }

    ##return
    rst
}
