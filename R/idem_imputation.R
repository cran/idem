##------------------------------------------------------
##
##           PREPARATION
##
##------------------------------------------------------


#' Check parameter specification
#'
#' Check if the \code{idem-parameters} are correctly specified and consistent
#' with the data
#'
#' @inheritParams imFitModel
#'
#' @param html logic indicator for the format of the error messages
#'
#' @return
#'
#'   \code{NULL} if the specification is correct.
#'
#'   Text messages if \code{html=TRUE}.
#'
#' @examples
#'
#' err.lst.var <- list(trt="TRT", outcome=c("Y1","Y2"),
#'                 y0=NULL, endfml="Y3", bounds=c(10,20),
#'                 duration=365);
#'
#' imChkPars(abc, err.lst.var);
#'
#' @export
#'
#'
#'
imChkPars <- function(data.all, lst.var, html = FALSE) {

    if (is.null(data.all)|is.null(lst.var))
        return(NULL);

    ep <- function(msg) {
        if (!html) {
            rst <- paste(err.msg, "\n----", msg, sep="");
        } else {
            rst <- paste(err.msg, "<li>", msg, "</li>", sep="");
        }
    }

    err.msg  <- NULL;
    if (0 == length(lst.var$trt))
        err.msg <- ep("No treatment specified");
    if (1 < length(lst.var$trt))
        err.msg <- ep("More than one treatment specified");
    if (0 == length(lst.var$surv))
        err.msg <- ep("No survival time specified");
    if (1 < length(lst.var$surv))
        err.msg <- ep("More than one survival time specified");
    if (0 == length(lst.var$outcome))
        err.msg <- ep("No outcome specified");
    if (0 == length(lst.var$endp))
        err.msg <- ep("Endpoint does not involve outcome");

    var.out <- c(lst.var$outcome, lst.var$y0);

    ##endpoints
    if (is.null(lst.var$endfml)) {
        err.msg <- ep("Please specify endpoint");
    } else if (0 == nchar(gsub("\\s","",lst.var$endfml))) {
        err.msg <- ep("Please specify endpoint");
    } else {
        chk.1 <- try({exp.end <- parse(text=lst.var$endfml);})
        if ("try-error" == class(chk.1)) {
            err.msg <- ep(paste("Endpoint error:", chk.1[1], sep=""));
        } else if (1 < length(exp.end)) {
            err.msg <- ep("Endpoint expression contains more than one line");
        } else {
            chk.end <- NULL;
            eval(parse(text=paste("chk.end <- try(with(data.all[, var.out],
                                  {",lst.var$endfml,"}))")));
            if ("try-error" == class(chk.end))
                err.msg <- ep( paste("Endpoint error:", chk.end[1], sep=""));
        }
    }

    ##duration
    if (is.null(lst.var$duration)) {
        err.msg <- ep("Study duration is not specified");
    } else if (is.na(lst.var$duration) | lst.var$duration <=0) {
        err.msg <- ep("Study duration is not a proper positive number");
    }

    ##boundary
    if (!is.null(lst.var$bounds)) {
        if (any(is.na(lst.var$bounds))) {
            err.msg <- ep("Lower or upper bound is not specified");
        } else {
            if (lst.var$bounds[1] > min(data.all[, var.out], na.rm=TRUE))
                err.msg <- ep("Lower bound is bigger than some observed outcomes");

            if (lst.var$bounds[2] < max(data.all[, var.out], na.rm=TRUE))
                err.msg <- ep("Upper bound is bigger than some observed outcomes");
        }
    }

    ##return
    if (is.null(err.msg)) {
        return(NULL);
    } else {
        rst <- "Model specification in not valid. Please check the following:";
        if (html) {
            rst <- paste(rst, "<ul>", err.msg, "</ul>");
        } else {
            rst <- paste(rst, err.msg, "\n");
        }
    }

    rst
}


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
#' @param data.all Original dataset
#' @param lst.var see \code{\link{idem-parameters}}
#'
#' @return A class \code{IDEM.FIT} list of modeling fitting results with the following items
#' \describe{
#' \item{lst.var}{List of parameters}
#' \item{rst.mdl}{A list of modeling fitting results for each model with
#'     \describe{
#'       \item{lm}{results from function \code{lm}}
#'       \item{formula}{model formula}
#'       \item{coef}{model coefficients}
#'       \item{res}{residuals}
#'       \item{h}{bandwidth of residuals for kernel density estimation}}
#' }}
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#'
#' rst.fit <- imFitModel(abc, lst.var);
#'
#' @export
#'
#'
imFitModel <- function(data.all=NULL, lst.var=NULL) {

    if (is.null(data.all)|is.null(lst.var))
        return(NULL);

    cpara <- imChkPars(data.all, lst.var);
    if (!is.null(cpara)) {
        cat(cpara);
        return(NULL);
    }

    vtrt     <- NULL;
    voutcome <- NULL;
    eoutcome <- NULL;
    vsurv    <- NULL;
    duration <- NULL;
    endfml   <- NULL;
    tmp.endp <- NULL;
    vy0      <- NULL;
    bounds   <- NULL;


    ##get parameters in current enviroment
    get.para(lst.var, environment());

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

    rtn.rst <- list(lst.var = lst.var,
                    rst.mdl = rst);

    class(rtn.rst) <- get.const("FIT.CLASS");;

    rtn.rst;
}


#' Get subjects that need imputation
#'
#' Get the index of subjects in a dataset that need imputation, i.e. survivors with
#' functional endpoint missing
#'
#' @inheritParams imFitModel
#'
#' @param endponly Logical variable indicating whether clinical outcomes not
#'     used in calculating the final clinical outcome will be imputed. The
#'     default is FALSE, indicating that all missing clinical outcomes will be
#'     imputed sequentially
#'
#' @return
#'
#' Vector of indices of subjects that need imputation
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' inx.imp <- imNeedImp(abc, lst.var);
#'
#' @export
#'
imNeedImp <- function(data.all, lst.var, endponly=TRUE) {

    cpara <- imChkPars(data.all, lst.var);
    if (!is.null(cpara)) {
        cat(cpara);
        return(NULL);
    }

    voutcome <- NULL;
    eoutcome <- NULL;
    vsurv    <- NULL;
    duration <- NULL;
    endfml   <- NULL;
    tmp.endp <- NULL;

    ##get parameters in current enviroment
    get.para(lst.var, environment());

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
}


#' Impute missing data under benchmark assumption
#'
#' Call STAN model to impute missing data for an individual subject under
#' benchmark assumption
#'
#' @inheritParams imFitModel
#' @inheritParams imImpAll
#'
#' @param dsub original individual subject data
#' @param chains STAN parameter. Number of Markov chainsm
#' @param iter STAN parameter. Number of iterations
#' @param warmup STAN parameter. Number of burnin.
#' @param control STAN parameter. See \code{rstan::stan} for details.
#' @param ... other options to call STAN sampling such as \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.

#' @return
#'
#' \code{NULL} if there is no missing data for the current subject.
#'
#' Otherwise, return a class \code{IDEM.IMPSUB} that contains a list with two components
#' \describe{
#'     \item{dsub}{original data of the subject}
#'     \item{rst.stan}{A \code{stan.fit} class result returned from \code{rstan::sampling}}
#'     \item{complete}{A dataframe with complete data for the selected subject}
#' }
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- imFitModel(abc, lst.var);
#' rst.imp <- imImpSingle(abc[1,], rst.fit, chains = 4, iter = 2000, warmup = 1000);
#' rstan::traceplot(rst.imp$rst.stan, "YMIS");
#'
#' @export
#'
imImpSingle <- function(dsub, fit.rst, normal=TRUE,
                        chains = 4, iter = 5000, warmup = 1000, control = list(adapt_delta=0.95),
                        ...) {

    stopifnot(class(fit.rst) == get.const("FIT.CLASS"));
    stopifnot(1 == nrow(dsub));

    lst.var <- fit.rst$lst.var;
    fit.all <- fit.rst$rst.mdl;

    vtrt     <- NULL;
    voutcome <- NULL;
    eoutcome <- NULL;
    vsurv    <- NULL;
    duration <- NULL;
    endfml   <- NULL;
    tmp.endp <- NULL;
    vy0      <- NULL;
    bounds   <- NULL;

    ##get parameters in current enviroment
    get.para(lst.var, environment());

    ## transfer data
    csub <- as.matrix(dsub);
    for (i in 1:length(voutcome)) {
        csub[, voutcome[i]] <- get.transfer(csub[, voutcome[i]], bounds);
    }

    outcome  <- csub[1, voutcome];
    y0       <- csub[1, vy0];
    trt      <- csub[1, vtrt];
    vx       <- csub[1, c(vy0, vcov)];

    if (any(is.na(vx)))
        return(NULL);

    fit.trt  <- fit.all[[as.character(trt)]];

    ##missing voutcome that needs imputation
    INX.MIS <- which(is.na(outcome));
    INX.OBS <- which(!is.na(outcome));
    if (0 == length(INX.MIS)) {
        return(NULL);
    }

    ##stan data
    IMIS         <- as.numeric(is.na(outcome));
    INX          <- rep(NA, length(outcome));
    INX[INX.MIS] <- 1:length(INX.MIS);
    INX[INX.OBS] <- 1:length(INX.OBS);

    YOBS <- 999999; ##dummy
    if (0 < length(INX.OBS)) {
        YOBS <- c(outcome[INX.OBS], YOBS);
    }

    COEF         <- NULL;
    RESIDUAL     <- NULL;
    H            <- NULL;
    for (i in 1:length(fit.trt)) {
        RESIDUAL <- cbind(RESIDUAL, fit.trt[[i]]$res);
        H        <- c(H, fit.trt[[i]]$h);
        cur.coef <- fit.trt[[i]]$coef;
        if (1 == i) {
            cur.coef <- c(cur.coef[1:2], 0, cur.coef[-(1:2)]);
        }
        COEF <- rbind(COEF, cur.coef);
    }

    list.stan <- list(NY   = length(outcome),
                      NOBS = length(INX.OBS),
                      YOBS = as.array(YOBS),
                      NX   = length(vx),
                      X    = as.array(vx),
                      IMIS = IMIS,
                      INX  = INX,
                      COEF = COEF,
                      ASSUMENORMAL = as.numeric(normal),
                      RESIDUAL     = RESIDUAL,
                      NRES         = nrow(RESIDUAL),
                      H            = H);
    ##stan sampling
    stan.rst <- sampling(stanmodels[["idem"]],
                         data=list.stan,
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         control = control,
                         ...);

    ##get samples
    ymis <- rstan::extract(stan.rst, "YMIS")$YMIS;
    ymis <- get.inv.transfer(ymis, bounds);

    ycomplete <- dsub[rep(1,nrow(ymis)),];
    for(j in 1:ncol(ymis)) {
        ycomplete[,voutcome[INX.MIS[j]]] <- ymis[,j];
    }

    ##compute endpoint
    eval(parse(text=paste("tmp.endp <- with(ycomplete, {", endfml,"})")));
    ycomplete[[get.const("TXT.ENDP")]] <- tmp.endp;


    ##return
    rst <- list(dsub      = dsub,
                rst.stan  = stan.rst,
                complete  = ycomplete);

    class(rst) <- get.const("BENCH.CLASS");

    rst
}



#' Impute missing data
#'
#' Impute missing data for all the subjects or a small sample of the subjects
#'
#' @inheritParams imFitModel
#' @inheritParams imNeedImp
#'
#' @param fit.rst A class \code{IDEM.FIT} results of linear regression.
#'     See \code{\link{imFitModel}}.
#' @param normal Logical variable indicating whether normality assumption should
#'     be made for the residuals
#' @param n.imp Number of complete datasets required
#' @param deltas Vector of imputation sensitivity parameters
#' @param update.progress Parameter reserved for run \code{idem} in GUI mode
#' @param imputeNone If \code{TRUE}, return subjects that do not need imputation
#' @param ... options to call STAN sampling. These options include
#'     \code{chains}, \code{iter}, \code{warmup}, \code{thin}, \code{algorithm}.
#'     See \code{rstan::sampling} for details.
#'
#'
#' @return
#'
#' If \code{imputeNone} is TRUE, return a dataset with the original data for the
#' subset of subjects who died at the end of the study or had no missing outcomes.
#'
#' Otherwise, return a class \code{IDEM.IMP} list with components
#' \describe{
#' \item{lst.var}{List of parameters}
#' \item{complete}{A dataset with  the original data for
#' the subset of subjects who died at the end of the study or had no missing
#' outcomes and the \code{n.imp} imputed missing outcomes for subjects who need
#' missing value imputation.
#' }
#' \item{n.imp}{Number of imputed complete datasets}
#' \item{deltas}{Imputation sensitivity parameters}
#' \item{org.data}{Original dataset}
#' \item{normal}{Normal assumption for the imputation}
#' \item{stan.par}{parameters in \code{...}}
#' }
#'
#' @examples
#'
#' \dontrun{
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
#'                 endp=c("Y2"), unitTime="days",
#'                 trt.label = c("UC+SBT", "SAT+SBT"),
#'                 cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));
#' rst.fit <- imFitModel(abc, lst.var);
#' rst.imp <- imImpAll(abc, rst.fit, deltas=c(-0.25,0,0.25),
#'                     normal=TRUE, chains = 2, iter = 2000, warmup = 1000);}
#'
#' @export
#'
imImpAll <- function(data.all,
                     fit.rst,
                     normal=TRUE,
                     n.imp=5,
                     endponly=TRUE,
                     deltas=0,
                     update.progress=NULL,
                     imputeNone=FALSE,
                     ...) {

    stopifnot(class(fit.rst) == get.const("FIT.CLASS"));

    f.addcols <- function(dset) {
        cbind('ID'=1:nrow(dset),
              'DELTA'=0,
              'IMP'=NA,
              dset);
    }

    lst.var <- fit.rst$lst.var;
    fit.all <- fit.rst$rst.mdl;

    if (is.null(deltas)) {
        deltas <- 0;
    } else {
        if (!(0 %in% deltas))
            deltas <- c(deltas, 0);
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

    ##subjects that need imputation
    need.imp <- imNeedImp(data.all, lst.var, endponly=endponly);

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
        class(rst) <- c(class(rst), get.const("IMP.CLASS"));
        return(rst);
    } else if (nrow(data.all) == length(need.imp)) {
        data.comp <- NULL;
    } else {
        data.comp <- data.all[-need.imp,];
        data.comp <- f.addcols(data.comp);
        data.comp <- data.frame(data.comp);
        eval(parse(text=paste("tmp.endp <- with(data.comp, {", endfml,"})")));
        data.comp[[get.const("TXT.ENDP")]] <- tmp.endp;
    }

    ## start with all subjects with complete information
    rst    <- data.comp;
    n.need <- length(need.imp);
    for (i in 1:n.need) {
        cur.bench <- imImpSingle(data.all[need.imp[i],],
                                 fit.rst,
                                 normal=normal, ...);
        cur.rst <- imp.exponential(cur.bench, deltas=deltas, n.imp=n.imp);
        cur.rst <- cbind('ID'=nrow(data.comp)+i, cur.rst);
        rst     <- rbind(rst, cur.rst);

        if ("PROGRESS" %in% toupper(class(update.progress))) {
            update.progress$set(value=i/n.need,
                                detail=paste(i, " out of ", n.need, sep=""));
        } else {
            print(i);
        }
    }

    ##return
    rtn.rst <- list(lst.var  = lst.var,
                    deltas   = deltas,
                    normal   = normal,
                    org.data = data.all,
                    n.imp    = n.imp,
                    stan.par = list(...),
                    complete = rst);

    class(rtn.rst) <- c(class(rtn.rst), get.const("IMP.CLASS"));
    rtn.rst
}
