
#' Check parameter specification
#'
#' Check if the \code{idem-parameters} are correctly specified and consistent
#' with the data
#'
#' @inheritParams plotCompleters
#'
#' @param html logic indicator for the format of the error messages
#'
#' @return
#'
#'   Text messages if \code{html=TRUE} .
#'
#' @examples
#'
#' err.lst.var <- list(trt="TRT", outcome=c("Y1","Y2"),
#'                 y0=NULL, endfml="Y3", bounds=c(10,20),
#'                 duration=365);
#'
#' chkPara(err.lst.var, abc);
#'
#' @export
#'
#'
#'
chkPara <- function(lst.var, data.all, html = FALSE) {

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
    if (is.null(lst.var$endfml) | 0 == nchar(gsub("\\s","",lst.var$endfml))) {
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
    if (is.na(lst.var$duration) |
        is.null(lst.var$duration) |
        lst.var$duration <=0)
        err.msg <- ep("Study duration is not a proper positive number");

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
        if (html) {
            rst <- NULL;
        } else {
            rst <- "No error found. \n";
        }
    } else {
        rst <- "Model specification in not valid. Please check the following:";
        if (html) {
            rst <- paste(rst, "<ul>", err.msg, "</ul>");
        } else {
            rst <- paste(rst, err.msg, "\n");
        }
    }

    if (html) {
        return(rst);
    } else {
        cat(rst);
    }
}

#' Plot data of completer
#'
#' Spaghetti plot for subjects alive at the end of the study without missing data
#'
#'
#' @param data.all Original dataset
#'
#' @param lst.var see \code{\link{idem-parameters}}
#'
#' @param fname File name of the result pdf file. If \code{fname} is null,
#'     result pdf file will not be generated
#'
#' @param ... Options for \code{pdf} function
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
#'                 y0=NULL, trt.label = c("UC+SBT", "SAT+SBT"),
#'                 duration=365);
#'
#' plotCompleters(abc, lst.var);
#'
#' @export
#'
#'
#'
plotCompleters <- function(data.all,
                           lst.var,
                           fname=NULL,
                           ...) {

    vsurv     <- NULL;
    duration  <- NULL;
    vy0       <- NULL;
    voutcome  <- NULL;
    vtrt      <- NULL;
    trt.len   <- NULL;
    get.para(lst.var, env=environment());

    ## change completers to survivors
    data.all <- as.matrix(data.all);
    data.all <- data.all[data.all[, vsurv] > duration,]

    all.y <- data.all[,c(vy0, voutcome)];
    ylims <- range(all.y, na.rm=TRUE);
    all.x <- seq(0, duration, length.out=length(c(vy0, voutcome)));
    a.trt <- sort(unique(data.all[,vtrt]));

    times <- as.numeric(substr(voutcome, 2, nchar(voutcome)));

    if (sum(!is.na(times)) == length(voutcome)) {
        if(is.null(vy0)) {
            all.x <- times;
        } else {
            all.x <- c(0, times);
        }
    }

    if (!is.null(fname))
        pdf(file=fname, ...);

    par(cex.lab=0.9,
        mfrow=c(1,length(a.trt)),
        mar=c(4,4,2,2),
        cex.axis=0.9);

    for (i in 1:length(a.trt)) {
        cur.data <- data.all[which(a.trt[i] == data.all[,vtrt]),
                             c(vy0,voutcome)];

        if (is.null(vy0)) {
            t0 <- NULL;
        } else {
            t0 <- 'y0';
        }

        plot(NULL,
             xlim=(range(all.x)), ylim=ylims,
             axes=F,
             xlab="Functional Outcome",
             ylab="Observed Value",
             main=trt.len[i]);

        axis(1, at=all.x, labels=c(t0, paste("y", times, sep="")));
        axis(2, at=round(seq(ylims[1],ylims[2],length=5)))

        for (j in 1:nrow(cur.data)) {
            cur.y <- cur.data[j,];
            cur.x <- all.x
            cols <- 'black'
            if(sum(is.na(cur.y)) == length(all.x)) next
            if(any(is.na(cur.y))){
                cols <- 'purple'
                cur.x <- all.x[-which(is.na(cur.y))]
                cur.y <- cur.y[-which(is.na(cur.y))]
            }
            lines(cur.x, cur.y,type='b',pch=16, col = cols)

        }
        lines(all.x, colMeans(cur.data,na.rm=TRUE),col='red',type='b',lwd=5.5)
    }

    if (!is.null(fname))
        dev.off();

}


#' Generate table of missingness pattern frequencies
#'
#' @inheritParams plotCompleters
#'
#' @return A matrix with frequencies of each missing pattern
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
#'                 y0=NULL, trt.label = c("UC+SBT", "SAT+SBT"),
#'                 duration=365);
#'
#' get.mis.table(abc, lst.var);
#'
#' @export
#'
get.mis.table <- function(data.all, lst.var) {

    vtrt     <- NULL;
    voutcome <- NULL;
    duration <- NULL;
    vsurv    <- NULL;
    trt.len  <- NULL;

    get.para(lst.var, environment());

    data.all <- as.matrix(data.all);
    a.trt    <- get.trt(data.all[,vtrt]);
    mis.pat  <- get.miss.pattern(length(voutcome));

    mis.pat[mis.pat == 1] <- 'Observed';
    mis.pat[mis.pat == 0] <- 'Missing';

    if (is.null(trt.len)) {
        trt.len <- paste(toupper(vtrt), "=", a.trt, sep="");
    }

    rst <- NULL;
    for (i in 1:length(a.trt)) {
        subg      <- data.all[which(a.trt[i] == data.all[,vtrt]),
                              c(vsurv, voutcome)];
        nsub      <- nrow(subg);
        cur.alive <- subg[which(subg[,vsurv] > duration),voutcome];
        n.dead    <- nsub-nrow(cur.alive);

        cur.y <- !is.na(cur.alive);
        inx.y <- table(1+get.compose(cur.y));

        n.p               <- rep(0, nrow(mis.pat));
        names(n.p)        <- rownames(mis.pat);
        n.p[names(inx.y)] <- inx.y;

        char.rst <- sapply( c(n.dead, n.p), function(x) { sprintf("%i (%.0f%%)", x, 100*x/nsub)});
        char.rst <- c(char.rst, nsub);
        rst      <- cbind(rst, char.rst);
    }

    rst <- cbind(rbind("", mis.pat, ""), rst);

    colnames(rst) <- c(voutcome, 'Control', 'Intervention');
    rownames(rst) <- c("Deaths on study",
                       paste("S=", 1:nrow(mis.pat), sep=""),
                       "Total");

    rst
}


#' Plot missing patterns
#'
#' Plot the missing patterns of the observed data
#'
#' @inheritParams plotCompleters
#'
#' @param cols Color of observed and missing values
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", outcome=c("Y1","Y2"),
#'                 trt.label = c("UC+SBT", "SAT+SBT"));
#'
#' plotMisPattern(abc, lst.var);
#'
#' @export
#'
plotMisPattern <- function(data.all,
                           lst.var,
                           cols=c("blue", "gray"),
                           fname=NULL, ...) {

    voutcome <- NULL;
    vtrt     <- NULL;

    get.para(lst.var, environment());

    n.time  <- length(voutcome);
    n.sub   <- max(table(data.all[,vtrt]));
    a.trt   <- sort(unique(data.all[,vtrt]));

    if (is.null(trt.len)) {
        trt.len <- paste(toupper(vtrt), "=", a.trt, sep="");
    }

    if (!is.null(fname))
        pdf(file=fname, ...);

    par(cex.lab=1, mfrow=c(1,length(a.trt)), mar=c(4,2,2,2), cex.axis=1);

    for (i in 1:length(a.trt)) {
        cur.data <- data.all[which(a.trt[i] == data.all[,vtrt]),
                             voutcome];

        cur.data <- cur.data[order(apply(cur.data,1,sum)),];
        plot(NULL, NULL, xlim=c(0.5, n.time+0.5), ylim=c(0, n.sub+8), axes=FALSE,
             xlab="Outcomes", ylab="Subjects", main=trt.len[i]);
        axis(1, at=1:n.time, labels=voutcome);
        box();
        for (j in 1:nrow(cur.data)) {
            for (k in 1:n.time) {
                if (is.na(cur.data[j,k])) {
                    mis <- 1;
                } else {
                    mis <- 0;
                }

                rect(k-0.48, j-1, k+0.48, j, col=cols[mis+1],
                     border=FALSE);
            }
        }
    }
    if (!is.null(fname))
        dev.off();
}


#' Plot survival curves
#'
#' Plot Kaplan-Meier survival curves
#'
#' @inheritParams plotCompleters
#'
#' @param cols Curve colors of the treatment and control arm
#'
#' @examples
#'
#' lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
#'                 y0=NULL, trt.label = c("UC+SBT", "SAT+SBT"),
#'                 duration=365);
#' plotSurv(abc, lst.var);
#'
#' @export
#'
plotSurv <- function(data.all,
                     lst.var,
                     cols=c("black", "blue"),
                     fname=NULL, ...) {

    vsurv    <- NULL;
    vtrt     <- NULL;
    duration <- NULL;
    unitTime <- NULL;

    get.para(lst.var, environment());
    vtime    <- data.all[, vsurv];
    grp      <- data.all[,vtrt];

    vevent   <- rep(1, nrow(data.all));
    vevent[which(duration < vtime)] <- 0;
    vtime[which(duration < vtime)]  <- duration;

    ##by group
    sfit   <- survfit(Surv(vtime, vevent) ~ grp);
    sdif   <- survdiff(Surv(vtime, vevent) ~ grp);
    p.val  <- 1 - pchisq(sdif$chisq, length(sdif$n) - 1);

    if (!is.null(fname))
        pdf(fname, ...);

    par(cex.lab=1.3,las=1)
    ylims <- c(min(sfit$lower),1)
    plot(sfit,
         xlab=paste("Time (", unitTime,")", sep=''),
    		 ylab="Survival Probability",
         ylim=ylims,yaxt='n',
         cex=1, conf.int=F, lty=c(1,1), lwd=2, col=cols,
         mark.time=FALSE,
         main = 'Survival Curves');

    axis(2, at=round(seq(ylims[1],ylims[2],len=5),2),
         labels=100*round(seq(ylims[1],ylims[2],len=5),2))
    box(lty='solid')
    grid(5,5);

    if (is.null(trt.len)) {
        trt.len <- paste(vtrt, "=", levels(as.factor(grp)), sep="");
    }

    legend("topright",
           legend=c('Control','Intervention', sprintf("p-value = %5.3f",p.val)),
           lty=c(1,1,0),
           col=c(cols,'black'),
           bty="n", cex=1.2);

    if (!is.null(fname))
        dev.off();
}




#' Plot density of imputed values
#'
#' Plot density of imputed values and the density of the observed outcomes
#'
#' @param imp.data A class \code{IDEM.IMP} dataframe containing complete data
#'     with relevant missing values imputed. See \code{\link{get.imp.all}}.
#'
#' @param deltas Imputation sensitivity parameter for which to generate the results
#'
#' @param endp If TRUE, plot the densities of the imputed functional outcomes.
#'     Otherwise, plot the densities of the imputed outcomes
#'
#' @param cols \code{plot} options
#'
#' @param ltys \code{plot} options
#'
#' @param xlim \code{plot} options
#'
#' @param ylim \code{plot} options
#'
#' @param mfrow \code{plot} options
#'
#' @param adj \code{density} estimation option
#'
#' @inheritParams plotCompleters
#'
#' @examples
#' \dontrun{
#' rst.imp <- get.imp.all(abc, rst.fit, lst.var, deltas=c(-0.25,0,0.25),
#'                        normal=TRUE, iter=300, n.imp=2, thin=10, p.scale=10);
#'
#' plotImputed(rst.imp, lst.var, deltas=c(-0.25,0,0.25),
#'             xlim=c(0,100), endp=FALSE);}
#'
#' @export
#'
plotImputed <- function(imp.data,
                        lst.var,
                        deltas=0,
                        endp=FALSE,
                        fname=NULL,
                        adj=1.5,
                        cols=c("red","cyan","blue","green","brown"),
                        ltys=rep(1, 6),
                        xlim=NULL,
                        ylim=NULL,
                        mfrow=NULL,
                        ...) {

    if (is.null(imp.data))
        return(NULL);

    stopifnot(any(class(imp.data) == get.const("IMP.CLASS")));

    TXT.ENDP   <- get.const("TXT.ENDP");
    ORG.PREFIX <- get.const("ORG.PREFIX");
    vtrt       <- NULL;
    voutcome   <- NULL;
    endfml     <- NULL;
    get.para(lst.var, environment());

    ##get trt
    a.trt   <- sort(unique(imp.data[,vtrt]));
    n.trt   <- length(a.trt);
    n.delta <- length(deltas);

    if (endp) {
        to.plot <- TXT.ENDP;
    } else {
        to.plot <- voutcome;
    }

    n.y <- length(to.plot);

    ##get densities
    lst.den <- rep(list(NULL), n.trt);
    lst.obs <- rep(list(NULL), n.trt);

    maxy <- 0;
    for (i in 1:n.trt) {
        ##observed
        lst.obs[[i]] <- list(NULL);
        inx          <- which(a.trt[i] == imp.data[,vtrt] & is.na(imp.data$IMP));
        lst.obs[[i]][[TXT.ENDP]] <- density(imp.data[inx, TXT.ENDP], adjust=adj, na.rm=TRUE);

        inx          <- which(a.trt[i] == imp.data[,vtrt] & 0 == imp.data$DELTA);

        for (k in 1:length(voutcome)) {
            cur.y <- imp.data[inx, paste(ORG.PREFIX, voutcome[k], sep="")];
            cur.d <- density(cur.y, adjust=adj, na.rm=TRUE);
            lst.obs[[i]][[voutcome[k]]] <- cur.d;
        }

        ##imputed
        lst.den[[i]] <- rep(list(NULL), n.delta);
        for (j in 1:n.delta) {
            lst.den[[i]][[j]] <- rep(list(NULL), n.y);
            inx               <- which(a.trt[i]  == imp.data[,vtrt] &
                                           !is.na(imp.data$IMP) &
                                           deltas[j] == imp.data$DELTA);

            for (k in 1:n.y) {
                cur.y <- imp.data[inx, to.plot[k]];

	            if(length(cur.y) == 0) break ## AL

                cur.d <- density(cur.y, adjust=adj, na.rm=TRUE);
                maxy  <- max(maxy, cur.d$y);
                lst.den[[i]][[j]][[k]] <- cur.d;
            }
        }
    }

    ##lims
    if (is.null(xlim)) {
        xlim <- range(imp.data[,to.plot], na.rm=TRUE);
    }

    if (is.null(ylim)) {
        ylim <- c(0, maxy*1.1);
    }

    if (is.null(trt.len)) {
        #trt.len <- paste(toupper(lst.var$trt), "=", a.trt, sep="");
        trt.len <- ifelse(a.trt == 0, 'Control','Intervention')  ## AL
    }

    ##outputfile
    if (!is.null(fname)) {
        pdf(fname, ...);
    }

    if (is.null(mfrow)) {
        if (n.y > 1) {
            mfrow <- c(n.trt, n.y);
        } else {
            mfrow <- c(1, n.trt);
        }
    }

    par(mfrow=mfrow, las = 1);

    if (endp) {
        ## xlabs <- "Z (Imputed)";
        xlabs <- endfml; ## AL
    } else {
        xlabs <- voutcome;
    }

    for (i in 1:n.trt) {
        for (k in 1:n.y) {
            plot(NULL,
                 xlab=xlabs[k], ylab="Density",
                 main=trt.len[i],
                 xlim=xlim, ylim=ylim);

            for (j in 1:n.delta) {
                lines(lst.den[[i]][[j]][[k]], col=cols[j], lty=ltys[j], lwd=2);
            }

            ##observed
            lines(lst.obs[[i]][[to.plot[k]]], col="gray", lty=2, lwd=2);

            ql <- as.expression(lapply(deltas,
                function(l) bquote(Delta==.(l))));

            legend("topleft",
                   c(ql, "Observed"),
                   lty=c(ltys[1:length(deltas)], 2),
                   col=c(cols[1:length(deltas)], "gray"),
                   bty="n");
        }
    }

    if (!is.null(fname))
        dev.off();

}


#' Cumulative Plot
#'
#' Generate cumulative plot of the composite survival and functional outcome
#'
#' @inheritParams plotImputed
#'
#' @param at.surv Sets the range of the survival times to plot in the cumulative distribution function.
#' By default the range is the range of survival values up to the duration of the study.
#'
#' @param at.z Sets the range of the functional outcome to plot in the cumulative distribution function.
#'  By defualt this is the range of the functional outcomes plus the buffer amount to improve visibility
#'  in the transition from survival to functional outcome.
#'
#' @param p.death Proportion of the plot width devoted to Survival. By default the
#' cumulative distribution will devote horizontal space to the survival portion that
#' is proportional to the number of subjects who die prior to duration.
#'
#' @param buffer Small horizontal gap used to better visually distinguish the
#' transition from survival to functional outcome.
#'
#' @param delta Imputation sensitivity parameter for which to generate the results
#'
#' @param seg.lab Labels for the two components of the composite outcome.
#'
#' @param main \code{plot} options
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
#'
#' plotComposite(rst.imp, lst.var, delta=0);}
#' @export
#'
plotComposite <- function(imp.data,
                          lst.var,
                          delta=0,
                          buffer=0.05,
                          at.surv=NULL,
                          at.z=NULL,
                          p.death=NULL,
                          seg.lab=c("Survival", "Functional"),
                          fname=NULL,
                          cols=rep(c("cyan", "red"),3),
                          ltys=rep(1, 6),
                          main="",
                          ...) {


    if (is.null(imp.data))
        return(NULL);

    stopifnot(any(class(imp.data) == get.const("IMP.CLASS")));

    f.y  <- function(x) {
        rst <- p.death+buffer + (1-p.death-buffer)*(x - range.y[1])/(range.y[2]-range.y[1]);
    }

    TXT.ENDP <- get.const("TXT.ENDP");
    vtrt     <- NULL;
    duration <- NULL;
    get.para(lst.var, environment());

    a.trt <- sort(unique(imp.data[,vtrt]));
    xlim  <- c(0,1);
    ylim  <- c(0,1);

    if (is.null(trt.len)) {
        trt.len <- paste("TRT=", a.trt, sep="");
    }

    ##convert xaxis
    delta.data <- subset(imp.data, imp.data$DELTA == delta | is.na(imp.data$DELTA));
    avg.data   <- sqldf(paste("select distinct ID, TRT, SURV, DELTA, avg(",
                              TXT.ENDP,
                              ")  as ENDP from 'delta.data' group by ID",
                              sep=""));

    inx.death  <- which(is.na(avg.data$ENDP));
    inx.alive  <- which(!is.na(avg.data$ENDP));
    if (is.null(p.death)) {
        p.death    <- length(inx.death)/nrow(avg.data);
    }
    ##add xaxis that combines surv and y
    avg.data$toplot <- NA;

    ##surv
    range.surv <- range(c(avg.data[inx.death,'SURV'], duration));
    f.surv <- function(x) {
        rst <- p.death * (x - range.surv[1])/(range.surv[2]-range.surv[1]);
    }
    avg.data$toplot[inx.death] <- f.surv(avg.data[inx.death,'SURV']);


    ##endp
    range.y <- range(avg.data[inx.alive, TXT.ENDP]);
    avg.data$toplot[inx.alive] <- f.y(avg.data[inx.alive,TXT.ENDP]);

    if (is.null(at.surv)) {
        at.surv <- range.surv;
    }
    if (is.null(at.z)) {
        at.z <- round(range.y,1);
    }

    ##outputfile
    if (!is.null(fname))
        pdf(fname, ...);

    par(mar=c(5.1,4.1,2.1,2.1),las=1)
    plot(NULL, xlab="", ylab="Percentile", xlim=xlim, ylim=ylim, main=main,
         axes=FALSE);

    box();
    ##axis(1, at=c(p.death/2, p.death+(1-p.death)/2), c("Survival","Functional"), tick=FALSE);
    axis(1, at=f.surv(at.surv), at.surv);
    axis(1, at=f.y(at.z), at.z);
    axis(2, at=seq(0,1,0.25));

    text(c(p.death/2, p.death+(1-p.death)/2),
         c(0.5,0.5),
         seg.lab,
         col="gray",
         cex=1.2);

    lines(c(p.death,p.death), c(-1,2), lwd=2, lty=2, col="gray");

    for (i in 1:length(a.trt)) {
        cur.d  <- subset(avg.data, avg.data$TRT==a.trt[i]);
        toplot <- c(sort(cur.d$toplot), 1);
        y      <- c(0, 1:nrow(cur.d)/nrow(cur.d), 1);
        lines(stepfun(toplot, y), lwd=2, col=cols[i], lty=ltys[i], do.points=FALSE);
    }

    legend("topleft",
           trt.len,
           lty=ltys[1:length(a.trt)],
           col=cols[1:length(a.trt)],
           bty="n");

    if (!is.null(fname))
        dev.off();
}



#' Contour plot of the sensitivity analysis results
#'
#' Generate contour plot of p-values for sensitivity analysis results
#'
#' @param rst.final A class \code{IDEM.TEST} list generated by
#'     \code{\link{get.overall.rst}}
#' @param con.v Levels of contour plot
#' @param nlevels Levels of color scale
#' @param ... Options for \code{filled.contour}
#'
#' @inheritParams plotCompleters
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
#'                        normal=TRUE, iter=300, n.imp=2, thin=10, p.scale=10);
#' rst.est <- get.theta.quantiles(rst.imp, lst.var,
#'                                quantiles=c(0.25,0.5,0.75));
#' rst.boot <- get.bs.all(n.boot = 10, n.cores = 5, data.all = abc,
#'                        lst.var = lst.var, deltas = c(-0.25, 0, 0.25),
#'                        quantiles = c(0.25,0.5,0.75), normal=TRUE,
#'                        iter=300, n.imp=2, thin=10, p.scale=10);
#'
#' rst.final <- get.overall.rst(rst.est, rst.boot);
#' plotContour(rst.final, lst.var, con.v=0.05, nlevels = 30);
#' }
#'
#' @export
#'
plotContour <- function(rst.final, lst.var, con.v=0.05, nlevels=30, ...) {

    stopifnot(class(rst.final) == get.const("TEST.CLASS"));

    trt.len <- NULL;
    get.para(lst.var, environment());

    cur.data <- rst.final$theta;
    alphas   <- sort(unique(cur.data$Delta0));
    ql       <- as.expression(lapply(trt.len, function(l) bquote(.(l)~Delta)));
    nalpha   <- length(alphas);
    rst      <- matrix(NA, nalpha, nalpha);
    for (i in 1:nalpha) {
        for (j in 1:nalpha) {
            c.d      <- subset(cur.data, cur.data$Delta0 == alphas[i] & cur.data$Delta1 == alphas[j]);
            rst[i,j] <- c.d[1,"PValue"];
        }
    }

    par(oma=c(1,0,0,0));
    filled.contour(alphas, alphas, rst, xlab = ql[1], ylab= ql[2],
                   ...,
                   col=grey(seq(1,0,length=nlevels)),
                   plot.axes={axis(side=1, at=alphas);
                              axis(side=2, at=alphas);
                              grid(nx=length(alphas)-1, ny=length(alphas)-1, col="black");
                              contour(alphas, alphas, rst, levels=con.v,
                                      labcex=1.2, lwd=2, add=T, drawlabels=TRUE)
                                 });
}


