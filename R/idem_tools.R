##------------------------------------------------------
##
##           FUNCTIONS
##
##------------------------------------------------------

##get coef from lm
get.coef <- function(reg.rst) {
    c(summary(reg.rst)$sigma, coef(reg.rst));
}

##get constants
get.const <- function(cname) {
    switch(cname,
           ORG.PREFIX = "ORG",
           TXT.ENDP   = "ENDP",
           FIT.CLASS  = "IDEM.FIT",
           IMP.CLASS  = "IDEM.IMP",
           RST.CLASS  = "IDEM.RST",
           BOOT.CLASS = "IDEM.BOOT",
           TEST.CLASS = "IDEM.TEST",
           "default"
           )
}

##parse endpoint formula
get.parsed.endfml <- function(endfml, data.name="d.frame") {
    if (is.null(endfml))
        return(NULL);

    parse(text=paste("with(", data.name, ", {",endfml,"})"));
}


##get parameters in a list to the current environment
get.para <- function(lst.var, env=environment()) {
    if (!is.environment(env))
        return();

    env$vtrt     <- lst.var$trt;
    env$voutcome <- lst.var$outcome;
    env$vy0      <- lst.var$y0;
    env$vcov     <- lst.var$cov;
    env$duration <- lst.var$duration;
    env$vsurv    <- lst.var$surv;
    env$endfml   <- lst.var$endfml;
    env$unitTime <- lst.var$unitTime;
    env$trt.len  <- lst.var$trt.label;
    env$bounds   <- lst.var$bounds;

    ##parsed endpoint formula
    env$parsed.endfml <- lst.var$parsed.endfml;
    if (is.null(env$parsed.endfml))
        env$parsed.endfml <- get.parsed.endfml(lst.var$endfml);

    if (is.null(env$unitTime))
        env$unitTime <- "Days";

    if (is.null(env$trt.len))
        env$trt.len <- c('Control', 'Intervention');

    if (!is.null(env$endfml)
        & !is.null(env$voutcome)) {
        tmp <- which(sapply(env$voutcome, grepl, env$endfml));
        if (length(tmp) > 0)
            env$eoutcome <- env$voutcome[tmp];
    }

}


##decompose numbers
get.decompose <- function(numbers, digits, base=2) {

    if (digits <=0 ) return(NULL);

    numbers <- as.matrix(numbers);

    if (1 == length(base)) {
        base <- base*array(1, digits-1);
    }

    ##bs order: [100 10 0]
    bs <- NULL;
    for (i in 1:(digits-1)) {
        bs <- c(bs, prod(base[1:(digits-i)]));
    }

    ngroup <- numbers;
    value  <- NULL;
    for (j in 1:dim(numbers)[1]) {
        cur.n      <- numbers[j];
        cur.value  <- NULL;

        if (digits >= 2) {
            for (i in 1:(digits-1)) {
                curw  <- bs[i];
                curg  <- floor(cur.n / curw);
                cur.n <- cur.n - curg * curw;
                cur.value <- c(cur.value, curg);
            }
        }
        cur.value <- c(cur.value, cur.n);
        value     <- rbind(value, cur.value);
    }

    value
}

##inverse of decompose
get.compose <- function(numbers, base=2) {

    numbers <- as.matrix(numbers);
    digits  <- ncol(numbers);

    ##bs order: [100 10 0]
    bs <- NULL;
    for (i in 1:digits) {
        bs <- c(2^(i-1), bs);
    }

    rst <- apply(numbers, 1, function(x){sum(x*bs)});
    rst;
}

##get trt groups
get.trt <- function(trts) {
    rst  <- sort(unique(trts));
}

##get missing pattern
##nt:number of time points
##mono.first: put monotone pattern at first
get.miss.pattern <- function(nt, mono.first=TRUE) {
    all.pattern <- get.decompose((2^nt-1):0, nt);
    rownames(all.pattern) <- 2^nt:1;

    if (mono.first & nt>1) {
        mono <- sapply(1:(nt-1), function(x){c(rep(1, nt-x), rep(0, x))})
        mono <- rbind(rep(1, nt), t(mono));
        inx.mono <- 2^nt - get.compose(mono);
        rst <- all.pattern[inx.mono,];
        rst <- rbind(rst, all.pattern[-inx.mono,]);
    } else {
        rst <- all.pattern;
    }

    rst
}

##read all data
get.all.data <- function(fname) {
    data.all        <- read.table(fname, header=TRUE);
    names(data.all) <- toupper(names(data.all));
    data.all$pid    <- 1:nrow(data.all);
    data.all
}


##transfer data to -inf to inf
get.transfer <- function(x, bounds) {
    ux  <- (x - bounds[1])/(bounds[2] - bounds[1]);
    rst <- log(ux/(1-ux));
    rst
}

##get jacobian for data transfermation

get.jacob <- function(x, bounds) {
    if (is.null(bounds)) {
        dphi <- 1/(x-bounds[1]) + 1/(bounds[2]-x);
    } else {
        dphi <- 1;
    }
    rst  <- prod(dphi);
}

##transfer back to original scale
get.inv.transfer <- function(x, bounds){
    if (is.null(bounds)) {
        rst <- x;
    } else {
        ux  <- exp(x)/(1+exp(x));
        rst <- bounds[1] + ux * (bounds[2] - bounds[1]);
    }
    rst
}

##transfer all data set
get.transfer.all <- function(data.all, lst.var, bounds=lst.var$bounds, inverse=FALSE) {
    if (is.null(bounds))
        return(data.all);

    ally <- lst.var$outcome;
    for (i in 1:length(ally)) {
        if (!inverse) {
            data.all[, ally[i]] <- get.transfer(data.all[, ally[i]], bounds);
        } else {
            data.all[, ally[i]] <- get.inv.transfer(data.all[, ally[i]], bounds);
        }
    }

    data.all
}


##------------------------------------------------------
##
##           Multivariate Normal
##
##------------------------------------------------------

# Returns conditional mean and variance of x[req.ind]
# Given x[given.ind] = x.given
# where X is multivariate Normal with
# mean = mu and covariance = sigma
#
get.cond.Normal <- function(mu, sigma, vec.y) {
    ##missing
    req.ind <- which(is.na(vec.y));

    if (length(mu) == length(req.ind)) {
        rst <- list(mu=mu, sigma=sigma);
    } else if (0 == length(req.ind)) {
        rst <- NULL;
    } else {
        given.ind <- which(!is.na(vec.y));
        B         <- sigma[req.ind, req.ind];
        C         <- sigma[req.ind, given.ind, drop=FALSE];
        D         <- sigma[given.ind, given.ind];
        CDinv     <- C %*% solve(D);
        cMu       <- c(mu[req.ind] + CDinv %*% (vec.y[given.ind] - mu[given.ind]));
        cVar      <- B - CDinv %*% t(C);
        rst       <- list(mu=cMu, sigma=cVar);
    }

    ##return
    rst
}


##get joint normal from the sequential model fitting results
get.joint.Normal.sigma <- function(fit.trt) {
    cur.trt <- fit.trt;
    ny      <- length(cur.trt);
    sigma   <- array(NA, dim=c(ny, ny));

    for (i in 1:ny) {
        cur.coef   <- cur.trt[[i]]$coef;
        sigma[i,i] <- cur.coef[1]^2;
        if (1 == i) next;

        py <- cur.coef[3:(3+i-1)];

        for (j in 1:(i-1)) {
            sigma[i,i] <- sigma[i,i] + py[j]^2*sigma[j,j];
            sigma[i,j] <- sigma[j,i] <- sum(py * sigma[1:(i-1), j]);
        }
    }

    sigma;
}

##get joint normal mus
get.joint.Normal.mu <- function(cov.x, fit.trt) {
    cur.trt <- fit.trt;
    ny      <- length(cur.trt);
    mu      <- rep(NA, ny);

    for (i in 1:ny) {
        cur.coef <- cur.trt[[i]]$coef;
        if (1 == i) {
            prev.y <- NULL;
        } else {
            prev.y <- mu[1:(i-1)];
        }
        mu[i] <- sum(cur.coef[-1] * c(1, prev.y, cov.x));
    }
    mu;
}

##compute kernel density band
get.band.h <- function(res) {
    rst <- 1.06*sd(res)*(length(res)^(-0.2));
}

##kernel density pdf
c.kdpdf <- function(err, res, h=NULL, log.v=TRUE) {

    if (is.null(h)) {
        h <- get.band.h(res);
    }

    tmp <- 0;
    rst.c <- .C("kdpdf",
                as.double(err), as.double(res),
                as.double(h),   as.integer(length(res)),
                as.double(tmp));
    rst <- rst.c[[5]];
    if (log.v) {
        rst <- log(rst);
    }

    rst
}

##empirical distributions
get.loglike <- function(org.all.y, cov.x, fit.trt, bounds, normal=TRUE) {
    ##out of boundary checking
    if (!is.null(bounds)) {
        if (any(org.all.y < bounds[1] | org.all.y > bounds[2] | is.na(org.all.y))) ## AL 06/08/2016
            return(-Inf);
    }

    ##transfer data
    all.y <- get.transfer(org.all.y, bounds);
    ny  <- length(fit.trt);
    rst <- 0;
    for (i in 1:ny) {
        cur.coef <- fit.trt[[i]]$coef;
        cur.res  <- fit.trt[[i]]$res;
        cur.h    <- fit.trt[[i]]$h;

        if (1 == i) {
            prev.y <- NULL;
        } else {
            prev.y <- all.y[1:(i-1)];
        }

        cur.mean <- sum(cur.coef[-1] * c(1, prev.y, cov.x));
        cur.err  <- all.y[i] - cur.mean;

        if (normal) {
            cur.l <- dnorm(cur.err, 0, cur.coef[1], log=TRUE);
        } else {

            ##cur.l <- predict(fit.trt[[i]]$kde, x=cur.res);
            cur.l <- c.kdpdf(cur.err, cur.res, cur.h);
        }

        rst <- rst + cur.l;
    }

    ##jabobian
    cur.jacob <- get.jacob(org.all.y, bounds);
    rst       <- rst + log(cur.jacob);

    ##return
    rst
}


##------------------------------------------------------
##
##           Metroplist Hasting
##
##------------------------------------------------------

##objective log likelihood function for sampling
smp.log.like <- function(miss.y, all.y, base.y, cov.x,
                         fit.trt,
                         bounds,
                         parsed.endfml,
                         delta=0,
                         normal=TRUE) {

    all.y[which(is.na(all.y))] <- miss.y;
    ll <- get.loglike(all.y, cov.x, fit.trt, bounds, normal);
    ##exponential tilting
    if (delta != 0 & !is.infinite(ll)) {
        ##name of d.frame can not be changed
        ##its used in endfml
        d.frame <- data.frame(rbind(c(all.y, base.y)));
        endp    <- eval(parsed.endfml);
        ll      <- ll + delta * endp;
    }


    ##return
    ll
}


##random walk metroplis hasting
##m is number of iteration
##scale is the size of random walk
rwmh.onestep <- function(log.fx, cur.par, scale, m=1, cur.pos=NULL, ...) {

    if (is.null(cur.pos)) {
        cur.pos <- log.fx(cur.par, ...);
    }

    accept  <- 0;
    rst.b0  <- array(NA, dim=c(m,length(cur.par)));

    for (i in 1:m) {
        next.b0          <- cur.par + scale*rnorm(length(scale));
        next.par         <- cur.par;
        next.par         <- next.b0;
        next.pos         <- log.fx(next.par, ...);
        prob             <- exp(next.pos - cur.pos);

        if (runif(1) < prob) {
            cur.par <- next.par;
            cur.pos <- next.pos;
            accept  <- accept + 1;
        }

        rst.b0[i,] <- cur.par;
    }

    list(par=rst.b0, naccept=accept, ll=cur.pos);
}

## random walk metroplis hasting
rwmh.sample <- function(log.fx,
                        init=NULL,
                        iteration=2000,
                        best.scale=NULL,
                        thin=1,
                        burnin=NULL,
                        p.scale=0.25,
                        p.m=5, ...){

    npars      <- length(init);
    cur.par    <- init;

    if (is.null(best.scale)) {
        best.scale <- array(p.scale, dim=c(1, npars));
    }

    if (is.null(burnin)) {
        burnin <- floor(iteration/2);
    }

    ##burin process
    rstall.mcmc <- array(NA, dim=c(iteration, length(init)));

    ##temporarilty keep the likelihood
    nxt.ll      <- NULL;
    for (i in 1:iteration) {
        rb0             <- rwmh.onestep(log.fx, cur.par, m=1, scale=best.scale, cur.pos=nxt.ll, ...);
        nxt.par         <- rb0$par;
        nxt.ll          <- rb0$ll;
        cur.par         <- nxt.par;
        rstall.mcmc[i,] <- nxt.par;
        ##}
    }

    ##burnin
    rstall.mcmc <- rstall.mcmc[-(1:burnin),,drop=FALSE];

    ##thin
    rstall.mcmc <- rstall.mcmc[seq(thin, by=thin,
                                   length.out=floor((iteration-burnin)/thin)),,drop=FALSE];

    rstall.mcmc;
}



##------------------------------------------------------
##
##           Imputation
##
##------------------------------------------------------

##impute individual
##p.scale: adjustment factor for metroplis random walking step. using sd(y) as the base
imp.single <- function(csub,
                       fit.all,
                       lst.var,
                       sd.y,
                       trace=FALSE,
                       normal=TRUE,
                       deltas=0, n.imp=5, iter=2000, thin=10, p.scale=0.25) {

    stopifnot(class(fit.all) == get.const("FIT.CLASS"));

    voutcome      <- NULL;
    vy0           <- NULL;
    vtrt          <- NULL;
    bounds        <- NULL;
    parsed.endfml <- NULL;

    get.para(lst.var, environment());

    outcome  <- csub[1, voutcome];
    y0       <- csub[1, vy0];
    trt      <- csub[1, vtrt];
    vx       <- csub[1, c(vy0, vcov)];
    fit.trt  <- fit.all[[as.character(trt)]];
    scales   <- sd.y[as.character(trt),]*p.scale;

    names(y0)      <- vy0;
    names(outcome) <- voutcome;

    ##---get trace plot---
    if (trace) {
        deltas <- deltas[1];
        n.imp  <- 1;
        burnin <- 1;
        thin   <- 1;
    } else {
        ##only keep numer of imp samples
        burnin <- iter;
    }

    ##missing voutcome that needs imputation
    inx.1 <- which(is.na(outcome));
    if (0 == length(inx.1)) {
        return(NULL);
    }

    ##initial values using the mean
    mean.y <- get.inv.transfer(get.joint.Normal.mu(vx, fit.trt),
                               bounds);

    ##iterations
    n.iter <- iter + n.imp*thin;
    ##imputation
    n.delta <- length(deltas);
    imp.all <- NULL;
    for (l in 1:n.delta) {
        cur.smp <- rwmh.sample(smp.log.like,
                               init=mean.y[inx.1],
                               iteration=n.iter,
                               best.scale=scales[inx.1],
                               burnin=burnin,
                               thin=thin,
                               all.y=outcome, base.y=y0, cov.x=vx,
                               fit.trt=fit.trt,
                               bounds=bounds,
                               parsed.endfml=parsed.endfml,
                               delta=deltas[l],
                               normal=normal);
        imp.all <- rbind(imp.all, cur.smp);
    }

    ##return
    if (trace) {
        rst <- coda::mcmc(imp.all);
    } else {
        rst  <- cbind("DELTA"=rep(deltas,  each=n.imp),
                      "IMP"=rep(1:n.imp, length(deltas)),
                      csub[rep(1, length(deltas)*n.imp), , drop=FALSE]);
        rst[, voutcome[inx.1]] <- imp.all;
    }
    rst
}



##------------------------------------------------------
##
##           RANK ANALYSIS
##
##------------------------------------------------------
c.rankij <- function(val.i, val.j, duration, cut.surv=0, cut.z=0) {
    tmp <- 0;
    rst.c <- .C("rankij",
                as.double(val.i[1]), as.double(val.i[2]),
                as.double(val.j[1]), as.double(val.j[2]),
                as.double(duration),
                as.double(cut.surv),
                as.double(cut.z),
                as.integer(tmp));
    rst <- rst.c[[8]];
    rst
}

##get the rank test statistic \theta
c.rankall <- function(val.trt1, val.trt2, duration, cut.surv=0, cut.z=0) {
    tmp <- 0;
    rst.c <- .C("rankall",
                as.double(t(val.trt1)), as.double(t(val.trt2)),
                as.integer(nrow(val.trt1)), as.integer(nrow(val.trt2)),
                as.double(duration), as.double(cut.surv), as.double(cut.z),
                as.double(tmp));
    rst <- rst.c[[8]];
    rst
}

##bubble sort for ordering subjects from one trt to get median
c.bubblesort <- function(val.trt, duration, cut.surv=0, cut.z=0) {

    if (1 == nrow(val.trt))
        return(c(val.trt,1));

    v <- cbind(val.trt, 1:nrow(val.trt));
    rst.c <- .C("bsort",
                as.double(t(v)),
                as.integer(nrow(v)),
                as.double(duration),
                as.double(cut.surv),
                as.double(cut.z));
    rst <- rst.c[[1]];
    dim(rst) <- c(3, nrow(val.trt));
    t(rst);
}


##get median
get.median <- function(val.trt, duration, quantiles=0.5, ...) {
    sort.val <- c.bubblesort(val.trt, duration, ...);
    med.inx  <- quantile(1:nrow(sort.val), probs=quantiles);

    cbind(quantiles,
          sort.val[ceiling(med.inx), 1:2])
}

##get rank and median from complete (after imputation) dataset
get.data.ready <- function(data.full, lst.var, atrt) {

    if (0 == nrow(data.full))
        return(NULL);

    vtrt     <- lst.var$trt;
    duration <- lst.var$duration;
    vsurv    <- lst.var$surv;

    endp <- get.const("TXT.ENDP");

    ##deceased
    inx.d <- which(data.full[,vsurv] <= duration);
    if (length(inx.d) > 0) {
        data.dead <- data.full[inx.d, c(vtrt, vsurv)];
        ##this -1 is for c program sorting
        ##should not be changed to na
        data.dead[,as.character(endp)] <- -1;
    } else {
        data.dead <- NULL;
    }

    ##alive
    inx.a <- which(data.full[,vsurv] > duration);
    if (length(inx.a) > 0) {
        data.alive <- data.full[inx.a, c(vtrt, vsurv, endp)];
    } else {
        data.alive <- NULL;
    }

    ##
    dat.red <- rbind(data.dead, data.alive);
    rst     <- list(NULL);
    for (i in 1:length(atrt)) {
        rst[[i]] <- dat.red[which(atrt[i] == dat.red[,1]), 2:3];
    }
    rst
}

##------------------------------------------------------
##
##           BOOTSTRAP
##
##------------------------------------------------------

##get bootstrap sample
get.bs.sample <- function(data.all=NULL, lst.var, a.trt=NULL, bs=TRUE) {

    ## treatment name
    vtrt <- lst.var$trt;

    ##return original dataset
    if (!bs)
        rst <- 1:nrow(data.all);

    ##
    if (is.null(a.trt))
        a.trt <- sort(unique(data.all[,vtrt]));

    rst <- NULL;
    for (i in 1:length(a.trt)) {
        cur.inx <- which(a.trt[i] == data.all[,vtrt]);
        cur.smp <- sample(cur.inx, length(cur.inx), TRUE);
        rst     <- c(rst, cur.smp);
    }

    rst;
}

get.sum.rank <- function(mat.rank) {
    n.imp    <- max(mat.rank[,1]);
    n.each   <- nrow(mat.rank)/n.imp;
    m.r      <- mat.rank[,ncol(mat.rank)];
    dim(m.r) <- c(n.each, n.imp);
    rst      <- cbind(mat.rank[1:n.each, 2:3],
                      apply(m.r, 1, mean));
}

get.sum.median <- function(mat.median, offset=10000) {
    n.imp    <- max(mat.median[,1]);
    n.each   <- nrow(mat.median)/n.imp;

    m.r      <- apply(mat.median, 1, function(x) {
                      if (!is.na(x[5])) {
                          rst <- x[5];
                      } else {
                          rst <- x[4] - offset;
                      }
                      rst
                      });
    dim(m.r) <- c(n.each, n.imp);
    rst      <- cbind(mat.median[1:n.each, 2:3],
                      apply(m.r, 1, median));
}


##combine bootstrap results if using cluster
combine.boot.all <- function(n.boot, prefix="boot_rst", rst.name="rst.bs") {
    rst <- rep(list(NULL), n.boot);
    for (i in 1:n.boot) {
        cur.f <- paste(prefix, i, ".Rdata", sep="");
        load(cur.f);
        cur.rst  <- get(rst.name);
        rst[[i]] <- cur.rst[[1]];
    }
    rst
}

##bootstrap single case
get.boot.single<- function(data.all,
                           lst.var,
                           deltas = NULL,
                           quantiles=NULL,
                           ...) {
    goodImp <- FALSE;
    while (!goodImp) {
        rst <- tryCatch({
            smp.inx  <- get.bs.sample(data.all, lst.var);
            cur.smp  <- data.all[smp.inx,];
            fit.rst  <- fit.model(cur.smp, lst.var);
            cur.full <- get.imp.all(cur.smp, fit.rst, lst.var, deltas=deltas, ...);
            rst      <- get.theta.quantiles(cur.full, lst.var, deltas=deltas, quantiles=quantiles);
            goodImp  <- TRUE;
            rst
        }, error = function(e) {
            goodImp <- FALSE;
			     	print('Error in Imputation! Resampling bootstrap sample');
			     	print(e);
        })
    }
    rst
}
