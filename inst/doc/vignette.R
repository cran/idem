## ---- eval=T, echo=FALSE-------------------------------------------------
require(idem);
set.seed(1000);

## ---- eval=T, echo=TRUE--------------------------------------------------
head(abc);

## ---- echo=TRUE, eval=T--------------------------------------------------
 err.lst.var <- list(trt="TRT", outcome=c("Y1","Y2"),
                 y0=NULL, endfml="Y2", bounds=c(10,20),
                 duration=365);

 chkPara(err.lst.var, abc);

 lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
                 y0=NULL, endp=c("Y2"),
                 unitTime="days",
                 trt.label = c("UC+SBT", "SAT+SBT"),
                 cov=c("AGE"), endfml="Y2",
                 duration=365, bounds=c(0,100));

 chkPara(lst.var, abc);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotCompleters(abc, lst.var);

## ---- echo=TRUE----------------------------------------------------------
get.mis.table(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotMisPattern(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotSurv(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
rst.fit <- fit.model(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
rst.mixing <- get.imp.all(abc, rst.fit, lst.var, deltas=0,
                          trace.n=1, normal=TRUE, iter=500,
                          p.scale=10);

rst.chain <- list();
for (i in 1:ncol(rst.mixing$mcmc[[1]])) {
    rst.chain[[i]] <- rst.mixing$mcmc[[1]][,i,drop=FALSE];
}
coda::traceplot(rst.chain);

## ---- echo=TRUE, results="hide"------------------------------------------
rst.imp <- get.imp.all(abc, rst.fit, lst.var, deltas=c(-0.25,0,0.25),
                       normal=TRUE, iter=20, n.imp=2, thin=1,
                       p.scale=10);

## ---- echo=TRUE----------------------------------------------------------
tail(rst.imp, n=10);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotImputed(rst.imp, lst.var, deltas=c(-0.25,0,0.25), xlim=c(0,100), endp=FALSE);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotImputed(rst.imp, lst.var, deltas=c(-0.25,0,0.25), xlim=c(0,100), endp=TRUE);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotComposite(rst.imp, lst.var, delta=0);

## ---- echo=TRUE----------------------------------------------------------
rst.est <- get.theta.quantiles(rst.imp, lst.var, quantiles=c(0.25,0.5,0.75));
print(rst.est$theta);
print(rst.est$quantiles);

## ---- echo=TRUE----------------------------------------------------------
rst.boot <- get.bs.all(n.boot = 2,
                       n.cores = 1,
                       data.all = abc,
                       lst.var = lst.var,
                       deltas = c(-0.25, 0, 0.25),
                       quantiles = c(0.25,0.5,0.75),
                       normal=TRUE, iter=20,
                       n.imp=2, thin=1,
                       p.scale=10);

## ---- echo=TRUE----------------------------------------------------------
rst.final <- get.overall.rst(rst.est, rst.boot);
print(rst.final);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
plotContour(rst.final, lst.var, nlevels = 30,
            con.v=0.05, zlim=c(0, 0.05));

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  run.idem();

