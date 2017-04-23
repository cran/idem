## ---- eval=T, echo=FALSE-------------------------------------------------
require(idem);
set.seed(1000);

## ---- eval=T, echo=TRUE--------------------------------------------------
head(abc);

## ---- echo=TRUE, eval=T--------------------------------------------------
 err.lst.var <- list(trt="TRT", outcome=c("Y1","Y2"),
                 y0=NULL, endfml="Y2", bounds=c(10,20),
                 duration=365);
 cat(imChkPars(abc, err.lst.var));

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
lst.var <- list(trt="TRT", surv="SURV", outcome=c("Y1","Y2"), y0=NULL,
                endp=c("Y2"), unitTime="days",
                trt.label = c("UC+SBT", "SAT+SBT"),
                cov=c("AGE"), endfml="Y2", duration=365, bounds=c(0,100));

imPlotCompleters(abc, lst.var);

## ---- echo=TRUE----------------------------------------------------------
imMisTable(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
imPlotMisPattern(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
imPlotSurv(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
rst.fit <- imFitModel(abc, lst.var);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
rst.mixing <- imImpSingle(abc[1,], rst.fit, chains = 4, iter = 2000, warmup = 1000);
rstan::traceplot(rst.mixing$rst.stan, "YMIS");

## ---- echo=TRUE, results="hide"------------------------------------------
rst.imp <- imImpAll(abc, rst.fit, deltas=c(-0.25,0,0.25),
                    normal=TRUE, chains = 4, iter = 300, warmup = 100);


## ---- echo=TRUE----------------------------------------------------------
tail(rst.imp$complete, n=10);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
imPlotImputed(rst.imp, deltas = c(-0.25,0,0.25), xlim=c(0,100), endp=FALSE);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
imPlotImputed(rst.imp, deltas=c(-0.25,0,0.25), xlim=c(0,100), endp=TRUE);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
imPlotComposite(rst.imp, delta=0);

## ---- echo=TRUE----------------------------------------------------------
rst.est <- imEstimate(rst.imp, quantiles=c(0.25,0.5,0.75));
print(rst.est$theta);
print(rst.est$quantiles);

## ---- echo=TRUE----------------------------------------------------------
rst.boot <- imBs(rst.imp,
                 n.boot = 2,
                 n.cores = 1,
                 quantiles = c(0.25,0.5,0.75));

## ---- echo=TRUE----------------------------------------------------------
rst.final <- imTest(rst.boot);
print(rst.final);

## ---- echo=TRUE, fig.width=6, fig.height=5-------------------------------
imPlotContour(rst.final, nlevels = 30, con.v=0.05, zlim=c(0, 0.05));

## ---- echo=TRUE, eval=FALSE----------------------------------------------
#  imShiny();

