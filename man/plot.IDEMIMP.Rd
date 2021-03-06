% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idem_imputation.R
\name{plot.IDEMIMP}
\alias{plot.IDEMIMP}
\title{Plot imputation results}
\usage{
\method{plot}{IDEMIMP}(x, opt = c("imputed", "composite"), fname = NULL, ...)
}
\arguments{
\item{x}{A class \code{IDEMIMP} object returned from \code{\link{imImpAll}}}

\item{opt}{Types of the plot
   \itemize{
    \item{\code{imputed}: }{Plot density of imputed values and the density of the observed outcomes}
    \item{\code{composite}: }{Generate cumulative plot of the composite survival and functional outcome}
}}

\item{fname}{File name of the result pdf file. If \code{fname} is null,
result pdf file will not be generated}

\item{...}{Options for generating the plots.

   \describe{
\item{type = imputed}{
\itemize{}
\itemize{
\item{\code{deltas}: }{Imputation sensitivity parameter for which to generate
    the results}

\item{\code{endp}: }{If \code{TRUE}, plot the densities of the imputed
    functional outcomes. Otherwise, plot the densities of the imputed
    outcomes}

\item{\code{adj}}{\code{density} estimation option}

\item{\code{cols}}{\code{plot} option for colors}

\item{\code{ltys}}{\code{plot} options for line types}

\item{\code{xlim}}{\code{plot} options}

\item{\code{ylim}}{\code{plot} options}

\item{\code{mfrow}}{\code{plot} options}

}}

\item{type = composite}{
    \itemize{}

\itemize{ \item{\code{at.surv}: }{Sets the range of the survival times to
plot in the cumulative distribution function. By default the range is the
range of survival values up to the duration of the study}

\item{\code{at.z}: }{Sets the range of the functional outcome to plot in the
    cumulative distribution function. By defualt this is the range of the
    functional outcomes plus the buffer amount to improve visibility in the
    transition from survival to functional outcome}

\item{\code{p.death}: }{Proportion
    of the plot width devoted to Survival. By default the cumulative
    distribution will devote horizontal space to the survival portion that is
    proportional to the number of subjects who die prior to duration}

\item{\code{buffer}: }{Small horizontal gap used to better visually distinguish
    the transition from survival to functional outcome}


\item{\code{delta}: }{Imputation sensitivity parameter for which to generate the
    results}

\item{\code{seg.lab}: }{Labels for the two components of the composite
    outcome}

\item{\code{main}: }{\code{plot} options}

}}
}}
}
\description{
Generate different types of plots for class \code{IDEMIMP} objects generated
by \code{\link{imImpAll}}
}
\examples{
\dontrun{
im.abc <- imData(abc, trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
                 y0=NULL, endfml="Y2",
                 trt.label = c("UC+SBT", "SAT+SBT"),
                 cov=c("AGE"), duration=365, bounds=c(0,100));
rst.fit  <- imFitModel(im.abc);
rst.imp <- imImpAll(rst.fit, deltas=c(-0.25,0,0.25),
                    normal=TRUE, chains = 2, iter = 2000, warmup = 1000);
plot(rst.imp, opt = "imputed"),
plot(rst.imp, opt = "composite")}
}
\seealso{
\code{\link{imImpAll}}
}
