% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/idem_data.R
\name{plot.IDEMDATA}
\alias{plot.IDEMDATA}
\title{Plot of IDEMDATA object}
\usage{
\method{plot}{IDEMDATA}(
  x,
  opt = c("survivor", "missing", "KM"),
  cols = c("black", "blue"),
  fname = NULL,
  ...
)
}
\arguments{
\item{x}{A class \code{IDEMDATA} object generated by \code{\link{imData}}}

\item{opt}{Types of the plot

\itemize{

\item{\code{survivor: }}{Spaghetti plot for subjects alive at the end of the
    study}

\item{\code{missing}: }{Plot the missing patterns of the observed data}

\item{\code{KM}: }{Plot Kaplan-Meier survival curves}
}}

\item{cols}{Curve colors of the treatment and control arm for survival plot
or colors of the observed and missing data for missingness plot.}

\item{fname}{File name of the result pdf file. If \code{fname} is null,
result pdf file will not be generated}

\item{...}{Extra arguments for \code{plot}}
}
\description{
Generate different types of plots for class \code{IDEMDATA} objects.
}
\examples{
rst.data <- imData(abc, trt="TRT", surv="SURV", outcome=c("Y1","Y2"),
                 y0=NULL, endfml="Y2",
                 trt.label = c("UC+SBT", "SAT+SBT"),
                 cov=c("AGE"), duration=365, bounds=c(0,100));
plot(rst.data, opt = "survivor");
plot(rst.data, opt = "missing", cols = c("blue", "gray"));
plot(rst.data, opt = "KM");

}
\seealso{
\code{\link{imData}}
}
