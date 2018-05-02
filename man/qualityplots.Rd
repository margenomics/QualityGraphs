\name{qualityplots}
\alias{qualityplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make plots for microarray quality assessment
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qualityplots(ds, picname, estimates_m = NULL, est_noctrls, labels = NULL, ds.N = NULL, plmTr = NULL, conditions = NULL, colors = NULL, groupn = NULL, all = FALSE, Density = TRUE, Boxplots = TRUE, Clusters = TRUE, PCA = TRUE, estimates = FALSE, noctrls = TRUE, resDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ds}{
%%     ~~Describe \code{ds} here~~
}
  \item{picname}{
%%     ~~Describe \code{picname} here~~
}
  \item{estimates_m}{
%%     ~~Describe \code{estimates_m} here~~
}
  \item{est_noctrls}{
%%     ~~Describe \code{est_noctrls} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{ds.N}{
%%     ~~Describe \code{ds.N} here~~
}
  \item{plmTr}{
%%     ~~Describe \code{plmTr} here~~
}
  \item{conditions}{
%%     ~~Describe \code{conditions} here~~
}
  \item{colors}{
%%     ~~Describe \code{colors} here~~
}
  \item{groupn}{
%%     ~~Describe \code{groupn} here~~
}
  \item{all}{
%%     ~~Describe \code{all} here~~
}
  \item{Density}{
%%     ~~Describe \code{Density} here~~
}
  \item{Boxplots}{
%%     ~~Describe \code{Boxplots} here~~
}
  \item{Clusters}{
%%     ~~Describe \code{Clusters} here~~
}
  \item{PCA}{
%%     ~~Describe \code{PCA} here~~
}
  \item{estimates}{
%%     ~~Describe \code{estimates} here~~
}
  \item{noctrls}{
%%     ~~Describe \code{noctrls} here~~
}
  \item{resDir}{
%%     ~~Describe \code{resDir} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ds, picname, estimates_m = NULL, est_noctrls, labels = NULL, 
    ds.N = NULL, plmTr = NULL, conditions = NULL, colors = NULL, 
    groupn = NULL, all = FALSE, Density = TRUE, Boxplots = TRUE, 
    Clusters = TRUE, PCA = TRUE, estimates = FALSE, noctrls = TRUE, 
    resDir = NULL) 
{
    library(RColorBrewer)
    require(gtools)
    if (!is.null(resDir)) {
        resultsDir = resDir
    }
    if (Density) {
        densityplot_all(ds, picname = picname, groupn = groupn, 
            all = all, resDir = resultsDir)
    }
    if (Boxplots) {
        boxplot_all(ds, ds.N = ds.N, plmTr = plmTr, labels = labels, 
            picname = picname, conditions = conditions, colors = colors, 
            resDir = resultsDir)
    }
    if (Clusters) {
        clusterdend(estimates_m, est_noctrls, picname = picname, 
            conditions = conditions, colors = colors, estimates = estimates, 
            noctrls = noctrls, resDir = resultsDir)
    }
    if (PCA) {
        makePCA(est_noctrls, picname = picname, conditions = conditions, 
            colors = colors, resDir = resultsDir)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
