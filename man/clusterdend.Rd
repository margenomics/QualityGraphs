\name{clusterdend}
\alias{clusterdend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clusters for normalized microarray data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clusterdend(estimates_m = NULL, est_noctrls = NULL, picname, conditions = NULL, colors = NULL, estimates = FALSE, noctrls = TRUE, resDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimates_m}{
%%     ~~Describe \code{estimates_m} here~~
}
  \item{est_noctrls}{
%%     ~~Describe \code{est_noctrls} here~~
}
  \item{picname}{
%%     ~~Describe \code{picname} here~~
}
  \item{conditions}{
%%     ~~Describe \code{conditions} here~~
}
  \item{colors}{
%%     ~~Describe \code{colors} here~~
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
function (estimates_m = NULL, est_noctrls = NULL, picname, conditions = NULL, 
    colors = NULL, estimates = FALSE, noctrls = TRUE, resDir = NULL) 
{
    if (!is.null(resDir)) {
        resultsDir = resDir
    }
    if (estimates) {
        labels <- colnames(estimates_m)
        parameters <- setparam(labels)
        many.clusters_V(estimates_m, resultsDir, paste(picname, 
            "ClustersWithCtrls", sep = "_"), "Clusters with ctrls", 
            parameters = parameters, conditions = conditions, 
            colors = colors)
    }
    if (noctrls) {
        labels <- colnames(est_noctrls)
        parameters <- setparam(labels)
        many.clusters_V(est_noctrls, resultsDir, paste(picname, 
            "Clusters", sep = "_"), "Clusters", parameters = parameters, 
            conditions = conditions, colors = colors)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
