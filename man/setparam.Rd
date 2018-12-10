\name{setparam}
\alias{setparam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set parameters to make plots.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
setparam(labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{labels}{
%%     ~~Describe \code{labels} here~~
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
function (labels) 
{
    long <- length(labels)
    wid <- 3500
    hei <- 3500
    res <- 400
    if (max(nchar(labels)) < 10) {
        ce <- 0.8
    }
    else if (max(nchar(labels)) < 15) {
        ce <- 0.6
    }
    else {
        ce <- 0.4
    }
    if (long < 35) {
        if (max(nchar(labels)) < 10) {
            ce <- 0.8
        }
        else if (max(nchar(labels)) < 15) {
            ce <- 0.6
        }
        else {
            ce <- 0.4
        }
    }
    else if (long < 90) {
        ce <- 0.4
    }
    else if (long < 140) {
        ce <- 0.3
    }
    else {
        ce <- 0.2
    }
    return(list(wid = wid, hei = hei, res = res, ce = ce))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
