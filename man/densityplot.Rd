\name{densityplot}
\alias{densityplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to make density plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
densityplot(ds, strt = NULL, nd = NULL, lab, parameters)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ds}{
%%     ~~Describe \code{ds} here~~
}
  \item{strt}{
%%     ~~Describe \code{strt} here~~
}
  \item{nd}{
%%     ~~Describe \code{nd} here~~
}
  \item{lab}{
%%     ~~Describe \code{lab} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
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
function (ds, strt = NULL, nd = NULL, lab, parameters) 
{
    if (is.null(strt) & is.null(nd)) {
        strt <- 1
        nd <- length(ds)
    }
    else if (is.null(strt)) {
        strt <- 1
    }
    else if (is.null(nd)) {
        nd <- length(ds)
    }
    yvalues <- vector()
    for (o in strt:nd) {
        fitxer <- getFile(ds, idx = o)
        dl <- log(aroma.affymetrix::getData(fitxer, fields = "intensities", 
            indices = 1:nbrOfCells(fitxer))$intensities, base = 2)
        ycoor <- ceiling(max(density(dl)$y))
        yvalues <- c(yvalues, ycoor)
    }
    colors <- rainbow(nd - strt + 1)
    fitxer <- getFile(ds, idx = strt)
    dl <- log(aroma.affymetrix::getData(fitxer, fields = "intensities", 
        indices = 1:nbrOfCells(fitxer))$intensities, base = 2)
    plot(density(dl), col = colors[1], main = paste("Samples", 
        strt, "to", nd), xlab = "log2(y)", ylim = c(0, max(yvalues)), 
        lty = strt, ylab = "")
    for (e in (strt + 1):nd) {
        fitxer <- getFile(ds, idx = e)
        dl <- log(aroma.affymetrix::getData(fitxer, fields = "intensities", 
            indices = 1:nbrOfCells(fitxer))$intensities, base = 2)
        lines(density(dl), col = colors[(e - strt) + 1], lty = e)
    }
    if (length(lab) < 15) {
        legend("topright", legend = lab, cex = 0.6, col = colors[1:nd], 
            lty = strt:nd, lwd = 1.6)
    }
    else {
        legend("topright", legend = lab, cex = parameters$ce, 
            col = colors[1:nd], lty = strt:nd, lwd = 1.6)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
