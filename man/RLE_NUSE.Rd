\name{RLE_NUSE}
\alias{RLE_NUSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make relative log expression values (RLE) and normalized unscaled standard errors (NUSE) boxplots.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RLE_NUSE(bpstats, kind, strt = NULL, nd = NULL, parameters, lab, conditions = NULL, colors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bpstats}{
%%     ~~Describe \code{bpstats} here~~
}
  \item{kind}{
%%     ~~Describe \code{kind} here~~
}
  \item{strt}{
%%     ~~Describe \code{strt} here~~
}
  \item{nd}{
%%     ~~Describe \code{nd} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{lab}{
%%     ~~Describe \code{lab} here~~
}
  \item{conditions}{
%%     ~~Describe \code{conditions} here~~
}
  \item{colors}{
%%     ~~Describe \code{colors} here~~
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
function (bpstats, kind, strt = NULL, nd = NULL, parameters, 
    lab, conditions = NULL, colors = NULL) 
{
    if (is.null(labels)) {
        labels <- gsub(paste("_(", xip, ")", sep = ""), "", ds$Names, 
            fixed = TRUE)
    }
    if (is.null(strt) & is.null(nd)) {
        strt <- 1
        nd <- length(bpstats)
    }
    else if (is.null(strt)) {
        strt <- 1
    }
    else if (is.null(nd)) {
        nd <- length(bpstats)
    }
    fitxer_int <- as.list(1:nd)
    for (i in (strt:nd)) {
        fitxer_int[[i]] <- bpstats[[i]]$stats
    }
    if (is.null(conditions)) {
        boxplot(fitxer_int[strt:nd], main = paste(kind, "plot"), 
            xaxt = "n")
        axis(1, at = 1:length(lab), labels = lab, cex.axis = parameters$ce, 
            las = 2)
    }
    else if (is.null(colors)) {
        CondNames <- mixedsort(lab)
        CondTable <- data.frame(CondNames, conditions)
        conditions.o <- CondTable[match(lab, CondTable$CondNames), 
            "conditions"]
        list1 <- unique(as.character(sort(conditions.o)))
        ColVect <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
        list2 <- ColVect[1:length(unique(conditions.o))]
        map = setNames(list2, list1)
        colors <- map[conditions.o]
        boxplot(fitxer_int[strt:nd], main = paste(kind, "plot"), 
            xaxt = "n", col = colors)
        axis(1, at = 1:length(lab), labels = lab, cex.axis = parameters$ce, 
            las = 2)
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
    }
    else {
        CondNames <- mixedsort(lab)
        CondTable <- data.frame(CondNames, conditions)
        conditions.o <- CondTable[match(lab, CondTable$CondNames), 
            "conditions"]
        list1 <- unique(conditions.o)
        list2 <- colors
        map = setNames(list2, list1)
        colors <- map[conditions.o]
        boxplot(fitxer_int[strt:nd], main = paste(kind, "plot"), 
            xaxt = "n", col = colors)
        axis(1, at = 1:length(lab), labels = lab, cex.axis = parameters$ce, 
            las = 2)
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
