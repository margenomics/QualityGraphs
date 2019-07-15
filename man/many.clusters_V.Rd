\name{many.clusters_V}
\alias{many.clusters_V}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make dendrogram with different clustering techniques
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
many.clusters_V(x, resultsDir, Filename, Title, parameters, toPDF = TRUE, toPNG=TRUE, conditions = NULL, colors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{resultsDir}{
%%     ~~Describe \code{resultsDir} here~~
}
  \item{Filename}{
%%     ~~Describe \code{Filename} here~~
}
  \item{Title}{
%%     ~~Describe \code{Title} here~~
}
  \item{parameters}{
%%     ~~Describe \code{parameters} here~~
}
  \item{toPDF}{
%%     ~~Describe \code{toPDF} here~~
}
  \item{toPNG}{
%%     ~~Describe \code{toPDF} here~~
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
function (x, resultsDir, Filename, Title, parameters, toPDF = TRUE, 
    conditions = NULL, colors = NULL) 
{
    labels <- colnames(x)
    use.cor = "pairwise.complete.obs"
    if (toPDF) {
        pdf(file = file.path(resultsDir, paste(Filename, "pdf", 
            sep = ".")))
    }
    if (is.null(conditions)) {
        opt <- par(cex.main = 2, cex = parameters$ce, cex.lab = 1.5, 
            cex.axis = 1.5)
        clust.cor.ward <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "ward.D2")
        plot(clust.cor.ward, main = Title, hang = -1)
        clust.cor.average <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "average")
        plot(clust.cor.average, main = Title, hang = -1)
        clust.cor.complete <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "complete")
        plot(clust.cor.complete, main = Title, hang = -1)
        clust.euclid.ward <- hclust(dist(t(x)), method = "ward.D2")
        plot(clust.euclid.ward, main = Title, hang = -1)
        clust.euclid.average <- hclust(dist(t(x)), method = "average")
        plot(clust.euclid.average, main = Title, hang = -1)
        clust.euclid.complete <- hclust(dist(t(x)), method = "complete")
        plot(clust.euclid.complete, main = Title, hang = -1)
        par(opt)
    }
    else if (is.null(colors)) {
        list1 <- unique(as.character(sort(conditions)))
        ColVect <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
        list2 <- ColVect[1:length(unique(conditions))]
        map = setNames(list2, list1)
        colors <- map[conditions]
        color_cluster <- function(hclus, condition, ce) {
            sampleDendrogram <- as.dendrogram(hclus)
            names(condition) <- hclus$labels
            sampleDendrogram <- dendrapply(sampleDendrogram, 
                function(x, batch) {
                  if (is.leaf(x)) {
                    label <- attr(x, "label")
                    attr(x, "nodePar") <- list(lab.col = as.vector(batch[label]), 
                      pch = "", cex = 0.8, lab.cex = ce)
                  }
                  x
                }, condition)
        }
        opt <- par(cex.main = 1, cex.axis = 0.8, cex = 0.8)
        clust.cor.ward <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "ward.D2")
        clust.cor.ward <- color_cluster(clust.cor.ward, colors, 
            parameters$ce)
        plot(clust.cor.ward, main = Title, xlab = "Ward")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.cor.average <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "average")
        clust.cor.average <- color_cluster(clust.cor.average, 
            colors, parameters$ce)
        plot(clust.cor.average, main = Title, xlab = "Average")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.cor.complete <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "complete")
        clust.cor.complete <- color_cluster(clust.cor.complete, 
            colors, parameters$ce)
        plot(clust.cor.complete, main = Title, xlab = "Complete")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.euclid.ward <- hclust(dist(t(x)), method = "ward.D2")
        clust.euclid.ward <- color_cluster(clust.euclid.ward, 
            colors, parameters$ce)
        plot(clust.euclid.ward, main = Title, xlab = "Euclidean Ward")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.euclid.average <- hclust(dist(t(x)), method = "average")
        clust.euclid.average <- color_cluster(clust.euclid.average, 
            colors, parameters$ce)
        plot(clust.euclid.average, main = Title, xlab = "Euclidean Average")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.euclid.complete <- hclust(dist(t(x)), method = "complete")
        clust.euclid.complete <- color_cluster(clust.euclid.complete, 
            colors, parameters$ce)
        plot(clust.euclid.complete, main = Title, xlab = "Euclidean Complete")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        par(opt)
    }
    else {
        list1 <- unique(as.character(conditions))
        list2 <- unique(colors)
        map = setNames(list2, list1)
        colors <- map[conditions]
        color_cluster <- function(hclus, condition, ce) {
            sampleDendrogram <- as.dendrogram(hclus)
            names(condition) <- hclus$labels
            sampleDendrogram <- dendrapply(sampleDendrogram, 
                function(x, batch) {
                  if (is.leaf(x)) {
                    label <- attr(x, "label")
                    attr(x, "nodePar") <- list(lab.col = as.vector(batch[label]), 
                      pch = "", cex = 0.8, lab.cex = ce)
                  }
                  x
                }, condition)
        }
        opt <- par(cex.main = 1, cex.axis = 0.8, cex = 0.8)
        clust.cor.ward <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "ward.D2")
        clust.cor.ward <- color_cluster(clust.cor.ward, colors, 
            parameters$ce)
        plot(clust.cor.ward, main = Title, xlab = "Ward")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.cor.average <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "average")
        clust.cor.average <- color_cluster(clust.cor.average, 
            colors, parameters$ce)
        plot(clust.cor.average, main = Title, xlab = "Average")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.cor.complete <- hclust(as.dist(1 - cor(x, use = use.cor)), 
            method = "complete")
        clust.cor.complete <- color_cluster(clust.cor.complete, 
            colors, parameters$ce)
        plot(clust.cor.complete, main = Title, xlab = "Complete")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.euclid.ward <- hclust(dist(t(x)), method = "ward.D2")
        clust.euclid.ward <- color_cluster(clust.euclid.ward, 
            colors, parameters$ce)
        plot(clust.euclid.ward, main = Title, xlab = "Euclidean Ward")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.euclid.average <- hclust(dist(t(x)), method = "average")
        clust.euclid.average <- color_cluster(clust.euclid.average, 
            colors, parameters$ce)
        plot(clust.euclid.average, main = Title, xlab = "Euclidean Average")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        clust.euclid.complete <- hclust(dist(t(x)), method = "complete")
        clust.euclid.complete <- color_cluster(clust.euclid.complete, 
            colors, parameters$ce)
        plot(clust.euclid.complete, main = Title, xlab = "Euclidean Complete")
        legend("topright", legend = list1, cex = parameters$ce + 
            0.2, fill = list2)
        par(opt)
    }
    if (toPDF) {
        dev.off()
    }
    return(list(Corr.ward = clust.cor.ward, Corr.avg = clust.cor.average, 
        Corr.compl = clust.cor.complete, Euclid.ward = clust.euclid.ward, 
        Euclid.avg = clust.euclid.average, Euclid.compl = clust.euclid.complete))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
