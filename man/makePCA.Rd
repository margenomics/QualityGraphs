\name{makePCA}
\alias{makePCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D Plot with Principal Component Analysis (PCA)
}
\description{
This function performs a principal component analysis with the prcomp function and uses the prcomp class object to plot the three first principal components.
}
\usage{
makePCA(est_noctrls, picname, conditions = NULL, colors = NULL, dist = 2, resDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est_noctrls}{
Numeric matrix where columns are samples and rows genes. est_noctrls corresponds to a matrix without controls.
}
  \item{picname}{
Character with the output file name.
}
  \item{conditions}{
Optional. Vector with sample conditions to colour sample names in the plot. conditions vector has to have the same order as the order of samples in the columns of the input expression matrix.
}
  \item{colors}{
Optional. Vector with colors corresponding to each sample. If this parameter is not defined, default colors will be used.
}
  \item{dist}{
Distance where the name of the sample is located with respect to the corresponing point in the plot. Optimal distance is different deppending on the sample number and PC distances. Default is 2. 
}
  \item{resDir}{
Output results directory. Default is ResultsDir.
}
}
\details{
The first three principal components obtained with the prcomp function capture part of the variability explaining the sample distribution. Each PC is represented in an axis of the plot, x.axis=PC1, y.axis=PC2 and z.axis=PC3. This plot is useful to check the overall similarity among samples. Column names in the expression matrix are used as sample names.
}
\value{
Output is a .png file with the results in resDir.
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
function (est_noctrls, picname, conditions = NULL, colors = NULL, 
    dist = 2, resDir = NULL) 
{
    require(scatterplot3d)
    if (!is.null(resDir)) {
        resultsDir = resDir
    }
    labels <- colnames(est_noctrls)
    parameters <- setparam(labels)
    summary(pca.filt <- prcomp(t(est_noctrls), center = TRUE, 
        scale. = TRUE))
    vars <- apply(pca.filt$x, 2, var)
    props <- vars/sum(vars)
    PCAvec <- round(props * 100, 2)
    png(file.path(resultsDir, paste(picname, "PCA.png", sep = "_")), 
        width = parameters$wid, height = parameters$hei, res = parameters$res)
    if (is.null(conditions)) {
        pca3d <- scatterplot3d(x = pca.filt$x[, 1], y = pca.filt$x[, 
            2], z = pca.filt$x[, 3], xlab = sprintf("PC1 \%.0f\%\%", 
            PCAvec[1]), ylab = sprintf("PC2 \%.0f\%\%", PCAvec[2]), 
            zlab = sprintf("PC3 \%.0f\%\%", PCAvec[3]), main = "PCA", 
            col.grid = "lightblue", cex.symbols = parameters$ce + 
                0.2)
        text(pca3d$xyz.convert(pca.filt$x + dist), labels = rownames(pca.filt$x), 
            cex = parameters$ce + 0.2)
    }
    else if (is.null(colors)) {
        list1 <- unique(as.character(sort(conditions)))
        ColVect <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))
        list2 <- ColVect[1:length(unique(conditions))]
        map = setNames(list2, list1)
        colors <- map[conditions]
        pca3d <- scatterplot3d(x = pca.filt$x[, 1], y = pca.filt$x[, 
            2], z = pca.filt$x[, 3], xlab = sprintf("PC1 \%.0f\%\%", 
            PCAvec[1]), ylab = sprintf("PC2 \%.0f\%\%", PCAvec[2]), 
            zlab = sprintf("PC3 \%.0f\%\%", PCAvec[3]), main = "PCA", 
            col.grid = "lightblue", cex.symbols = parameters$ce + 
                0.2, color = colors, pch = 16)
        text(pca3d$xyz.convert(pca.filt$x + dist), labels = rownames(pca.filt$x), 
            cex = parameters$ce + 0.2, col = colors)
        legend("bottomright", legend = list1, col = list2, pch = 16, 
            ncol = 1, cex = 0.9)
    }
    else {
        list1 <- unique(conditions)
        list2 <- unique(colors)
        map = setNames(list2, list1)
        colors <- map[conditions]
        pca3d <- scatterplot3d(x = pca.filt$x[, 1], y = pca.filt$x[, 
            2], z = pca.filt$x[, 3], xlab = sprintf("PC1 \%.0f\%\%", 
            PCAvec[1]), ylab = sprintf("PC2 \%.0f\%\%", PCAvec[2]), 
            zlab = sprintf("PC3 \%.0f\%\%", PCAvec[3]), main = "PCA", 
            col.grid = "lightblue", cex.symbols = parameters$ce + 
                0.2, color = colors, pch = 16)
        text(pca3d$xyz.convert(pca.filt$x + dist), labels = rownames(pca.filt$x), 
            cex = parameters$ce + 0.2, col = colors)
        legend("bottomright", legend = list1, col = list2, pch = 16, 
            ncol = 1, cex = 0.9)
    }
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
