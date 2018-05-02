\name{qualityplots}
\alias{qualityplots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make different plots for microarray quality assessment.
}
\description{
This function can be used as a shortcut to make all quality plots included in QualityGraphs package with a single call. 
}
\usage{
qualityplots(ds, picname, estimates_m = NULL, est_noctrls, labels = NULL, ds.N = NULL, plmTr = NULL, conditions = NULL, colors = NULL, groupn = NULL, all = FALSE, Density = TRUE, Boxplots = TRUE, Clusters = TRUE, PCA = TRUE, estimates = FALSE, noctrls = TRUE, resDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ds}{
Object obtained with aroma.affymetrix and the function AffymetrixCelSet in the aroma.affymetrix package.
}
  \item{picname}{
Name of output files.
}
  \item{estimates_m}{
Numeric matrix where columns are samples and rows genes. estimates_m corresponds to the expression matrix with controls.
}
  \item{est_noctrls}{
Numeric matrix where columns are samples and rows genes. est_noctrls corresponds to the expression matrix without controls.
}
  \item{labels}{
Optional. Vector with sample names. Samples are in the same order as in the ds object. Default get sample names from the ds object. It only works for boxplots and density plots
}
  \item{ds.N}{
Aroma affymetrix ds.N object with the normalized data from .CEL files using RMA and quantile normalization.
}
  \item{plmTr}{
plmTr object obtained with the ExonRmaPlm function in aroma.affymetrix() package.
}
  \item{conditions}{
Optional, vector with sample condition in the same order as samples appear in the ds object.
}
  \item{colors}{
Optional, vector with the colous corresponding to each condition. Same order as samples in ds object in the case of boxplots, and same order as columns in the case of clusters and PCA. If the colors are set, clusters, PCA and boxplots have to be generated separately, as samples have a different order.
}
  \item{groupn}{
Optional. Number of samples to include in each density plot. If the number is lower than total number of samples, more than one plot is going to be shown. Default is groupn = all samples.
}
  \item{all}{
Optional. If TRUE, an additional plot with all density plots made with a groupn is made.
}
  \item{Density}{
Call densityplot_all function to make density plots. Default is TRUE.
}
  \item{Boxplots}{
Call boxplot_all function to make boxplots. Default is TRUE.
}
  \item{Clusters}{
Call clusterdend function to make cluster dendrograms. Default is TRUE.
}
  \item{PCA}{
Call makePCA function to make PCA plots. Default is TRUE
}
  \item{estimates}{
If TRUE, clusters for estimates_m are made. Default is FALSE.
}
  \item{noctrls}{
If TRUE, clusters for est_noctrls are made. Default is TRUE.
}
  \item{resDir}{
Output results directory. Default is ResultsDir.
}
}
\details{
This function is a shortcut to call functions: makePCA, clusterdend, boxplot_all and densityplot_all. It makes all the plots for quality assessment in microarray data. \strong{Density plots} show the distribution of log2 raw intensities to check possible technical biases. \strong{Boxplots} show the log2 raw intensities in boxplots to check possible technical biases. \strong{BoxplotsRMA} show normalized log2 intensities with quantile normalization and Robust Multiarray average (RMA). It is used to check if normalization is made properly. \strong{Boxplot NUSE} plot with normalized unscaled standard errors (NUSE). The standard error from the probe-level model are visualized as boxplots. The deviant arrays can be identified by not being centered at 1, or being more spread out than the other arrays.\strong{Boxplot RLE} plot with relative log expression values (RLE). Assuming that most probesets on arrays are not changing, most of the RLE values are close 0. When examining the boxplot, the deviant arrays can be identified by not being centered at 0, or being more spread out than the other arrays.\strong{Clusters} show the sample distribution using different methods. Samples are expected to group together based on their similarity. \strong{PCA} show the sample distribution using a principal component analysis (PCA) in a 3D plot. }
\value{
Files are created in the resDir directory with the picname and .png or .pdf extensions.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Magdalena Arnal Segura <marnal@imim.es>.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}