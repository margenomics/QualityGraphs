\name{boxplot_all}
\alias{boxplot_all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BoxPlots for quality control
}
\description{
Function to make all boxplots required for quality control.
}
\usage{
boxplot_all(ds, ds.N = NULL, plmTr = NULL, picname, labels = labels, conditions = NULL, colors = NULL, raw.bp = TRUE, RMA.bp = TRUE, RLE = TRUE, NUSE = TRUE, resDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ds}{
ds object created with aroma.affymetrix and the function AffymetrixCelSet() from the .CEl files and the cdf file.
}
  \item{ds.N}{
object with normalized data from .CEL files with RMA and quantile normalization.
}
  \item{plmTr}{
plmTr object obtained with the ExonRmaPlm function in aroma.affymetrix() package.
}
  \item{picname}{
character corresponding to the name for the output graphs.
}
  \item{labels}{
vector with sample names with the same order as in the expression matrix.
}
  \item{conditions}{
Optional, vector with sample condition in the same order as in the expression matrix.
}
  \item{colors}{
Optional, vector with the colours corresponding to each condition.
}
  \item{raw.bp}{
Plot boxplots with raw data at log2 scale. Default TRUE.
}
  \item{RMA.bp}{
Plot boxplots with normalized log2 data. Default TRUE.
}
  \item{RLE}{
Plot boxplots with relative log expression values (RLE). Default TRUE.
}
  \item{NUSE}{
Plot boxplots with normalized unscaled standard errors (NUSE). Default TRUE.
}
  \item{resDir}{
Results directory. Default ResultsDir.
}
}
\details{
This function outputs boxplots for quality assessment in microarray data using objects obtained with aroma.affymetrix package pre-processing. A boxplot with raw data, and a boxplot with Normalized data. Additionally reports relative log expression values (RLE), and the second the normalized unscaled standard errors (NUSE).
}
\value{
This function outputs plots in .png format in the resDir directory.
}
\references{
H. Bengtsson, K. Simpson, J. Bullard & K. Hansen. aroma.affymetrix: A generic framework in R for analyzing small to very large Affymetrix data sets in bounded memory, Tech Report #745, Department of Statistics, University of California, Berkeley, February 2008.
}
\author{
Magdalena Arnal 
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
