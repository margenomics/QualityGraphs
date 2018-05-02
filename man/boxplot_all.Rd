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
Aroma affymetrix ds.N object with the normalized data from .CEL files using RMA and quantile normalization.
}
  \item{plmTr}{
plmTr object obtained with the ExonRmaPlm function in aroma.affymetrix() package.
}
  \item{picname}{
Character with the output file name.
}
  \item{labels}{
Optional. Vector with sample names. Samples are in the same order as in the ds object. Default get sample names from the ds object.
}
  \item{conditions}{
Optional, vector with sample condition in the same order as samples appear in the ds object.
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
Output results directory. Default is ResultsDir.
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
Magdalena Arnal Segura <marnal@imim.es>
}

\examples{


}

