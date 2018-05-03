\name{clusterdend}
\alias{clusterdend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clusters for normalized microarray data
}
\description{
Function to make plots with dendrograms using different clustering methods. Methods are Ward.D2, average and complete with correlations and euclidean distances.
}
\usage{
clusterdend(estimates_m = NULL, est_noctrls = NULL, picname, conditions = NULL, colors = NULL, estimates = FALSE, noctrls = TRUE, resDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{estimates_m}{
Numeric matrix where columns are samples and rows genes. estimates_m corresponds to the expression matrix with controls.
}
  \item{est_noctrls}{
Numeric matrix where columns are samples and rows genes. est_noctrls corresponds to the expression matrix without controls.
}
  \item{picname}{
Character with the output file name.
}
  \item{conditions}{
Optional. Vector with sample conditions to colour sample names in the plot. conditions vector has to have the same order as the order of samples in the columns of the input expression matrix. IMPORTANT Sample columns in the expression matrix and sample conditions have to be in order mixedsort(labels).
}
  \item{colors}{
Optional. Vector with colors corresponding to each sample with the order mixedsort(samplename). If this parameter is not defined, default colors will be used.
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
This function makes a cluster analysis for a numeric matrix. Optionally two cluster analysis can be made, one for a estimates_m (expression matrix with controls control AffyID) and est_noctrls (expression matrix without controls). Methods used for clustering are: Ward.D2, average and complete with correlations and euclidean distances. Column names in the expression matrix are used as sample names. 
}
\value{
This function outputs a .pdf file with several dendrograms for each method and distance.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Magdalena Arnal Segura <marnal@imim.es>
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

