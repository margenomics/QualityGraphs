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
Optional. Vector with sample conditions to colour sample names in the plot. conditions vector has to have the same order as the order of samples in the columns of the input expression matrix. IMPORTANT Sample columns in the expression matrix and sample conditions have to be in order mixedsort(labels).
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

}
