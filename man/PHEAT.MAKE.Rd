% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Pheat.make.R
\name{PHEAT.MAKE}
\alias{PHEAT.MAKE}
\title{Heatmaps using pheatmap package}
\usage{
PHEAT.MAKE(
  matrix,
  cols.s = NULL,
  filename,
  filedir,
  rownames.vect = "",
  cex.col = 0.8,
  scale = "none",
  toPNG = FALSE,
  annotation_col = NULL,
  showRownames = F
)
}
\arguments{
\item{matrix}{Numeric matrix with the intensity values to plot.}

\item{cols.s}{List for specifying annotation_row and annotation_col track 
colors manually. It is possible to define the colors for only some of the features.}

\item{filename}{Output file name.}

\item{filedir}{Output directory, it can be a path.}

\item{rownames.vect}{Vector with the names of the genes in the heatmap.}

\item{cex.col}{Size of the column labels}

\item{scale}{Character indicating if the values should be centered and 
scaled in either the row direction or the column direction, or none. 
Corresponding values are "row", "column" and "none"}

\item{toPNG}{TRUE if we want .png output, FALSE if we want .pdf output.}

\item{annotation_col}{Similar to annotation_row, but for columns.}

\item{showRownames}{Indicates wether you want to show rownames or not. Each row defines the features for a specific row.}
}
\value{
The plot is created in the "filedir" directory with the name "filename"
in the especified format (default pdf).
}
\description{
Function to make heatmaps from a numeric matrix using pheatmap package
}
\details{
This function makes a plot in png format using the results from the differential 
expression analysis with "limma" package. It uses the columns "adj.P.Val" or 
"P.Value" depending on the value "adj" and the column logFC. Cut-off is used in 
order to indicate the DE genes in the scatterplot.
}
\author{
Ariadna Acedo Terrades <aacedo@imim.es>
}
