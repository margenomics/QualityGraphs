\name{Heat.make}
\alias{Heat.make}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make heatmaps from a numeric matrix.
}
\description{
  This function makes a heatmap from a numeric matrix using the package "gplots". The output is an image in pdf format. 
}
\usage{
  Heat.make(data.clus, cols.s, filename, filedir, rownames=TRUE, rownames.vect=NULL, cex.col=0.8, cex.row=0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.clus}{
    Matrix to plot in the heatmap.
  }
  \item{cols.s}{
    Vector of colors corresponding to the columns in the matrix data.clus.
  }
  \item{filename}{
    Name of the output file in .pdf format.
  }
  \item{filedir}{
    Name of the resultsDir.
  }
  \item{rownames}{
    If TRUE row names are shown in the heatmap.
  }
  \item{rownames.vect}{
    Vector with row names to be shown in the heatmap.
  }
  \item{cex.col}{
    Size of column names in the plot.
  }
  \item{cex.row}{
    Size of row names in the plot, (if rownames=TRUE).
  }
}
\details{
  This function uses the package gplots to make a heatmap. The algorithm used for the dendrogram above (columns) is "ward.D2" and the algorithm used for the rows (dendrogram not shown) is "average". Column names are always in the plot, and with the option rownames=TRUE you chose weather to show or not the rownames, not recommended if you have a large number of rows. options cex.col and cex.row are used in order to change the size of row and column names respectively.
}
\value{
  The plot is created in the "filedir" directory with the name "filename" in ".pdf" format.
}
\author{
  Magdalena Arnal Segura <marnal@imim.es>.
}

