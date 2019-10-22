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
  Heat.make(data.clus, cols.s=NULL, filename, filedir, rownames.vect="", cex.col=0.8, clustinrows=T,clustincols=T, toPNG=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.clus}{
    Numeric matrix with the intensity values to plot.
  }
  \item{cols.s}{
    Vector with colors corresponding to the conditions in the heatmap with the same order as the columns in the matrix.
  }
  \item{filename}{
    Name of the output file.
  }
  \item{filedir}{
    Name of the resultsDir.
  }
  \item{rownames.vect}{
    Vector with the names of the genes in the heatmap. If rownames.vect="" no rownames are shown.
  }
  \item{cex.col}{
    Size of the column labels.
  }
  \item{clustinrows}{
    TRUE if we want to order rows based on hierarchical clustering average+correlation distances.
  }
  \item{clustincols}{
    TRUE if we want to order columns based on hierarchical clustering Ward D2 + euclidean distances.
  }
  \item{toPNG}{
    TRUE if we want .png output, FALSE if we want .pdf output
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

