\name{MakeVolcano}
\alias{MakeVolcano}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make volcano plots from the results of DE expression
}
\description{
  This function makes a volcano plot using the statistics in the matrix obtained with topTable() in limma. DE genes are coloured (positive or negative expressed) using p-value and logFC cut offs. 
}
\usage{
  MakeVolcano(TabDiff, adj=TRUE, p.cf, FC.cf, filePath, fileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TabDiff}{
    Matrix with the statistics obtained with limma and the function topTable()
  }
  \item{adj}{
    If TRUE adjusted p.Value in the TabDiff table is used
  }
  \item{p.cf}{
    P-value cut-off to colour the points in the scatterplot
  }
  \item{FC.cf}{
    LogFC cut-off to colour the points in the scatterplot
  }
  \item{filePath}{
    Output directory
  }
  \item{fileName}{
    Name of the output file
  }
  
}
\details{
  This function makes a plot in png format using the results from the differential expression analysis with "limma" package. It uses the columns "adj.P.Val" or "P.Value" deppending on the value "adj" and the column logFC. Cut-off is used in order to indicate the DE genes in the scatterplot.
}
\value{
  The plot is created in the "filePath" directory with the name "fileName" in ".png" format.
}
\author{
  Magdalena Arnal Segura <marnal@imim.es>.
}

