##########################
### PHEAT.MAKE FUNCTION ###
##########################

##' Heatmaps using pheatmap package
##' 
##' Function to make heatmaps from a numeric matrix using pheatmap package
##'
##' This function makes a plot in png format using the results from the differential 
##' expression analysis with "limma" package. It uses the columns "adj.P.Val" or 
##' "P.Value" depending on the value "adj" and the column logFC. Cut-off is used in 
##' order to indicate the DE genes in the scatterplot.
##'
##' @param matrix Numeric matrix with the intensity values to plot.
##' @param cols.s List for specifying annotation_row and annotation_col track 
##' colors manually. It is possible to define the colors for only some of the features.
##' @param filename Output file name.
##' @param filedir Output directory, it can be a path.
##' @param rownames.vect Vector with the names of the genes in the heatmap.
##' @param cex.col Size of the column labels
##' @param scale Character indicating if the values should be centered and 
##' scaled in either the row direction or the column direction, or none. 
##' Corresponding values are "row", "column" and "none"
##' @param toPNG TRUE if we want .png output, FALSE if we want .pdf output.
##' @param annotation_col Similar to annotation_row, but for columns.
##' @param showRownames Indicates wether you want to show rownames or not. Each row defines the features for a specific row.
##' @return The plot is created in the "filedir" directory with the name "filename"
##' in the especified format (default pdf).
##' @author Ariadna Acedo Terrades <aacedo@imim.es>
##' @export
##' @importFrom pheatmap pheatmap
##' @import gplots 

PHEAT.MAKE <- function(matrix, cols.s=NULL, filename, filedir,
         rownames.vect="", cex.col=0.8, scale="none", toPNG=FALSE, annotation_col=NULL, showRownames=F) {
  #Matrix: Numeric matrix with the intensity values to plot
  #cols.s: list for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features. 
  # !! use names of colors in R !!
  
  # exemple of cols.s  :
  # ann_colors = list(
  # Time = c("white", "firebrick"),
  # CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  # GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
  # )
  
  #filename: name of the output file
  #filedir: name of the resultsdir
  #rownames.vect: vector with the names of the genes in the heatmap. 
  #cex.col: size of the column labels
  #toPNG: TRUE if we want .png output, FALSE if we want .pdf output
  #annotation_col : data frame that specifies the annotations shown on left side of the heatmap
  #annotation_row : similar to annotation_row, but for columns.
  
  require(pheatmap)
  require(gplots)
  heatcol<-colorRampPalette(c("blue", "white","red"), space = "rgb")
  
  ##################################################################################
  #In case you want to set the colors
  # library(RColorBrewer)
  # cond <- gsub("_.*","",colnames(num.mtrx), perl=T)
  # pal <- brewer.pal(length(unique(cond)),name="Dark2")
  # map <- setNames(pal,unique(cond))
  # colors <- map[cond]
  ##################################################################################
  
  #################################################################################
  #Set the cex values for rows
  if (length(rownames.vect) == 1) {
    cex.row <- 0.1
  } else if (length(rownames.vect) < 40) {
    cex.row <- 0.8
  } else if (length(rownames.vect) < 70) {
    cex.row <- 0.5
  } else if (length(rownames.vect) < 90) {
    cex.row <- 0.4
  } else if (length(rownames.vect) < 120) {
    cex.row <- 0.3
  } else if (length(rownames.vect) < 190) {
    cex.row <- 0.2
  } else {
    warning(paste(length(rownames.vect),"rows are too many rows in the heatmap for cex.row to be defined. Maximum is 190 rows"), call. = FALSE)
    cex.row <- 0.1
  }
  ################################################################################
  if (toPNG) {
    parameters <- setparam(colnames(data.clus))
    png(file=file.path(filedir, paste(filename,"png", sep=".")),width=parameters$wid,height=parameters$hei,res=parameters$res)
  } else {
    pdf(file=file.path(filedir, paste(filename,"pdf", sep=".")))
  }
  ################################################################################
  
  
  if(is.null(cols.s)) { #no hi ha colors
    if(is.null(annotation_col)){ #no hi ha col annotation
      if(showRownames){ #no hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T, cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = T, treeheight_row = 0)
      }else{ # si hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0)
      }
    }else{ #si hi ha col annotation
      if(showRownames){ #no hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = T, treeheight_row = 0, annotation_col = annotation_col)
      }else{ # si hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0, annotation_col = annotation_col)
      }
    }
  }else{ # si hi ha cols.s
    if(is.null(annotation_col)){ #no hi ha col annotation
      if(showRownames){ #no hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = T, treeheight_row = 0,annotation_colors=cols.s)
      }else{ # si hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0, annotation_colors=cols.s)
      }
    }else{ #si hi ha col annotation
      if(showRownames){ #no hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = T, treeheight_row = 0, annotation_col = annotation_col,
                         annotation_colors=cols.s)
      }else{ # si hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="", cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0, annotation_col = annotation_col,
                         annotation_colors=cols.s)
      }
    }
  }
  
  dev.off()
  
}
