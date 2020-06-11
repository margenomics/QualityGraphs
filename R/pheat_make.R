PHEAT.MAKE <- function(matrix, cols.s=NULL, filename, filedir,
         rownames.vect="", cex.col=0.8, scale="none", toPNG=FALSE, annotation_col=NULL, annotation_row=NULL) {
  #Matrix: Numeric matrix with the intensity values to plot
  #cols.s: Vector with colors corresponding to the conditions in the heatmap
  #filename: name of the output file
  #filedir: name of the resultsdir
  #rownames.vect: vector with the names of the genes in the heatmap. 
  #cex.col: size of the column labels
  #toPNG: TRUE if we want .png output, FALSE if we want .pdf output
  
  
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
  

if(!is.null(cols.s)) { #no hi ha colors
  if(!is.null(annotation_col)){ #no hi ha col annotation
    if(!is.null(annotation_row)){ #no hi ha row annotation
      pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T, cluster_cols = T ,
                        clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0)
    }else{ # si hi ha row annotation
       pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                        clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = heatcol,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0,
                        annotation_row = annotation_row)
    }
  }else{ #si hi ha col annotation
    if(!is.null(annotation_row)){ #no hi ha row annotation
      pheatm<-pheatmap(matrix, col = heatcol(256),
                       cluster_row = T,cluster_cols = T ,
                       clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                       scale=scale, colCol = heatcol,
                       main="",cexRow=cex.row, cexCol=cex.col,
                       show_rownames = F, treeheight_row = 0, annotation_col = annotation_col)
    }else{ # si hi ha row annotation
      pheatm<-pheatmap(matrix, col = heatcol(256),
                       cluster_row = T,cluster_cols = T ,
                       clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                       scale=scale, colCol = heatcol,
                       main="",cexRow=cex.row, cexCol=cex.col,
                       show_rownames = F, treeheight_row = 0,
                       annotation_row = annotation_row, annotation_col = annotation_col)
    }
  }
}else{ # si hi ha cols.s
    if(!is.null(annotation_col)){ #no hi ha col annotation
      if(!is.null(annotation_row)){ #no hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = cols.s,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0)
      }else{ # si hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = cols.s,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0,
                         annotation_row = annotation_row)
      }
    }else{ #si hi ha col annotation
      if(!is.null(annotation_row)){ #no hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = cols.s,
                         main="",cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0, annotation_col = annotation_col)
      }else{ # si hi ha row annotation
        pheatm<-pheatmap(matrix, col = heatcol(256),
                         cluster_row = T,cluster_cols = T ,
                         clustering_method = "ward.D2",clustering_distance_cols = "correlation",
                         scale=scale, colCol = cols.s,
                         main="", cexRow=cex.row, cexCol=cex.col,
                         show_rownames = F, treeheight_row = 0,
                         annotation_row = annotation_row, annotation_col = annotation_col)
      }
    }
}
  
dev.off()

}
