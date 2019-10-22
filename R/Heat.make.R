Heat.make <- function(data.clus, cols.s=NULL, filename, filedir,
                        rownames.vect="", cex.col=0.8, clustinrows=T,
                        clustincols=T, toPNG=FALSE) {
  #data.clus: Numeric matrix with the intensity values to plot
  #cols.s: Vector with colors corresponding to the conditions in the heatmap
  #filename: name of the output file
  #filedir: name of the resultsdir
  #rownames.vect: vector with the names of the genes in the heatmap. If rownames.vect="" no rownames are shown
  #cex.col: size of the column labels
  #clustinrows:TRUE if we want to order rows based on hierarchical clustering average+correlation distances
  #clustincols: TRUE if we want to order columns based on hierarchical clustering Ward D2 + euclidean distances
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
  #Set dendrogram for rows
  if (clustinrows){
    clust.rows <- as.dendrogram(hclust(as.dist(1-cor(t(data.clus))),method="average"))
  } else {
    clust.rows <- F
  }
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
  
  if(!is.null(cols.s)) {
    if(clustincols) {
      #Dendrograma de les columnes
      clust.cols <- as.dendrogram(hclust(dist(t(data.clus)),method="ward.D2"))
      #Fem el plot
      heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                       dendrogram="column", Colv=clust.cols,
                       Rowv=clust.rows,labRow=rownames.vect,
                       scale="row",cexRow=cex.row, cexCol=cex.col, ColSideColors= cols.s,
                       main="",key=TRUE,keysize=1,density.info="none",trace="none")
    } else {
      #Fem el plot
      heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                       dendrogram="none", Colv=F,
                       Rowv=clust.rows,labRow=rownames.vect,
                       scale="row",cexRow=cex.row, cexCol=cex.col, ColSideColors= cols.s,
                       main="",key=TRUE,keysize=1,density.info="none",trace="none")
    }
    
  } else {
    
    if(clustincols) {
      #Dendrograma de les columnes
      clust.cols <- as.dendrogram(hclust(dist(t(data.clus)),method="ward.D2"))
      #Fem el plot
      heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                       dendrogram="column", Colv=clust.cols,
                       Rowv=clust.rows,labRow=rownames.vect,
                       scale="row",cexRow=cex.row, cexCol=cex.col,
                       main="",key=TRUE,keysize=1,density.info="none",trace="none")
    } else {
      #Fem el plot
      heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                       dendrogram="none", Colv=F,
                       Rowv=clust.rows,labRow=rownames.vect,
                       scale="row",cexRow=cex.row, cexCol=cex.col,
                       main="",key=TRUE,keysize=1,density.info="none",trace="none")
    }
  }
  dev.off()
}