Heat.make <- function(data.clus, cols.s, filename, filedir, rownames=TRUE, rownames.vect=NULL, cex.col=0.8, cex.row=0.3, toPNG=FALSE) {
  #data.clus: matrix to plot in the heatmap
  #cols.s: vector of colors corresponding to the columns in the matrix data.clus 
  #filename: Name of the output file in .pdf format
  #filedir: Name of the resultsDir
  #rownames: If TRUE row names are shown in the heatmap
  #rownames.vect: Vector with row names to be shown in the heatmap
  #cex.col: size of colnames
  #cex.row: size of rownames, if rownames=TRUE
  #toPNG: If TRUE heatmap is generated in png format
  require(gplots)
  
  heatcol<-colorRampPalette(c("blue", "white","red"), space = "rgb")
  #Dendrograma de les columnes
  clust.cols <- hclust(dist(t(data.clus)),method="ward.D2")
  #Dendrograma de les rows
  clust.rows <- hclust(as.dist(1-cor(t(data.clus))),method="average")
  if(rownames) {
    #Fem el plot
    if (toPNG) {
       parameters <- setparam(colnames(data.clus))
       png(file=file.path(filedir, paste(filename,"pdf", sep=".")),width=parameters$wid,height=parameters$hei,res=parameters$res)
    } else {
       pdf(file=file.path(filedir, paste(filename,"pdf", sep=".")))
    }
    heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                     dendrogram="column", Colv=as.dendrogram(clust.cols),
                     Rowv=as.dendrogram(clust.rows),labRow=rownames.vect,
                     scale="row",cexRow=cex.row, cexCol=cex.col, ColSideColors= cols.s,
                     main="",key=TRUE,keysize=1,density.info="none",trace="none")
    dev.off()
  } else {
    #Fem el plot
    if (toPNG) {
       parameters <- setparam(colnames(data.clus))
       png(file=file.path(filedir, paste(filename,"pdf", sep=".")),width=parameters$wid,height=parameters$hei,res=parameters$res)
    } else {
       pdf(file=file.path(filedir, paste(filename,"pdf", sep=".")))
    }
    heatm<-heatmap.2(as.matrix(data.clus), col = heatcol(256),
                     dendrogram="column", Colv=as.dendrogram(clust.cols),
                     Rowv=as.dendrogram(clust.rows),labRow="",
                     scale="row",cexRow=cex.row, cexCol=cex.col, ColSideColors= cols.s,
                     main="",key=TRUE,keysize=1,density.info="none",trace="none")
    dev.off()
  }
  
  
}
