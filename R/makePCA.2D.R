##' Make 2D PCA from a numeric matrix
##'
##' Function that returns a 2D PCA plot with PC1 in x axis and PC2 in y axis.
##' @param est_noctrls numeric matrix of intensities or normalized counts (eg. TMM log2CPM)
##' @param conditions vector with the conditions (same length as colnames(est_noctrls))
##' @param colors vector with the colors assigned to each condition (in order of the unique(conditions))
##' @param dist distance between the dots of the PCa and the labels
##' @param resDir character vector with output results directory. Default is working directory.
##' @param picname character vector with the name of output files.
##' @return files are created in the resDir directory with the picname and .png extension.
##' @author Julia Perera Bel <jperera@imim.es>
##' @export

makePCA.2D <-
  function(est_noctrls, picname, conditions=NULL, colors=NULL, dist=2, 
           resDir=NULL){

    if (!is.null(resDir)){
      resultsDir=resDir
    }
    labels <- colnames(est_noctrls)
    parameters <- setparam(labels)
    summary(pca.filt <- prcomp(t(est_noctrls), center = TRUE, scale.= TRUE )) #38% 
    #summary(pca.filt <- prcomp(t(est_noctrls)))
    vars <- apply(pca.filt$x, 2, var)  
    props <- vars / sum(vars)
    PCAvec <- round(props*100, 2)
    png(file.path(resultsDir,paste(picname,"PCA.2D.png", sep="_")),width=parameters$wid,height=parameters$hei,res=parameters$res)
    
    if (is.null(conditions)) {
      pca2d<-plot(x=pca.filt$x[,1],y=pca.filt$x[,2],  
                           xlab= sprintf("PC1 %.0f%%", PCAvec[1]), 
                           ylab=sprintf("PC2 %.0f%%", PCAvec[2]), 
                           main='PCA', pch=20)
      ##per canviar la dist?ncia del text amb els punts s'ha de sumar a les coordenades pca3d$xyz.convert(pca.filt$x)
      text(pca.filt$x[,1]+dist,pca.filt$x[,2]+dist, labels=rownames(pca.filt$x), cex=parameters$ce+0.2)
      
    }else if(is.null(colors)){
      
      list1 <- unique(as.character(sort(conditions)))
      ColVect <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired")) #20 colors en total
      list2 <- ColVect[1:length(unique(conditions))]
      map = setNames(list2, list1)
      colors <- map[conditions]
      pca2d<-plot(x=pca.filt$x[,1],y=pca.filt$x[,2],  
                  xlab= sprintf("PC1 %.0f%%", PCAvec[1]), 
                  ylab=sprintf("PC2 %.0f%%", PCAvec[2]), 
                  main='PCA', col= colors, pch=20)
      ##per canviar la dist?ncia del text amb els punts s'ha de sumar a les coordenades pca3d$xyz.convert(pca.filt$x)
      text(pca.filt$x[,1]+dist,pca.filt$x[,2]+dist, labels=rownames(pca.filt$x), cex=parameters$ce+0.2,col=colors)
      legend("bottomright",legend=list1,col=list2,pch=16,ncol=1,cex=0.9)
      
    } else {
      list1 <- unique(conditions)
      list2 <- unique(colors)
      map = setNames(list2, list1)
      colors <- map[conditions]
      pca2d<-plot(x=pca.filt$x[,1],y=pca.filt$x[,2],  
                  xlab= sprintf("PC1 %.0f%%", PCAvec[1]), 
                  ylab=sprintf("PC2 %.0f%%", PCAvec[2]), 
                  main='PCA', col= colors, pch=20)
      ##per canviar la dist?ncia del text amb els punts s'ha de sumar a les coordenades pca3d$xyz.convert(pca.filt$x)
      text(pca.filt$x[,1]+dist,pca.filt$x[,2]+dist, labels=rownames(pca.filt$x), cex=parameters$ce+0.2,col=colors)
      legend("bottomright",legend=list1,col=list2,pch=16,ncol=1,cex=0.9)
    }
    dev.off()
  }
