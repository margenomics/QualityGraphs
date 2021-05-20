makePCA <-
function(est_noctrls, picname, conditions=NULL, colors=NULL, dist=2, 
                    resDir=NULL,palette=c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired"))){
    #est_noctrls: Matriu amb les intensitats sense controls
    #labels: Array names
    #picname: Name of the picture
    #conditions: Vector with the different conditions
    #colors: Vector with the colors assigned to each condition (in order of the unique(conditions))
    #dist: Distancia de les labels amb els punts del plot PCA
    #resDir: resultsDir per defecte, es el directori de resultats
    require(scatterplot3d)
  
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
    png(file.path(resultsDir,paste(picname,"PCA.png", sep="_")),width=parameters$wid,height=parameters$hei,res=parameters$res)
    
    if (is.null(conditions)) {
        pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],  
                             xlab= sprintf("PC1 %.0f%%", PCAvec[1]), 
                             ylab=sprintf("PC2 %.0f%%", PCAvec[2]), 
                             zlab=sprintf("PC3 %.0f%%", PCAvec[3]), 
                             main='PCA',col.grid="lightblue", 
                             cex.symbols=parameters$ce+0.2)
        ##per canviar la dist?ncia del text amb els punts s'ha de sumar a les coordenades pca3d$xyz.convert(pca.filt$x)
        text(pca3d$xyz.convert(pca.filt$x+dist), labels=rownames(pca.filt$x), cex=parameters$ce+0.2)
    
    }else if(is.null(colors)){
        
        list1 <- unique(as.character(sort(conditions)))
        ColVect <-palette #20 colors en total
        list2 <- ColVect[1:length(unique(conditions))]
        map = setNames(list2, list1)
        colors <- map[conditions]
        pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],
                             xlab=sprintf("PC1 %.0f%%", PCAvec[1]),
                             ylab=sprintf("PC2 %.0f%%", PCAvec[2]),
                             zlab=sprintf("PC3 %.0f%%", PCAvec[3]),
                             main='PCA',col.grid="lightblue",
                             cex.symbols=parameters$ce+0.2, color= colors, pch=16)
        ##per canviar la dist?ncia del text amb els punts s'ha de sumar a les coordenades pca3d$xyz.convert(pca.filt$x)
        text(pca3d$xyz.convert(pca.filt$x+dist), labels=rownames(pca.filt$x), cex=parameters$ce+0.2, col= colors)
        legend("bottomright",legend=list1,col=list2,pch=16,ncol=1,cex=0.9)
        
    } else {
        list1 <- unique(conditions)
        list2 <- unique(colors)
        map = setNames(list2, list1)
        colors <- map[conditions]
        pca3d<-scatterplot3d(x=pca.filt$x[,1],y=pca.filt$x[,2],z=pca.filt$x[,3],  
                             xlab=sprintf("PC1 %.0f%%", PCAvec[1]),
                             ylab=sprintf("PC2 %.0f%%", PCAvec[2]),
                             zlab=sprintf("PC3 %.0f%%", PCAvec[3]),
                             main='PCA',col.grid="lightblue", 
                             cex.symbols=parameters$ce+0.2, color= colors, pch=16)
        ##per canviar la dist?ncia del text amb els punts s'ha de sumar a les coordenades pca3d$xyz.convert(pca.filt$x)
        text(pca3d$xyz.convert(pca.filt$x+dist), labels=rownames(pca.filt$x), cex=parameters$ce+0.2, col= colors)
        legend("bottomright",legend=list1,col=list2,pch=16,ncol=1,cex=0.9)
    }
    dev.off()
}
