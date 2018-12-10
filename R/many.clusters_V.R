many.clusters_V <-
function(x, resultsDir,Filename, Title, parameters, toPDF=TRUE, 
                           conditions= NULL, colors=NULL){
    #as.dist ho converteix en un objecte dist mentre que dist calcula l'Euclidiana per defecte 
    #dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
    #050214: afegim a la correlaci? la possibilitat d'acceptar NAs
    #x: matriu amb les intensitats(rows) per cada mostra (columns)
    #parameters: param list with the ce values
    #   Per canviar la mida dels labels:
    #       -No colors: S'ha de canviar el par?metre cex general i ajustar els altres en funci? d'aquest (normalment m?s grans de 1)
    #       -Colors: ?nicament s'ha de canviar el par?metre lab.cex de dins de la funci? color_cluster()
    #conditions: Vector amb el nom de les condicions
    labels <- colnames(x)
    
    use.cor="pairwise.complete.obs"
    #Es un parche per poder fer servir aquesta funciÃ³ en linux
    if (toPDF) 
    {pdf(file=file.path(resultsDir, paste(Filename,"pdf", sep=".")))}
    
    if (is.null(conditions)) {
        
        opt<-par(cex.main=2, cex = parameters$ce, cex.lab=1.5, cex.axis=1.5)
        clust.cor.ward<- hclust(as.dist(1-cor(x,use=use.cor)),method="ward.D2")
        plot(clust.cor.ward, main=Title, hang=-1)
        clust.cor.average <- hclust(as.dist(1-cor(x,use=use.cor)),method="average")
        plot(clust.cor.average, main=Title, hang=-1)
        clust.cor.complete <- hclust(as.dist(1-cor(x,use=use.cor)),method="complete")
        plot(clust.cor.complete, main=Title, hang=-1)
        clust.euclid.ward <- hclust(dist(t(x)),method="ward.D2")
        plot(clust.euclid.ward, main=Title, hang=-1)
        clust.euclid.average <- hclust(dist(t(x)),method="average")
        plot(clust.euclid.average, main=Title, hang=-1)    
        clust.euclid.complete <- hclust(dist(t(x)),method="complete")
        plot(clust.euclid.complete, main=Title, hang=-1)           
        par(opt)
        
        
        
    }else if(is.null(colors)){
        
        list1 <- unique(as.character(sort(conditions)))
        ColVect <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired")) #20 colors en total
        list2 <- ColVect[1:length(unique(conditions))]
        map = setNames(list2, list1)
        colors <- map[conditions]
        
        color_cluster <- function(hclus, condition, ce) {
            ##hclus: is the hclust object obtained with the hclust function
            ##condition: a numeric vector of the length=numb of samples, each number is a condition(ex: condition=c(1,1,1,1,2,2,2,2))
            ##Ce: cex parameter that is set in function of the length of the characters of the vector
            
            sampleDendrogram <- as.dendrogram(hclus)
            names(condition) <- hclus$labels
            
            sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch) {
                if(is.leaf(x)){
                    label <- attr(x, "label")
                    attr(x, "nodePar") <- list(lab.col = as.vector(batch[label]), pch="", cex=0.8, lab.cex=ce)
                    
                }
                x
                
            }, condition)
            
        } 
        
        opt<-par(cex.main=1,cex.axis=0.8, cex=0.8)                                            
        clust.cor.ward<- hclust(as.dist(1-cor(x,use=use.cor)),method="ward.D2")
        clust.cor.ward <- color_cluster(clust.cor.ward, colors, parameters$ce)
        plot(clust.cor.ward, main=Title, xlab="Ward")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.cor.average <- hclust(as.dist(1-cor(x,use=use.cor)),method="average")
        clust.cor.average <- color_cluster(clust.cor.average, colors, parameters$ce)
        plot(clust.cor.average, main=Title, xlab="Average")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.cor.complete <- hclust(as.dist(1-cor(x,use=use.cor)),method="complete")
        clust.cor.complete <- color_cluster(clust.cor.complete, colors, parameters$ce)
        plot(clust.cor.complete, main=Title, xlab="Complete")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.euclid.ward <- hclust(dist(t(x)),method="ward.D2")
        clust.euclid.ward <- color_cluster(clust.euclid.ward, colors, parameters$ce)
        plot(clust.euclid.ward, main=Title, xlab="Euclidean Ward")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.euclid.average <- hclust(dist(t(x)),method="average")
        clust.euclid.average <- color_cluster(clust.euclid.average, colors, parameters$ce)
        plot(clust.euclid.average, main=Title, xlab="Euclidean Average")
        legend("topright",legend=list1, cex=parameters$ce+0.2, fill=list2)
        clust.euclid.complete <- hclust(dist(t(x)),method="complete")
        clust.euclid.complete <- color_cluster(clust.euclid.complete, colors, parameters$ce)
        plot(clust.euclid.complete, main=Title, xlab= "Euclidean Complete")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        par(opt)
        
    } else {
        list1 <- unique(as.character(conditions))
        list2 <- unique(colors)
        map = setNames(list2, list1)
        colors <- map[conditions]
        
        color_cluster <- function(hclus, condition, ce) {
            ##hclus: is the hclust object obtained with the hclust function
            ##condition: a numeric vector of the length=numb of samples, each number is a condition(ex: condition=c(1,1,1,1,2,2,2,2))
            ##Ce: cex parameter that is set in function of the length of the characters of the vector
            
            sampleDendrogram <- as.dendrogram(hclus)
            names(condition) <- hclus$labels
            
            sampleDendrogram <- dendrapply(sampleDendrogram, function(x, batch) {
                if(is.leaf(x)){
                    label <- attr(x, "label")
                    attr(x, "nodePar") <- list(lab.col = as.vector(batch[label]), pch="", cex=0.8, lab.cex=ce)
                    
                }
                x
                
            }, condition)
            
        } 
        
        opt<-par(cex.main=1,cex.axis=0.8, cex=0.8)                                            
        clust.cor.ward<- hclust(as.dist(1-cor(x,use=use.cor)),method="ward.D2")
        clust.cor.ward <- color_cluster(clust.cor.ward, colors, parameters$ce)
        plot(clust.cor.ward, main=Title, xlab="Ward")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.cor.average <- hclust(as.dist(1-cor(x,use=use.cor)),method="average")
        clust.cor.average <- color_cluster(clust.cor.average, colors, parameters$ce)
        plot(clust.cor.average, main=Title, xlab="Average")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.cor.complete <- hclust(as.dist(1-cor(x,use=use.cor)),method="complete")
        clust.cor.complete <- color_cluster(clust.cor.complete, colors, parameters$ce)
        plot(clust.cor.complete, main=Title, xlab="Complete")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.euclid.ward <- hclust(dist(t(x)),method="ward.D2")
        clust.euclid.ward <- color_cluster(clust.euclid.ward, colors, parameters$ce)
        plot(clust.euclid.ward, main=Title, xlab="Euclidean Ward")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        clust.euclid.average <- hclust(dist(t(x)),method="average")
        clust.euclid.average <- color_cluster(clust.euclid.average, colors, parameters$ce)
        plot(clust.euclid.average, main=Title, xlab="Euclidean Average")
        legend("topright",legend=list1, cex=parameters$ce+0.2, fill=list2)
        clust.euclid.complete <- hclust(dist(t(x)),method="complete")
        clust.euclid.complete <- color_cluster(clust.euclid.complete, colors, parameters$ce)
        plot(clust.euclid.complete, main=Title, xlab= "Euclidean Complete")
        legend("topright",legend=list1, cex=parameters$ce+0.2,fill=list2)
        par(opt)
        
    }
    if (toPDF) 
    {dev.off()}
    return(list(Corr.ward=clust.cor.ward,Corr.avg=clust.cor.average,Corr.compl=clust.cor.complete,

                Euclid.ward=clust.euclid.ward,Euclid.avg=clust.euclid.average,Euclid.compl=clust.euclid.complete))
}
