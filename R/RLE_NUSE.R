RLE_NUSE <-
function (bpstats, kind, strt=NULL,nd=NULL, parameters, lab, conditions=NULL, colors=NULL) {
    #bpstats: NUSE or RLE statistics obtained by boxplotstats
    #strt: first sample to be analyzed
    #nd: last sample to be analyzed
    #conditions: Vector with the different conditions
    #colors: Vector with the colors assigned to each condition (in order of the unique(conditions))
  
    if (is.null(labels)){
        labels <- gsub(paste("_(",xip,")",sep=""),"",ds$Names,fixed=TRUE)
    }
    
    if (is.null(strt) & is.null(nd)) {
        strt <- 1
        nd <- length(bpstats)
    } else if (is.null(strt)){
        strt <- 1
    } else if (is.null(nd)) {
        nd <- length(bpstats)
    }
    
    fitxer_int<-as.list(1:nd)
    for (i in (strt:nd)){
        
        fitxer_int[[i]]<- bpstats[[i]]$stats
    }
    if (is.null(conditions)) {
        boxplot(fitxer_int[strt:nd],main=paste(kind, "plot"),xaxt="n")
        axis(1,at=1:length(lab),labels=lab,cex.axis=parameters$ce,las=2)
    }else if(is.null(colors)){
        CondNames <- mixedsort(lab)
        CondTable <- data.frame(CondNames,conditions)
        conditions.o <- CondTable[match(lab, CondTable$CondNames),"conditions"]
        list1 <- unique(as.character(sort(conditions.o)))
        ColVect <- c(brewer.pal(8, "Dark2"), brewer.pal(12, "Paired")) #20 colors en total
        list2 <- ColVect[1:length(unique(conditions.o))]
        map = setNames(list2, list1)
        colors <- map[conditions.o]
        boxplot(fitxer_int[strt:nd],main=paste(kind, "plot"),xaxt="n", col=colors)
        axis(1,at=1:length(lab),labels=lab,cex.axis=parameters$ce,las=2)
        legend("topright",legend=list1, cex=parameters$ce+0.2, fill=list2)
    } else {
        CondNames <- mixedsort(lab)
        #CondNames <- lab
        CondTable <- data.frame(CondNames,conditions)
        conditions.o <- CondTable[match(lab, CondTable$CondNames),"conditions"]
        list1 <- unique(conditions.o)
        list2 <- colors
        map = setNames(list2, list1)
        colors <- map[conditions.o]
        boxplot(fitxer_int[strt:nd],main=paste(kind, "plot"),xaxt="n", col=colors)
        axis(1,at=1:length(lab),labels=lab,cex.axis=parameters$ce,las=2)
        legend("topright",legend=list1, cex=parameters$ce+0.2, fill=list2)
    }
}
