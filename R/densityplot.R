densityplot <-
function(ds,strt=NULL,nd=NULL, lab, parameters){
    #NOTA: si no es posa nd i strt sha de posar lab= i parameters=
    #ds: Aroma affimetrix object
    #strt: first sample to be analyzed
    #nd: last sample to be analyzed
    #lab: labels of the sample
    if (is.null(strt) & is.null(nd)) {
        strt <- 1
        nd <- length(ds)
    } else if (is.null(strt)){
        strt <- 1
    } else if (is.null(nd)) {
        nd <- length(ds)
    }
    #############################Set the ylim in the graph #################
    #Agafem el valor d'intensitat m?s elevat i l'arrodonim a l'al?a
    yvalues <- vector()
    for (o in strt:nd){
          fitxer <- getFile(ds,idx=o)
          dl <- log(aroma.affymetrix::getData(fitxer,fields="intensities",
                            indices=1:nbrOfCells(fitxer))$intensities,base=2)
          ycoor <- ceiling(max(density(dl)$y))
          yvalues <- c(yvalues,ycoor )
    }

    ###########################################################
    colors <-rainbow(nd-strt+1)
    #Dibuixar la primera linia
    fitxer <- getFile(ds,idx=strt)
    dl <- log(aroma.affymetrix::getData(fitxer,fields="intensities",
                          indices=1:nbrOfCells(fitxer))$intensities,base=2)

    
    plot(density(dl), col=colors[1], main=paste("Samples", strt, "to", nd), xlab= "log2(y)", ylim=c(0,max(yvalues)), lty=strt, ylab="")
    
    for (e in (strt+1):nd){
          fitxer <- getFile(ds,idx=e)
          dl <- log(aroma.affymetrix::getData(fitxer,fields="intensities",
                            indices=1:nbrOfCells(fitxer))$intensities,base=2)
          lines(density(dl), col=colors[(e-strt)+1], lty = e)
    }
    
    if (length(lab) < 15){
        legend("topright",legend=lab, cex=0.6,col=colors[1:nd], lty= strt:nd, lwd=1.6)
    } else {
        legend("topright",legend=lab, cex=parameters$ce,col=colors[1:nd], lty= strt:nd, lwd=1.6)
    }
    
    
}
