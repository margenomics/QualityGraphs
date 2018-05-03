boxplot_all <-
function(ds, ds.N=NULL, plmTr=NULL, picname, labels=labels, conditions=NULL,
                         colors=NULL, raw.bp=TRUE, RMA.bp=TRUE, RLE=TRUE, NUSE=TRUE,
                         resDir=NULL){
    #ds: Aroma affimetrix object
    #ds.N: Aroma affimetrix object normalized
    #plmTr: Prove level model object to generate NUSE and RLE
    #labels: Array names
    #picname: Name of the picture
    #conditions: Vector with the different conditions
    #raw.bp: If FALSE doesn't plot the raw intensities boxplot
    #RMA.bp: If false doesn't plot normalized boxplot
    #RLE: If false doesn't plot the RLE plot
    #NUSE: If false doesn't plot the NUSE plot
    library(aroma.affymetrix)
    if (is.null(labels)){
        labels <- gsub(paste("_(",xip,")",sep=""),"",ds$Names,fixed=TRUE)
    }
    if (!is.null(resDir)){
      resultsDir=resDir
    }
    parameters <- setparam(labels)
    
    #Boxplot rawdata
    if (raw.bp){
        png(file.path(resultsDir, paste(picname,"Boxplots.png", sep="_")),width=parameters$wid,height=parameters$hei,res=parameters$res)
        boxplot_raw(ds, parameters=parameters, lab =labels, conditions=conditions, colors=colors)
        dev.off()
    }
    
    #Boxplot RMA
    if (RMA.bp){
        if (is.null(ds.N)){
            bc <- RmaBackgroundCorrection(ds)
            ds.BC<-process(bc,verbose=verbose)
            qn <- QuantileNormalization(ds.BC, typesToUpdate="pm")
            ds.N <- process(qn, verbose=verbose)
        }
        png(file.path(resultsDir,paste(picname, "BoxplotsRMA.png", sep="_")),width=parameters$wid,height=parameters$hei,res=parameters$res)
        boxplot_norm(ds.N, parameters=parameters, lab =labels, conditions=conditions, colors=colors)
        dev.off()
    }
    #RLE plot
    if (RLE) {
        if (is.null(plmTr)){
            plmTr <- ExonRmaPlm(ds.N, mergeGroups=TRUE)
            fit(plmTr,verbose=verbose)
        }
        
        png(file.path(resultsDir,paste(picname, "RLE.png", sep="_")),width=parameters$wid,height=parameters$hei,res=parameters$res)
        plmTr.QA<-QualityAssessmentModel(plmTr) #per a generar el RLE i NUSE
        bpstats1 <- boxplotStats( getChipEffectSet(plmTr.QA), type="RLE" )
        RLE_NUSE(bpstats1, "RLE", parameters=parameters, lab=labels, conditions=conditions, colors=colors)
        dev.off()
    }
    
    #NUSE plot
    if(NUSE) {
        png(file.path(resultsDir,paste(picname, "NUSE.png", sep="_")),width=parameters$wid,height=parameters$hei,res=parameters$res)
        bpstats2 <- boxplotStats( getChipEffectSet(plmTr.QA), type="NUSE" )
        RLE_NUSE(bpstats2, "NUSE", parameters=parameters, lab=labels, conditions=conditions, colors=colors)
        dev.off()
    }

}
