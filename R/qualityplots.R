qualityplots <-
function(ds, picname, estimates_m=NULL, est_noctrls, labels=NULL, ds.N=NULL, 
                         plmTr=NULL, conditions=NULL,colors=NULL, groupn=NULL, all=FALSE, 
                         Density=TRUE, Boxplots=TRUE, Clusters=TRUE, PCA=TRUE, estimates=FALSE, 
                         noctrls=TRUE, resDir=NULL) {
    library(RColorBrewer)
    require(gtools)
    #ds: Aroma affimetrix object
    #picname: Name of the picture
    #estimates_m: matriu amb les intensitats
    #est_noctrls: Matriu amb les intensitats sense controls
    #ds.N: Aroma affimetrix object normalized
    #plmTr: Prove level model object to generate NUSE and RLE
    #conditions: Vector with the different conditions. It has to be a vector of characters!!
    #colors: Vector with the colors assigned to each condition (in order of the unique(conditions))
    #groupn: subset to represent in the density plot if all=TRUE
    #all: If TRUE a summary with all the density plots in the same picture is generated
    #Density: If TRUE the density plots are generated
    #Boxplots: If TRUE the boxplots are generated
    #Clusters: If TRUE clusters are generated
    #PCA: If TRUE PCA are generated, only with est_noctrls
    #estimates: If TRUE clusters with all the estimates are generated
    #noctrls: If TRUE clusters with no controls are generated. Use this matrix for other analyses
    #resDir: Directori de resultats per defecte es el resultsDir
    if (!is.null(resDir)){
      resultsDir=resDir
    }
    #Density plots per separat (subsets of 16 samples) i tots en una mateixa imatge
    if(Density){
        densityplot_all(ds, picname=picname, groupn= groupn, all=all, resDir=resultsDir)
    }
    #Boxplots
    if(Boxplots){
        boxplot_all(ds, ds.N=ds.N, plmTr=plmTr, labels=labels, picname=picname, 
                    conditions=conditions, colors=colors, resDir=resultsDir)
    }
    #Clusters with and without controls
    if(Clusters){
        clusterdend(estimates_m, est_noctrls, picname=picname, conditions=conditions, 
                    colors=colors, estimates=estimates, noctrls=noctrls, resDir=resultsDir)
    }
    
    # PCA 
    if(PCA){
      makePCA(est_noctrls, picname=picname, conditions=conditions, colors=colors, 
              resDir=resultsDir)
    }
}
