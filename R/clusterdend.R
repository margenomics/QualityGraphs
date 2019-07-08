clusterdend <-
function(estimates_m=NULL, est_noctrls=NULL, picname, conditions=NULL, 
                        colors=NULL, estimates=FALSE, noctrls=TRUE, resDir=NULL, toPNG=TRUE) {
    #picname: Name of the picture
    #conditions: Vector with the different conditions
    #estimates_m: matriu amb les intensitats
    #est_noctrls: Matriu amb les intensitats sense controls
    #colors: Vector with the colors assigned to each condition (in order of the unique(conditions))
    
  #Afegim aquest loop per tal de poder determinar el directori de resultats
    if (!is.null(resDir)){
      resultsDir=resDir
    }
    if(estimates) {
      labels <- colnames(estimates_m)
      parameters <- setparam(labels)
      many.clusters_V(estimates_m,resultsDir,paste(picname, "ClustersWithCtrls", sep="_"),"Clusters with ctrls", parameters=parameters, conditions=conditions, colors=colors)
    } 
    
    if(noctrls) {
      labels <- colnames(est_noctrls)
      parameters <- setparam(labels)
      many.clusters_V(est_noctrls,resultsDir,paste(picname, "Clusters", sep="_"),"Clusters", parameters=parameters, toPNG=TRUE, conditions=conditions, colors=colors)
    }
}
