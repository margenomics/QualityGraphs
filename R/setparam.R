setparam <-
function (labels) {
    #labels: nom dels arrays
    #Par?metres:
    #   -res: ?s el gruix de les l?nies del plot. Si volem incrementar la resoluci? tamb? em d'incrementar la mida, sin? les linies seran massa gruixudes
    #   -pointsize: ?s la mida relativa del ce dins de la imatge, en aquest cas es pot canviar quan es crida png() si la mida de les lletres 
    #       a la llegenda no quadra
    
    long <- length(labels)
    
    wid <- 3500
    hei <- 3500
    res <- 400
    
    if (max(nchar(labels)) < 10) {
        ce <- 0.8
    }else if (max(nchar(labels)) < 15){
        ce <- 0.6
    }else {
        ce <- 0.4
    }

    if (long < 35) {
        if (max(nchar(labels)) < 10) {
            ce <- 0.8
        }else if (max(nchar(labels)) < 15){
            ce <- 0.6
        }else {
            ce <- 0.4
        }
    } else if (long < 90) {
        ce <- 0.4
    } else if (long < 140){
        ce <- 0.3
    } else {
        ce <- 0.2
    }
    return(list(wid=wid,hei=hei,res=res,ce=ce))
}
