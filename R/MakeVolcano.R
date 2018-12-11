MakeVolcano <- function(TabDiff, adj=TRUE, p.cf, FC.cf, filePath, fileName) {
  #TabDiff: Data frame with the statistics obtained with limma and the function topTable()
  #adj: If TRUE adjusted p.Value in the TabDiff table is used
  #p.cf: p-value cut-off to colour the points in the scatterplot
  #FC.cf: logFC cut-off to colour the points in the scatterplot
  #filePath: output directory
  #fileName: Name of the output file
  
  png(file= paste(file.path(filePath, fileName), "png", sep="."),width=1000,height=1000,res=100)
  
  if (adj) {
    
    plot(TabDiff$logFC, -log10(TabDiff$adj.P.Val), pch=".", 
         cex=4, col=grey(0.75),cex.axis=1.2, las=1, cex.lab=1.5,xlab=expression(paste(log[2], " Fold change")), 
         ylab=expression(paste(-log[10], "adj P.value")))
    points(TabDiff[TabDiff$adj.P.Val < p.cf & TabDiff$logFC < -FC.cf ,"logFC"],
           -log10(TabDiff[TabDiff$adj.P.Val < p.cf & TabDiff$logFC < -FC.cf , "adj.P.Val"]),pch=".", cex=4, col="blue")
    points(TabDiff[TabDiff$adj.P.Val < p.cf & TabDiff$logFC > FC.cf,"logFC"],
           -log10(TabDiff[TabDiff$adj.P.Val < p.cf & TabDiff$logFC > FC.cf, "adj.P.Val"]),pch=".", cex=4, col="red")
    
    abline(h=-log10(max(TabDiff[TabDiff$adj.P.Val < p.cf, "adj.P.Val"])), col= grey(0.1), lty=2)
    abline(v=FC.cf, col= grey(0.1), lty=2 )
    abline(v=-FC.cf, col= grey(0.1), lty=2 )
    
  } else {
    
    plot(TabDiff$logFC, -log10(TabDiff$P.Value), pch=".", 
         cex=4, col=grey(0.75),cex.axis=1.2, las=1, cex.lab=1.5,xlab=expression(paste(log[2], " Fold change")), 
         ylab=expression(paste(-log[10], "P.Value")))
    points(TabDiff[TabDiff$P.Value < p.cf & TabDiff$logFC < -FC.cf ,"logFC"],
           -log10(TabDiff[TabDiff$P.Value < p.cf & TabDiff$logFC < -FC.cf , "P.Value"]),pch=".", cex=4, col="blue")
    points(TabDiff[TabDiff$P.Value < p.cf & TabDiff$logFC > FC.cf,"logFC"],
           -log10(TabDiff[TabDiff$P.Value < p.cf & TabDiff$logFC > FC.cf, "P.Value"]),pch=".", cex=4, col="red")
    
    abline(h=-log10(max(TabDiff[TabDiff$P.Value < p.cf, "P.Value"])), col= grey(0.1), lty=2)
    abline(v=FC.cf, col= grey(0.1), lty=2 )
    abline(v=-FC.cf, col= grey(0.1), lty=2 )
    
  }
  
  dev.off()
  
}