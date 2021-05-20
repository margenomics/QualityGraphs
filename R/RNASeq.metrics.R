##' Make RNASeq metrics plots for alignment assessment
##'
##' Function that returns up to three plots for the assessment of 
##' RNASeq alignment quality. Possible plots are: 1. read mapping stats provided 
##' by STAR (eg. uniquely, multimapped), 2. mapping by genomic regions (eg. coding, 
##' intergenic, intronic) provided by Picard tools CollectRNASeqMetrics and 
##' 3. counts by chromosome (X, Y and M are shown).
##' @param starDir Character vector with the directory where STAR 
##' *Log.final.out files are stored. If more than one run, provide a list with
##' all directories eg. list("/dir/run1","/dir/run2","/dir/run3").
##' @param rnametricsDir Character vector with the directory where 
##'  CollectRnaSeqMetrics output files are stored. Files must be named *.RNA_Metrics.
##' @param counts Data frame with annotations (first 7 columns) and counts 
##' (one column per sample), as read from CountsTable.txt from FeatureCounts. 
##' Only needed if "Counts by Chromosome" plot is desired.
##' @param resDir Character vector with output results directory. Default is working directory.
##' @param picname Character vector with the name of output files.
##' @param toPNG Generate additional plots in png format  Default FALSE.
##' @return Files are created in the resDir directory with the picname and .png or .pdf extensions.
##' @author Julia Perera Bel <jperera@imim.es>
##' @export
##' @importFrom reshape2 melt
##' @import ggplot2 
RNASeq.metrics.plots <- function(starDir=NULL,rnametricsDir=NULL,counts=NULL,
                                 resDir=getwd(),picname="RNASeq.metrics",toPNG=F){
  require(reshape2)
  require(ggplot2)
  
  if(is.null(starDir)&is.null(rnametricsDir)&is.null(counts)){stop("You need to provide at least one of: starDir, rnametricsDir or counts")}

  pdf(file=file.path(resDir, paste(picname,"pdf", sep=".")))
  
  ######################## STAR ALignment
  if(!is.null(starDir)){
    cat("--- STAR plot for",length(starDir),"runs ----------------------- \n")
    
    for (i in 1:length(starDir)){
      cat("Reading files run",i)
      df=read.STAR.Logs(starDir[[i]])
      
      if (i==1){df.1=df}
      if (i>1){df.1=df+df.1}
      
    }
    # Clean variables to use in the plot
    df2 <- df.1[,2:ncol(df.1)]/df.1$`TOTAL READS`
    rownames(df2) <- rownames(df.1)
    df.m=melt(t(df2))
    
    # PLOT STAR METRICS
    print(ggplot(data=df.m, aes(x=Var2, y=value, fill=as.factor(Var1))) +
      geom_bar(stat="identity", position="stack")+theme_minimal()+
      theme(axis.text.x = element_text( angle = 90, hjust = 1),
            plot.title = element_text(color="black", size=14, face="bold",hjust = 0.5))+
      scale_x_discrete(name="Samples") +
      scale_y_continuous(name="% mappings")+
      labs(fill = "")+
      coord_flip()+
      scale_fill_manual(values=rev(rainbow(7)) )+ggtitle ("Read Mapping Stats"))
    
    # Save png
    if (toPNG){ ggsave(file=file.path(resDir, paste(picname,"MappingStats","png", sep=".")))}
    
    }
  
  
  ######################## RNA biotypes 
  if(!is.null(rnametricsDir)){
    cat("--- RNA Metrics plot ---------------------------\n")
    # Read RNA_Metrics files
    cat("Reading files")
    df=read.RNA_Metrics(rnametricsDir)
    
    
    df.m=melt(t(df[,c("RIBOSOMAL_BASES","INTRONIC_BASES", 
                      "INTERGENIC_BASES", "MRNA_BASES")]))
    
    # PLOT GENOMIC REGIONS
    print(ggplot(data=df.m, aes(x=Var2, y=value, fill=as.factor(Var1))) +
      geom_bar(stat="identity", position="stack")+theme_minimal()+
      theme(axis.text.x = element_text( angle = 90, hjust = 1),
            plot.title = element_text(color="black", size=14, face="bold",hjust = 0.5))+
      scale_x_discrete(name="Samples") +
      scale_y_continuous(name="% mappings")+
      labs(fill = "")+
      coord_flip()+
      scale_fill_manual(values=rev(rainbow(5)) )+ggtitle ("Mapping by Genomic Regions"))
    
    # Save png
    if (toPNG){ ggsave(file=file.path(resDir, paste(picname,"GenomicRegions","png", sep=".")))}
  } 
  
  
  
  ########################## COUNTS BY CHROMOSOME
  if(!is.null(counts)){
    cat("--- Chromosome counts plot ---------------------\n")
    cat("Doing plot")
    cat(".")
    # Clean chr annotations (many separated by ;)
    for (i in 1:nrow(counts)){
      chr=paste(unique(strsplit(as.character(counts$Chr[i]),";")[[1]]),collapse = ";")
      counts$Chr_unique[i]=chr
      
    }
    cat(".")
    # Sum gene counts from chrX 
    counts_x=counts[counts$Chr_unique=="chrX",7:(ncol(counts)-1)]
    # Sum gene counts from chrY
    counts_y=counts[counts$Chr_unique=="chrY",7:(ncol(counts)-1)]
    # Sum gene counts from chrM
    counts_m=counts[counts$Chr_unique=="chrM",7:(ncol(counts)-1)]
    
    df=rbind(colSums(counts_x),colSums(counts_y),colSums(counts_m))
    rownames(df)=c("X", "Y","M")
    df.m=melt(df)
    df.m$value=log2(df.m$value+1)
    cat(".")
    
    # PLOT CHR COUNTS
    print(ggplot(data=df.m, aes(x=Var2, y=value, fill=as.factor(Var1))) +
      geom_bar(stat="identity", position=position_dodge())+
      theme(axis.text.x = element_text( size = 10, angle = 45, hjust = 1),
            plot.title = element_text(color="black", size=14, face="bold",hjust = 0.5))+
      scale_x_discrete(name="Samples") +
      scale_y_continuous(name="log2 counts")+
      coord_flip()+
      labs(fill = "Chromosome")+ggtitle("Counts by Chromosomes (X, Y, M)"))
    cat(".\n")
    
    # Save png
    if (toPNG){ ggsave(file=file.path(resDir, paste(picname,"Chromosomes","png", sep=".")))}
  }
  dev.off()
  cat("all plots done!\n")
  
  

}

##' Read log files from STAR alignment
##'
##' @param starDir Character vector with the directory where STAR *Log.final.out files are stored
##' @return Numeric dataframe with one row per file and columns:  number of total reads, uniquely mapped, 
##' multimapped, multmapped (too manany), mapped (too maney missmatches), 
##' unmapped (too short), unmapped (other) and chimeric. 
##' @author Julia Perera Bel <jperera@imim.es>
##' @importFrom dplyr mutate_all
read.STAR.Logs <- function(starDir){
  
  require(dplyr)
  # Prepare matrix to fill in loop
  files=list.files(starDir,pattern = "*Log.final.out")
  df=as.data.frame(matrix(NA,nrow=length(files),ncol=8))
  colnames(df)=c("TOTAL READS","UNIQUELY MAPPED",
                 "MULTIMAPPED","MULTIMAPPED (too many)",
                 "UNMAPPED (too many mm)","UNMAPPED (too short)",
                 "UNMAPPED (other)","CHIMERIC")
  rownames(df)=gsub("Log.final.out","",files)
  
  # Read Log.final.out files
  telliter <- 5
  iter=0
  for (file in gsub("Log.final.out","",files)){
    stats=read.delim(paste0(starDir,"/",file,"Log.final.out"),sep="\t")
    df[file,]=stats[c(4,7,22,24,27,29,31,34),2] # Total reads, uniquely, muli, multi too many, unmapped too many missm, unmapped too short, unmapped other
    #Print progress
    iter=iter+1
    if( iter %% telliter == 0 ) cat(".")
  }
  cat("done! \n")
  df.n <- mutate_all(df, function(x) as.numeric(as.character(x)))
  rownames(df.n)=rownames(df)
  return(df.n)
}


##' Read log files from Picard CollectRnaSeqMetrics tool
##'
##' @param rnametricsDir Character vector with the directory where CollectRnaSeqMetrics 
##' output files are stored
##' @return Numeric dataframe with one row per file and columns:  percentage of ribosomal,  
##' coding, UTR, intronic, intergenic and mRNA bases. 
##' @author Julia Perera Bel <jperera@imim.es>
read.RNA_Metrics <- function(rnametricsDir){
  files=list.files(rnametricsDir,pattern = "*RNA_Metrics")
  df=as.data.frame(matrix(NA,nrow=length(files),ncol=6))
  colnames(df)=c("RIBOSOMAL_BASES","CODING_BASES",
                 "UTR_BASES","INTRONIC_BASES",
                 "INTERGENIC_BASES","MRNA_BASES")
  rownames(df)=gsub(".RNA_Metrics","",files)
  telliter <- 10
  iter=0
  for (file in gsub(".RNA_Metrics","",files)){
    map=read.delim(paste0(rnametricsDir,"/",file,".RNA_Metrics"),
                   nrows = 1,skip = 6)
    map=map[,11:16]
    df[file,]=map
    #Print progress
    iter=iter+1
    if( iter %% telliter == 0 ) cat(".")
  }
  cat("done! \n")
  return(df)
}
