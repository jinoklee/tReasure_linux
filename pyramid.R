#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("--input"), type="character", default="isoacceptor.output.txt",
              help="The output matrix of statistical analysis [default= %default]", metavar="character"),
  make_option(c("--pvalue"), type="character", default=NULL, 
              help="p.value threshold for statistical significance", metavar="numeric"),
  make_option(c("--foldchange"), type="character", default=NULL, 
              help="foldchange threshold for statistical significance", metavar="numeric"),
  make_option(c("-o","--output"), type="character", default="output.png", 
              help="The output file name [default= %default]", metavar="character"),
  make_option(c("--width"), type="numeric",default=800,
              help = "The value of width(px) [default= %default]", metavar = "numeric"),
  make_option(c("--height"), type="numeric",default=500,
              help = "The value of height(px) [default= %default]", metavar = "numeric")
); 


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Argument must be supplied input file", call.=FALSE)
}

if (is.null(opt$pvalue)){
  print_help(opt_parser)
  stop("Argument must be supplied the threshold of p.value", call.=FALSE)
}

if (is.null(opt$foldchange)){
  print_help(opt_parser)
  stop("Argument must be supplied the threshold of foldchange ", call.=FALSE)
}

input <- opt$input
pval <- as.numeric(opt$pvalue)
fc <- as.numeric(opt$foldchange)
outfile <- opt$output
width <- opt$width
height <- opt$height

library(dplyr)
library(ggplot2)
library(plotrix)
# install.packages("plotrix")

  tRNA_aa <- function(name, t){
    name <- gsub("tRNA-","", name)
    out <- data.frame(aa= substr(name, 1, 3))
    out <- data.frame(table(out$aa))
    colnames(out) <- c("Var1",t)
    return(out)
  }
  detRNA <- read.delim(input)
  detRNA <- mutate(detRNA, sig=ifelse(detRNA $FDR<=pval & detRNA $logFC<=-fc , "Down.sig", ifelse(detRNA$FDR<=pval & detRNA$logFC>=fc,"Up.sig", "Not.sig")))
  dw <- subset(detRNA, sig == "Down.sig")
  up <- subset(detRNA, sig == "Up.sig")
  
  if(nrow(up)==0){
    invisible()}else{up <- tRNA_aa(as.character(rownames(up)), "Up_DEtRNA")
    }
  if(nrow(dw)==0){
    invisible()}else{dw <- tRNA_aa(as.character(rownames(dw)), "Down_DEtRNA")
    }

  geneplot <- full_join(dw, up)
  png(outfile, width=width, height=height)
  pyramid.plot(geneplot$Down_DEtRNA, geneplot$Up_DEtRNA,labels= geneplot$Var1,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Freqency",gap=0.3, space=0.15, top.labels = c("Down_DEtRNAs", "tRNA-AA","Up_DEtRNAs"),laxlab=c(0,1,2,3), raxlab=c(0,1,2,3))
  dev.off()
  
  
  
  
