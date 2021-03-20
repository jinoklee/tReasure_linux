#!/usr/bin/env Rscript

library(optparse)


option_list <- list(
  make_option(c("--input"), type="character", default="individual.tRNA.output.txt",
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

detRNA <- read.delim(input)
  
detRNA <- mutate(detRNA, sig=ifelse(detRNA$FDR<=pval & detRNA$logFC<=-fc , "Down.sig", ifelse(detRNA$FDR<=pval & detRNA$logFC>=fc,"Up.sig", "Not.sig")))
  png(outfile, width=width, height=height)
  ggplot(detRNA, aes(x = logFC, y = -log10(FDR)))+
    geom_point(aes(col=sig)) +
    xlab(" log2 Fold Change") +
    ylab("-log10 Adjusted P value ") +
    geom_vline(xintercept = c(-1.5,1.5),col = "red",linetype = "dotted",size = 0.5) +
    geom_hline(yintercept = c(-log10(0.05)),col = "red", linetype = "dotted",size = 0.5) +
    theme_bw() +
    theme(legend.position = "none")+
    scale_colour_manual(values = c("Not.sig"="grey", "Up.sig"="tomato","Down.sig"="#67A9CF"))
  dev.off()

