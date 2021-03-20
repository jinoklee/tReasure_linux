#!/usr/bin/env Rscript

library(optparse)


option_list <- list(
  make_option(c("--bed"), type="character", default="Hsapi38.tRNAscan_pre-tRNAs.bed12",
	     help="bed file of tRNAs [default = %default]", metavar="character"),
  make_option(c("--bam"), type="character", default=NULL, 
              help="premapped bam file", metavar="character"),
  make_option(c("fq","--fastq"), type="character", default=NULL, 
             help="trimmed fastq file", metavar="character"),
  make_option(c("-o","--output"), type="character", default="mature.fastq", 
              help="output file name [default= %default]", metavar="character")
); 


opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$bed)){
	          print_help(opt_parser)
  stop("Argument must be supplied bed file", call.=FALSE)
}

if (is.null(opt$bam)){
		print_help(opt_parser)
  stop("Argument must be supplied bam file", call.=FALSE)
}

if (is.null(opt$fastq)){
  print_help(opt_parser)
  stop("Argument must be supplied trimmed fastq file", call.=FALSE)
}



bed <- opt$bed
bam <- opt$bam
fq <- opt$fastq
out <- opt$output

library(dplyr)
library(stringr)
library(Rsamtools)
library(seqinr)
library(ShortRead)
library(tidyverse)

bed <- read.table(bed, sep = "\t")
len <- data.frame(refname = bed$V4, reflen= bed$V11)
len %>% mutate_if(is.factor, as.character) -> len

for ( i in 1:nrow(len)){
  if(length(grep("[,]",len$reflen[i])) != 0){
    d <- data.frame(do.call('rbind',strsplit(as.character(len$reflen[i]), split=",")), stringsAsFactors = F)
    d <- lapply(d,as.numeric)
    len$reflen[i]<- as.character(d[[1]] + d[[2]])
  }
}
len$reflen <- as.numeric(len$reflen)

# read bam
  print(paste("load :", bam))
  b1 <- scanBam(bam)
  b1 <- data.frame(b1[[1]], stringsAsFactors = F)
  b1$rname <- as.character(b1$rname)
  tb1 <- b1[grep(".tRNA", b1$rname),]
  gb1 <- b1[!grepl(".tRNA", b1$rname),]
  
  tb1 <- filter(tb1, tb1$pos > 50 )
  tb1$refname <- gsub("\\::.*", "", tb1$rname)
  
  print(paste("loading", fq))

  mfun <- function(x){
    tid <- gsub("\\s.*","" ,data.frame(id(x))[,1])
    x[tid%in%tb1$qname]
  }
  
# filtering 
  # function ---------------
  matcher <- function(pattern, x) {
    ind = gregexpr(pattern, x)[[1]]
    start = as.numeric(ind)
    end = start + attr(ind, "match.length")- 2
    apply(cbind(start,end), 1, function(y) substr(x, start=y[1], stop=y[2]));
  }
  doone <- function(c, cigar) {
    pat <- paste("\\d+", c , sep="")
    sum(as.numeric(matcher(pat, cigar)), na.rm=T)
  }
  ## takes a cigar string and parses it, not very fast but...
  cigarsums <- function(cigar, chars=c("D","I")) {
    sapply (chars, doone, cigar)
  }
  
  
  # cigar ------------
  con <- unique(c(grep("D",tb1$cigar),grep("I", tb1$cigar)))
  if(length(con) == 0 ){
    tb1$end <- tb1$pos + tb1$qwidth -1
  }else{
    tb1t <- tb1[con,]
    tb1o <- tb1[-con,]
    tb1tt <- sapply(tb1t$cigar, cigarsums)
    tb1tt <- data.frame(t(tb1tt))
    tb1t$end <- tb1t$pos + tb1t$qwidth -1 + tb1tt$D - tb1tt$I
    tb1o$end <- tb1o$pos + tb1o$qwidth -1 
    
    tb1 <- rbind(tb1t, tb1o)
  }
  df <- left_join(tb1, len)
  df <- filter(df, end <= (df$reflen-44))
  
  
  # filter function
  fun <- function(x){
    tid <- gsub("\\s.*","" ,data.frame(id(x))[,1])
    x[tid%in%df$qname]
  }
  
 # samve matrue fastq 
  filterFastq(fq, destinations = out, filter = fun , compress= FALSE)  
  print(paste("compelete :", out))

