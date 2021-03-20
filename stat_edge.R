#!/usr/bin/env Rscript

library(optparse)


option_list <- list(
  make_option(opt_str="--control", type="character", default=NULL, 
              help="The directory contained the readcount table of the control group", metavar="character"),
  make_option("--test", type="character", default=NULL,
              help="The directory contained the readcount table of the test group", metavar="character"),
  make_option("--stat", type="character", default=NULL,
              help="The method of statistical analysis.
               This must be one of the stringss (Exact, Quasi-likelihood, likelihood) ", metavar="character"),
  make_option("--adj", type="character", default="BH",
              help="The adjust method of statistical analysis.
               This must be one of the strings (holm, hochberg, hommel, bonferroni, BH, BY, fdr, none) [default= %default]",
              metavar="character"),
  make_option(c("-o","--output"), type="character", default="output.txt",
              help="The output file name [default= %default]", metavar="character")
);


#
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$control)){
  print_help(opt_parser)
  stop("Argument must be supplied the input path of control group", call.=FALSE)
}

if (is.null(opt$test)){
  print_help(opt_parser)
  stop("Argument must be supplied the output path of test group", call.=FALSE)
}

if (is.null(opt$stat)){
  print_help(opt_parser)
  stop("Argument must be supplied the stat method", call.=FALSE)
}

if (is.null(opt$adj)){
  print_help(opt_parser)
  stop("Argument must be supplied the adj.pval method", call.=FALSE)
}

  control <- as.character(opt$control)
  test <- as.character(opt$test)
  stat <- opt$stat
  adj <- opt$adj
  out <- opt$output

  library(dplyr)
  library(ggplot2)
  library(edgeR)

  ### combined count matrix 
  combine.df <- function(path){
    list <- list.files(path, recursive = F, full.names = F)
    c <- read.table(file.path(path,list[1]), sep = "\t", stringsAsFactors = F)
    c <- c[,1, drop = FALSE]
    for(i in list){
      df <- read.table(file.path(path,i), sep = "\t", stringsAsFactors = F)
      df <- df[,c(1,3)]
      colnames(df)[2] <- gsub(".expression.txt","",i)
      c <- left_join(c, df)
    }
    c <- c[grep("high", c$V1),]
    return(c)
  }
  
  c <- combine.df(control)
  t <- combine.df(test)
  
  name <- data.frame(do.call('rbind', strsplit(as.character(c$V1), split = "[.]")))[,2]
  trna <- data.frame(do.call('rbind', strsplit(as.character(name), split = "[_]")))[,2]
  iso <- data.frame(do.call('rbind', strsplit(as.character(name), split = "[_]")))[,3]
  c <- c[,-1]
  t <- t[,-1]
  
  ### sample list
  saminfo <- data.frame(samplename=c(colnames(c), colnames(t)) , group=c(rep("control", ncol(c)), rep("test", ncol(t))))
  
  ### individual tRNA count matrix
  t.count <- cbind(c,t)
  rownames(t.count) <- trna
  
  ### isodecoder count matrix
  count <- t.count
  count$iso <- iso

  sum.df <- function(count){
    mat <- c()
    for(i in unique(count$iso)){
      df <- count[grep(i, count$iso),]
      if(nrow(df) > 1){
        df <- colSums(subset(df, select = -c(iso)))
        df <- data.frame(t(df))
        df$iso <- i
      }else{
        rownames(df)<- NULL
      }
      mat <- rbind(mat, df) 
    }
    rownames(mat) <- mat$iso
    mat <- subset(mat, select=-c(iso))
    return(mat)
  }

  d.count <- sum.df(count)
  
  ### isoaccepotr count matrix
  count <- d.count
  count$iso <- substr(rownames(count), 1, 12)
  a.count <- sum.df(count)
  
  ### make design
  case <- factor(saminfo$group)
  case <- relevel(case, ref="control")
  design <- model.matrix(~case)
  coef = 2
  rownames(design) <- saminfo$samplename

  
  ## stat 
  
  stat.out <- function(count, stat, adj){
    y <- DGEList(count, group = case)
    keep <- filterByExpr(y)
    y <- y[keep,,keep.lib.sizes=FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust = TRUE)
    
    if(stat == "Exact"){
      test <- exactTest(y)
    }else if(stat == "Quasi-likelihood"){
      fit <- glmQLFit(y, design)
      test <- glmQLFTest(fit, coef=coef)
    }else{
      fit <- glmFit(y, design)
      test <- glmLRT(fit, coef=coef)
    }
    topTags(test, n = Inf, adjust.method=adj)
  }
  
  tstat <- stat.out(t.count, stat,adj)
  dstat <- stat.out(d.count,stat,adj)
  astat <- stat.out(a.count,stat,adj) 
  
  
  write.table(tstat$table, paste0("individual.tRNA.", out), sep = "\t", quote = F)
  write.table(dstat$table, paste0("isodecoders.", out), sep = "\t", quote = F)
  write.table(astat$table, paste0("isoacceptor.", out), sep = "\t", quote = F)
