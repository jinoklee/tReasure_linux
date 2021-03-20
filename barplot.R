#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("--input"), type="character", default="isodecoders.output.txt",
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

  tRNA_aa_codon <- function(name, t){
    name <- gsub("tRNA-","", name)
    out <- data.frame(aa_codon = substr(name, 1, 7))
    out <- data.frame(table(out$aa_codon))
    out <- transform(out, aa = substr(Var1, 1, 3), codon = substr(Var1, 5,7))
    colnames(out) <- c("Var1","Freq","aa","codon")
    out$aa <- as.factor(out$aa)
    out$Group <- rep(t,nrow(out))
    return(out)
  }
  
  detRNA <- read.delim(input)
  
  detRNA <- mutate(detRNA, sig=ifelse(detRNA$FDR<=pval & detRNA$logFC<=-fc , "Down.sig", ifelse(detRNA$FDR<=pval & detRNA$logFC>=fc,"Up.sig", "Not.sig")))
  dw <- subset(detRNA, sig == "Down.sig")
  up <- subset(detRNA, sig == "Up.sig")
  no <- subset(detRNA, sig =="Not.sig")
  
  if(nrow(up)== 0){
    invisible()}else{up <- tRNA_aa_codon(as.character(rownames(up)), "Up_DEtRNA")
    }
  
  if(nrow(dw)== 0){
    invisible()}else{dw <- tRNA_aa_codon(as.character(rownames(dw)), "Down_DEtRNA")
    }
  
  if(nrow(no)== 0){
    invisible()}else{no <- tRNA_aa_codon(as.character(rownames(no)), "Non_DEtRNA")
    }
  
  c <- rbind(up,dw, no) #, rm
  c<- c[order(c$aa),]
  
  empty_bar <- 1
  
  to_add <- data.frame(matrix(NA, empty_bar*nlevels(as.factor(c$aa)), ncol(c)))
  colnames(to_add) <- colnames(c)
  to_add$aa <- rep(levels(as.factor(c$aa)), each=empty_bar)
  to_add$codon <- paste("a", seq(1, nrow(to_add)))
  to_add$Var1 <- paste("a", seq(1, nrow(to_add)))
  c <- rbind(c, to_add)
  
  
  c$aa <- as.character(c$aa)
  c$codon <- as.character(c$codon)
  c$Var1 <- gsub("iMet-CA","iMet-iCAT", c$Var1)
  c$codon <- gsub("-CA","iCAT", c$codon)
  
  c <- c %>% arrange(aa,codon)
  
  c$Freq[is.na(c$Freq)] <- 0
  
  
  base <- data.frame(Var1=c$Var1)
  base$Var1 <- as.character(base$Var1)
  n <- length(grep("iMet", base$Var1))
  base$Var1[grep("iMet", base$Var1)] <- rep("iMet-iCAT", n)
  base_data <- data.frame(Var1= base[!duplicated(base$Var1),])
  out <- data.frame(do.call('rbind',strsplit(as.character(base_data$Var1), split="-")))
  
  base_data$aa <- as.character(out$X1)
  base_data$codon <- as.character(out$X2)
  base_data$id <- as.numeric(seq(1, nrow(base_data)))
  
  
  b1_data<- base_data %>%
    group_by(aa) %>%
    dplyr::summarize(start=min(id), end=max(id)) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  
  b1_data<- b1_data[!grepl("^a", b1_data$aa),]
  
  b2_data <- base_data%>%
    group_by(codon) %>%
    dplyr::summarise(start=min(id), end=max(id)) %>%
    rowwise() %>%
    mutate(title=mean(c(start, end)))
  
  b2_data<- b2_data[!grepl("^a", b2_data$codon),]
  
  png(outfile, width=width, height=height)
  p <- ggplot(c, aes(x=codon, y=Freq, fill=Group))+
    geom_bar( stat="identity")+
    scale_fill_manual(values=c("Up_DEtRNA"="tomato", "Down_DEtRNA"="#67A9CF", "Non_DEtRNA"="grey89", filter="white"), breaks = c("Up_DEtRNA","Down_DEtRNA","Non_DEtRNA"))+
    scale_x_discrete(limits=unique(c$codon))+
    theme_minimal()+
    scale_y_continuous(breaks = seq(0,max(c$Freq),2))+
    theme(axis.text = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "top")+
    ylab("Frequency")+
    xlab("Aminoacid : Anticodon")+
    geom_text(data=b1_data, aes(x = title, y = -2, label=aa), colour ="black", size=4, inherit.aes = FALSE)+
    geom_segment(data=b1_data, aes(x=start, y = -1.5, xend=end, yend=-1.5), colour="black", alpha=1, size=1,inherit.aes = FALSE)+
    geom_text(data=b2_data, aes(x = title, y= -0.65, label=codon), size=3, colour = "black", angle=90, inherit.aes = FALSE)
  print(p)
  dev.off()
