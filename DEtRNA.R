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
  make_option(c("--pvalue"), type="character", default=NULL, 
              help="p.value threshold for statistical significance", metavar="numeric"),
  make_option(c("--foldchange"), type="character", default=NULL, 
              help="foldchange threshold for statistical significance", metavar="numeric"),
  make_option(c("-pre","--prefix"), type="character", default="output.png", 
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


control <- as.character(opt$control)
test <- as.character(opt$test)
stat <- opt$stat
adj <- opt$adj
pval <- as.numeric(opt$pvalue)
fc <- as.numeric(opt$foldchange)
prefix <- as.character(opt$prefix)
width <- as.character(opt$width)
height <- as.character(opt$height)

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


write.table(tstat$table, paste0(prefix, ".individual.tRNA.txt"), sep = "\t", quote = F)
write.table(dstat$table, paste0(prefix, ".isodecoders.txt"), sep = "\t", quote = F)
write.table(astat$table, paste0(prefix, ".isoacceptor.txt"), sep = "\t", quote = F)

#1.Volcano

detRNA <- tstat$table

detRNA <- mutate(detRNA, sig=ifelse(detRNA$FDR<=pval & detRNA$logFC<=-fc , "Down.sig", ifelse(detRNA$FDR<=pval & detRNA$logFC>=fc,"Up.sig", "Not.sig")))
png(outfile)
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
rm(detRNA)

#2.Barplot

detRNA <- dstat$table

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

png(outfile)
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

rm(detRNA)

#3.Pyramid

detRNA <- astat$table

tRNA_aa <- function(name, t){
  name <- gsub("tRNA-","", name)
  out <- data.frame(aa= substr(name, 1, 3))
  out <- data.frame(table(out$aa))
  colnames(out) <- c("Var1",t)
  return(out)
}

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
png(outfile)
pyramid.plot(geneplot$Down_DEtRNA, geneplot$Up_DEtRNA,labels= geneplot$Var1,lxcol="#67A9CF", rxcol="#EF8A62",unit = "Freqency",gap=0.3, space=0.15, top.labels = c("Down_DEtRNAs", "tRNA-AA","Up_DEtRNAs"),laxlab=c(0,1,2,3), raxlab=c(0,1,2,3))
dev.off()




