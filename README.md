# CLOSHA  

### 새롭게 등록해야할 프로그램  
  * removeNontRNA.R  
  * DEtRNA_edge.R

### R 코드 실행시 필요 library  
~~~   
  library(optparse)  
  library(dplyr) 
  library(tidyverse)  
  library(stringr)  
  library(Rsamtools)  
  library(seqinr)  
  library(ShortRead)  
  library(edgeR)  
  library(ggplot2)  
  library(plotrix)   
~~~

### tRNA expression analysis pipeline
![Pipeline](./bioexpress_pipeline.png)
  
  
*****************  
### USAGE  
### removeRead.R   
~~~
Usage: Rscript removeNontRNA.R [options]
	Rscript removeNontRNA.R --bed Hsapi38.tRNAscan_pre-tRNAs.bed12 \
			     --bamdir ${path}/firstmapping \
			     --fastqdir ${path}/trimmed_fastq \
			     --outdir ${path}/removed_read 


Options:
	--bed=CHARACTER
		bed file of tRNAs

	--bamdir=CHARACTER
		directory of 1st round mapped bam files

	--fastqdir=CHARACTER
		directory trimmed fastq files

	--outdir=CHARACTER
		directory output files

	-h, --help
		Show this help message and exit
~~~   
 
### DEtRNA.R
~~~
Usage: DEtRNA_edge.R [options]
	Rscript DEtRNA_edge.R --control ./controlDir/ --test ./testDir/ \
			--stat Exact \
			--adj BH \
			--pvalue 0.05 \
			--foldchange 1 \
			--prefix BH


Options:
	--control=CHARACTER
		directory contained the readcount table of the control group

	--test=CHARACTER
		directory contained the readcount table of the test group

	--stat=CHARACTER
		The method of statistical analysis.
               This must be one of the stringss (Exact, Quasi-likelihood, likelihood) 

	--adj=CHARACTER
		The adjust method of statistical analysis.
               This must be one of the strings (holm, hochberg, hommel, bonferroni, BH, BY, fdr, none)

	--pvalue=NUMERIC
		p.value threshold for statistical significance (0.01, 0.05)

	--foldchange=NUMERIC
		foldchange threshold for statistical significance (1, 1.5, 2)

	--prefix=CHARACTER
		The prefix output file name

	-h, --help
		Show this help message and exit   
~~~
