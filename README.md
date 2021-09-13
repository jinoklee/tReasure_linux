# CLOSHA  

### 새롭게 등록해야할 프로그램  
  * removeRead.R  
  * DEtRNA.R

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
Usage: Rscript removeRead.R [options]
	Rscript removeRead.R --bed Hsapi38.tRNAscan_pre-tRNAs.bed12 \
			     --bamdir ${path}/firstmapping \
			     --fastqdir ${path}/trimmed_fastq \
			     --outdir ${path}/removed_read 


Options:
	--bed=CHARACTER
		bed file of tRNAs

	--bamdir=CHARACTER
		directory of premapped bam files

	--fastqdir=CHARACTER
		directory trimmed fastq files

	--outdir=CHARACTER
		directory output files

	-h, --help
		Show this help message and exit
~~~   
 
### DEtRNA.R
~~~
Usage: DEtRNA.R [options]
	Rscript DEtRNA.R --control ./controlDir/ --test ./testDir/ \
			--stat Exact \
			--adj BH \
			--pvalue 0.05 \
			--foldchange 1 \
			--prefix stat


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
		p.value threshold for statistical significance

	--foldchange=NUMERIC
		foldchange threshold for statistical significance

	--prefix=CHARACTER
		The output file name

	-h, --help
		Show this help message and exit   
~~~
