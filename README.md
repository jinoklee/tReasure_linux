
### tRNA expression analysis pipeline
![Pipeline](./bioexpress_pipeline.png)
## Required tools
~~~  
  cutadapt
  bowtie
  samtools
~~~
### Required R pacakges  
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

### Run script after change variables and paths
~~~
  sh tReasure_linux.sh
~~~

### Statisic analysis  for differentially expressed tRNA using EdgeR  
* before run "DEtRNA_dege.R", move readcount files to each groups.
   * move **control sample.expression.txt files** to **[resultsPATH]/control** directory
   * move **test sample.expression.txt" files** to **[resultsPATH]/test** directory

~~~
  cd [resultsPATH]/statistics
  Rscript [scriptPATH}/DEtRNA_edge.R --control [resultsPATH]/control --test  [resultsPATH]/test  --stat Exact --adj BH -pvalue 0.05  --foldchange 1 --prefix BH
~~~
