
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
