#!/bin/bash

#####################################################################################################
#change here variable and paths
#####################################################################################################
genomeName=Hsapi38 
# https://treasure.pmrc.re.kr/data/genome/Hsapi38/Hsapi38.tar.gz
# https://treasure.pmrc.re.kr/data/genome/Mmusc10/Mmusc10.tar.gz

# reference path
genomeDir=/refer/${genomeName}

# tool path
tool="/home/tools"
bowtie="${tool}/bowtie-1.3.0-linux-x86_64/bowtie"
cutadapt="${tool}/cutadapt"
samtools="${tool}/samtools"		
script="/tRNA/script"  # R 

# workdir
ngsDir=/fastq # raw fastq
outDir=/results # results (from trim to rc.txt)

###  directory of outputs ############################################################################
# trim  - outputs of trimming ( _trimmed.fastq)
# fmapping - outputs of first mapping (.bam, .bai)
# remove - outputs of remove the premature reads (.mature.fastq)
# smapping - outputs of second mapping (.bam, .bai)
# rc - outputs of read counting (.expression.txt)
######################################################################################################

# before the start, fastqc running ==> check quality and adapter 
# 3' small adapter : TGGAATTCTCGGGTGCCAAGG
# universal adapter : AGATCGGAAGAG
# fastq have 3' small adapter 

# 1. Trimming using cutadapt

	mkdir -p ${outDir}/trim

	for i in $(ls ${ngsDir}/*.fastq)
	do
		bi=$(basename $i .fastq)
		echo ${bi}
		${cutadapt} -m 10 -M 50 -q 25 -a TGGAATTCTCGGGTGCCAAGG -o ${outDir}/trim/${bi}_trimmed.fastq ${ngsDir}/${bi}.fastq 

	done


# 2. 1st round mapping using bowtie

	mkdir -p ${outDir}/fmapping/log
	log=${outDir}/fmapping/log

	cd ${outDir}/fmapping
	
	for i in $(ls ${outDir}/trim/*_trimmed.fastq)
	do
		bi=$(basename ${i} _trimmed.fastq)
		echo ${bi}
		(${bowtie}  -v 3 --best -p 15 ${genomeDir}/bowtie/${genomeName}  ${outDir}/trim/${bi}_trimmed.fastq -S ${outDir}/fmapping/${bi}.sam) 2> ${log}/${bi}.fmapping.log
		(${samtools} view -bS ${outDir}/fmapping/${bi}.sam | samtools sort -@ 2 > ${bi}.bam) 2>> ${log}/${bi}.samtools.log
		${samtools} index ${bi}.bam
		rm ${bi}.sam
	done

# 3. remove the premature reads 
	
	mkdir ${outDir}/remove

 	Rscript ${script}/removeNontRNA.R --bed ${genomeDir}/${genoneName}.tRNAscan_pre-tRNAs.bed12 --bamdir ${outDir}/fmapping --fastqdir ${outDir}/trim  --outdir ${outDir}/remove


# 4. 2nd round mapping using bowtie

	mkdir -p ${outDir}/smapping/log
	log=${outDir}/smapping/log

	cd ${outDir}/smapping

	for i in $(ls ${outDir}/remove/*.mature.fastq)
	do
		bi=$(basename ${i} .mature.fastq)
		echo ${bi}
		(${bowtie} -v 3 --best -p 15 ${genomeDir}/bowtie/${genomeName}.tRNAscan_mature ${outDir}/remove/${bi}.mature.fastq -S ${outDir}/smapping/${bi}.sam) 2> ${log}/${bi}.smapping.log
		(${samtools} view -bS ${outDir}/smapping/${bi}.sam | samtools sort -@ 2 > ${bi}.bam) 2>> ${log}/${bi}.samtools.log
	       ${samtools} index ${bi}.bam
	       rm ${bi}.sam
       done

# 5. count

	mkdir -p ${outDir}/rc
	
	cd ${outDir}/smapping

	for i in $(ls ${outDir}/smapping/*.bam)
	do
		bi=$(basename ${i} .bam)
		echo ${bi}
		${samtools} idxstats ${bi}.bam > ${outDir}/rc/${bi}.expression.txt
	done

# 6. statistics using edgeR 
 mkdir -p ${outDir}/control
 mkdir -p ${outDir}/test
 mkdir -p ${outDir}/statistics


######################################################################################################
# before run "DEtRNA_dege.R", move readcount files to each groups. 
# mv  ${outDir}/rc/*.expression.txt  ${outDir}/control
# mv  ${outDir}/rc/*.expression.txt  ${outDir}/test
######################################################################################################
 
# cd ${outDir}/statistics
# Rscript ${script}/DEtRNA_edge.R --control ${outDir}/control --test  ${outDir}/test  --stat Exact --adj BH -pvalue 0.05  --foldchange 1 --prefix BH
