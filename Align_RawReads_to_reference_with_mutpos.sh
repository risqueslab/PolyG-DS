#!/bin/bash

#  Align_RawReads_to_reference_with_mutpos.sh
#  This script aligns raw CRISPR-DS_Poly-G data
#  and produces the number of reads on target as
#  well as a pileup and mutpos to give more specific 
#  information about the positions of those reads.

set -e
set -o pipefail
set -u

# 1. Set tag/adapter length
tagLen=10
spacerLen=1

# 2. SET FILE LOCATIONS AND PATHS
DS_PATH=/Users/RRisques/Desktop/Duplex_Sequencing
ALIGN_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/Human_Genome/human_g1k_v37.fasta

# 3. SET SAMPLES TO ANALYZE:
# Folder names separated by spaces containing R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and samplename.seq2.fastq.gz
# NOTE: This script can use compressed Fastq files, no need to unzip files before running this script

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and copy of this script
# >> bash Align_RawReads_to_reference_with_mutpos.sh > Align_RawReads_Poly-G.txt 

# Automated work flow follows:

folderList='7'

for elmt in $folderList
do	
	cd ${elmt}
	
	# If raw reads directory doesn't exist then make the directory
	if [ ! -d "Raw_Reads" ]; then
  	mkdir Raw_Reads
	fi
	
	# Remove tag from header
   	python $DS_PATH/Programs/tag_to_headertemp.py --infile1 ${elmt}.seq1.fastq.gz --infile2 ${elmt}.seq2.fastq.gz --taglen $tagLen --spacerlen $spacerLen --outprefix ${elmt}.raw --tagstats
    
    # Align raw reads using BWA algorithm MEM
    bwa mem $ALIGN_REF ${elmt}.raw.seq1.smi.fq.gz ${elmt}.raw.seq2.smi.fq.gz|samtools sort -o ${elmt}.raw_mem.sort.bam -
	    
    # Index raw reads 
    samtools index ${elmt}.raw_mem.sort.bam
    
	# Move all of the 'raw reads' analysis files to another folder
	for file in ${elmt}.raw*
	do	
		mv "$file" Raw_Reads/
	done
	
	cd Raw_reads
    
    # Print the number of reads total and the number of reads for each Poly-G region
    echo "Sample: ${elmt}" > ${elmt}.RawReads_OnTarget.txt
    echo "" >> ${elmt}.RawReads_OnTarget.txt
    echo "$(samtools view ${elmt}.raw_mem.sort.bam | wc -l) Reads total" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
    echo "$(samtools view ${elmt}.raw_mem.sort.bam 1:46990396-46990976 | wc -l) PolyG14 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 10:99726450-99727030 | wc -l) PolyG70 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 2:187271278-187271858 | wc -l) PolyG20 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 7:3092447-3093027 | wc -l) PolyG53 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 5:158562946-158563526 | wc -l) PolyG43 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 6:47496647-47497227 | wc -l) PolyG50 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 16:24756847-24757427 | wc -l) PolyG91 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 17:35451209-35451789 | wc -l) PolyG92 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 22:44895793-44896373 | wc -l) PolyG106 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 10:44878929-44879170 | wc -l) PolyG108 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 16:1321061-1321387 | wc -l) PolyG116 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 16:585761-586088 | wc -l) PolyG117 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 12:45801203-45801428 | wc -l) PolyG112 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 19:49154979-49155217 | wc -l) PolyG100 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 19:32450693-32450966 | wc -l) PolyG99 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 15:61181311-61181638 | wc -l) PolyG87 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 3:17292658-17292991 | wc -l) PolyG33 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 5:72511715-72511932 | wc -l) PolyG45 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 1:244081897-244082217 | wc -l) PolyG13 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	echo "$(samtools view ${elmt}.raw_mem.sort.bam 8:102121274-102121545 | wc -l) PolyG55 Reads" | sed -e 's/^[ \t]*//' >> ${elmt}.RawReads_OnTarget.txt
	
	# Make a pileup file	
	samtools mpileup -Q0 -B -A -d 500000 -f $ALIGN_REF ${elmt}.raw_mem.sort.bam > ${elmt}.raw_mem.pileup
	
	# Make a mutpos file
	cat ${elmt}.raw_mem.pileup | python $DS_PATH/Programs/mut-position.1.33.py -C 1 -d 0 > ${elmt}.raw.mutpos
	cd ..
	cd ..
done