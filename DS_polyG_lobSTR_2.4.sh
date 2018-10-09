#!/bin/bash


# DS_polyG_lobSTR.sh
# 
# By Dana Nachmanson (1)
# (1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
# July 2017
# Modified by Brendan Kohrn (1)
# 
# Creating single stranded and double stranded consensus
# based on lobSTR alignment of polyG/polyC regions. Can be used 
# for other repeats and repeat units. 
# 
# FastQ --> TagtoHeader.py -> lobSTR alignment -> consensus_by_alignment_lobSTR2.py
# 
# ---------------------------------------------------------------------------------
# Written for Python 2.7
# Required modules: Pysam, samtools, lobSTR
# 
# Inputs:
#     1: Raw paired-end fastQ files with tag information in the first section of the read
#     2: A bed file with the positions of loci of interest
#     3: A fasta reference genome



# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u
set -x

# 1. SET RUN VARIABLES
minMem=3            # Minimum number of reads to reach consensus
minDiff=0			# The frequency difference between first and second most abundant allele   
                    # must be GREATER than this for a consensus to be called for that tag.
refDiff=10			# Largest unit difference between reference allele and actual allele.
motif='G'			# STR repeating unit
tagLen=10           # Adapter sequence length
spacerLen=1         # Spacer sequence length (adaptor end, begin target gene)

# 2. SET FILE LOCATIONS AND PATHS
DS_PATH=/Users/RRisques/Desktop/Duplex_Sequencing/Programs/PolyG
ALIGN_REF=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/PolyG/lobSTR_ref/hg19_v3.0.2/lobstr_v3.0.2_hg19_ref/lobSTR_
REGION_BED=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/PolyG/PolyG_20_sites.bed # Bed file with genomic region

# 3. SET SAMPLES TO ANALYZE:
# Folder names separated by spaces containing R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and samplename.seq2.fastq.gz
# NOTE: This script can use compressed Fastq files, no need to unzip files before running this script

folderList='s11'

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and a copy of this script
# >> bash -x DS_polyG_lobSTR.sh > DS_polyG_lobSTR.txt
for elmt in ${folderList}
do	
	cd ${elmt}
	
	#Remove tag from read and put into header
    python /Users/RRisques/Desktop/Duplex_Sequencing/Programs/tag_to_header.py --infile1 ${elmt}.seq1.fastq.gz --infile2 ${elmt}.seq2.fastq.gz \
    --taglen $tagLen --spacerlen $spacerLen --outprefix ${elmt} --tagstats 
  	
	#Align fastQ file with lobSTR
	#Extra parameters can be added in this command and usage found here: http://lobstr.teamerlich.org/usage.html 
 	lobSTR --p1 ${elmt}.seq1.smi.fq.gz --p2 ${elmt}.seq2.smi.fq.gz -q --rg-sample ${elmt} --index-prefix $ALIGN_REF -o ${elmt}.smi --rg-lib spike_in --gzip --multi --mismatch 3
 	
	#Sort bam file
 	samtools sort ${elmt}.smi.aligned.bam -o ${elmt}.smi.aligned.sorted.bam

	#Index bam file
 	samtools index ${elmt}.smi.aligned.sorted.bam 
	
	# Perform consensus of lobSTR allele calls
	python ${DS_PATH}/consensus_by_alignment_lobSTR.2.3.py --input ${elmt}.smi.aligned.sorted.bam --bed ${REGION_BED} --prefix ${elmt} \
	--taglen ${tagLen} --minmem ${minMem} --mindiff ${minDiff} --motif ${motif} --rawreads
	
	python ~/Desktop/Duplex_Sequencing/Programs/PolyG/filterPolyGCalls_v0.1.0.py \
	-i ${elmt}.lobSTR_DCSReads.txt -p ${elmt}.lobSTR_DCSReads -b -m ${motif} -d 1
	cd ..
done
