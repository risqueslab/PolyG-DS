#!/bin/bash


# DS_polyG_lobSTR_troubleshooting.sh
# 
# By Dana Nachmanson (1)
# (1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
# July 2017
# 
# Can be used AFTER DS_polyG_lobSTR.sh has been run. 
# 
# 
# FastQ --> TagtoHeader.py -> lobSTR alignment -> consensus_by_alignment_lobSTR2.py
# 
# ---------------------------------------------------------------------------------
# Written for Python 2.7
# Required modules: Pysam, samtools
# 
# Inputs:
#    1: A bam file that has been aligned with lobSTR.
# 	2: A bed file with the positions of loci of interest

clear

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u

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
REGION_BED=/Users/RRisques/Desktop/Duplex_Sequencing/Reference/PolyG/PolyG_20_sites.bed # Bed file with genomic region


# 3. SET SAMPLES TO ANALYZE:
# Folder names separated by spaces containing R1 and R2 Fastq.gz files named: samplename.seq1.fastq.gz and samplename.seq2.fastq.gz
# NOTE: This script can use compressed Fastq files, no need to unzip files before running this script

folderList='s1 s2 s3 s4 s5 s6 s11 s66 7 8 9 '
 #7e 8 8e 9 9e s1 s1seperate s2 s2e s3 s3e s4 s4e s5 s5e s6 s6e s10min s11 s11e s66'

# 4. RUN SCRIPT:
# From the terminal-
# >> cd into the directory containing your SAMPLE FOLDERS and a copy of this script
# >> bash DS_polyG_lobSTR_troubleshooting.sh > DS_polyG_lobSTR_troubleshooting.txt   

for elmt in $folderList
do	
	cd ${elmt}
	# If troubleshooting directory doesn't exist then make the directory
	if [ ! -d "Troubleshooting" ]; then
  	mkdir Troubleshooting
	fi
	cp ${elmt}.smi.aligned.sorted.bam Troubleshooting/${elmt}.smi.aligned.sorted.bam 
	cp ${elmt}.smi.aligned.sorted.bam.bai Troubleshooting/${elmt}.smi.aligned.sorted.bam.bai 
	cd Troubleshooting
	
		# Prints read information by tag
		python $DS_PATH/consensus_by_alignment_lobSTR_troubleshoot1.2.0.py --input ${elmt}.smi.aligned.sorted.bam --bed $REGION_BED --prefix ${elmt} \
		--taglen $tagLen --minmem $minMem --mindiff $minDiff --refdiff $refDiff --motif $motif
    	for file in *.all_tags.txt
		do	
			python $DS_PATH/sum_tags.py "$file" False
		done
		python $DS_PATH/sum_tags.py ${elmt}.ALL.all_tags.txt True
		
		# Prints read and frequency information by tag ***NOTE TAGS ARE REPRESENTED TWICE IN THIS OUTPUT***
		python $DS_PATH/consensus_by_alignment_lobSTR_troubleshoot3.2.0.py --input ${elmt}.smi.aligned.sorted.bam --bed $REGION_BED --prefix ${elmt} \
		--taglen $tagLen --minmem $minMem --mindiff $minDiff --refdiff $refDiff --motif $motif

	cd ..
	cd ..
done