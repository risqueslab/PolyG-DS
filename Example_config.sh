#!/bin/bash

# This is the configuration file template for the CRISPR-DS-Poly_G pipeline
# Version 3
#
# Copy this into the directory where you want your output files, 
# Then fill in the variables as appropriate.

# 1: Set up paths to input files:

SEQ1=seq1.fq
SEQ2=seq2.fq
REF_PATH=/path/to/lobSTR_ref/hg19_v3.0.2/lobstr_v3.0.2_hg19_ref/lobSTR_
BED_PATH=/path/to/myBed.bed
# 2: Set up run variables
RUN_ID=MyRunID
minMem=3
minDiff=0
motif=G
tagLen=10
spacerLen=1
rawOut=TRUE
cleanup=FALSE
startStep=0         # Specifies which step of the process to start on; 
                    #  0: Run all steps, or start at last step run
                    #  1: Start at tag_to_header. 
                    #  2: Start at alignment with lobSTR
                    #  3: Start at sorting bam
                    #  4: Start at indexing bam
                    #  5: Start at consensus making
                    #  6: Start at filtering
                    #  7: Start at statistics compilation
