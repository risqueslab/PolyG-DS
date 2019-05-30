#!/bin/bash


# DS_polyG_lobSTR.sh
# Version 0.3.0
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
# FastQ --> TagtoHeader.py -> lobSTR alignment -> consensus_by_alignment_lobSTR.py
# 
# ---------------------------------------------------------------------------------
# Written for Python 3
# Required modules: Pysam, samtools, lobSTR
# 
# Inputs:
#     1: Raw paired-end fastQ files with tag information in the first section of the read
#     2: A bed file with the positions of loci of interest
#     3: A fasta reference genome indexed for lobSTR 
#       (can be downloaded from http://lobstr.teamerlich.org/download.html)

# Stop on any error inside or outside pipeline or on an unassigned variable.
set -e
set -o pipefail
set -u
set -x

source $1

# 1. SET FILE LOCATIONS AND PATHS
DS_PATH=/Users/RRisques/Desktop/BK/risques/CRISPR-DS-Poly_G
SAMTOOLS_PATH=samtools
lobSTR_PATH=lobSTR

source $1

if [ ${startStep} -eq 0 ]; then
    if [ -e ".prog_step" ]; then 
        startStep=$(cat .prog_step)
    else
        echo 1 > .prog_step
        startStep=1
    fi
fi

progStep=1

if [ "${startStep}" -le "${progStep}" ]; then
# Remove tag from read and put into header
python3 ${DS_PATH}/tag_to_header.py \
    --infile1 ${SEQ1} \
    --infile2 ${SEQ2} \
    --taglen ${tagLen} \
    --spacerlen ${spacerLen} \
    --outprefix ${RUN_ID} \
    --tagstats 
# Outputs:
#     ${elmt}.seq1.smi.fq.gz
#     ${elmt}.seq2.smi.fq.gz
#     ${elmt}.tagstats.txt
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${startStep}" -le "${progStep}" ]; then
# Align fastQ file with lobSTR
# Extra parameters can be added in this command and usage found 
#   here: http://lobstr.teamerlich.org/usage.html 
${lobSTR_PATH} \
    --p1 ${RUN_ID}.seq1.smi.fq.gz \
    --p2 ${RUN_ID}.seq2.smi.fq.gz \
    -q \
    --rg-sample ${RUN_ID} \
    --index-prefix ${REF_PATH} \
    -o ${RUN_ID}.smi \
    --rg-lib spike_in \
    --gzip \
    --multi \
    --mismatch 3
# Output: ${elmt}.smi.aligned.bam
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${startStep}" -le "${progStep}" ]; then
# Sort bam file
${SAMTOOLS_PATH} sort ${RUN_ID}.smi.aligned.bam \
    -o ${RUN_ID}.smi.aligned.sorted.bam
# Output: ${elmt}.smi.aligned.sorted.bam
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${startStep}" -le "${progStep}" ]; then
# Index bam file
${SAMTOOLS_PATH} index ${RUN_ID}.smi.aligned.sorted.bam 
# Output: ${elmt}.smi.aligned.sorted.bai
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${startStep}" -le "${progStep}" ]; then
# Perform consensus of buffered lobSTR allele calls
python3 ${DS_PATH}/consensus_by_alignment_lobSTR.py \
    --input ${RUN_ID}.smi.aligned.sorted.bam \
    --bed ${BED_PATH} \
    --prefix ${RUN_ID} \
    --taglen ${tagLen} \
    --minmem ${minMem} \
    --mindiff ${minDiff} \
    --motif ${motif} \
    --rawreads
# Outputs:
#     ${elmt}.lobSTR_RawReads.txt
#     ${elmt}.lobSTR_SSCSReads.txt
#     ${elmt}.lobSTR_DCSReads.txt
#     ${elmt}_lobSTR.tagstats.txt
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${startStep}" -le "${progStep}" ]; then
# Filter polyG calls
# Raw
python3 ${DS_PATH}/filterPolyGCalls.py \
    -i ${RUN_ID}.lobSTR_RawReads.txt \
    -p ${RUN_ID}.lobSTR_RawReads \
    -b \
    -m ${motif} \
    -d 2 \
    -D 2
# SSCS
python3 ${DS_PATH}/filterPolyGCalls.py \
    -i ${RUN_ID}.lobSTR_SSCSReads.txt \
    -p ${RUN_ID}.lobSTR_SSCSReads \
    -b \
    -m ${motif} \
    -d 2 \
    -D 2
# DCS
python3 ${DS_PATH}/filterPolyGCalls.py \
    -i ${RUN_ID}.lobSTR_DCSReads.txt \
    -p ${RUN_ID}.lobSTR_DCSReads \
    -b \
    -m ${motif} \
    -d 2 \
    -D 2
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${startStep}" -le "${progStep}" ]; then
# Combine statistics from each filter
echo -e "#Call Type\tPolyG\tGood Calls\tTotal Calls\t% Good Calls\tGood Alleles\tTotal Alleles\t% Good Alleles" > ${RUN_ID}.lobSTR_summary_stats.txt
firstLine=1
while IFS="" read -r myLine || [ -n "$myLine" ]; do
    if (( "${firstLine}" == "1" )); then
        firstLine=0
    else
        echo -e "Raw\t${myLine}" >> ${RUN_ID}.lobSTR_summary_stats.txt
    fi
done < ${RUN_ID}.lobSTR_RawReads_stats.txt

firstLine=1
while IFS="" read -r myLine || [ -n "$myLine" ]; do
    if (( "${firstLine}" == "1" )); then
        firstLine=0
    else
        echo -e "SSCS\t${myLine}" >> ${RUN_ID}.lobSTR_summary_stats.txt
    fi
done < ${RUN_ID}.lobSTR_SSCSReads_stats.txt

firstLine=1
while IFS="" read -r myLine || [ -n "$myLine" ]; do
    if (( "${firstLine}" == "1" )); then
        firstLine=0
    else
        echo -e "DCS\t${myLine}" >> ${RUN_ID}.lobSTR_summary_stats.txt
    fi
done < ${RUN_ID}.lobSTR_DCSReads_stats.txt
cd ..
echo $((progStep+1)) > .prog_step
fi

progStep=$((progStep+1))

if [ "${cleanup}" = "TRUE" ]; then
rm ${RUN_ID}*.smi.fq.gz
rm ${RUN_ID}.smi.aligned.ba*
rm ${RUN_ID}.smi.aligned.sorted.ba*
fi

