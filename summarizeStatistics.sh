#!/bin/bash
# Summarize multiple summary statistics files in a single file
# Run as:
# bash summarizeStatistics.sh label fileList.txt

set -e
set -o pipefail
set -u
set -x

currDir=$(pwd)
label=$1
folderList=$(cat $2)
baseDir=$3
# Set up header
echo -e "#Sample\tCall Type\tPolyG\tGood Calls\tTotal Calls\t% Good Calls\tGood Alleles\tTotal Alleles\t% Good Alleles" > ${label}.lobSTR_summary_stats.txt

# 4. RUN SCRIPT:
for elmt in ${folderList}
do    
    echo ${elmt}
    cd ${baseDir}/${elmt}
    if [ -f ${elmt}.lobSTR_summary_stats.txt ]; then
        firstLine=1
        while IFS="" read -r myLine || [ -n "$myLine" ]; do
            if (( "${firstLine}" == "1" )); then
                firstLine=0
            else
                echo -e "${elmt}\t${myLine}" >> ${currDir}/${label}.lobSTR_summary_stats.txt
            fi
        done < ${elmt}.lobSTR_summary_stats.txt
    else
        echo "${elmt}.lobSTR_summary_stats.txt not found!"
    fi
    cd ${currDir}
done