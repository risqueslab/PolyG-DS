#!/bin/bash

set -e
set -o pipefail
set -u
set -x

fileList=$(cat $1)

for sampName in ${fileList}; do
    cd ${sampName}
    echo "Running sample ${sampName}"
    
    time bash /bioinformatics/programs/CRISPR-DS-Poly_G/DS_polyG_lobSTR.sh ${sampName}_config.sh
    
    echo "Done running sample ${sampName}"
    cd ../
done