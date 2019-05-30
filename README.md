
# CRISPR-DS-Poly_G

**CRISPR-DS with Poly_G data processing pipeline**
By Dana Nachmanson (1)
July 2017
Modified by Brendan Kohrn (1)
May 2019
(1) Department of Pathology, University of Washington School of Medicine,
 Seattle, WA 98195 

## Dependencies:
This pipeline can be run on a computer running either MacOS, Linux, or 
Windows 10.  In the case of Windows 10, all operation of this program 
must happen within the context of Windows Subsystem for Linux (WSL).  
In addition, the following programs and modules must be installed for 
full functionality of the pipeline.  

|Program                            |Written with version 
|-----------------------------------|------------
|Samtools                           |>=1.3.1
|Python                             |3.7.0
|Python module Pysam                |>=0.15.1
|Python module MatPlotLib (optional)|>=1.5.1
|lobSTR                             |4.0.6

## SETUP OF PIPELINE
When you first download this set of scripts, go to DS_polyG_lobSTR.sh 
and change the following variables.  Make sure there are no spaces 
around the equals sign after editing.  

|Variable     |Line #|Value
|-------------|------|------------
|DS_PATH      |37    |Path to your directory with these scripts.
|SAMTOOLS_PATH|38    |Path to your installation of samtools.
|lobSTR_PATH  |39    |Path to your installation of lobSTR.

## BED FILE PREPERATION
The bed file used in this pipeline has the following columns:

|Column    |Description
|----------|------------------
|chrom     |The chromosome containing the polyG
|start     |The start position of the polyG
|stop      |The end position of the polyG
|annotation|An anotaition, in the form "PolyG_name:Motif:ReferenceLength"

## BASIC USAGE
### CONFIGURING A RUN
In order to configure a run, set up a pre-config file with all the 
relevant information about each sample to be run.  A excel template has 
been provided to help with this operation.  In the event that the 
template does not work, or if you don't want to use it, a pre-config 
file is a CSV file with the following headers (note: capitalization 
matters): 

|Header    |default|description
|----------|-------|----------
|RUN_ID    |       |A label for this sample.  Required.  
|ref       |       |The lobSTR-formated reference to use for this sample.  Required.  
|target_bed|       |The bed file of target PolyG loci to use for this sample.  Required.  
|seq1      |seq1.fq|The input read 1 file for this sample.  May be gziped.
|seq2      |seq2.fq|The input read 2 file for this sample.  May be gziped.
|minMem    |3      |The minimum number of reads passing lobSTR required for consensus making.
|minDiff   |0      |The stringency of consensus making.  The most common allele call for a family must occur more that the second most common allele by more than this.
|motif     |G      |The expected motif; note that the reverse compliment of this motif is also expected.
|tagLen    |10     |The length of duplex tags used in library preparation.
|spacerLen |1      |The length of the spacer sequence used in library preparation.
|cleanup   |FALSE  |Sets whether to clean up intermediate files at end of run.

Make sure that there are no spaces in this file; if spaces are present, this could mess up processing later.  

After saving this file as a .CSV, run `python3 makeConfigs.py -I myPreconfigFile.csv` where 'myPreconfigFile.csv' is the file you created.  Move the resulting config files into the appropriate folders for their samples.  

Alternatively, you can take copy the config example (Example_config.sh) and edit it for a particular sample.  

### Running the pipeline
The pipeline can be run by calling `bash DS_polyG_lobSTR.sh RUN_ID_config.sh` where RUN_ID is the run id for a particular sample (as specified in the preconfig file.  You can also use a loop script that looks something like this:

    #!/bin/bash
    
    set -e
    set -o pipefail
    set -u
        
    folderList=$(cat $1)
    
    # 4. RUN SCRIPT:
    for fileIter in ${folderList}
    do    
        cd ${fileIter}
        bash /path/to/DS_PolyG/programs/DS_polyG_lobSTR.sh \
        ${fileIter}_config.sh
        cd ../
    done

which would take a list of directories you want to process as input.

### FINAL OUTPUTS:
For each sample, this pipeline generates a number of output files, described below. 

|Output                     |Description
|---------------------------|-----------
|*.tagstats.txt             |A description of family sizes found in the raw data.
|*.png                      |A graphical representation of the tagstats file, above.
|*.zoom.png                 |Family size plot zoomed in to families of size 1-40.
|*.smi.aligned.stats        |Statistics from lobSTR.
|*_lobSTR_tagstats.txt      |A description of family sizes found in post-lobSTR data.
|*.lobSTR_DCSReads_bad.txt  |DCS calls that failed filtering.
|*.lobSTR_DCSReads_good.txt |DCS calls that passed filtering.
|*.lobSTR_SSCSReads_bad.txt |SSCS calls that failed filtering.
|*.lobSTR_SSCSReads_good.txt|SSCS calls that passed filtering.
|*.lobSTR_RawReads_bad.txt  |Raw calls that failed filtering.  This file contains an example read for each allele.
|*.lobSTR_RawReads_good.txt |Raw calls that passed filtering.  
|*.logSTR_summary_stats.txt |Combined filtering statistics

## ADVANCED USE
### OVERVIEW OF PROCESSING STEPS
 1. Remove tags from reads
 2. Align reads with lobSTR 
 3. Sort aligned reads
 4. Index sorted reads
 5. Consensus of genotypes by alignment
 6. Allele filtering

### 1. Remove tags from reads
**Command:** 

    python3 ${DS_PATH}/tag_to_header.py \
        --infile1 read1.fq.gz \
        --infile2 read2.fq.gz \
        --taglen 10 \
        --spacerlen 1 \
        --outprefix myRunID \
        --tagstats

**Inputs:** 

|Input         |Required / default|Description
|--------------|------------------|------
|`--infile1`   |Required          |The read 1 file for this sample.  May be gziped.
|`--infile2`   |Required          |The read 2 file for this sample.  May be gziped.
|`--outprefix` |Required          |Prefix for output files. Will prepend onto file name .smi.fq. Will be gziped if inputs are gziped.  
|`--taglen`    |12                |The length of duplex tags used in library preparation.
|`--spacerLen` |5                 |The length of the spacer sequence used in library preparation.
|`--filtSpacer`|`None`            |Optional: Filter out sequences lacking the inputed spacer sequence. Not recommended due to significant base calling issues with the invariant spacer  
|`--tagstats`  |False             |Optional Flag: Output tagstats file and make distribution plot of tag family sizes.  Requires matplotlib to be installed.
|`--readout`   |1000000           |How many reads are processed before progress is reported.|
|`--reduce`    |False             |Optional Flag: Only output reads that will make a final DCS read.  Will only work when the --tagstats option is invoked.   Not recommended

**Outputs:**

|Output|Description
|--------------------------|--
|outPrefix.seq1.smi.fq(.gz)|Read 1 file with duplex tags moved to the read name.
|outPrefix.seq2.smi.fq(.gz)|Read 2 file with duplex tags moved to the read name.
|outPrefix.tagstats.txt    |Family size data as proportion of reads with each family size.
|outPrefix.png             |Full family size plot.
|outPrefix.zoom.png        |Family size plot zoomed in to families of size 1-40.

### 2. Align reads with lobSTR 
See lobSTR documentation ([http://lobstr.teamerlich.org/documentation.html](http://lobstr.teamerlich.org/documentation.html)) for description of lobSTR arguments and command structure.  
**Inputs:**   `outPrefix.seq1.smi.fq(.gz)` and `outPrefix.seq2.smi.fq(.gz)` from step 1.  
**Outputs:** 
 - outPrefix.smi.aligned.bam: aligned and genotyped raw reads.
 - outPrefix.smi.aligned.stats: statistics generated by lobSTR during its run.

### 3. Sort aligned reads
Read sorting by position is done using samtools sort; see samtools website ([http://www.htslib.org/](http://www.htslib.org/)) for full documentation.  

### 4. Index sorted reads
Read indexing is done using samtools index; see samtools website ([http://www.htslib.org/](http://www.htslib.org/)) for full documentation.  

### 5. Consensus of genotypes by alignment
**Command:** 

        python3 ${DS_PATH}/consensus_by_alignment_lobSTR.py \
        --input ${RUN_ID}.smi.aligned.sorted.bam \
        --bed ${BED_PATH} \
        --prefix ${RUN_ID} \
        --taglen ${tagLen} \
        --minmem ${minMem} \
        --mindiff ${minDiff} \
        --motif ${motif} \
        --rawreads

**Inputs:**

|Input       |Required / default|Description
|------------|------------------|---------
|`--input`   |Required          |Path to lobSTR-aligned, paired-end, sorted bam file.
|`--bed`     |Required          |Path to bed file which has the loci identifiers.
|`--motif`   |Required          |The expected motif; note that the reverse compliment of this motif is also expected.
|`--minmem`  |3                 |The minimum number of reads passing lobSTR required for consensus making.
|`--mindiff` |0                 |The frequency difference between first and second most abundant allele must be GREATER than this for a consensus to be called for that tag.
|`--rawreads`|False             |Option to have a text file with raw reads printed out.
|`--taglen`  |10                |The length of duplex tags used in library preparation. 
|`--prefix`  |Required          |Sample name to uniquely identify samples.
|`--debug`   |False             |Use to print out families that don't make consensus to stdout.

**Outputs**

|Output                |Description
|----------------------|-----------
|*.lobSTR_RawReads.txt |Raw calls.  This file contains an example read for each allele.
|*.lobSTR_SSCSReads.txt|SSCS calls.
|*.lobSTR_DCSReads.txt |DCS calls.
|*_lobSTR_tagstats.txt |A description of family sizes found in post-lobSTR data.


### 6. Allele filtering
Filter out alleles that diverge too much from a pure polyG.  
This step is run once for each of RawReads, SSCSReads, and DCSReads.  The example shown below is for DCSReads
**Command:** 

    python3 ${DS_PATH}/filterPolyGCalls.py \
        -i ${RUN_ID}.lobSTR_DCSReads.txt \
        -p ${RUN_ID}.lobSTR_DCSReads \
        -b \
        -m ${motif} \
        -d 2 \
        -D 2

**Inputs:**

|Input         |Required / default|Description
|--------------|------------------|-----------
|-i, --input   |Required          |A file containing tabulated Raw, SSCS, or DCS polyG calls.
|-p, --prefix  |Required          |A prefix for output files.
|-b, --badOut  |False             |Flag to output a file containing bad calls.
|-m, --motif   |Required          |The motif (A, T, C, or G) that you are expecting.
|-d, --badDiff |2                 |The number of bases that are neither the motif nor the reverse compliment of the motif.
|-D, --goodDiff|2                 |The number of reverse compliment of motif bases you are willing to tolerate.


**Outputs**

|Output                     |Description
|---------------------------|-----------
|*.lobSTR_DCSReads_bad.txt  |DCS calls that failed filtering.
|*.lobSTR_DCSReads_good.txt |DCS calls that passed filtering.
|*.lobSTR_DCSReads_stats.txt|DCS call filtering statistics.

