# CRISPR-DS-Poly_G
---------------------
CRISPR-DS with Poly_G data processing pipeline



see below for dependencies.

# OVERVIEW OF PROCESS
---------------------
[Remove tags from reads]
[Align reads with lobSTR]
[consensus of genotypes by alignment]

# USAGE
---------------------
python consensus_by_alignment_lobSTR.py --input aligned_bam_file

Default values are indicated with brackets.

Arguments:

-h Show this help message and exit

--input IN_BAM Path to paired-end bam file aligned with lobSTR.

--taglen TAG_LEN Length in bases of the duplex tag sequence.[10]

--spacerlen SPCR_LEN Length in bases of the spacer sequence between duplex tag and the start of target DNA. [1]

--minmem MINMEM Minimum number of reads allowed to comprise a consensus. [3]

--cutoff CUTOFF Minimum difference in frequency between the most abundant allele and  [0.7]

--write-sscs Print the SSCS reads to file in FASTQ format [False]

--prefix PREFIX Sample name to uniquely identify samples that will be appended as a prefix to the output files [None]


# LIST OF DEPENDENCIES 
---------------------
The following software and packages must be installed on your computer.

Package	Written with version
Samtools	>=1.3.1
Python	2.7.X
Pysam	>=0.9.0
MatPlotLib	>=1.5.1 (optional)
lobSTR  v4.0.6
