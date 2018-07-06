#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
***********

consensus_by_alignment_lobSTR.2.0.py
Version 2.0

By Dana Nachmanson (1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
July 2017

Creating single stranded and double stranded consensus
based on lobSTR alignment of polyG/polyC regions.

FastQ --> TagtoHeader.py -> lobSTR alignment -> consensus_by_alignment_lobSTR2.py

---------------------------------------------------------------------------------
Written for Python 2.7
Required modules: Pysam

Inputs:
    1: A position-sorted paired-end BAM file containing reads aligned
    with lobSTR 4.0.6
    2: A bed file with the positions of loci of interest
Outputs: 
    1: A text file containing DCS allele calls for each region
    2: A text file containing SSCS allele calls for each region
    3: A text file containing raw allele calls for each region

usage: consensus_by_alignment_lobSTR.2.0.py [-h] [--input in_bam] [--bed in_bed]
                      [--minmem minmem]  --mindiff min_diff --refdiff ref_diff
                      [--rawreads raw_reads] [--taglen tag_len]
                      [--prefix prefix] [--motif motif]

optional arguments:
  -h, --help            show this help message and exit

  --taglen tag_len      Length of the duplex tag sequence.[10]
  --rawreads raw_reads  Option to have a text file with raw reads printed out. [True]
  --minmem minmem       Min. number of reads allowed to comprise a consensus. [3]
  --mindiff min_diff    The frequency difference between first and second most abundant allele 
                        must be GREATER than this for a consensus to be called for that tag. [0]
  --refdiff ref_diff    Largest unit difference between reference allele and actual allele. [10]
  --rawreads raw_reads  Option to have a text file with raw reads printed out. [True]
  --motif motif         Motif of repeat unit. [G]

---------------------------------------------------------------------------------
"""
import pysam
from argparse import ArgumentParser
from collections import Counter
from copy import deepcopy

def createTagDict(bam, bed, raw_reads, ref_dict_raw, ref_lens, motif):
    '''
    Parses through a bam file and retains information on reads 
    with tags that have usable polyG info in the regions 
    specified in the bed file.
    '''
    print("Going through tags...")
    print('\n')
    tag_dict = {}
    with open(bed) as region_file:
        # For each loci in the bed file, extract those regions from the bam file
        for i, line in enumerate(region_file):
            chrom, start, end, polyG, ref_length = line.split()
            for read in bam.fetch(chrom, int(start), int(end)):
                read_id, tag_info = read.query_name.split('|')
                tag = str(tag_info.split('/')[0])
                # This checks whether the read spanned the PolyG/PolyC 
                # and has a genotype.
                if read.has_tag('XD') and (read.get_tag('XR') == motif or read.get_tag('XR') == revCompl(motif)):
                    delta_ref = str(read.get_tag('XD'))
                    allele_len = int(ref_lens[polyG]) + int(delta_ref)
                    if tag not in tag_dict.keys():
                        tag_dict[tag] = {}
                        tag_dict[tag]['READ_ID'] = []
                        tag_dict[tag]['XD'] = []
                        tag_dict[tag]['POLY_G'] = []
                    # Make sure reads that are unpaired by lobSTR are not genotyped
                    # and counted separately.
                    if read_id not in tag_dict[tag]['READ_ID']:
                        tag_dict[tag]['READ_ID'].append(read_id)
                        tag_dict[tag]['XD'].append(allele_len)
                        tag_dict[tag]['POLY_G'].append(polyG)
                    # If you want raw read information, record info here.
                    if raw_reads:
                        if polyG in ref_dict_raw.keys() and allele_len in ref_dict_raw[polyG].keys():
                            ref_dict_raw[polyG][allele_len] += 1
                            ref_dict_raw[polyG]['TOTAL_READS'] += 1
    region_file.close()
    return tag_dict, ref_dict_raw
                                                                
def createSSCSDict(tag_dict, minmem, min_diff):
    '''
    Input is a dictionary of tags and raw reads. 
    Returns a dictionary of consensus tags.
    '''
    print("Creating SSCS Dict...")
    print('\n')
    sscs_dict = {}
    too_few_reads = 0
    no_consensus = 0
    for tag in tag_dict.keys():
        # Check that tag has minimum amount of reads
        num_reads = len(tag_dict[tag]['READ_ID'])
        if num_reads >= minmem: 
            sscs_dict[tag] = determineMajorAllele(tag_dict[tag], min_diff)
            if sscs_dict[tag] == None:
                no_consensus += 1
        else:
            too_few_reads += 1 
    print(str(len(tag_dict.keys())) + ' total tags processed.')
    print(str(too_few_reads) + ' tags were lost because they had too few reads')
    print(str(no_consensus) + ' tags were lost because they could not find a major allele')
    print(str(len(sscs_dict.keys())) + ' SSCS tags left.') 
    print('Now only ' + str(len(sscs_dict.keys())/2) + ' DCS tags are possible.')   
    #print(str(below_ratio) + ' tags were lost because they were below the minimum frequency difference')
    print('\n')
    return sscs_dict
    
def createDCSDict(sscs_dict_input, tag_len):
    '''
    Input is a dictionary of tags that have been assigned and allele based on 
    alignment consensus. Searches for the complementary tag and tried to make 
    a consensus based on both tags. 
    Returns a dictionary of the consensus of complementary tags.
    '''
    print("Creating DCS Dict...")
    print('\n')
    dcs_dict = {}
    num_mismatch = 0
    rev_tag_not_there = 0
    how_many_nones = 0
    for tag in sscs_dict_input.keys():
        rev_tag = reverse(tag, tag_len)
        if str(rev_tag) in sscs_dict_input.keys():
            # Makes sure that both tags have a major allele
            if sscs_dict_input[tag] and sscs_dict_input[rev_tag]:
                # Check if both tags have same major allele and mapped to sample PolyG
                if (sscs_dict_input[tag][0] == sscs_dict_input[rev_tag][0]) and (sscs_dict_input[tag][1] == sscs_dict_input[rev_tag][1]):
                    if (str(rev_tag) + '|' + str(tag)) not in dcs_dict.keys():
                        dcs_dict[str(tag) + '|' + str(rev_tag)] = [sscs_dict_input[tag][0], sscs_dict_input[tag][1]] 
                else:
                    num_mismatch += 1
            else:
                how_many_nones += 1
        else:
            rev_tag_not_there += 1
    print(str(num_mismatch) + ' tags did not agree.')
    print(str(how_many_nones) + ' tags equal None.')
    print(str(rev_tag_not_there) + ' tags could not find their partners.')
    print(str(len(dcs_dict.keys())) + ' tags left in the end.')
    print('\n')
    return dcs_dict

def determineMajorAllele(one_tag_dict, min_diff):
    '''
    Input is tag information and returns an 'SSCS' read which
    takes the most abundant allele for that tag if the difference
    between it and the second most abundant are greater than indicated difference. 
    '''
    alleles = Counter(one_tag_dict['XD'])
    most_common_alleles = alleles.most_common() 
    top_allele = most_common_alleles[0][0]                   
    top_reads = most_common_alleles[0][1]
    if len(most_common_alleles) == 1:
        second_reads = 0
    else:
        second_reads = most_common_alleles[1][1]                   
    total_reads = top_reads + second_reads                        
    top_allele_freq = float(top_reads)/total_reads 
    second_allele_freq = float(second_reads)/total_reads  
    diff_freq = top_allele_freq - second_allele_freq
    polyG = one_tag_dict['POLY_G'][one_tag_dict['XD'].index(top_allele)]
    if diff_freq > min_diff:
        return [polyG, top_allele]
    else:
        return None 

def createRefDicts(in_bed_file, ref_diff):
    '''
    Using a minimum allele length, maximum allele length and a csv file 
    of loci, returns a dictionary with that reference information.
    '''
    bed_file = open(in_bed_file, 'r')
    allele_dict = {}
    ref_len_dict = {}
    for line in bed_file:
        line_lst = line.split()
        if len(line_lst) > 1:
            poly_g = str(line_lst[3].strip())
            ref_len = int(line_lst[4].strip())
            ref_len_dict[poly_g] = ref_len
            # Adds a key which is a polyG
            allele_dict[poly_g] = {}
            allele_dict[poly_g]['REF_LEN'] = ref_len
            allele_dict[poly_g]['TOTAL_READS'] = 0
            for rep_length in range(ref_len - int(ref_diff), ref_len + int(ref_diff + 1)):
                allele_dict[poly_g][int(rep_length)] = 0
    bed_file.close()
    return allele_dict, ref_len_dict
    
def writeOutputFile(ref_dict, prefix, read_type):
    '''
    Takes in a polyG reference dictionary 
    and write all of the info to an output file.
    '''
    out_file = file(prefix + '.' + read_type + 'Reads.txt', 'w+')
    # Write header for file
    out_file.write('PolyG' + '\t' + 'Allele Length' + '\t' + str(read_type) + 'Reads' + '\t' + 'Allele Freq' + '\n')
    for polyG, alleles in sorted(ref_dict.items()):
        for allele in alleles:
            if str(allele).isdigit():
                if ref_dict[polyG]['TOTAL_READS'] > 0:
                    allele_freq = float(ref_dict[polyG][allele])/int(ref_dict[polyG]['TOTAL_READS'])
                    rounded_allele_freq = "{0:.4f}".format(allele_freq)
                else:
                    rounded_allele_freq = '0'
                out_file.write(str(polyG) + '\t' + str(allele) + '\t' + str(ref_dict[polyG][allele]) + '\t' + str(rounded_allele_freq) + '\n')
    out_file.close()    

def fillInRefDict(ref_dict, tag_dict):
    '''
    Uses a template dictionary with all the regions and possible
    alleles and fills in with actual read information.
    '''
    for tag in tag_dict.keys():
        if tag_dict[tag]:
            polyG = tag_dict[tag][0]
            allele = int(tag_dict[tag][1])
            if polyG in ref_dict.keys() and allele in ref_dict[polyG].keys():
                ref_dict[polyG][allele] += 1
                ref_dict[polyG]['TOTAL_READS'] += 1
            else:
                print("This allele is out of range as specified by 'refdiff':")
                print(str(polyG) + '\t' + str(allele))
    return ref_dict
    
def reverse(tag, tag_len):
    ''' 
    Input is duplex tag
    ouput is the complementary tag.
    '''
    return tag[tag_len:] + tag[:tag_len]

def revCompl(dna):
    '''
    Takes a dna sequence and returns
    the reverse complement of sequence.
    '''
    comp_bases = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    rev_comp_dna = []
    rev_dna = (list(dna))[::-1]
    for base in rev_dna:
        rev_comp_dna.append(comp_bases[base])
    return ''.join(rev_comp_dna)        
    
def main():
    parser = ArgumentParser()
    parser.add_argument('--input', dest='in_bam', required=True,
                        help='Path to aligned, paired-end, bam file.')
    parser.add_argument('--bed', dest='in_bed', required=True,
                        help='Path to bed file which has the loci identifiers.')
    parser.add_argument('--motif', dest='motif', type=str, default='G', required=True,
                        help='Motif of repeat unit. [G]')
    parser.add_argument('--minmem', dest='minmem', type=int,  default=3,
                        help='Min. number of reads allowed to comprise a consensus. [3]')
    parser.add_argument('--mindiff', dest='min_diff', type=float, default=0,
                        help='The frequency difference between first and second most abundant allele' +  
                        'must be GREATER than this for a consensus to be called for that tag. [0]')
    parser.add_argument('--refdiff', dest='ref_diff', type=float, default=10,
                        help='Largest unit difference between reference allele and actual allele. [10]')
    parser.add_argument('--rawreads', dest='raw_reads', type=bool, default=True,
                        help='Option to have a text file with raw reads printed out. [True]')
    parser.add_argument('--taglen', dest='tag_len', type=int,  default=10,
                        help='Length of the duplex tag sequence. [10]')
    parser.add_argument('--prefix', dest='prefix', type=str, required=True,
                        help='Sample name to uniquely identify samples')
    o = parser.parse_args()

    # Read in reference alleles
    polyG_ref_dict, ref_lens = createRefDicts(o.in_bed, o.ref_diff)
    
    # Create a dictionary of tags in our areas of interest and where they align to
    bam = pysam.AlignmentFile(o.in_bam, 'rb')
    tag_dict, raw_ref_dict = createTagDict(bam, o.in_bed, o.raw_reads, deepcopy(polyG_ref_dict), ref_lens, o.motif)
    
    sscs_dict = createSSCSDict(tag_dict, o.minmem, o.min_diff)
    dcs_dict = createDCSDict(sscs_dict, o.tag_len)
    dcs_ref_dict = fillInRefDict(deepcopy(polyG_ref_dict), dcs_dict)
    sscs_ref_dict = fillInRefDict(deepcopy(polyG_ref_dict), sscs_dict)
    # Create output files
    writeOutputFile(dcs_ref_dict, o.prefix, 'lobSTR_DCS')
    writeOutputFile(sscs_ref_dict, o.prefix, 'lobSTR_SSCS')
    # Create output files
    if o.raw_reads:
        writeOutputFile(raw_ref_dict, o.prefix, 'lobSTR_Raw')
    
if __name__ == "__main__":
    main()