#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
***********

consensus_by_alignment_lobSTR.2.3.py
Version 2.3

By Dana Nachmanson (1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
July 2017
Modified August 2018 by Brendan Kohrn (1)

Creating single stranded and double stranded consensus
based on lobSTR alignment of polyG/polyC regions.

FastQ --> TagtoHeader.py -> lobSTR alignment -> consensus_by_alignment_lobSTR.2.3.py

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
    4: A text file containing tagstats information for the processed reads

usage: consensus_by_alignment_lobSTR.2.3.py [-h] --input IN_BAM --bed IN_BED
                                            [--motif MOTIF] [--minmem MINMEM]
                                            [--mindiff MIN_DIFF] [--rawreads]
                                            [--taglen TAG_LEN] --prefix PREFIX

optional arguments:
  -h, --help          show this help message and exit
  --input IN_BAM      Path to aligned, paired-end, bam file.
  --bed IN_BED        Path to bed file which has the loci identifiers.
  --motif MOTIF       Motif of repeat unit. [G]
  --minmem MINMEM     Min. number of reads allowed to comprise a consensus.
                      [3]
  --mindiff MIN_DIFF  The frequency difference between first and second most
                      abundant allelemust be GREATER than this for a consensus
                      to be called for that tag. [0]
  --rawreads          Option to have a text file with raw reads printed out.
  --taglen TAG_LEN    Length of the duplex tag sequence. [10]
  --prefix PREFIX     Sample name to uniquely identify samples
---------------------------------------------------------------------------------
"""
import pysam
import sys

from argparse import ArgumentParser
from collections import Counter
from copy import deepcopy


class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class BedInputError(Error):
    """Exception raised for errors in the input.
    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    
    def __init__(self, lineNum, line):
        self.lineNum = lineNum
        self.line = line
        self.message = "ERROR: line %s in the bed file is malformed!\n" % lineNum
        sys.stderr.write(self.message)

class TagError(Error):
    """Exception raised if a tag is the wrong length"""
    
    def __init__(self, bad_tag, tag_len):
        self.bad_tag = bad_tag
        self.tag_len = tag_len
        if len(bad_tag) > tag_len:
            self.message = "ERROR: Tag %s is too long for declared tag length of %s.  \n" % (bad_tag, tag_len)
        else:
            self.message = "ERROR: Tag %s is too long for declared tag length of %s.  \n" % (bad_tag, tag_len)
        sys.stderr.write(self.message)



class PolyG:
    """This class holds information about a particular PolyG in a manner
    that makes it easy to use later.  It takes a line from a bed file
    as input.  The line number is only used for error reporting.  
    """
    
    def __init__(self, inLine, inLineNum):
        linebins = inLine.strip().split()
        if len(linebins) >= 5:
            self.name = linebins[3]
            self.ref_len = int(linebins[4])
            self.start   = int(linebins[1])
            self.end     = int(linebins[2])
            self.chrom   = linebins[0]
        else:
            raise BedInputError(inLineNum, inLine)

def get_polyG_info (inBedFile):
    bed_file = open(inBedFile, 'rb')
    polyG_info = []
    lineNum = 0
    for line in bed_file:
        lineNum += 1
        polyG_info.append(PolyG(line, lineNum))
    return(polyG_info)

def createTagDict(bam, polyG_info, motif, tagLen):
    '''
    Parses through a bam file and retains information on reads 
    with tags that have usable polyNt info in the regions 
    specified in the bed file.
    '''
    sys.stderr.write("Going through tags...\n")
    tag_dict = {}
    tag_ctr = 0
    badReads = 0
    for polyG in polyG_info:
        # For each polyG, extract those regions from 
        # the bam file
        for read in bam.fetch(polyG.chrom,  
                              polyG.start, 
                              polyG.end 
                              ):
            tag_ctr += 1
            if tag_ctr % 10000 == 0:
                sys.stderr.write("%s reads processed.....\n" % tag_ctr)
            read_id, tag_info = read.query_name.split('|')
            tag = str(tag_info.split('/')[0])
            if len(tag) != 2 * tagLen:
                raise TagError(tag, tagLen)
            # This checks whether the read spanned the PolyNT 
            # and has a genotype.
            if (read.has_tag('XD') and 
                    (read.get_tag('XR') == motif or 
                    read.get_tag('XR') == revCompl(motif))
                    ):
                if tag not in tag_dict.keys():
                    tag_dict[tag] = {}
                    tag_dict[tag]['READ_ID'] = []
                    tag_dict[tag]['XD'] = []
                    tag_dict[tag]['XG'] = []
                    tag_dict[tag]['POLY_G'] = []
                if read_id not in tag_dict[tag]['READ_ID']:
                    tag_dict[tag]['READ_ID'].append(read_id)
                    tag_dict[tag]['XD'].append(
                        polyG.ref_len + read.get_tag('XD')
                        )
                    tag_dict[tag]['XG'].append(read.get_tag('XG'))
                    tag_dict[tag]['POLY_G'].append(polyG.name)
            else:
                badReads += 1
    sys.stderr.write("\n%s reads processed\n%s bad reads\n%s good reads remaining\n\n" % 
                     (tag_ctr, badReads, tag_ctr - badReads))
    return(tag_dict)

def createSSCSDict(tag_dict, minmem, min_diff, debug):
    '''
    Input is a dictionary of tags and raw reads. 
    Returns a dictionary of consensus tags.
    '''
    sys.stderr.write("Creating SSCS Dict...\n")
    sscs_dict = {}
    too_few_reads = 0
    no_consensus = 0
    familySizeDict = {}
    for tag in tag_dict.keys():
        # Check that tag has minimum amount of reads
        num_reads = len(tag_dict[tag]['READ_ID'])
        if num_reads in familySizeDict:
            familySizeDict[num_reads] += 1
        else:
            familySizeDict[num_reads] = 1.
        if num_reads >= minmem: 
            sscs_dict[tag] = determineMajorAllele(tag_dict[tag], min_diff)
            if sscs_dict[tag] == None:
                no_consensus += 1
                if debug:
                    printTag = ["%s,%s,%s,%s" % (tag_dict[tag]['POLY_G'][x],tag,tag_dict[tag]['XD'][x], tag_dict[tag]['XG'][x]) for x in xrange(len(tag_dict[tag]['READ_ID']))]
                    for entry in printTag:
                        print(entry)
        else:
            too_few_reads += 1 
    sys.stderr.write(str(len(tag_dict.keys())) + ' total tags processed.\n')
    sys.stderr.write(str(too_few_reads) + ' tags were lost because they had too few reads\n')
    sys.stderr.write(str(no_consensus) + ' tags were lost because they could not find a major allele\n')
    sys.stderr.write(str(len(sscs_dict.keys())) + ' SSCS tags left.\n') 
    sys.stderr.write('Now only ' + str(len(sscs_dict.keys())/2) + ' DCS tags are possible.\n\n')   
    #print(str(below_ratio) + ' tags were lost because they were below the minimum frequency difference')
    return(sscs_dict, familySizeDict)
    
def createDCSDict(sscs_dict_input, tag_len):
    '''
    Input is a dictionary of tags that have been assigned and allele based on 
    alignment consensus. Searches for the complementary tag and tried to make 
    a consensus based on both tags. 
    Returns a dictionary of the consensus of complementary tags.
    '''
    sys.stderr.write("Creating DCS Dict...\n")
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
    sys.stderr.write(str(num_mismatch) + ' tags did not agree.\n')
    sys.stderr.write(str(how_many_nones) + ' tags equal None.\n')
    sys.stderr.write(str(rev_tag_not_there) + ' tags could not find their partners.\n')
    sys.stderr.write(str(len(dcs_dict.keys())) + ' tags left in the end.\n\n')
    return(dcs_dict)
    
def determineMajorAllele(one_tag_dict, min_diff):
    '''
    Input is tag information and returns an 'SSCS' read which
    takes the most abundant allele for that tag if the difference
    between it and the second most abundant are greater than indicated difference. 
    '''
    alleles = Counter(one_tag_dict['XG'])
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
    polyG = one_tag_dict['POLY_G'][one_tag_dict['XG'].index(top_allele)]
    if diff_freq > min_diff:
        return([polyG, top_allele])
    else:
        return(None)

def makeAlleleDict(polyG_info, tag_dict, dictType = "consensus"):
    '''Counts the number of reads with each allele for each polyG'''
    allele_dict = {x.name:{"TOTAL_READS": 0} for x in polyG_info}
    if dictType == "raw":
        for tag in tag_dict:
            for readIter in xrange(len(tag_dict[tag]["POLY_G"])):
                polyG = tag_dict[tag]["POLY_G"][readIter]
                allele = tag_dict[tag]["XG"][readIter]
                if polyG in allele_dict.keys():
                    if allele in allele_dict[polyG]:
                        allele_dict[polyG][allele]['count'] += 1
                    else:
                        allele_dict[polyG][allele] = {"len": len(allele) + 1, 
                                                      "count": 1.
                                                      }
                    allele_dict[polyG]['TOTAL_READS'] += 1
    else:
        for tag in tag_dict:
            if tag_dict[tag] != None:
                polyG = tag_dict[tag][0]
                allele = tag_dict[tag][1]
                if polyG in allele_dict.keys():
                    if allele in allele_dict[polyG]:
                        allele_dict[polyG][allele]['count'] += 1
                    else:
                        allele_dict[polyG][allele] = {"len": len(allele) + 1, 
                                                        "count": 1.
                                                      }
                    allele_dict[polyG]['TOTAL_READS'] += 1
    return(allele_dict)

def writeOutputFile(ref_dict, prefix, read_type):
    '''
    Takes in a polyG reference dictionary 
    and write all of the info to an output file.
    '''
    out_file = file(prefix + '.' + read_type + 'Reads.txt', 'wb')
    # Write header for file
    out_file.write('#PolyG\tAllele\tAllele_Length\t%s_Reads\tAllele_Freq\n' % 
                   (read_type)
                   )
    for polyG, alleles in sorted(ref_dict.items()):
        #print(polyG)
        #print(alleles)
        if ref_dict[polyG]['TOTAL_READS'] > 0:
            for allele in sorted(alleles, key=lambda x: len(x) + 1):
                if allele != "TOTAL_READS":
                    allele_freq = ref_dict[polyG][allele]["count"] / ref_dict[polyG]['TOTAL_READS']
                    rounded_allele_freq = "{0:.4f}".format(allele_freq)
                    out_file.write("%s\t%s\t%s\t%s\t%s\n" %
                            (polyG, allele, ref_dict[polyG][allele]["len"], 
                             ref_dict[polyG][allele]["count"], 
                             rounded_allele_freq
                             ))
    out_file.close() 

def reverse(tag, tag_len):
    ''' 
    Input is duplex tag
    ouput is the complementary tag.
    '''
    return(tag[tag_len:] + tag[:tag_len])

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
    return(''.join(rev_comp_dna))
    
def main():
    parser = ArgumentParser()
    parser.add_argument('--input', 
                        dest='in_bam', 
                        required=True,
                        help='Path to lobSTR-aligned, paired-end, bam file.'
                        )
    parser.add_argument('--bed', 
                        dest='in_bed', 
                        required=True,
                        help='Path to bed file which has the loci identifiers.'
                        )
    parser.add_argument('--motif', 
                        dest='motif', 
                        default='G', 
                        help='Motif of repeat unit. [G]'
                        )
    parser.add_argument('--minmem', 
                        dest='minmem', 
                        type=int, 
                        default=3,
                        help='Min. number of reads allowed to comprise a consensus. [3]'
                        )
    parser.add_argument('--mindiff', 
                        dest='min_diff', 
                        type=float, 
                        default=0,
                        help='The frequency difference between first and second most abundant allele' +  
                        'must be GREATER than this for a consensus to be called for that tag. [0]'
                        )
    parser.add_argument('--rawreads', 
                        dest='raw_reads', 
                        action="store_true",
                        help='Option to have a text file with raw reads printed out. '
                        )
    parser.add_argument('--taglen', 
                        dest='tag_len', 
                        type=int,  
                        default=10,
                        help='Length of the duplex tag sequence. [10]'
                        )
    parser.add_argument('--prefix', dest='prefix',
                        required=True,
                        help='Sample name to uniquely identify samples'
                        )
    parser.add_argument('--debug', 
                        dest='debug', 
                        action="store_true",
                        help="Use to print out families that don't make consensus to stdout"
                        )
    o = parser.parse_args()
    
    polyG_info = get_polyG_info(o.in_bed)
    tag_dict = createTagDict(pysam.AlignmentFile(o.in_bam, 'rb'), 
                             polyG_info, 
                             o.motif,
                             o.tag_len
                             )
    sscs_dict, tagstats = createSSCSDict(tag_dict, o.minmem, o.min_diff, o.debug)
    dcs_dict = createDCSDict(sscs_dict, o.tag_len)
    # write output files
    if o.raw_reads:
        raw_allele_dict = makeAlleleDict(polyG_info, tag_dict, dictType = "raw")
        writeOutputFile(raw_allele_dict, o.prefix, "lobSTR_Raw")
    sscs_allele_dict = makeAlleleDict(polyG_info, sscs_dict)
    writeOutputFile(sscs_allele_dict, o.prefix, "lobSTR_SSCS")
    dcs_allele_dict = makeAlleleDict(polyG_info, dcs_dict)
    writeOutputFile(dcs_allele_dict, o.prefix, "lobSTR_DCS")
    # write tagstats file
    tagstatsOut = open("%s_lobSTR_tagstats.txt" % o.prefix, 'wb')
    totalReads = sum([x * tagstats[x] for x in tagstats])
    for x in sorted(tagstats.keys()):
        tagstatsOut.write("%s\t%s\t%s\n" % (x, x * tagstats[x], x * tagstats[x]/totalReads))
    tagstatsOut.close()
    
if __name__ == "__main__":
    main()
