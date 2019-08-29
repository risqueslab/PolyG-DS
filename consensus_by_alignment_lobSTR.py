#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
***********

consensus_by_alignment_lobSTR.2.3.py
Version 0.3.0

By Dana Nachmanson (1) 
(1) Department of Pathology, University of Washington School of 
Medicine, Seattle, WA 98195 
July 2017
Modified November 2018 by Brendan Kohrn (1)

Creating single stranded and double stranded consensus
based on lobSTR alignment of polyG/polyC regions.

FastQ 
--> TagtoHeader.py 
-> lobSTR alignment 
-> PolyG_ID_forDCS_fromLobSTR.py 
-> consensus_by_alignment_lobSTR.2.4_get_seqs.py

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
    4: A text file containing tagstats information for the processed 
       reads

usage: consensus_by_alignment_lobSTR.2.3.py [-h] --input IN_BAM 
                                            --bed IN_BED [--motif MOTIF]
                                            [--minmem MINMEM] 
                                            [--mindiff MIN_DIFF]
                                            [--rawcalls] 
                                            [--taglen TAG_LEN] 
                                            --prefix PREFIX

optional arguments:
  -h, --help          show this help message and exit
  --input IN_BAM      Path to aligned, paired-end, bam file.
  --bed IN_BED        Path to bed file which has the loci identifiers.
  --motif MOTIF       Motif of repeat unit. [G]
  --minmem MINMEM     Min. number of reads allowed to comprise a 
                      consensus.  [3]
  --mindiff MIN_DIFF  The frequency difference between first and second 
                      most abundant allelemust be GREATER than this for 
                      a consensus to be called for that tag. [0]
  --rawcalls          Option to have a text file with raw calls printed 
                      out.
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
        self.message = (f"ERROR: line {lineNum} in the bed file is "
                        f"malformed!\n"
                        )
        sys.stderr.write(self.message)

class TagError(Error):
    """Exception raised if a tag is the wrong length"""
    
    def __init__(self, bad_tag, tag_len):
        self.bad_tag = bad_tag
        self.tag_len = tag_len
        if len(bad_tag) > tag_len:
            self.message = (f"ERROR: Tag {bad_tag} is too long for "
                            f"declared tag length of {tag_len}.  \n"
                            )
        else:
            self.message = (f"ERROR: Tag bad_tag is too long for "
                            f"declared tag length of tag_len.  \n"
                            )
        sys.stderr.write(self.message)



class PolyG:
    """This class holds information about a particular PolyG in a manner
    that makes it easy to use later.  It takes a line from a bed file
    as input.  The line number is only used for error reporting.  
    """
    
    def __init__(self, inLine, inLineNum):
        linebins = inLine.strip().split()
        if len(linebins) >= 4:
            self.chrom = linebins[0]
            self.start = int(linebins[1])
            self.end = int(linebins[2])
            notebins = linebins[3].split(':')
            if len(notebins) >= 3:
                self.name = notebins[0]
                self.motif = notebins[1]
                self.ref_len = int(notebins[2])
            else:
                raise BedInputError(inLineNum, inLine)
        else:
            raise BedInputError(inLineNum, inLine)


def get_polyG_info (inBedFile):
    bed_file = open(inBedFile, 'r')
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
    badRalls = 0
    for polyG in polyG_info:
        # For each polyG, extract those regions from 
        # the bam file
        for read in bam.fetch(polyG.chrom,  
                              polyG.start, 
                              polyG.end 
                              ):
            tag_ctr += 1
            if tag_ctr % 10000 == 0:
                sys.stderr.write(f"{tag_ctr} reads processed.....\n")
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
                    tag_dict[tag]['read'] = []
                if read_id not in tag_dict[tag]['READ_ID']:
                    tag_dict[tag]['READ_ID'].append(read_id)
                    tag_dict[tag]['XD'].append(
                        polyG.ref_len + read.get_tag('XD')
                        )
                    tag_dict[tag]['XG'].append(read.get_tag('XG'))
                    tag_dict[tag]['POLY_G'].append(polyG.name)
                    tag_dict[tag]['read'].append(read.query_sequence)
            else:
                badRalls += 1
    sys.stderr.write(f"\n{tag_ctr} reads processed\n"
                     f"\t{badCalls} bad reads\n"
                     f"\t{tag_ctr - badRalls} good reads remaining\n\n"
                     )
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
                    printTag = [(f"{tag_dict[tag]['POLY_G'][x]},"
                                 f"{tag},{tag_dict[tag]['XD'][x]},"
                                 f"{tag_dict[tag]['XG'][x]}" 
                                 )
                                for x in range(
                                    len(tag_dict[tag]['READ_ID'])
                                    )
                                ]
                    for entry in printTag:
                        print(entry)
        else:
            too_few_reads += 1 
    sys.stderr.write(f"{len(tag_dict.keys())} total tags processed.\n"
                     f"\t{too_few_reads} tags were lost because they "
                     f"had too few reads\n"
                     f"\t{no_consensus} tags were lost because they "
                     f"could not find a major allele\n"
                     f"\t{len(sscs_dict.keys())} SSCS tags left.\n"
                     f"\tNow only {len(sscs_dict.keys())/2} DCS tags "
                     f"are possible.\n\n"
                     )   
    return(sscs_dict, familySizeDict)
    
def createDCSDict(sscs_dict_input, tag_len):
    '''
    Input is a dictionary of tags that have been assigned and allele 
    based on alignment consensus. Searches for the complementary tag and
    tried to make a consensus based on both tags.  Returns a dictionary 
    of the consensus of complementary tags.
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
                if (
                  (sscs_dict_input[tag][0] 
                   == sscs_dict_input[rev_tag][0]
                   ) 
                  and (sscs_dict_input[tag][1] 
                       == sscs_dict_input[rev_tag][1])
                  ):
                    if f"{rev_tag}|{tag}" not in dcs_dict.keys():
                        dcs_dict[f"{tag}|{rev_tag}"] = [
                            sscs_dict_input[tag][0], 
                            sscs_dict_input[tag][1]
                            ] 
                else:
                    num_mismatch += 1
            else:
                how_many_nones += 1
        else:
            rev_tag_not_there += 1
    sys.stderr.write(f"{num_mismatch} tags did not agree.\n"
                     f"\t{how_many_nones} tags equal None.\n"
                     f"\t{rev_tag_not_there} tags could not find their "
                     f"partners.\n"
                     f"\t{len(dcs_dict.keys())} tags left in the end."
                     f"\n\n"
                     )
    return(dcs_dict)
    
def determineMajorAllele(one_tag_dict, min_diff):
    '''
    Input is tag information and returns an 'SSCS' call which
    takes the most abundant allele for that tag if the difference
    between it and the second most abundant are greater than indicated 
    difference. 
    '''
    alleles = Counter(one_tag_dict['XG'])
    most_common_alleles = alleles.most_common() 
    top_allele = most_common_alleles[0][0]                   
    top_calls = most_common_alleles[0][1]
    if len(most_common_alleles) == 1:
        second_calls = 0
    else:
        second_calls = most_common_alleles[1][1]                   
    total_calls = top_calls + second_calls                        
    top_allele_freq = float(top_calls)/total_calls 
    second_allele_freq = float(second_calls)/total_calls  
    diff_freq = top_allele_freq - second_allele_freq
    polyG = one_tag_dict['POLY_G'][one_tag_dict['XG'].index(top_allele)]
    if diff_freq > min_diff:
        return([polyG, top_allele])
    else:
        return(None)

def makeAlleleDict(polyG_info, tag_dict, dictType="consensus"):
    '''Counts the number of reads with each allele for each polyG'''
    allele_dict = {x.name:{"TOTAL_READS": 0} for x in polyG_info}
    if dictType == "raw":
        for tag in tag_dict:
            for readIter in range(len(tag_dict[tag]["POLY_G"])):
                polyG = tag_dict[tag]["POLY_G"][readIter]
                allele = tag_dict[tag]["XG"][readIter]
                if polyG in allele_dict.keys():
                    if allele in allele_dict[polyG]:
                        allele_dict[polyG][allele]['count'] += 1
                    else:
                        allele_dict[polyG][allele] = {
                            "len": len(allele), 
                            "count": 1.,
                            "XD": tag_dict[tag]["XD"][readIter],  
                            "read": tag_dict[tag]["read"][readIter]
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
                        allele_dict[polyG][allele] = {
                            "len": len(allele), 
                            "count": 1.
                            }
                    allele_dict[polyG]['TOTAL_READS'] += 1
    return(allele_dict)

def writeOutputFile(ref_dict, prefix, call_type):
    '''
    Takes in a polyG reference dictionary 
    and write all of the info to an output file.
    '''
    out_file = open(f"{prefix}.{call_type}Calls.txt", 'w')
    # Write header for file
    if call_type == "lobSTR_Raw":
        out_file.write(f"#PolyG\tAllele\tAllele_Length\t"
                       f"{call_type}_Calls\tAllele_Freq\tXD\tSequence\n"
                       )
    else:
        out_file.write(f"#PolyG\tAllele\tAllele_Length\t"
                       f"{call_type}_Calls\tAllele_Freq\n"
                       )
    for polyG, alleles in sorted(ref_dict.items()):
        #print(polyG)
        #print(alleles)
        if ref_dict[polyG]['TOTAL_READS'] > 0:
            for allele in sorted(alleles, key=lambda x: len(x) + 1):
                if allele != "TOTAL_READS":
                    allele_freq = (ref_dict[polyG][allele]["count"] 
                                   / ref_dict[polyG]['TOTAL_READS']
                                   )
                    rounded_allele_freq = f"{allele_freq:.4f}"
                    if call_type == "lobSTR_Raw":
                        out_file.write(
                            f"{polyG}\t"
                            f"{allele}\t"
                            f"{ref_dict[polyG][allele]['len']}\t"
                            f"{ref_dict[polyG][allele]['count']}\t"
                            f"{rounded_allele_freq}\t"
                            f"{ref_dict[polyG][allele]['XD']}\t"
                            f"{ref_dict[polyG][allele]['call']}\n"
                            )
                    else:
                        out_file.write(
                            f"{polyG}\t"
                            f"{allele}\t"
                            f"{ref_dict[polyG][allele]['len']}\t"
                            f"{ref_dict[polyG][allele]['count']}\t"
                            f"{rounded_allele_freq}\n"
                            )
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
    parser.add_argument(
        '--input', 
        dest = 'in_bam', 
        required = True,
        help = 'Path to lobSTR-aligned, paired-end, sorted bam file.'
        )
    parser.add_argument(
        '--bed', 
        dest = 'in_bed', 
        required = True,
        help = 'Path to bed file which has the loci identifiers.'
        )
    parser.add_argument(
        '--motif', 
        dest = 'motif', 
        default = 'G', 
        help = 'Motif of repeat unit. [G]'
        )
    parser.add_argument(
        '--minmem', 
        dest = 'minmem', 
        type = int, 
        default = 3,
        help = (f'Min. number of reads allowed to comprise a '
                f'consensus. [3]'
                )
        )
    parser.add_argument(
        '--mindiff', 
        dest = 'min_diff', 
        type = float, 
        default = 0,
        help = (f'The frequency difference between first and second '
                f'most abundant allele must be GREATER than this value for a '
                f'consensus to be called for that family. [0]'
                )
        )
    parser.add_argument(
        '--rawcalls', 
        dest = 'raw_calls', 
        action = "store_true",
        help = 'Option to have a text file with raw calls printed out. '
        )
    parser.add_argument(
        '--taglen', 
        dest = 'tag_len', 
        type = int,  
        default = 10,
        help = 'Length of the duplex tag sequence. [10]'
        )
    parser.add_argument(
        '--prefix', 
        dest = 'prefix',
        required = True,
        help = 'Sample name to uniquely identify samples'
        )
    parser.add_argument(
        '--debug', 
        dest = 'debug', 
        action = "store_true",
        help = (f"Use to print out families that don't make consensus "
                f"to stdout"
                )
        )
    o = parser.parse_args()
    
    polyG_info = get_polyG_info(o.in_bed)
    tag_dict = createTagDict(pysam.AlignmentFile(o.in_bam, 'rb'), 
                             polyG_info, 
                             o.motif,
                             o.tag_len
                             )
    sscs_dict, tagstats = createSSCSDict(tag_dict, 
                                         o.minmem, 
                                         o.min_diff, 
                                         o.debug
                                         )
    dcs_dict = createDCSDict(sscs_dict, o.tag_len)
    # write output files
    if o.raw_calls:
        raw_allele_dict = makeAlleleDict(polyG_info, 
                                         tag_dict, 
                                         dictType = "raw"
                                         )
        writeOutputFile(raw_allele_dict, o.prefix, "lobSTR_Raw")
    sscs_allele_dict = makeAlleleDict(polyG_info, sscs_dict)
    writeOutputFile(sscs_allele_dict, o.prefix, "lobSTR_SSCS")
    dcs_allele_dict = makeAlleleDict(polyG_info, dcs_dict)
    writeOutputFile(dcs_allele_dict, o.prefix, "lobSTR_DCS")
    # write tagstats file
    tagstatsOut = open(f"{o.prefix}_lobSTR_tagstats.txt", 'w')
    totalCalls = sum([x * tagstats[x] for x in tagstats])
    for x in sorted(tagstats.keys()):
        tagstatsOut.write(f"{x}\t"
                          f"{x * tagstats[x]}\t"
                          f"{x * tagstats[x]/totalCalls}\n"
                          )
    tagstatsOut.close()
    
if __name__ == "__main__":
    main()
