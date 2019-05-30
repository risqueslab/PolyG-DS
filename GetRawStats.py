#!/bin/python

import pysam
import sys
import subprocess

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
        '-N','--Name', 
        dest = 'samp_name', 
        required = True,
        help = 'Sample Name'
        )
    parser.add_argument(
        '-R','--RAW', 
        dest = 'raw_bam', 
        required = True,
        help = 'Path to bwa-aligned, paired-end, bam file.'
        )
    parser.add_argument(
        '-L','--lobSTR', 
        dest = 'lobstr_bam', 
        required = True,
        help = 'Path to lobSTR-aligned, paired-end, bam file.'
        )
    parser.add_argument(
        '-b', '--bed', 
        dest = 'in_bed', 
        required = True,
        help = 'Path to bed file which has the loci identifiers.'
        )
    o = parser.parse_args()
    polyG_info = get_polyG_info(o.in_bed)
    #inRawBam = pysam.AlignmentFile(o.raw_bam, 'rb')
    #inLobSTR_Bam = pysam.AlignmentFile(o.lobSTR_bam, 'rb')
    # Get counts of total raw reads
    
    totRawReads = subprocess.check_output(
        ["samtools", "view", "-c", o.raw_bam]
        ).decode().strip()
    for polyG in polyG_info:
        inRegion = f"{polyG.chrom}:{polyG.start}-{polyG.end}"
        # Get counts of raw reads at polyG
        polyG_raw_reads = subprocess.check_output(
            ["samtools", "view", "-c", o.raw_bam, inRegion]
            ).decode().strip()
        # Get counts of lobSTR reads at polyG
        polyG_lobSTR_reads = subprocess.check_output(
            ["samtools", "view", "-c", o.lobstr_bam, inRegion]
            ).decode().strip()
        # output line
        sys.stdout.write(
            f"{o.samp_name}\t{polyG.name}\t{totRawReads}\t"
            f"{polyG_raw_reads}\t{polyG_lobSTR_reads}\n"
            )
    #inRawBam.close()
    #inLobSTR_Bam.close()

if __name__ == "__main__":
    main()