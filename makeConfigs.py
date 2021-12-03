#!/usr/bin/env python3
import logging
from argparse import ArgumentParser
import sys

class Error(Exception):
    """Base class for exceptions in this module.
    """
    pass

class InputError(Error):
    """Exception raised for errors in the input
    """
    
    def __init__(self):
        self.message = (f"ERROR: Input file must contain all of: \n"
                        f"\tRUN_ID\n"
                        f"\tref\n"
                        f"\ttarget_bed\n"
                        )
        sys.stderr.write(self.message)

parser = ArgumentParser()
parser.add_argument(
    '-I',
    '--in',
    action="store", 
    dest="inFile", 
    required=True, 
    help="Input file containing information about the different samples."
    )
o = parser.parse_args()


optionDefaults = {"minMem": 3,  
                  "minDiff": 0, 
                  "motif": 'G', 
                  "tagLen": 10, 
                  "spacerLen": 1,
                  "locLen": 0, 
                  "rawOut": "TRUE", 
                  "cleanup": "FALSE"
                  }

inFile = open(o.inFile, 'r')
firstLine = True
columnLabels = []
for line in inFile:
    if firstLine:
        firstLine = False
        columnLabels = line.strip().split(',')
        if (
                "RUN_ID" not in columnLabels 
                or "ref" not in columnLabels 
                or "target_bed" not in columnLabels
                ):
            raise InputError()
    else:
        linebins = line.strip().split(',')
        print(f"Writing {linebins[columnLabels.index('RUN_ID')]}_config.sh")
        outfile=open(f"{linebins[columnLabels.index('RUN_ID')]}_config.sh", 'w')
        outfile.write(
            f"#!/bin/bash\n"
            f"\n"
            f"# This is the configuration file template for the "
            f"CRISPR-DS-Poly_G pipeline\n"
            f"# Version 3\n"
            f"#\n"
            f"# Copy this into the directory where you want your output "
            f"files, \n"
            f"# Then fill in the variables as appropriate.\n"
            f"\n"
            f"# 1: Set up paths to input files:\n"
            f"\n"
            )
        
        if "seq1" in columnLabels:
            outfile.write(
                f"SEQ1={linebins[columnLabels.index('seq1')]}\n"
                )
        else:
            outfile.write(f"SEQ1=seq1.fq\n")
        if "seq2" in columnLabels:
            outfile.write(
                f"SEQ2={linebins[columnLabels.index('seq2')]}\n"
                )
        else:
            outfile.write(f"SEQ2=seq2.fq\n")
        
        outfile.write(f"REF_PATH={linebins[columnLabels.index('ref')]}\n")
        outfile.write(
            f"BED_PATH={linebins[columnLabels.index('target_bed')]}\n"
            )
        
        outfile.write("# 2: Set up run variables\n")
        
        outfile.write(
            f"RUN_ID={linebins[columnLabels.index('RUN_ID')]}\n"
            )
        for columnIter in optionDefaults:
            if columnIter in columnLabels:
                outfile.write(
                    f"{columnIter}="
                    f"{linebins[columnLabels.index(columnIter)]}\n"
                    )
            else:
                outfile.write(
                    f"{columnIter}={optionDefaults[columnIter]}\n"
                    )
        outfile.write(
            f"startStep=0         # Specifies which step of the process to "
            f"start on; \n"
            f"                    #  0: Run all steps, or start at last step run\n"
            f"                    #  1: Start at tag_to_header. \n"
            f"                    #  2: Start at alignment with lobSTR\n"
            f"                    #  3: Start at sorting bam\n"
            f"                    #  4: Start at indexing bam\n"
            f"                    #  5: Start at consensus making\n"
            f"                    #  6: Start at filtering\n"
            f"                    #  7: Start at statistics compilation\n"
            )
        outfile.close()
print("Done")