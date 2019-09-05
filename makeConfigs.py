#!/local/bin/python
"""
BERKELEY SOFTWARE DISTRIBUTION LICENSE
Copyright (c) 2019, Brendan Kohrn, University of Washington
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation and/or 
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors 
may be used to endorse or promote products derived from this software without 
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
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
