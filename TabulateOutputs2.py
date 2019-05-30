import gzip
from argparse import ArgumentParser
import sys

parser = ArgumentParser()
parser.add_argument('-L', '--fileList', action="store", dest="inFileList", help = "Tab-delimiated file with sample names and sample file paths.  ")
parser.add_argument('-o', '--outputPrefix', action="store", dest="prefix", help = "A prefix to use for output files.  ")
parser.add_argument('-t', '--type', action="store", dest="dataType", help="What type of data is this? One of Raw, SSCS, or DCS.  ")
o = parser.parse_args()

# Open fileList
fileList = open(o.inFileList, 'r')

# Read each VCF file
OutLines = []
OutLines.append(f"#Sample\tPolyG\tAllele\tAllele_Length\t"
                f"{o.dataType}_Reads\tAllele_Freq\n"
                )

for inFileLine in fileList:
    linebins = inFileLine.strip().split('\t')
    sampName = linebins[0]
    sampFile = linebins[1]
    print(sampName)
    inFile = open(sampFile, 'r')

    for line in inFile:
        if line[0] == "#":
            pass
        else:
            OutLines.append(f"{sampName}\t{line}")
    inFile.close()

outFile = open(f"{o.prefix}.txt", 'w')

for outIter in range(len(OutLines)):
    outFile.write(OutLines[outIter])
outFile.close()
