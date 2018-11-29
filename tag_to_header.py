#!/usr/bin/env python3
#
# Tag To Header
# Version 2.0.1
# By Joe Hiatt, Scott Kennedy(1), Brendan Kohrn (1) and Mike Schmitt(1)
# (1) Department of Pathology, University of Washington School of 
#     Medicine, Seattle, WA 98195
# March 24, 2014
#
# Isolate duplex tags, move them from within the sequenced read to the 
# header region, and remove the spacer region.
#
# usage: tag_to_header.py [-h] [--infile1 INFILE1] [--infile2 INFILE2]
#                        [--outfile1 OUTFILE1] [--outfile2 OUTFILE2]
#                        [--taglen BLENGTH] [--spacerlen SLENGTH]
#                        [--read_out ROUT] [--filt_spacer ADAPTERSEQ] 
#                        --tagstats
#
# Optional arguments:
#  -h, --help                   show this help message and exit
#  --infile1 INFILE1            First input raw fastq file.
#  --infile2 INFILE2            Second input raw fastq file.
#  --outfile1 OUTFILE1          Output file for first fastq reads.
#  --outfile2 OUTFILE2          Output file for second fastq reads.
#  --taglen BLENGTH             Length of the duplex tag sequence. [12]
#  --spacerlen SLENGTH          Length of the spacer sequences used. [5]
#  --read_out ROUT              How often you want to be told what the 
#                               program is doing. [1000000]
#  --filt_spacer ADAPTERSEQ     Optional: Spacer sequence for filtering 
#                               on the presence of the spacer. This 
#                               could be thrown off by low quality 
#                               scores.
#  --tagstats                   Optional: Output tagstats file and make 
#                               distribution plot of tag family sizes.
#                               Requires matplotlib to be installed


import sys
import gzip
from argparse import ArgumentParser
from collections import defaultdict


def fastq_general_iterator(read1_fastq, read2_fastq):
    read1_readline = read1_fastq.readline
    read2_readline = read2_fastq.readline

    while True:
        read1_line = read1_readline()
        read2_line = read2_readline()

        if not read1_line and read2_line:
            return()
        if read1_line[0] == '@' and read2_line[0] == '@':
            break
        if (
                isinstance(read1_line[0], int) 
                or isinstance(read2_line[0], int)
                ):
            raise ValueError((f"FASTQ files may contain binary "
                              f"information or are compressed"
                              ))

    while read1_line and read2_line:

        if read1_line[0] != '@' or read2_line[0] != '@':
            print(f"{read1_line}, {read2_line}")
            raise ValueError((f"Records in FASTQ files should start "
                              f"with a '@' character. Files may be "
                              f"malformed or out of synch."
                              ))

        title_read1_line = read1_line[1:].rstrip()
        title_read2_line = read2_line[1:].rstrip()

        read1_seq_string = read1_readline().rstrip()
        read2_seq_string = read2_readline().rstrip()

        while True:
            read1_line = read1_readline()
            read2_line = read2_readline()

            if not read1_line and read2_line:
                raise ValueError((f"End of file without quality "
                                  f"information. Files may be malformed"
                                  f" or out of synch"
                                  ))
            if read1_line[0] == '+' and read2_line[0] == '+':
                break

            read1_seq_string += read1_line.rstrip()
            read2_seq_string += read2_line.rstrip()

        read1_quality_string = read1_readline().rstrip()
        read2_quality_string = read2_readline().rstrip()

        while True:
            read1_line = read1_readline()
            read2_line = read2_readline()

            if not read1_line or read2_line:
                break  # end of file
            if (
                    read1_line[0] == '@' 
                    and read2_line[0] == '@' 
                    and read1_line.isalpha() is not True 
                    and read2_line.isalpha() is not True
                    ):
                break

            read1_quality_string += read1_line.rstrip()
            read2_quality_string += read2_line.rstrip()

        yield (title_read1_line, 
               title_read2_line, 
               read1_seq_string, 
               read2_seq_string, 
               read1_quality_string, 
               read2_quality_string
               )
    #return()
    #raise StopIteration


def tag_extract_fxn(read_seq, blen):
    # This is the function that extracts the UID tags from both the
    # forward and reverse read.  Assigns read1 the sequence from some
    # position to the end, then read2 from some position to the end,
    # then assigns tag1 from the 5'-end to length of the UID tag for
    # read1 and then read 2.
    return(read_seq[0][:blen], read_seq[1][:blen])


def hdr_rename_fxn(read_title, read1_tag, read2_tag):
    # This function renames the header with the formatting of
    # *header coordinates,etc*, *tag from read1*, *tag from read2*,
    # *read designation from original header (for paired reads)*

    illumina = read_title.split(" ")[0].split(":")

    if len(illumina) == 7:
        #Illumina CASAVA >=1.8
        #e.g. @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
        readnum = read_title.split(" ")[1].split(":")[0]
        return(f"{read_title.split(' ')[0]}"
               f"|{read1_tag}{read2_tag}"
               f"/{readnum}"
               )
    elif len(illumina) == 5:
        #Illumina CASAVA >=1.4?
        #e.g. @HWUSI-EAS100R:6:73:941:1973#ATCGAT/1
        read_title = read_title.replace(' ', '_')
        return(f"{read_title.split('/')[0]}"
               f"|{read1_tag}{read2_tag}"
               f"/{read_title.split('/')[1]}"
               )
    else :
        raise ValueError("Unknown read name format: %s" % read_title)


def tag_stats(barcode_counts, outfile):
    family_size_dict = defaultdict(lambda:0)
    tagstat_file =  open(outfile + '.tagstats', 'w')
    total_tags = 0

    for value in barcode_counts:
        family_size_dict[value] += 1

    for family_size in family_size_dict.keys():
        family_size_dict[family_size] *= int(family_size)
        total_tags += int(family_size_dict[family_size])

    for family_size in sorted(family_size_dict.keys()):
        tagstat_file.write(
            f"{family_size}\t"
            f"{float(family_size_dict[family_size])/float(total_tags)}"
            f"\n"
            )

    tagstat_file.close()
    return(family_size_dict, total_tags)

def open_fastq(infile, outfile):
    if infile.endswith(".gz"):
        in_fh = gzip.open(infile, 'rt')
        out_fh = gzip.open(outfile + ".gz", 'wt')
    else:
        in_fh = open(infile, 'r')
        out_fh = open(outfile, 'w')
    return(in_fh, out_fh)

def main():
    parser =  ArgumentParser()
    parser.add_argument(
        '--infile1', 
        dest = 'infile1', 
        help = 'Path to FASTQ file for Read 1.', 
        required = True
        )
    parser.add_argument(
        '--infile2', 
        dest = 'infile2', 
        help = 'Path to FASTQ file for Read 2.', 
        required = True
        )
    parser.add_argument(
        '--outprefix', 
        dest = 'outfile', 
        help = (f'Prefix for output files. Will prepend onto file name '
                f'of ".fq.smi"'
                ), 
        required=True
        )
    parser.add_argument(
        '--taglen', 
        dest = 'taglen', 
        type = int, 
        default = 12,
        help = 'Length in bases of the duplex tag sequence.[12]'
        )
    parser.add_argument(
        '--spacerlen', 
        dest = 'spclen', 
        type = int, 
        default = 5,
        help = (f'Length in bases of the spacer sequence between duplex'
                f' tag and the start of target DNA. [5]'
                )
        )
    parser.add_argument(
        '--readout', 
        dest = 'readout', 
        type = int, 
        default = 1000000,
        help = (f'How many reads are processed before progress is '
                f'reported. [1000000]'
                )
        )
    parser.add_argument(
        '--filtspacer', 
        dest = 'spacer_seq', 
        type = str, 
        default = None,
        help = (f'Optional: Filter out sequences lacking the inputed '
                f'spacer sequence. Not recommended due to significant '
                f'base calling issues with the invariant spacer '
                f'sequence'
                )
        )
    parser.add_argument(
        '--tagstats', 
        dest = 'tagstats', 
        action = "store_true",
        help = (f'Optional: Output tagstats file and make distribution '
                f'plot of tag family sizes.  Requires matplotlib to be '
                f'installed.'
                )
        )
    parser.add_argument(
        '--reduce', 
        dest = 'reduce', 
        action = "store_true", 
        help = (f'Optional: Only output reads that will make a final '
                f'DCS read.  Will only work when the --tagstats option '
                f'is invoked.'
                )
        )
    o = parser.parse_args()

    if o.reduce and not o.tagstats:
        raise ValueError(f"--reduce option must be invoked with the "
                         f"--tagstats option."
                         )

    (read1_fastq, read1_output) = open_fastq(o.infile1, 
                                            o.outfile + '.seq1.smi.fq'
                                            )
    (read2_fastq, read2_output) = open_fastq(o.infile2, 
                                            o.outfile + '.seq2.smi.fq'
                                            )

    readctr = 0
    nospacer = 0
    goodreads = 0
    badtag = 0
    oldBad = 0
    barcode_dict = defaultdict(lambda:0)

    for readParts  in fastq_general_iterator(read1_fastq, read2_fastq):
        (read1_title, read2_title, read1_seq, 
         read2_seq, read1_qual, read2_qual) = readParts
        readctr += 1

        if (
                o.spacer_seq != None 
                and (
                    read1_seq[o.taglen:o.taglen + o.spclen] 
                    != o.spacer_seq 
                    or read2_seq[o.taglen:o.taglen + o.spclen] 
                    != o.spacer_seq
                    )
                ):
            nospacer += 1
        else:
            tag1, tag2 = tag_extract_fxn((read1_seq, read2_seq), 
                                         o.taglen
                                         )

            if (
                    (tag1.isalpha() 
                     and tag1.count('N') == 0
                     ) 
                    and (tag2.isalpha() 
                         and tag2.count('N') == 0
                         )
                    ):
                renamed_read1_title =  hdr_rename_fxn(read1_title,
                                                      tag1, 
                                                      tag2
                                                      )
                renamed_read2_title =  hdr_rename_fxn(read2_title, 
                                                      tag1, 
                                                      tag2
                                                      )
                read1_output.write((
                    f'@{renamed_read1_title}\n'
                    f'{read1_seq[o.taglen+o.spclen:]}\n'
                    f'+\n'
                    f'{read1_qual[o.taglen + o.spclen:]}\n'
                    ))
                read2_output.write((
                    f'@{renamed_read2_title}\n'
                    f'{read2_seq[o.taglen+o.spclen:]}\n'
                    f'+\n'
                    f'{read2_qual[o.taglen + o.spclen:]}\n'
                    ))
                goodreads += 1

                if o.tagstats:
                    barcode_dict[tag1 + tag2] += 1

            else:
                badtag += 1

        if readctr % o.readout is 0:
            sys.stderr.write((
                f"Total sequences processed: {readctr}\n"
                f"Sequences with passing tags: {goodreads}\n"
                f"Missing spacers: {nospacer}\n"
                f"Bad tags: {badtag}\n" 
                ))
            if badtag == oldBad + o.readout:
                sys.stderr.write((
                    f"Warning! Potential file error between lines "
                    f"{(readctr - o.readout) * 4} and {readctr * 4}."
                    ))
                oldBad = badtag

    read1_fastq.close()
    read2_fastq.close()
    read1_output.close()
    read2_output.close()

    sys.stderr.write((
        f"Total sequences processed: {readctr}\n"
        f"Sequences with passing tags: {goodreads}\n"
        f"Missing spacers: {nospacer}\n"
        f"Bad tags: {badtag}\n" 
        ))

    if o.tagstats:
        read_data_file = open(o.outfile + '_data.txt', 'w')
        sscs_count = 0
        dcs_count = 0
        dcs_tags_list = []
        family_size_dict, total_tags = tag_stats(barcode_dict.values(), 
                                                 o.outfile
                                                 )

        for tag in barcode_dict.keys():

            if barcode_dict[tag] >= 3:
                sscs_count += 1

                if (
                  tag[o.taglen:] + tag[:o.taglen] in barcode_dict 
                  and barcode_dict[tag[o.taglen:] + tag[:o.taglen]] >= 3
                  ):
                    dcs_count += 1

                    if o.reduce and tag not in dcs_tags_list:
                        dcs_tags_list.append(tag)
                        dcs_tags_list.append(tag[o.taglen:] 
                                             + tag[:o.taglen]
                                             )

        read_data_file.write(f'# Passing Reads\t'
                             f'# SSCS Reads\t'
                             f'# DCS Reads\t'
                             f'SSCS:DCS\n'
                             f'{goodreads}\t'
                             f'{sscs_count}\t'
                             f'{dcs_count}\t'
                             f'{float(sscs_count)/float(dcs_count)}\n'
                             )
        read_data_file.close()

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            x_value = []
            y_value = []

            for family_size in sorted(family_size_dict.keys()):
                x_value.append(family_size)
                y_value.append(float(family_size_dict[family_size]) 
                               / float(total_tags)
                               )

            plt.bar(x_value, y_value)
            plt.xlabel('Family Size')
            plt.ylabel('Proportion of Total Reads')
            plt.savefig(o.outfile + '.png', bbox_inches='tight')
                        
            plt.bar(x_value, y_value)
            plt.xlabel('Family Size')
            plt.ylabel('Proportion of Total Reads')
            plt.xlim([0,40])
            plt.savefig(o.outfile + '.zoom.png', bbox_inches='tight')

        except ImportError:
            sys.stderr.write(f'matplotlib not present. Only tagstats '
                             f'file will be generated.'
                             )

        if o.reduce:

            read1_fastq = open(o.outfile + '.seq1.smi.fq', 'r')
            read2_fastq = open(o.outfile + '.seq2.smi.fq', 'r')
            read1_output = open(o.outfile + '.seq1.reduced.fq', 'w')
            read2_output = open(o.outfile + '.seq2.reduced.fq', 'w')
            
            for readParts  in fastq_general_iterator(read1_fastq, 
                                                     read2_fastq
                                                     ):
                (read1_title, read2_title, read1_seq, 
                 read2_seq, read1_qual, read2_qual) = readParts
                testTag1 = read1_title.split('|')[1].split('/')[0]
                testTag2 = read2_title.split('|')[1].split('/')[0]
                if (
                        testTag1 in dcs_tags_list 
                        and testTag2 in dcs_tags_list
                        ):
                    read1_output.write(f'@{read1_title}\n'
                                       f'{read1_seq}\n'
                                       f'+\n'
                                       f'{read1_qual}\n'
                                       )
                    read2_output.write(f'@{read2_title}\n'
                                       f'{read2_seq}\n'
                                       f'+\n'
                                       f'{read2_qual}\n'
                                       )

            read1_fastq.close()
            read2_fastq.close()
            read1_output.close()
            read2_output.close()

if __name__ == "__main__":
    main()