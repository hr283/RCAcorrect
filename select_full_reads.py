from __future__ import print_function
import pysam
import argparse
import csv

"""
A function to select reads for which there is a >= X length primary mapping 
More specifically, the script measures the length of the reference section that is aligned to and compares this to the 
 cut-off
"""

def primary_read_positions(read_align_file):
    """
    Takes the tab-delimited text file with the mapped position(s) of a read
    and then returns the chromosome and position of the primary read.
    :param read_align_file: path and filename of the text file
    :return name, flag, chromosome, position: lists containing these details for all mapped reads
    """
    flag = []
    chrom = []
    pos = []
    name =[]
    with open(read_align_file) as tabdelim:
        filereader = csv.reader(tabdelim, delimiter='\t')
        for row in filereader:
            this_flag = int(row[1])
            if this_flag & 256 | this_flag & 2048:
                pass
            else:
                if this_flag & 4:
                    print("Read {} unmapped".format(row[0]))
                else:
                    chrom.append(row[2])
                    pos.append(int(row[3]))
                    flag.append(this_flag)
                    name.append(row[0])

    return name, flag, chrom, pos


def get_sequence_anchor(read_ref, rc_flag, rc_func):
    if rc_flag:
        # then this read is reverse complemented - select the end of the reference sequence and
        # reverse complement it to get the anchor
        anchor = rc_func(read_ref[-100:])
    else:
        # Select the first x bases of the reference sequence as the anchor.
        anchor = read_ref[:100]

    return ''.join(anchor)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=str, help="path to bamfile")
    parser.add_argument('bam_pos', type=str, help="bam positions file")
    parser.add_argument('out', type=str, help="path to output file")
    args = parser.parse_args()

    # Cut-off for what counts as 'full length'
    hbv_length = 3200
    revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}[B] for B in x][::-1])


    with pysam.AlignmentFile(args.bam, 'rb') as bamfile:
        with open(args.out, 'w') as outfile:
            # For each primary alignment, use pysam to get aligned length of the read on the reference genome
            # Print details of full length reads to file
            bam_iter = bamfile.fetch()
            for read in bam_iter:
                if not read.is_secondary: # this should also have been checked in the previous step of the pipeline
                    length = read.reference_length
                    if length >= hbv_length:
                        if read.is_reverse:
                            revcompl_flag=True
                        else:
                            revcompl_flag=False

                        print(read.query_name, "\t", length, "\t", revcompl_flag, file=outfile)





