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

    # Get info on primary alignments from alignment positions
    read_names, bam_flags, chromosomes, positions = primary_read_positions(args.bam_pos)
    print(len(read_names))

    with pysam.AlignmentFile(args.bam, 'rb') as bamfile:
        with open(args.out, 'w') as outfile:
            # For each primary alignment, use pysam to get aligned length of the read on the reference genome
            # Print details of full length reads to file
            # This method is not particularly efficient - better suited to large genomes with low coverage.
            for i, read_name in enumerate(read_names):
                bam_iter = bamfile.fetch(chromosomes[i], positions[i], positions[i] + 1)
                for read in bam_iter:
                    if read.query_name == read_name:
                        length = read.reference_length
                        if length >= hbv_length:
                            print("Read {} full length".format(read_name))
                            if bam_flags[i] & 16:
                                revcompl_flag=True
                            else:
                                revcompl_flag=False
                            # This anchor is no longer used to split concatamers - could be removed
                            seq_anchor = get_sequence_anchor(read.get_reference_sequence().upper(),
                                                             revcompl_flag, revcompl)

                            print(read_name, "\t", seq_anchor, "\t", length, "\t", revcompl_flag, file=outfile)
                        else:
                            print("Read {} too short: {}".format(read_name, length))





