from __future__ import print_function
import pysam
import argparse
import csv

"""
A function to select reads for which there is a >= X length primary mapping 
More specifically, the script measures the length of the reference section that is aligned to and compares this to the 
 cut-off
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=str, help="path to bamfile")
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





