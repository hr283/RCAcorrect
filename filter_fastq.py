from __future__ import print_function
import csv
import argparse
from Bio import SeqIO


"""
These functions filter a fastq file based on a list of read names
The output is a new fastq containing only those reads
"""

def get_read_names(read_name_handle):
    filereader = csv.reader(read_name_handle, delimiter='\t')
    read_names = []
    for row in filereader:
        read_names.append(row[0].strip(" "))

    return read_names


def filter_fastq(fastq_fname, read_names, output_handle):
    record_dict = SeqIO.index(fastq_fname, "fastq")
    for read_name in read_names:
        SeqIO.write(record_dict[read_name], output_handle, "fastq")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('read_file', type=str, help="path file with list of read names")
    parser.add_argument('fastq_in', type=str, help="large fastq to be filtered")
    parser.add_argument('fastq_out', type=str, help="output filtered fastq")
    args = parser.parse_args()

    with open(args.read_file) as tabdelim:
        read_names_list = get_read_names(tabdelim)

    with open(args.fastq_out, 'w') as outfile:
        filter_fastq(args.fastq_in, read_names_list, outfile)