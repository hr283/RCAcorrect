from __future__ import print_function
import argparse
from Bio import SeqIO
import collections

def select_variant_sites(fasta_name, sites, fasta_out):
    combinations = []
    for record in SeqIO.parse(fasta_name, "fasta"):

        # Some of the sequences do not cover the extreme ends of the reference - treat as missing rather than deletion
        raw_seq = list(record.seq)
        select_sites = "".join([raw_seq[site] for site in sites])
        combinations.append(select_sites)

    seq_set = set(combinations)
    seq_uniq_list = list(seq_set)
    counter = collections.Counter(combinations)
    for i, seq in enumerate(seq_uniq_list):
        # let's say we're interested in everything that appears at > 1% frequency...
        count = counter[seq]
        if float(count)/len(combinations) > 0.01:
            print(">seq_" + seq + "_" + str(count), file=fasta_out)
            print(seq, file=fasta_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_in', type=str, help="Path to corrected concatemers")
    parser.add_argument('fasta_out', type=str, help="Path to output fasta")
    parser.add_argument('sample_num', type=int, help="1331, 1332 or 1348")
    args = parser.parse_args()

    if args.sample_num == 1331:
        variants = [1041 - 1, 1054 - 1, 1936 - 1, 2134 - 1]  # account for zero-indexing
    elif args.sample_num == 1332:
        print("No variant sites for this patient")
        exit(0)
    elif args.sample_num == 1348:
        variants = [400 - 1, 841 - 1, 915 - 1, 1425 - 1, 2189 - 1]
    else:
        print("Incorrect sample number")
        exit(1)
    with open(args.fasta_out, 'w') as outfile:
        select_variant_sites(args.fasta_in, variants, outfile)