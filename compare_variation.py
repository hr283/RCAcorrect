"""
A script for calculating the % non-consensus sites in Illumina and Nanopore data for the same patient
Used to get data for plotting Figure 6A.
"""
from __future__ import print_function
import pysam
import argparse
from Bio import SeqIO
import numpy as np
import collections


def get_illumina_variation(bam_file_name, consensus):
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    # Naming could be confusing - These are proportions not percentages
    percentages = []
    for pileupcolumn in samfile.pileup(stepper="all"):
        position = pileupcolumn.pos
        coverage = pileupcolumn.n
        bases = []
        consensus_base = list(consensus)[position]
        pileupcolumn.set_min_base_quality(0)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
            if pileupread.is_del:
                bases.append("-")

        if len(bases) != coverage:
            # if set_min_base_quality(0) is not used above, this warning may be produced
            print("Warning, coverage at position {} not equal to number of bases in Illumina pileup"
                  "{} vs {}".format(position, coverage, len(bases)))

        counter = collections.Counter(bases)
        num_consensus = counter[consensus_base]
        non_consensus = 1 - float(num_consensus) / float(coverage)
        percentages.append(non_consensus)

    return percentages


def get_nanopore_variation(fasta_name, consensus):
    all_seqs = []
    for record in SeqIO.parse(fasta_name, "fasta"):

        # Some of the sequences do not cover the extreme ends of the reference - treat as missing rather than deletion
        raw_seq = list(record.seq)
        for i, pos in enumerate(raw_seq):
            if pos == "-":
                raw_seq[i] = "N"
            else:
                break
        for i, pos in enumerate(raw_seq[::-1]):
            if pos =="-":
                raw_seq[-(i+1)] = "N"
            else:
                break
        all_seqs.append(raw_seq)

    # The first sequence is the consensus sequence
    # Since the sequences are aligned they should be nicely converted into array form
    seq_array = np.array(all_seqs[1:])
    num_seqs = len(all_seqs) - 1
    percentages = []
    for position, base in enumerate(list(consensus)):
        counter = collections.Counter(seq_array[:,position])
        num_consensus = counter[base]
        if num_seqs - counter["N"] > 2:
            non_consensus = 1 - float(num_consensus)/float(num_seqs - counter["N"])
        else:
            non_consensus = "NA"
        percentages.append(non_consensus)

    return percentages


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('consensuses', type=str, help="Path to corrected consensuses (fasta file)")
    parser.add_argument('consensus_name', type=str, help="Name of consensus to use as truth (as given in fasta)")
    parser.add_argument('--Illumina_bam', type=str, help="path to bamfile with Illumina data", default=None)
    parser.add_argument('--Nanopore_fasta', type=str, help="path to fasta of aligned corrected reads "
                                                           "(must be aligned to the given consensus)", default=None)
    args = parser.parse_args()

    # get the 'truth' consensus sequence to use for comparison
    record_dict = SeqIO.to_dict(SeqIO.parse(args.consensuses, "fasta"))
    consensus_seq = record_dict[args.consensus_name].seq

    if args.Illumina_bam:
        illumina_proportions = get_illumina_variation(args.Illumina_bam, consensus_seq)
        print(illumina_proportions)

    if args.Nanopore_fasta:
        nanopore_proportions = get_nanopore_variation(args.Nanopore_fasta, consensus_seq)
        print(nanopore_proportions)
