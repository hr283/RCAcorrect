"""
A script for calculating the % non-consensus sties in Illumina and Nanopore data for the same patient
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
    total_aligned_bases = 0
    total_match = 0
    variable_positions = []
    for pileupcolumn in samfile.pileup(stepper="all"):
        position = pileupcolumn.pos
        coverage = pileupcolumn.n
        filt_coverage = 0
        bases = []
        consensus_base = list(consensus)[position]
        pileupcolumn.set_min_base_quality(20)
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.mapping_quality > 0:
                filt_coverage += 1
                if not pileupread.is_del and not pileupread.is_refskip:
                    bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
                if pileupread.is_del:
                    bases.append("-")

        if len(bases) != filt_coverage:
            # This was useful for debugging at one point but I think redundant now...
            print("Warning, coverage at position {} not equal to number of bases in Illumina pileup"
                  "{} vs {}".format(position, filt_coverage, len(bases)))

        counter = collections.Counter(bases)
        num_consensus = counter[consensus_base]
        non_consensus = 1 - float(num_consensus) / float(len(bases))
        percentages.append(non_consensus)
        if non_consensus < 0.01 or len(bases) - num_consensus < 2:
            total_aligned_bases += len(bases)
            total_match += num_consensus
        else:
            variable_positions.append(position)

    return percentages, variable_positions


def get_nanopore_variation(fasta_name, consensus, variable_positions=None):
    all_seqs = []
    total_aligned_bases = 0
    total_match = 0
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
        if num_seqs - counter["N"] - counter["-"] > 2:
            non_consensus = 1 - float(num_consensus)/float(num_seqs - counter["N"])
        else:
            non_consensus = "NA"
        percentages.append(non_consensus)
        if position not in variable_positions:
            total_aligned_bases += float(num_seqs)
            total_match += float(num_consensus)

    accuracy = float(total_match) / float(total_aligned_bases)
    return percentages, accuracy


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('consensuses', type=str, help="Path to corrected consensuses")
    parser.add_argument('consensus_name', type=str, help="Name of consensus to use as truth")
    parser.add_argument('--Illumina_bam', type=str, help="path to bamfile with Illumina", default=None)
    parser.add_argument('--Nanopore_fasta', type=str, help="path to fasta of aligned corrected reads", default=None)
    args = parser.parse_args()

    # get the 'truth' consensus sequence to use for comparison
    record_dict = SeqIO.to_dict(SeqIO.parse(args.consensuses, "fasta"))
    consensus_seq = record_dict[args.consensus_name].seq

    var_positions = None

    if args.Illumina_bam:
        illumina_proportions, var_positions = get_illumina_variation(args.Illumina_bam, consensus_seq)
        print(illumina_proportions)
        print(var_positions)
        print("Number of variable positions:", len(var_positions))

    if args.Nanopore_fasta:
        nanopore_proportions, nanopore_acc = get_nanopore_variation(args.Nanopore_fasta, consensus_seq, var_positions)
        print(nanopore_proportions)
        print("Nanopore accuracy:", nanopore_acc)

# python compare_variation.py ~/Dropbox/Nanopore/HBV/final_consensus_sequences.fasta p1331_RCA_Illumina_consensus1 --Illumina_bam ~/Dropbox/Nanopore/HBV/bams/P31_rep1_bwa.bam --Nanopore_fasta ~/Dropbox/Nanopore/HBV/corrected_fastas/p1331_corrected_filtered.fasta >~/Dropbox/Nanopore/HBV/p1331_variant_comparison_new.txt
# change in error rate with reps per read cutoff for 1331. Above command gives 98.0% (5 reps) accuracy
# 4reps: 97.9 %
# 5reps: 98.0 %
# 6reps: 98.2 %
# 8reps: 98.3 %
# 10reps: 98.5 %
# Number of reads: 276, 208, 158, 84, 41