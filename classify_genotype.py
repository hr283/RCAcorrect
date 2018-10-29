from __future__ import print_function
import pysam
import argparse
import os
import collections
from Bio import SeqIO
from filter_corrected_reads import detect_dual_subtype_chimeras
import numpy as np

def takeFirst(elem):
    return elem[0]


def takeSecond(elem):
    return elem[1]


def get_discordant_positions(aligned_seqs):
    record_dict = SeqIO.to_dict(SeqIO.parse(aligned_seqs, "fasta"))
    p1331 = record_dict["p1331_RCA_Illumina_consensus1"].seq
    p1332 = record_dict["p1332_RCA_Illumina_consensus1"].seq

    # record the discordant positions relative to numbering in C and E
    discordant_positions_relC = []
    pos = 0
    for i, base in enumerate(p1331):
        if base != "-":
            if base != p1332[i] and p1332[i] != "-":
                discordant_positions_relC.append((pos, base, p1332[i]))
            pos += 1

    discordant_positions_relE = []
    pos = 0
    for i, base in enumerate(p1332):
        if base != "-":
            if base != p1331[i] and p1331[i] != "-":
                discordant_positions_relE.append((pos, base, p1331[i]))
            pos += 1

    return discordant_positions_relC, discordant_positions_relE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('consensuses', type=str, help="Path to aligned consensuses")
    parser.add_argument('bam_file', type=str, help="Path to bamfile of chopped reads")
    parser.add_argument('out_file', type=str, help="output file")
    args = parser.parse_args()

    chimera_names = detect_dual_subtype_chimeras(args.bam_file)

    diffs_relC, diffs_relE = get_discordant_positions(args.consensuses)

    print(len(diffs_relC))

    samfile = pysam.AlignmentFile(args.bam_file, "rb")
    # now, I think we want to fetch reads aligned to both C and E and then later order the returned rows according to read name and number?
    RGs = []
    section_nums = []
    aligned_to = []
    classification_vectors = []
    classification_pc = []
    # where 0 = no data/neither C nor E, 1 = C, 2 = E
    for read in samfile.fetch("C_cons_seq"):
        classification_vec = []
        if not read.is_secondary and not read.is_supplementary:
            read_group = read.get_tag("RG")
            RGs.append(read_group)
            read_num = read.query_name.split("-")[-1]
            section_nums.append(int(read_num))
            seq = read.query_sequence
            pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
            # create a dictionary where the keys are the reference sequence positions
            dict_pairs = dict(zip([l[1] for l in pairs], [l[0] for l in pairs]))
            for position, C_base, E_base in diffs_relC:
                try:
                    read_base_pos = dict_pairs[position]
                    if read_base_pos:
                        read_base = seq[read_base_pos]
                    else:
                        read_base = "-"
                except KeyError:
                    read_base_pos = None
                    read_base = "-"

                if read_base == C_base:
                    classification_vec.append(1)
                elif read_base == E_base:
                    classification_vec.append(2)
                else:
                    classification_vec.append(0)

            counter = collections.Counter(classification_vec)
            pc = float(counter[1])/(counter[1] + counter[2])

            classification_vectors.append(classification_vec)
            classification_pc.append(pc)
            aligned_to.append("C")

    for read in samfile.fetch("E_cons_seq"):
        classification_vec = []
        if not read.is_secondary and not read.is_supplementary:
            read_group = read.get_tag("RG")
            RGs.append(read_group)
            read_num = read.query_name.split("-")[-1]
            section_nums.append(int(read_num))
            seq = read.query_sequence
            pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
            # create a dictionary where the keys are the reference sequence positions
            dict_pairs = dict(zip([l[1] for l in pairs], [l[0] for l in pairs]))
            for position, E_base, C_base in diffs_relE:
                try:
                    read_base_pos = dict_pairs[position]
                    if read_base_pos:
                        read_base = seq[read_base_pos]
                    else:
                        read_base = "-"
                except KeyError:
                    read_base_pos = None
                    read_base = "-"

                if read_base == C_base:
                    classification_vec.append(1)
                elif read_base == E_base:
                    classification_vec.append(2)
                else:
                    classification_vec.append(0)
            counter = collections.Counter(classification_vec)
            pc = float(counter[1]) / (counter[1] + counter[2])

            classification_vectors.append(classification_vec)
            classification_pc.append(pc)
            aligned_to.append("E")

    all_out_data = zip(RGs, section_nums, aligned_to, classification_pc, classification_vectors)

    # sort list with key
    out_data_sort1 = sorted(all_out_data, key=takeSecond)
    out_data_sort2 = sorted(out_data_sort1, key=takeFirst)

    with open(args.out_file, 'w') as f:
        for row in out_data_sort2:
            print(row[0], ", ", row[1], ", ", row[2], ", ", row[3], ", ", row[4], file=f)




