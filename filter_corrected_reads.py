from __future__ import print_function
import pysam
import argparse
from Bio import SeqIO
import collections

"""
Scripts to filter out read groups (=concatamers) that:
a) had (separate) sections mapping to the plus and minus strands
b) had (separate) sections mapping to both C and E genotype references
c) were not truly full length alignments
Read groups in category (a) and (c) were routinely filtered at the end of chop_and_correct_RCA.sh
(b) was used just to examine the 1331/1332 mix sample and look for potential recombinants.
"""


def detect_dual_strand_chimeras(bam_file_name, ref_name):
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    rv_RGs=[]
    fw_RGs=[]
    for read in samfile.fetch(ref_name):
        if not read.is_secondary and not read.is_supplementary:
            read_group = read.get_tag("RG")
            if read.is_reverse:
                rv_RGs.append(read_group)
            else:
                fw_RGs.append(read_group)

    unique_fw = set(fw_RGs)
    unique_rv = set(rv_RGs)
    dual_strand = unique_fw.intersection(unique_rv)
    return list(dual_strand)


def detect_dual_subtype_chimeras(bam_file_name):
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    C_RGs = []
    E_RGs = []
    for read in samfile.fetch("C_cons_seq"):
        if not read.is_secondary and not read.is_supplementary:
            read_group = read.get_tag("RG")
            C_RGs.append(read_group)

    for read in samfile.fetch("E_cons_seq"):
        if not read.is_secondary and not read.is_supplementary:
            read_group = read.get_tag("RG")
            E_RGs.append(read_group)

    unique_C = set(C_RGs)
    unique_E = set(E_RGs)
    dual_subtype = unique_C.intersection(unique_E)
    return list(dual_subtype)


def detect_non_fulllen(fasta_corrected):
    non_fulllen = []
    for record in SeqIO.parse(fasta_corrected, "fasta"):
        read_name = record.id
        num_blanks = collections.Counter(record.seq)["-"]
        if num_blanks > 300:
            non_fulllen.append(read_name)

    return non_fulllen


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('allbam', type=str, help="path to bamfile with all reads")
    parser.add_argument('fasta_corrected', type=str, help="Corrected fasta")
    parser.add_argument('fasta_out', type=str, help="Corrected, filtered fasta output")
    parser.add_argument('reference_name', type=str, help="Should be E_cons_seq or C_cons_seq")
    args = parser.parse_args()

    dual_strand_reads = detect_dual_strand_chimeras(args.allbam, args.reference_name)

    non_fulllen_reads = detect_non_fulllen(args.fasta_corrected)

    orig_reads = 0
    out_seqs=[]
    for record in SeqIO.parse(args.fasta_corrected, "fasta"):
        orig_reads+=1
        if record.id not in dual_strand_reads and record.id not in non_fulllen_reads:
            out_seqs.append(record)

    SeqIO.write(out_seqs, args.fasta_out, "fasta")

    print('\n{}/{} sequences have dual strand mappings:'.format(len(dual_strand_reads), orig_reads))
    print(dual_strand_reads)
    print('\n{}/{} sequences have excessive gaps:'.format(len(non_fulllen_reads), orig_reads))
    print(non_fulllen_reads)
    print('\n{} sequences written to {}'.format(len(out_seqs), args.fasta_out))
