from __future__ import print_function
import pysam
import argparse
import collections
from FisherExact import fisher_exact


def correct_from_pileup(bam_file_name, thresh_all, thresh_read, ref_len):
    """
    more sophisticated consensus correction
    inputs : bamfile of all reads, grouped by original concatamer with the RG tag
    Use bamfile of all reads to get consensus at each site (output N if no base >thresh_all present).
    Split reads based on whether they are +ve or -ve strand
    For each strand, go through site by site.
      Is there any read that seems to be non-consensus at this site (threshold of thresh_read to call confidently)
      If so then construct table with read group and base at that site. 
      Is there evidence of association between base and RG? (use Fisher's exact test, with 4 base categories)
    """
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    cons_seq = ["-"] * ref_len

    # make a dictionary of all read groups
    all_RGs = []
    for read in samfile.fetch():
        all_RGs.append(read.get_tag("RG"))

    unique_RGs = list(set(all_RGs))
    RG_dict = {RG: index for (index, RG) in enumerate(unique_RGs)}
    corrected_seqs = [["-"] * ref_len for _ in range(len(unique_RGs))]


    for pileupcolumn in samfile.pileup(stepper="all"):
        pos = pileupcolumn.pos
        # record base and read group. Want to be able to split reads by plus/minus strand
        bases_fw = []
        bases_rv = []
        RG_fw = []
        RG_rv = []
        # Don't filter on base quality, filter only on mapping quality
        pileupcolumn.set_min_base_quality(0)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                if pileupread.alignment.mapping_quality > 30:
                    if pileupread.alignment.is_reverse:
                        bases_rv.append(pileupread.alignment.query_sequence[pileupread.query_position])
                        RG_rv.append(pileupread.alignment.get_tag("RG"))
                    else:
                        bases_fw.append(pileupread.alignment.query_sequence[pileupread.query_position])
                        RG_fw.append(pileupread.alignment.get_tag("RG"))

        # First, work out whole sample consensus
        if len(bases_fw) == 0 or len(bases_rv) == 0:
            cons_base="-"
        else:
            counter_fw = collections.Counter(bases_fw)
            most_common_fw, freq_fw = counter_fw.most_common(1)[0]
            counter_rv = collections.Counter(bases_rv)
            most_common_rv, freq_rv = counter_rv.most_common(1)[0]
            if float(freq_fw)/len(bases_fw) > thresh_all and float(freq_rv)/len(bases_rv) > thresh_all \
                    and most_common_fw == most_common_rv:
                cons_base=most_common_fw
            else:
                cons_base="N"

        cons_seq[pos] = cons_base
        print(pos, cons_base)

        # Next, see if there is any evidence for polymorphism in reads - any non-consensus base at > thresh_read
        # in a read group with read depth > 1 at that position
        RG_fw_set = list(set(RG_fw))
        RG_rv_set = list(set(RG_rv))
        polymorphic_potential_fw = False
        polymorphic_potential_rv = False
        if cons_base!="-":
            for read_group in RG_fw_set:
                RG_bases = [base for i, base in enumerate(bases_fw) if RG_fw[i]==read_group]
                counter = collections.Counter(RG_bases)
                most_common, freq = counter.most_common(1)[0]
                if most_common!=cons_base and float(freq)/len(RG_bases) > thresh_read and freq > 1:
                    polymorphic_potential_fw = True
                    break

            for read_group in RG_rv_set:
                RG_bases = [base for i, base in enumerate(bases_rv) if RG_rv[i] == read_group]
                counter = collections.Counter(RG_bases)
                most_common, freq = counter.most_common(1)[0]
                if most_common != cons_base and float(freq) / len(RG_bases) > thresh_read and freq > 1:
                    polymorphic_potential_rv = True
                    break

        # Now, on both fw and rev strand data (separately) conduct fisher's exact test
        if polymorphic_potential_fw:
            print("fw polymorphism test")
            freq_mx = []
            for read_group in RG_fw_set:
                RG_bases = [base for i, base in enumerate(bases_fw) if RG_fw[i] == read_group]
                counter = collections.Counter(RG_bases)
                RG_freq = [counter[l] for l in "ACTG"]
                freq_mx.append(RG_freq)
            try:
                p_val = fisher_exact(freq_mx, simulate_pval=True)
            except ValueError:
                # ValueError will be raised if less than 2 non-zero columns or rows
                # i.e if all reads have the same base call, or if there is only one read group with coverage here
                # take the consensus within each read group in this case.
                p_val = 0.01

            if p_val > 0.02: # then just assign consensus base to all fw sequences
                for read_group in RG_fw_set:
                    corrected_seqs[RG_dict[read_group]][pos] = cons_base

            else:
                for seqi, read_group in enumerate(RG_fw_set):
                    base_distrib = freq_mx[seqi]
                    if float(max(base_distrib))/sum(base_distrib) > thresh_read:
                        base_index = [i for i, val in enumerate(base_distrib) if val==max(base_distrib)][0]
                        corrected_seqs[RG_dict[read_group]][pos] = 'ACTG'[base_index]
                    else:
                        corrected_seqs[RG_dict[read_group]][pos] = cons_base

        else: # no polymorphic potential, assign consensus base to all fw sequences
            for read_group in RG_fw_set:
                corrected_seqs[RG_dict[read_group]][pos] = cons_base

        # repeat all over again for reverse (should be a different set of read groups)
        if polymorphic_potential_rv:
            print("rv polymorphism test")
            freq_mx = []
            for read_group in RG_rv_set:
                RG_bases = [base for i, base in enumerate(bases_rv) if RG_rv[i] == read_group]
                counter = collections.Counter(RG_bases)
                RG_freq = [counter[l] for l in "ACTG"]
                freq_mx.append(RG_freq)
            try:
                p_val = fisher_exact(freq_mx, simulate_pval=True)
            except ValueError:
                p_val = 0.01

            if p_val > 0.02:  # then just assign consensus base to all sequences
                for read_group in RG_rv_set:
                    corrected_seqs[RG_dict[read_group]][pos] = cons_base

            else:
                for seqi, read_group in enumerate(RG_rv_set):
                    base_distrib = freq_mx[seqi]
                    if float(max(base_distrib)) / sum(base_distrib) > thresh_read:
                        base_index = [i for i, val in enumerate(base_distrib) if val == max(base_distrib)][0]
                        corrected_seqs[RG_dict[read_group]][pos] = 'ACTG'[base_index]
                    else:
                        corrected_seqs[RG_dict[read_group]][pos] = cons_base


        else:
            for read_group in RG_rv_set:
                corrected_seqs[RG_dict[read_group]][pos] = cons_base


    return "".join(cons_seq), corrected_seqs, unique_RGs


def correct_by_consensus_only(bam_file_name, thresh, ref_len):
    """
    correct reads simply by taking the consensus within a read group
    :param bam_file_name: bam file with all reads, grouped by cancatamer using RG tag
    :param thresh: threshold for base frequency to be considered consensus (otherwise record "N")
    :param ref_len: length of genotype reference that has been aligned to
    :return: 
    """
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    cons_seq = ["-"] * ref_len

    # make a dictionary of all read groups
    all_RGs = []
    for read in samfile.fetch():
        all_RGs.append(read.get_tag("RG"))

    unique_RGs = list(set(all_RGs))
    RG_dict = {RG: index for (index, RG) in enumerate(unique_RGs)}
    corrected_seqs = [["-"] * ref_len for _ in range(len(unique_RGs))]

    for pileupcolumn in samfile.pileup(stepper="all"):
        pos = pileupcolumn.pos
        # Don't filter on base quality, filter only on mapping quality
        pileupcolumn.set_min_base_quality(0)
        # record base and read group. Want to be able to split reads by plus/minus strand
        bases = []
        RG = []
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # query position is None if is_del or is_refskip is set.
                if pileupread.alignment.mapping_quality > 30:
                    bases.append(pileupread.alignment.query_sequence[pileupread.query_position])
                    RG.append(pileupread.alignment.get_tag("RG"))

        # First, work out whole sample consensus
        if len(bases) == 0:
            cons_base="-"
        else:
            counter = collections.Counter(bases)
            most_common, freq = counter.most_common(1)[0]
            if float(freq)/len(bases) > thresh:
                cons_base=most_common
            else:
                cons_base="N"

        cons_seq[pos] = cons_base

        # Get read group consensus, or set as whole sample consensus if base frequency is not above threshold
        RG_set = list(set(RG))
        for read_group in RG_set:
            RG_bases = [base for i, base in enumerate(bases) if RG[i] == read_group]
            counter = collections.Counter(RG_bases)
            base_distrib = [counter[l] for l in "ACTG"]

            if float(max(base_distrib))/sum(base_distrib) > thresh:
                base_index = [i for i, val in enumerate(base_distrib) if val==max(base_distrib)][0]
                corrected_seqs[RG_dict[read_group]][pos] = 'ACTG'[base_index]
            else:
                corrected_seqs[RG_dict[read_group]][pos] = cons_base

    print(''.join(cons_seq))

    return "".join(cons_seq), corrected_seqs, unique_RGs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('allbam', type=str, help="path to bamfile with all reads")
    parser.add_argument('fasta_out', type=str, help="output corrected fasta")
    parser.add_argument('--cons_type', type=str, default="error_correct",
                        help="Use error_correct or consensus_only")
    args = parser.parse_args()

    # I have over estimated the reference length, but the extra "-"s can be trimmed later
    if args.cons_type=="error_correct":
        patient_consensus, corrected_reads, unique_RGs = correct_from_pileup(args.allbam,
                                                                             thresh_all=0.4,
                                                                             thresh_read=0.6,
                                                                             ref_len=3400)
    else:
        patient_consensus, corrected_reads, unique_RGs = correct_by_consensus_only(args.allbam,
                                                                                   thresh=0.4,
                                                                                   ref_len=3400)
    # print results as a fasta file
    with open(args.fasta_out, 'w') as outfile:
        print(">patient_consensus",file=outfile)
        print(patient_consensus, file=outfile)
        for i, read_seq in enumerate(corrected_reads):
            print(">{}".format(unique_RGs[i]), file=outfile)
            print("".join(read_seq), file=outfile)

