from __future__ import print_function
import pysam
import argparse
import collections
import numpy as np
import datetime
import scipy
from FisherExact import fisher_exact
from Bio import SeqIO


def print_vcf_header(ref_name, ref_len, bam_path, vcf_handle):
    header_lines=["##fileformat=VCFv4.2",
                  "##source=correct_by_consensus.py",
                  "##fileDate={:%d-%m-%Y}".format(datetime.date.today()),
                  "##contig=<ID={},length={}>".format(ref_name, ref_len),
                  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth of coverage at this position\">",
                  "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Observations of the alternate allele on individual reads as a proportion of total coverage\">",
                  "##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of forward strand concatamers without alternate allele observations\">",
                  "##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reverse strand concatamers without alternate allele observations\">",
                  "##INFO=<ID=SAF,Number=A,Type=Integer,Description=\"Number of alternate observations in forward strand concatamers\">",
                  "##INFO=<ID=SAR,Number=A,Type=Integer,Description=\"Number of alternate observations in reverse strand concatamers\">",
                  "##INFO=<ID=SBP,Number=A,Type=Float,Description=\"Strand bias p value, calculated using chisq test on [[SRF, SRR], [SAF, SAR]]\">",
                  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"]
    [print(line, file=vcf_handle) for line in header_lines]
    vcf_header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", bam_path]
    print("\t".join(vcf_header), file=vcf_handle)


def add_to_vcf(chrom, pos, id, ref, alt, qual, info, strand_counts, p_val, gt, filehandle):

    info_field = 'DP={};AF={};SRF={};SRR={};SAF={};SAR={};SBP={}'.format(info[0], info[1],
                                                                  strand_counts[0], strand_counts[1],
                                                                  strand_counts[2], strand_counts[3],
                                                                  p_val)
    format_field = '{}'.format(gt)
    vcf_line=[chrom, str(pos + 1), str(id), ref, alt, str(qual), ".", info_field, "GT", format_field]
    print("\t".join(vcf_line), file=filehandle)


def correct_from_pileup(bam_file_name, thresh_all, thresh_read, ref_seq, ref_name, vcf_handle):
    """
    more sophisticated consensus correction
    inputs : bamfile of all reads, grouped by original concatamer with the RG tag
    Use bamfile of all reads to get consensus at each site (output N if no base >thresh_all present).
    Split reads based on whether they are +ve or -ve strand
    For each strand, go through site by site.
      Is there any read that seems to be non-consensus at this site (threshold of thresh_read to call confidently)
      If so then construct table with read group and base at that site. 
      Is there evidence of association between base and RG? (use Fisher's exact test, with 4 base categories)
    Additionally, is there any significant evidence of strand bias at this site?
      
    """
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    cons_seq = ["-"] * len(ref_seq)
    next_id = 1

    base_mapping = {key:value for (key, value) in zip(['A','C','T','G'], range(4))}

    # make a dictionary of all read groups
    all_RGs = []
    for read in samfile.fetch():
        all_RGs.append(read.get_tag("RG"))

    unique_RGs = list(set(all_RGs))
    RG_dict = {RG: index for (index, RG) in enumerate(unique_RGs)}
    corrected_seqs = [["-"] * len(ref_seq) for _ in range(len(unique_RGs))]


    for pileupcolumn in samfile.pileup(stepper="all"):
        pos = pileupcolumn.pos
        ref_base = ref_seq[pos]
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
            if (float(freq_fw) + float(freq_rv))/(len(bases_fw) + len(bases_rv)) > thresh_all \
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
        fw_p_val = None
        # this freq mx will be used later to record strandedness of potential variants
        freq_mx_fw = []
        for read_group in RG_fw_set:
            RG_bases = [base for i, base in enumerate(bases_fw) if RG_fw[i] == read_group]
            counter = collections.Counter(RG_bases)
            RG_freq = [counter[l] for l in "ACTG"]
            freq_mx_fw.append(RG_freq)

        # repeat all over again for reverse (should be a different set of read groups)
        rv_p_val = None
        freq_mx_rv = []
        for read_group in RG_rv_set:
            RG_bases = [base for i, base in enumerate(bases_rv) if RG_rv[i] == read_group]
            counter = collections.Counter(RG_bases)
            RG_freq = [counter[l] for l in "ACTG"]
            freq_mx_rv.append(RG_freq)

        if polymorphic_potential_fw and polymorphic_potential_rv:
            print("polymorphism tests")
            try:
                fw_p_val = fisher_exact(freq_mx_fw, simulate_pval=True)
            except ValueError:
                # ValueError will be raised if less than 2 non-zero columns or rows
                # i.e if all reads have the same base call, or if there is only one read group with coverage here
                # take the consensus across all sites here.
                fw_p_val = None

            try:
                rv_p_val = fisher_exact(freq_mx_rv, simulate_pval=True)
            except ValueError:
                rv_p_val = None

            # Filter for strand specificity
            counter_all = collections.Counter(bases_fw + bases_rv)
            top_two = counter_all.most_common(2)
            if cons_base == "N":
                working_ref = ref_base
            else:
                working_ref = cons_base
            if working_ref == top_two[0][0]:
                alt_base, freq_alt = top_two[1]
            else:
                alt_base, freq_alt = top_two[0]

            srf = sum(np.array(freq_mx_fw)[:, base_mapping[alt_base]] == 0)
            srr = sum(np.array(freq_mx_rv)[:, base_mapping[alt_base]] == 0)
            saf = sum(np.array(freq_mx_fw)[:, base_mapping[alt_base]] > 0)
            sar = sum(np.array(freq_mx_rv)[:, base_mapping[alt_base]] > 0)

            try:
                # calculate chisq statistic for strand/allele counts, with Yates' correction applied:
                chi2, p_val, dof, ex = scipy.stats.chi2_contingency([[srf, srr],
                                                                     [saf, sar]])
                if p_val > 0.01:
                    strand_bias_pass = True
                else:
                    strand_bias_pass = False

            except ValueError:  # this is due to zero expected frequency
                strand_bias_pass = False

            # correct sites as 'variants' if there is strong concatamer-association in both plus and minus strand data,
            # and if there is NOT strand specificity at a significance level of p <0.01
            if fw_p_val < 0.01 and rv_p_val < 0.01 and strand_bias_pass:
                print("variant at site ", pos)
                for seqi, read_group in enumerate(RG_rv_set):
                    base_distrib = freq_mx_rv[seqi]
                    if float(max(base_distrib)) / sum(base_distrib) > thresh_read:
                        base_index = [i for i, val in enumerate(base_distrib) if val == max(base_distrib)][0]
                        corrected_seqs[RG_dict[read_group]][pos] = 'ACTG'[base_index]
                    else:
                        corrected_seqs[RG_dict[read_group]][pos] = cons_base

                for seqi, read_group in enumerate(RG_fw_set):
                    base_distrib = freq_mx_fw[seqi]
                    if float(max(base_distrib)) / sum(base_distrib) > thresh_read:
                        base_index = [i for i, val in enumerate(base_distrib) if val == max(base_distrib)][0]
                        corrected_seqs[RG_dict[read_group]][pos] = 'ACTG'[base_index]
                    else:
                        corrected_seqs[RG_dict[read_group]][pos] = cons_base

            else: # Not enough evidence for polymorphism, assign consensus base at this site to all concatamers
                for read_group in RG_fw_set:
                    corrected_seqs[RG_dict[read_group]][pos] = cons_base
                for read_group in RG_rv_set:
                    corrected_seqs[RG_dict[read_group]][pos] = cons_base

        else: # no polymorphic potential, assign consensus base to all sequences
            for read_group in RG_fw_set:
                corrected_seqs[RG_dict[read_group]][pos] = cons_base
            for read_group in RG_rv_set:
                corrected_seqs[RG_dict[read_group]][pos] = cons_base

        # option to be lenient in initial vcf output (note the 'or') - can filter later
        if polymorphic_potential_fw or polymorphic_potential_rv:
            # note that polymorphic_potential can only be true if cons_base!="-"
            # we will define the alt base as the most common base that is not the reference base
            # these next few lines should be identical to those used to detect strand specificity for correction
            counter_all = collections.Counter(bases_fw + bases_rv)
            top_two = counter_all.most_common(2)
            if cons_base=="N":
                working_ref = ref_base
            else:
                working_ref = cons_base
            if working_ref == top_two[0][0]:
                alt_base, freq_alt = top_two[1]
            else:
                alt_base, freq_alt = top_two[0]

            corrected_bases = [seq[pos] for seq in corrected_seqs]
            corrected_counter = collections.Counter(corrected_bases)
            AO = corrected_counter[alt_base]
            DP = corrected_counter["A"] + corrected_counter["C"] + corrected_counter["G"] + corrected_counter["T"]
            if DP > 0:
                prop_alt = float(AO)/float(DP)
            else:
                prop_alt = 0.0

            srf = sum(np.array(freq_mx_fw)[:, base_mapping[alt_base]] == 0)
            srr = sum(np.array(freq_mx_rv)[:, base_mapping[alt_base]] == 0)
            saf = sum(np.array(freq_mx_fw)[:, base_mapping[alt_base]] > 0)
            sar = sum(np.array(freq_mx_rv)[:, base_mapping[alt_base]] > 0)
            try:
                chi2, p_val, dof, ex = scipy.stats.chi2_contingency([[srf, srr],
                                                                     [saf, sar]])
            except ValueError:  # this is due to zero expected frequency
                p_val = ""

            # calculate phred based qual score, based on -10log_10 prob(call in ALT is wrong)
            if fw_p_val > 0 and rv_p_val > 0:
                qual = -10 * np.log10(fw_p_val * rv_p_val)
            elif fw_p_val > 0:
                qual = -10 * np.log10(fw_p_val)
            elif rv_p_val > 0:
                qual = -10 * np.log10(rv_p_val)
            elif fw_p_val == 0 or rv_p_val == 0:
                # set an (arbirary) high qual score. NB: I don't think this ever happens.
                qual = 100
            else:
                # not enough data to assess polymorphic potential based on distribution across read groups
                qual = "."

            # This if clause makes sure only the variants in Supp Table 3 are output. Remove if want to ouput more
            # note that not necessary to add strand bias to if clause assuming that sequences have been corrected
            # with strand bias as a criterion - only concatemer-corrected sites will have prop_alt > 0.
            if fw_p_val < 0.01 and rv_p_val < 0.01 and prop_alt > 0.1:
                add_to_vcf(chrom=ref_name, pos=pos, id=next_id, ref=working_ref, alt=alt_base, qual=qual,
                           info=[len(bases_fw + bases_rv), prop_alt],
                           strand_counts=[srf, srr, saf, sar], p_val=p_val, gt="0/1", filehandle=vcf_handle)
                next_id += 1


    return "".join(cons_seq), corrected_seqs, unique_RGs


def correct_by_consensus_only(bam_file_name, thresh, ref_seq):
    """
    correct reads simply by taking the consensus within a read group
    :param bam_file_name: bam file with all reads, grouped by cancatamer using RG tag
    :param thresh: threshold for base frequency to be considered consensus (otherwise record "N")
    :param ref_seq: sequence of genotype reference that has been aligned to
    :return: 
    """
    samfile = pysam.AlignmentFile(bam_file_name, "rb")
    cons_seq = ["-"] * len(ref_seq)

    # make a dictionary of all read groups
    all_RGs = []
    for read in samfile.fetch():
        all_RGs.append(read.get_tag("RG"))

    unique_RGs = list(set(all_RGs))
    RG_dict = {RG: index for (index, RG) in enumerate(unique_RGs)}
    corrected_seqs = [["-"] * len(ref_seq) for _ in range(len(unique_RGs))]

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
    parser.add_argument('ref_seq', type=str, help="path to fasta file containing (single) reference sequence")
    parser.add_argument('fasta_out', type=str, help="output corrected fasta")
    parser.add_argument('vcf_out', type=str, help="output vcf of variant calls")
    parser.add_argument('--cons_type', type=str, default="error_correct",
                        help="Use error_correct or consensus_only")
    args = parser.parse_args()


    try:
        record = SeqIO.read(args.ref_seq, "fasta")
    except ValueError:
        print("error reading reference fasta. Does it contain exactly one record?")
        exit(1)

    ref = list(record.seq)
    ref_id = record.id

    with open(args.vcf_out, 'w') as vcf_out:
        print_vcf_header(ref_id, len(ref), args.allbam, vcf_out)

        if args.cons_type=="error_correct":
            patient_consensus, corrected_reads, unique_RGs = correct_from_pileup(args.allbam,
                                                                                 thresh_all=0.4,
                                                                                 thresh_read=0.6,
                                                                                 ref_seq=ref,
                                                                                 ref_name=ref_id,
                                                                                 vcf_handle=vcf_out)

        else:
            patient_consensus, corrected_reads, unique_RGs = correct_by_consensus_only(args.allbam,
                                                                                       thresh=0.4,
                                                                                       ref_seq=ref)
    # print results as a fasta file
    with open(args.fasta_out, 'w') as outfile:
        print(">patient_consensus",file=outfile)
        print(patient_consensus, file=outfile)
        for i, read_seq in enumerate(corrected_reads):
            print(">{}".format(unique_RGs[i]), file=outfile)
            print("".join(read_seq), file=outfile)

