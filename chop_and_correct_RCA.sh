#!/bin/bash
#$ -P lunter.prjc
#$ -S /bin/bash
#$ -q short.qc
#$ -cwd
# Hannah Roberts, Feb 2018

# (1) Select 100bp anchor at start of reference genome
# (2) For each full length read: Find this anchor within the read using blast
# (3)   Chop up read based on anchor locations -> create fastq for this read
# (4)   Align these reads to (single) reference genome
# (5)   Correct sequences and write corrected consensus to fasta file
# (6) Filter corrected reads to remove those with dual stand mappings or excessive gaps

# Requires blast, bwa, samtools
BWA_PATH=/apps/well/bwa/0.7.12/bwa
SAMTOOLS_PATH=/apps/well/samtools/1.4.1/bin/samtools

set -e
set -x

# Arguments to command line -
# (1) workspace, where output directory will be created. Also where the python scripts are assumed to be located.
# must be the same as that given to select_full_reads.sh
# (2) path to genotype reference fasta (single copy)
# (3) patient/sample ID. e.g. 31, 32, 1348RCA
# must be the same as that given to select_full_reads.sh
# (4) correction/consensus type. Should be one of 'consensus_only' (simple within-concatamer consensus) and 'error_correct'
# (involves detection and correction of repeated Nanopore errors).
WORKSPACE=$( echo $@ | awk '{print $1}' )
REFPATH=$( echo $@ | awk '{print $2}' )
pID=$( echo $@ | awk '{print $3}')
CORRECT_TYPE=$( echo $@ | awk '{print $4}' )

NEW_FASTQ_DIR=${WORKSPACE}/p${pID}/full_len_fastqs_3200
NEW_FASTA_DIR=${WORKSPACE}/p${pID}/full_len_fastas_3200

cd ${WORKSPACE}
mkdir -p ${NEW_FASTQ_DIR}
mkdir -p ${NEW_FASTA_DIR}

# These files are produced by select_full_reads.sh
FASTQ_FILTERED=${WORKSPACE}/p${pID}/${pID}_pass-trimmed-3200.fastq
READ_INFO=${WORKSPACE}/p${pID}/full_length_reads_3200.txt

# Step (1) Select first 100 bases of reference sequence (assumes reference sequence is single-line not multi-line fasta)
ANCHOR=$( awk 'NR==2{print $1}' $REFPATH | cut -c1-100 )
ANCHOR_FASTA=${WORKSPACE}/p${pID}/splitting_anchor.fa
printf "%s\n" ">anchor" $ANCHOR > $ANCHOR_FASTA
# and reverse complement version...
RC_ANCHOR=$( echo $ANCHOR | rev | tr "ATGC" "TACG" )
RC_ANCHOR_FASTA=${WORKSPACE}/p${pID}/splitting_anchor_rc.fa
printf "%s\n" ">rc_anchor" $RC_ANCHOR > $RC_ANCHOR_FASTA

TMP_BLAST_DIR=${WORKSPACE}/p${pID}/tmp_blast4correct
mkdir -p $TMP_BLAST_DIR

# Steps (2) and (3)
BLAST_RESULTS=${WORKSPACE}/p${pID}/blast_chop_results_3200.txt
if [[ ! -e ${BLAST_RESULTS} ]]; then
  while read -r firstline; do
    read_name=$( echo $firstline | awk '{print $1}' | cut -c2- )
    read -r sequence
    read -r p
    read -r qscores
    # make fasta file for read, for use with blase
    read_fasta=$TMP_BLAST_DIR/$read_name.fa
    printf "%s\n" ">"$read_name $sequence > $read_fasta

    # work out whether sequence is rc and choose correct anchor
    # need to add tab to read_name for cases where reads have been split by porechop
    rc_flag=$( grep "${read_name} " $READ_INFO | awk '{print $4}' )
    if [ $rc_flag == "True" ]; then
      anchor_fa=$RC_ANCHOR_FASTA
    else
      anchor_fa=$ANCHOR_FASTA
    fi   

    # make blast db and blast anchor against read sequence
    /apps/well/ncbi-blast+/2.6.0/bin/makeblastdb -dbtype nucl -in $read_fasta -out $TMP_BLAST_DIR/${read_name}_db
    /apps/well/ncbi-blast+/2.6.0/bin/blastn -outfmt "6 sseqid pident evalue length sstart send bitscore" \
      -query $anchor_fa -db $TMP_BLAST_DIR/${read_name}_db -out $TMP_BLAST_DIR/${read_name}-blastout.txt \
      -reward 5 -penalty -4 -gapopen 8 -gapextend 6 -dust no
    sort -k5 -n $TMP_BLAST_DIR/${read_name}-blastout.txt > $TMP_BLAST_DIR/${read_name}-blastsort.txt

    prev_start=1 # the start position of the previous section in the read
    # Reads with > 1 short sections (additional to the last section) will be discarded
    short_section_count=0
    while read blast_line; do
      this_start=$( echo $blast_line | awk '{print $5}' )
      gap=$(( ${this_start} - ${prev_start} ))
      if (( $gap < 2900 )); then short_section_count=$(( $short_section_count + 1 )); fi
      # If there has been a big gap since the last chop, further chop this section into equal length intervals
      # This logic assumes that big gaps are of approximate size 6000, 9000, 12000, etc. (i.e. whole genome multiples)
      if (( $gap > 6000 )); then 
        rep=$(( $gap / 2900 )) 
        interval=$(( $gap / $rep ))
        new_start=${prev_start}
        while ((${new_start} < $((${this_start} - 5800)))); do
          new_end=$(($new_start + $interval))
          this_sequence=$( echo $sequence | cut -c${new_start}-${new_end} )
          this_qscore=$( echo $qscores | cut -c${new_start}-${new_end} )
          printf "%s\n" "@"${read_name}-$new_start $this_sequence "+" $this_qscore >> $NEW_FASTQ_DIR/${read_name}-chop.fastq
          printf "%s\n" ">"${read_name}-$new_start $this_sequence >> $NEW_FASTA_DIR/${read_name}-chop.fasta
          new_start=$new_end
        done
        prev_start=$new_start
      fi
      this_sequence=$( echo $sequence | cut -c${prev_start}-${this_start} )
      this_qscore=$( echo $qscores | cut -c${prev_start}-${this_start} )
      printf "%s\n" "@"${read_name}-$prev_start $this_sequence "+" $this_qscore >> $NEW_FASTQ_DIR/${read_name}-chop.fastq
      printf "%s\n" ">"${read_name}-$prev_start $this_sequence >> $NEW_FASTA_DIR/${read_name}-chop.fasta
      prev_start=$this_start
    done < $TMP_BLAST_DIR/${read_name}-blastsort.txt

    # print the last sections to file
    this_sequence=$( echo $sequence | cut -c${prev_start}- )
    this_qscore=$( echo $qscores | cut -c${prev_start}- )
    printf "%s\n" "@"${read_name}-$prev_start $this_sequence "+" $this_qscore >> $NEW_FASTQ_DIR/${read_name}-chop.fastq
    printf "%s\n" ">"${read_name}-$prev_start $this_sequence >> $NEW_FASTA_DIR/${read_name}-chop.fasta
    cat $TMP_BLAST_DIR/${read_name}-blastsort.txt >> $BLAST_RESULTS
    if ((${short_section_count} > 1)); then rm $NEW_FASTQ_DIR/${read_name}-chop.fastq; rm $NEW_FASTA_DIR/${read_name}-chop.fasta; fi
    rm $TMP_BLAST_DIR/${read_name}*
  done < $FASTQ_FILTERED
  rm -r $TMP_BLAST_DIR
else
  echo "Not chopping reads; already done"
fi

# Step (4) Re-align fastqs to single-copy reference
FASTA_CORRECTED=$WORKSPACE/p${pID}/${pID}_pass-corrected.fasta
BAMS_PER_READ=$WORKSPACE/p${pID}/bams_per_read
VCF_OUT=$WORKSPACE/p${pID}/variants.vcf
mkdir -p $BAMS_PER_READ

module load python/2.7.11
source ~/hbv_py/bin/activate

for fastq in $( ls $NEW_FASTQ_DIR ); do
  fastq_len=$(wc -l $NEW_FASTQ_DIR/$fastq | awk '{print $1}')
  # want to just correct reads for which there are at least 5 concatenated (partial + full) copies of the genome
  # This corresponds to at least 3 full (or near-full) copies of the genome
  if [ "$fastq_len" -ge "20" ]; then
    $BWA_PATH mem -x ont2d -R "@RG\tID:1\tSM:HBV" -M \
      $REFPATH \
      $NEW_FASTQ_DIR/$fastq \
      | $SAMTOOLS_PATH view -bS - \
      | $SAMTOOLS_PATH sort -o $BAMS_PER_READ/$( basename $fastq .fastq ).bam -
    $SAMTOOLS_PATH index $BAMS_PER_READ/$( basename $fastq .fastq ).bam
  fi
done    

# Step (5) Correct reads
# merge all the individual bams and add read group tag, for use in correcting by concatamer
for r in $( ls $BAMS_PER_READ/*.bam ); do echo $r >> $BAMS_PER_READ/all_bam_names.txt ; done
$SAMTOOLS_PATH merge -b $BAMS_PER_READ/all_bam_names.txt -r $BAMS_PER_READ/all.bam
$SAMTOOLS_PATH index $BAMS_PER_READ/all.bam

python correct_by_consensus.py $BAMS_PER_READ/all.bam $REFPATH $FASTA_CORRECTED $VCF_OUT --cons_type $CORRECT_TYPE

# Step (6) Filter corrected reads to remove those with dual stand mappings or excessive gaps
REF_NAME=$( grep -m 1 '>' $REFPATH | cut -c2- )
FASTA_FILTERED=$WORKSPACE/p${pID}/${pID}_pass-corrected-filtered.fasta
python filter_corrected_reads.py $BAMS_PER_READ/all.bam $FASTA_CORRECTED $FASTA_FILTERED $REF_NAME

deactivate
