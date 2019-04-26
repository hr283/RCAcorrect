# RCAcorrect

The starting point for this pipeline is fastqs containing basecalled ONT reads.
It is assumed that the genotype of the sample has already been determined and a concatenated (two-copy) reference sequence has been created for that genotype
as well as a single-copy reference sequence.

First, run `select_full_reads.sh `\
This looks for fastqs in $RUN_ID_PATH/fastq/pass, concatenates them into a single fastq and then trims adapters using porechop. 
Trimmed reads are then mapped with bwa-mem, and pysam is used to select reads with full-length mappings (of 3.2kb or more).
The 3.2kb threshold can be changed by modifying the value of 'hbv_length' in select_full_reads.py

Inputs are:\
`./select_full_reads.sh $RUN_ID_PATH $OUTPUT_DIRECTORY $REFERENCE_FASTA $SAMPLE_ID`

where the reference sequence should be two concatenated copies of the relevant genotype reference sequence.

Then, run `chop_and_correct_RCA.sh` \
The steps carried out in this script are as follows:

(1) Select 100bp anchor at start of reference genome \
(2) For each full length read: Find this anchor within the read using blast \
(3)   Chop up read based on anchor locations -> create fastq for this read \
(4)   Align these reads to (single) reference genome \
(5)   Correct sequences and write corrected consensus to fasta file \
(6) Filter corrected reads to remove those with dual stand mappings or excessive gaps \

Inputs are:\
`./chop_and_correct_RCA.sh $OUTPUT_DIRECTORY $REFERENCE_FASTA $SAMPLE_ID $CORRECTION_TYPE`

The output directory and sample ID should be the same as that given to select_full_reads.sh. The reference should in this case be a single copy genotype reference.
Correction type should be one of consensus_only or error_correct.

Paths to relevant software are found near the top of each bash script and need to be appropriately modified before use.
The scripts also use two python virtualenvs, `$PORECHOP_VENV`, in which porechop 0.2.3 is installed, and `$PYSAM_VENV`, the requirements for which can be found in `pysam_venv_requirements.txt`.

