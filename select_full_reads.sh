#!/bin/bash
#$ -P bsg.prjc
#$ -S /bin/bash
#$ -q short.qc
#$ -cwd
# Hannah Roberts, Feb 2018

# (1) Use Porechop to trim fastq
# (2) Map to concatenated reference
# (3) Select reads with primary mapping >= full length

# Requires bwa mem, samtools, python with pysam. 

set -e
set -x

# Arguments to command line -
# (1) minion run id path in form "path_to_run/<MNname>_<date>_<FC>_<MN_ID>_sequencing_run_<runID>_<sample>_<num>"
# (2) workspace, where output directory will be created. Also where the python scripts are assumed to be located.
# (3) path to genotype reference fasta
# (4) patient/sample ID. e.g. 31, 32, 1348RCA
RUN_ID_PATH=$( echo $@ | awk '{print $1}' )
WORKSPACE=$( echo $@ | awk '{print $2}' )
# Reference needs to have all the correct indexes created (e.g. with bwa index)
REFPATH=$( echo $@ | awk '{print $3}' )
pID=$(echo $@ | awk '{print $4}')

RUN_ID=$(basename $RUN_ID_PATH)


BWA_PATH=/apps/well/bwa/0.7.12/bwa
SAMTOOLS_PATH=/apps/well/samtools/1.4.1/bin/samtools
PORECHOP_VENV=/well/ont/apps/porechop/virtualenvs/python3.5.2-gcc5.4.0/porechop-0.2.3/bin/activate
PYSAM_VENV=~/hbv_py/bin/activate # the scripts used in this venv are compatible with python 2.7.11

cd $WORKSPACE
mkdir -p p$pID

# Make concatenated Fastq from those in 'pass' directory. Note assumed minion run directory structure
FASTQ=$WORKSPACE/p${pID}/${pID}_pass.fastq
cat $( ls $RUN_ID_PATH/fastq/pass/* ) > $FASTQ

# Step (1) Trim reads with Porechop.
# Substitute with own porechop virtualenv.
source $PORECHOP_VENV

FASTQ_TRIMMED=$WORKSPACE/p${pID}/${pID}_pass-trimmed.fastq
LOG_TRIMMED=$WORKSPACE/p${pID}/log-trimmed.txt
porechop -i $FASTQ -o $FASTQ_TRIMMED --format fastq --verbosity 1 -t 1 > $LOG_TRIMMED

deactivate

# Step (2) Map with bwa mem
BAM_FILE_PATH=${WORKSPACE}/p${pID}/${pID}_pass-trimmed.bam

$BWA_PATH mem -x ont2d -R "@RG\tID:1\tSM:HBV" -M \
  $REFPATH \
  $FASTQ_TRIMMED \
  | $SAMTOOLS_PATH view -bS - \
  | $SAMTOOLS_PATH sort -o $BAM_FILE_PATH -

$SAMTOOLS_PATH index $BAM_FILE_PATH

# Step (3) Use pysam to detect reads with full length mappings
# Extract these from the trimmed fastq and make a new fastq with just these reads
source $PYSAM_VENV

BAM_PRIMARY=${WORKSPACE}/p${pID}/${pID}_pass-trimmed-primary.bam
$SAMTOOLS_PATH view -b $BAM_FILE_PATH -F 0x904 > $BAM_PRIMARY
$SAMTOOLS_PATH index $BAM_PRIMARY
FULL_LEN_READS=${WORKSPACE}/p${pID}/full_length_reads_3200.txt
FASTQ_FILTERED=${WORKSPACE}/p${pID}/${pID}_pass-trimmed-3200.fastq

python select_full_reads.py $BAM_PRIMARY $FULL_LEN_READS
python filter_fastq.py $FULL_LEN_READS $FASTQ_TRIMMED $FASTQ_FILTERED

deactivate

