MPRA_DIR="/pasteur/helix/projects/evo_immuno_pop/MPRA/"
DATA_DIR="${MPRA_DIR}/data/01_barcode_lib_2/seq1_B7505-4/fastq"

# source /pasteur/helix/projects/evo_immuno_pop/MPRA/MPRAflow/conda/envs/mpraflow_py36/bin/activate
#
# source deactivate

MPRAflow_DIR=${MPRA_DIR}/MPRAflow
module load java/13.0.2

### note 1:
# MPRAflow has been modified to allow running on the pasteur HPC cluster maestro.
# compare ${MPRAflow_DIR}/association_modified.nf and ${MPRAflow_DIR}/association.nf to see the exact changes
# mainly changes include

mkdir ${MPRAflow_DIR}/work/assoc_B7505-4_LIB1_noMapQ_270M_exactMatch
nextflow run -c ${MPRAflow_DIR}/nextflow.config \
  -w ${MPRAflow_DIR}/work/assoc_B7505-4_LIB1_noMapQ_270M_exactMatch/\
 ${MPRAflow_DIR}/association_exactMatch.nf \
 --name assoc_B7505-4_LIB1_noMapQ_270M_exactMatch\
 --fastq-insert "${DATA_DIR}/LIB1_S1_L001_R1_001.fastq.gz"\
 --fastq-insertPE "${DATA_DIR}/LIB1_S1_L001_R3_001.fastq.gz"\
 --fastq-bc "${DATA_DIR}/LIB1_S1_L001_R2_001.fastq.gz"\
 --design "${MPRA_DIR}/data/Final_oligoList_withControls_270nt_160pos_MPRA.fasta"\
 --mapq 0 --cigar '270M' -resume

mkdir ${MPRAflow_DIR}/work/assoc_B7505-4_LIB2_noMapQ_270M_exactMatch
nextflow run -c ${MPRAflow_DIR}/nextflow.config \
  -w ${MPRAflow_DIR}/work/assoc_B7505-4_LIB2_noMapQ_270M_exactMatch \
  ${MPRAflow_DIR}/association_exactMatch.nf \
  --name assoc_B7505-4_LIB2_noMapQ_270M_exactMatch \
  --fastq-insert "${DATA_DIR}/LIB2_S2_L001_R1_001.fastq.gz" \
  --fastq-insertPE "${DATA_DIR}/LIB2_S2_L001_R3_001.fastq.gz" \
  --fastq-bc "${DATA_DIR}/LIB2_S2_L001_R2_001.fastq.gz" \
  --design "${MPRA_DIR}/data/Final_oligoList_withControls_270nt_160pos_MPRA.fasta" \
  --mapq 0 --cigar '270M' -resume
