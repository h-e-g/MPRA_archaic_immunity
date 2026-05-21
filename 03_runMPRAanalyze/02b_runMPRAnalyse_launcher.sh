#!/bin/bash
################################################################################
################################################################################
# File name: 02b_runMPRAnalyse_launcher.sh
# Author: M. Rotival
################################################################################
################################################################################

#SBATCH --job-name=MPRA_emVars_launcher
#SBATCH --output=/pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/logfile_run_MPRA_chunks_%A_%a.log
#SBATCH --array=1-312
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem=4G                # Memory per task
#SBATCH --time=01:00:00           # Time limit
#SBATCH --qos=fast
#SBATCH --partition=common,dedicated,geh
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END

############### @ START SCRIPT
#TODO: adapt for running emVars only or calling a specific subtype (eg. celltype emVars) 

# local ()
DRIVE=/pasteur/helix/projects/evo_immuno_pop/
MPRA_DIR="$DRIVE/MPRA/"
ANALYSIS_DIR="MPRA_count_exp6_analysis3"
SCRIPT_DIR="$MPRA_DIR/scripts/$ANALYSIS_DIR/"
LOG_DIR=${DRIVE}/MPRA/logs/MPRA_emVars_launcher_logs/
OUT_DIR="${MPRA_DIR}/data/${ANALYSIS_DIR}/02b_runMPRAnalyze/"
# call :
# sbatch --array=1-312 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN1_Z3_nBC10 --nBC_min 10 --outlier_threshold 3

# run tests to find the optima0 settings
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z3_nBC10_adjLib --nBC_min 10 --outlier_threshold 3 --adj_lib TRUE -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z3_nBC10_adjRep --nBC_min 10 --outlier_threshold 3 --adj_rep TRUE -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z2_nBC10 --nBC_min 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z3_nBC10 --nBC_min 10 --outlier_threshold 3 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z4_nBC10 --nBC_min 10 --outlier_threshold 4 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z5_nBC10 --nBC_min 10 --outlier_threshold 5 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z7_nBC10 --nBC_min 10 --outlier_threshold 7 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z10_nBC10 --nBC_min 10 --outlier_threshold 10 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN0_Z100_nBC10 --nBC_min 10 --outlier_threshold 100 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt

## run exp4 only
# sbatch --array=1-25 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN1_Z3_nBC10 --nBC_min 10 --outlier_threshold 3 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_exp4.txt
# sbatch --array=38,37,41,40,43,50,54,70,66,65,77,78,116,117,119,120,122,129,133,132,120,144,145,149,157,156 02b_runMPRAnalyse_launcher.sh --output_dir $OUT_DIR/RUN1_Z3_nBC10 --nBC_min 10 --outlier_threshold 3 

# Default parameter values
VERSION_IN="1"
VERSION_OUT="1"
ADJ_LIB=FALSE
ADJ_REP=FALSE
N_BC_MAX=2000
N_BC_MIN=10
OUTLIER_Z_THRESHOLD=3
INSTRUCTION_FILE=${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted.txt
CHUNK_FILE_DIR=${MPRA_DIR}/chunk_files/${ANALYSIS_DIR}/
CHUNK_SIZE=1
# Usage function for help
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --version_in|-vi <value>        Set input version"
    echo "  --version_out|-vo <value>       Set output version"
    echo "  --instruction_file|-i <file>    Set instruction file"
    echo "  --nBC_max|-a <value>            Set max number of barcodes per sample (over this number, barcodes will be subsampled)"
    echo "  --nBC_min|-n <value>            Set min number of barcodes per sample (under this number, barcodes will be tested)"
    echo "  --outlier_threshold|-z <value>  Set outlier Z threshold"
    echo "  --output_dir|-o <dir>           Set output directory"
    echo "  --adj_lib|-l <value>            Adjust for plamsid library and lentivirus preparation (TRUE/FALSE). This is to preferred"
    echo "  --adj_rep|-r <value>            Adjust for replicates (TRUE/FALSE)"
    echo "  --chunk_file_dir|-c <dir>       Set chunk file directory"  
    echo "  --help, -h                      Display this help message"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --version_in|-vi)
            VERSION_IN="$2"
            shift 2
            ;;
        --version_out|-vo)
            VERSION_OUT="$2"
            shift 2
            ;;
        --nBC_max|-a)
            N_BC_MAX="$2"
            shift 2
            ;;
        --nBC_min|-n)
            N_BC_MIN="$2"
            shift 2
            ;;            
        --outlier_threshold|-z)
            OUTLIER_Z_THRESHOLD="$2"
            shift 2
            ;;
        --instruction_file|-i)
            INSTRUCTION_FILE="$2" 
            shift 2
            ;;
        --output_dir|-o)
            OUT_DIR="$2"
            shift 2
            ;;
        --chunk_file_dir|-c)
            CHUNK_FILE_DIR="$2"
            shift 2
            ;;
        --adj_lib|-l)
            ADJ_LIB="$2"
            shift 2
            ;;
        --adj_rep|-r)
            ADJ_REP="$2"
            shift 2
            ;;
        *)  
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
  esac
done

# loop through INSTRUCTION_FILE
# for each line run an array of Rscripts each pointing to a given chunk
# 1 chunk = 1 hour
# ~20-50 chunks per instruction_line


INSTRUCTION_LINE=$(echo $SLURM_ARRAY_TASK_ID | bc)
EXTRACT_LINE=$(echo $INSTRUCTION_LINE+1 | bc)
line=$(head -n $EXTRACT_LINE $INSTRUCTION_FILE | tail -n1)
echo $line
analysis_type=$(echo "$line" | cut -f1)
analysis_subtype=$(echo "$line" | cut -f2)
analysis_name=$(echo "$line" | cut -f3)
group1_samples=$(echo "$line" | cut -f4)
group2_samples=$(echo "$line" | cut -f5)
RUN_ID=$(basename $OUT_DIR)
CHUNK_DIR="${CHUNK_FILE_DIR}/${RUN_ID}/${analysis_type}/${analysis_subtype}/${analysis_name}/"
CHUNK_FILES=($(ls $CHUNK_DIR))

# estimate memory requirement based on number of replicates used
GROUP_COUNT=$(echo $group1_samples,$group2_samples | tr "," "\n" | wc -l)
MEM_REQUIRED=$(echo $GROUP_COUNT*2+15 | bc)G

mkdir -p -m 770 $OUT_DIR
mkdir -p -m 770 ${CHUNK_FILE_DIR}/launcher_chunks/${RUN_ID}/${analysis_type}/${analysis_subtype}/${analysis_name}/

if [ "$analysis_type" == "emVars" ]; then 
  if [ "$group2_samples" == "" ]; then
    PARAM_LIST=("" "--perm_allele 1" "--boot 1 --power 1" "--perm_allele 1 --boot 1 --power 1")
  else
    PARAM_LIST=("" "--perm_group 1" "--boot 1 --power 1" "--perm_group 1 --boot 1 --power 1")
  fi
fi

if [ "$analysis_type" == "oligo_activity" ]; then
  if [ "$group2_samples" == "" ]; then
    PARAM_LIST=("" "--perm_barcode 1" "--boot 1 --power 1" "--perm_barcode 1 --boot 1 --power 1")
  else 
    PARAM_LIST=("" "--perm_group 1" "--boot 1 --power 1" "--perm_group 1 --boot 1 --power 1")
  fi
fi

# list chunk_file & parameters to run

# Initialize combinations array
combinations=()
# Loop over chunks and parameters
for CHUNK in "${CHUNK_FILES[@]}"; do
  for PARAMS in "${PARAM_LIST[@]}"; do
    CHUNK_FILE="${CHUNK_DIR}/${CHUNK}"
    echo "$CHUNK_FILE:$PARAMS"
    combinations+=("$CHUNK_FILE:$PARAMS")
  done
done

# Function to generate SLURM job script for a chunk
generate_chunk_file() {
    local chunk_id=$1
    local chunk_combinations=("${!2}")  # Array of combinations
    local chunk_path="${CHUNK_FILE_DIR}/launcher_chunks/${RUN_ID}/${analysis_type}/${analysis_subtype}/${analysis_name}/chunk_${chunk_id}.txt"
    # Write combinations to a chunk file
    printf "%s\n" "${chunk_combinations[@]}" > "$chunk_path"
}

# Generate chunk files and keep track of chunk IDs
total_combinations=${#combinations[@]}
chunk_id=0
chunk_ids=()
chunk_id_list=''
mkdir -p -m 770 ${CHUNK_FILE_DIR}/launcher_chunks/${RUN_ID}
for ((i = 0; i < total_combinations; i += CHUNK_SIZE)); do
    chunk=("${combinations[@]:i:CHUNK_SIZE}")
    generate_chunk_file "$chunk_id" chunk[@]
    chunk_ids+=("$chunk_id")
    ## commented lines below are used to run only a subset of instructions
    #if [[ ${chunk} == *_1.txt* ]]; then 
    #if [[ ${chunk} == *perm_group* ]]; then
      chunk_id_list=$chunk_id_list,$chunk_id
      echo $chunk
    #fi
    chunk_id=$((chunk_id + 1))
    chunk_count=${#chunk_ids[@]}
done

# Define array script
ARRAY_SCRIPT="${CHUNK_FILE_DIR}/launcher_chunks/${RUN_ID}/${analysis_type}/${analysis_subtype}/${analysis_name}/run_chunks.sh"

# UNCOMMENT TO RUN ONLY 2 CHUNKS 
# chunk_count=2 

# Define array range #SBATCH --array=0-${chunk_count-1} 
cat <<EOL > "$ARRAY_SCRIPT"
#!/bin/bash
#SBATCH --job-name=MPRA_emVars_chunk_array
#SBATCH --output=${LOG_DIR}/MPRA_emVars_chunk_array_%A_chunk_%a.log
#SBATCH --error=${LOG_DIR}/MPRA_emVars_chunk_array_%A_chunk_%a.log
#SBATCH --time=04:00:00  # Adjust time limit as needed
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MEM_REQUIRED}  # Adjust memory per task as needed
#SBATCH --partition=common,depgg # Adjust SLURM partition as needed
#SBATCH --array ${chunk_id_list#,} # Define array range 

SCRIPT_DIR=$SCRIPT_DIR
# Read the chunk file corresponding to the task ID

module purge
module load gcc/9.2.0
module load lapack/3.10.0 
module load R/4.1.0

PARAM_FILE="${CHUNK_FILE_DIR}/launcher_chunks/${RUN_ID}/${analysis_type}/${analysis_subtype}/${analysis_name}/chunk_\${SLURM_ARRAY_TASK_ID}.txt"

while IFS= read -r line; do
    CHUNK_FILE=\$(echo "\$line" | cut -d':' -f1)
    PARAM=\$(echo "\$line" | cut -d':' -f2)
    CHUNK_FILE= echo \$CHUNK_FILE 
    echo "running instruction: 
          Rscript $SCRIPT_DIR/02b_runMPRAanalyse_core.R \
            --version_in $VERSION_IN \
            --version_out $VERSION_OUT \
            --sub_sample $N_BC_MAX \
            --min_nBC $N_BC_MIN \
            --outlier_threshold $OUTLIER_Z_THRESHOLD \
            --instruction_file $INSTRUCTION_FILE \
            --instruction_line $INSTRUCTION_LINE \
            --adj_lib $ADJ_LIB \
            --adj_rep $ADJ_REP \
            --output_dir $OUT_DIR \
            --chunk_file \$CHUNK_FILE\
            \$PARAM" 

      Rscript $SCRIPT_DIR/02b_runMPRAanalyse_core.R \
            --version_in $VERSION_IN \
            --version_out $VERSION_OUT \
            --sub_sample $N_BC_MAX \
            --min_nBC $N_BC_MIN \
            --outlier_threshold $OUTLIER_Z_THRESHOLD \
            --instruction_file $INSTRUCTION_FILE \
            --instruction_line $INSTRUCTION_LINE \
            --adj_lib $ADJ_LIB \
            --adj_rep $ADJ_REP \
            --output_dir $OUT_DIR \
            --chunk_file \$CHUNK_FILE\
            \$PARAM

    echo ""
    echo "${analysis_type}/${analysis_subtype}/${analysis_name} finished running \$(basename \$CHUNK_FILE) with parameters \$PARAM done"
    echo ""
    
done < "\$PARAM_FILE"
EOL

chmod +x "$ARRAY_SCRIPT"

# Submit the array job
JOB_ID=$(sbatch --parsable "$ARRAY_SCRIPT")
echo "Array job "$ARRAY_SCRIPT" submitted "
echo "Job ID: $JOB_ID"
echo " tasks ${chunk_id_list#,}"
echo "Check logs at ${LOG_DIR}/MPRA_emVars_chunk_array_${JOB_ID}_chunk_XXX.log"
echo "parameters in ${CHUNK_FILE_DIR}/launcher_chunks/${RUN_ID}/${analysis_type}/${analysis_subtype}/${analysis_name}/chunk_\${SLURM_ARRAY_TASK_ID}.txt"

echo "All done"