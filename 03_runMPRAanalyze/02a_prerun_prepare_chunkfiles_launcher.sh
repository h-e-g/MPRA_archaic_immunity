#!/bin/bash
################################################################################
################################################################################
# File name: 02a_prepare_chunkfiles_launcher.sh
# Author: M. Rotival
################################################################################
################################################################################

#SBATCH --job-name=prepare_chunk_emVars
#SBATCH --output=/pasteur/helix/projects/evo_immuno_pop/MPRA/logs/%x_logs/logfile_prepare_chunks_%A_%a.log
#SBATCH --array=1-313
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1         # Number of CPUs per task
#SBATCH --mem=100G                # Memory per task
#SBATCH --time=02:00:00           # Time limit
#SBATCH --qos=fast
#SBATCH --partition=common,dedicated
#SBATCH --mail-user=mrotival@pasteur.fr
#SBATCH --mail-type=END

# local ()
DRIVE=/pasteur/helix/projects/evo_immuno_pop/
ANALYSIS_DIR="MPRA_count_exp6_analysis3"
MPRA_DIR="$DRIVE/MPRA/"
SCRIPT_DIR="$MPRA_DIR/scripts/$ANALYSIS_DIR/"

# call :
# sbatch --array=1-313 02a_prerun_prepare_chunkfiles_launcher.sh 

# run tests to find the optimal settings
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z3_nBC10_adjLib --min_nBC 10 --outlier_threshold 3 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z3_nBC10_adjRep --min_nBC 10 --outlier_threshold 3 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z2_nBC10 --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z3_nBC10 --min_nBC 10 --outlier_threshold 3 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z4_nBC10 --min_nBC 10 --outlier_threshold 4 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z5_nBC10 --min_nBC 10 --outlier_threshold 5 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z7_nBC10 --min_nBC 10 --outlier_threshold 7 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z10_nBC10 --min_nBC 10 --outlier_threshold 10 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-30 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN0_Z100_nBC10 --min_nBC 10 --outlier_threshold 100 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_tests.txt
# sbatch --array=1-25 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN1_Z3_nBC10 --min_nBC 10 --outlier_threshold 3  -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_exp4.txt

#### run with optimal settings 
# step1 - oligo activity - Quanti only (no adj)
# sbatch --array=1-15,36-77 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN3_Z2_nBC10 --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted_noA549.txt

# step2 - oligo activity - Diff only (lib adj)
# sbatch --array=16-35 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN3_Z2_nBC10 --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted_noA549.txt

# step3 - emVars - Quanti_only (rep adj)
# sbatch --array=78-92,113-154 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN3_Z2_nBC10 --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted_noA549.txt

# step4 - emVars - Diff_only (lib adj)
# sbatch --array=93-112 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN3_Z2_nBC10 --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted_noA549.txt

# step5 - oligo activity - Diff only (no adj)
# sbatch --array=16-35 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN2_Z2_nBC10_noLibAdj --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted_noA549.txt

# step6 - emVars - Diff_only (no adj)
# sbatch --array=93-112 02a_prerun_prepare_chunkfiles_launcher.sh --output_dir RUN2_Z2_nBC10_noLibAdj --min_nBC 10 --outlier_threshold 2 -i ${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted_noA549.txt


# Default parameter values
VERSION_IN=1
VERSION_OUT=1
instruction_file="${SCRIPT_DIR}/02a_MPRA_analyse_instructions_full_sorted.txt"
sub_sample=2000
min_nBC=10 # not used so far
Outlier_Z_threshold=3
target_time_s=3600 # target of 1 hour per job (set limit between 2 and 4 h)
cst_time=6 # values based on empirical test + some margin
linear_complexity=0.004 # values based on empirical test + some margin
squared_complexity=0 # values based on empirical test + some margin

usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  --version_in <version>         Set input version"
    echo "  --version_out <version>        Set output version"
    echo "  --instruction_file <file>      Set instruction file"
    echo "  --subsample <number>           Set sub sample size"
    echo "  --min_nBC <number>             Set minimum number of barcode"
    echo "  --outlier_threshold <value>    Set outlier Z threshold"
    echo "  --target_time_s <seconds>      Set target time in seconds"
    echo "  --cst_time <value>             Set constant time"
    echo "  --linear_complexity <value>    Set linear complexity"
    echo "  --squared_complexity <value>   Set squared complexity"
    echo "  --output_dir <value>           Set subdirectory where chunk files should be written"
    exit 1
}

# Parse command line arguments
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
        --instruction_file|-i)
            instruction_file="$2"
            shift 2
            ;;
        --subsample|-u)
            sub_sample="$2"
            shift 2
            ;;
        --min_nBC|-n)
            min_nBC="$2"
            shift 2
            ;;
        --outlier_threshold|-z)
            Outlier_Z_threshold="$2"
            shift 2
            ;;
        --target_time_s| -s)
            target_time_s="$2"
            shift 2
            ;;
        --cst_time| -c)
            cst_time="$2"
            shift 2
            ;;
        --linear_complexity | -lc)
            linear_complexity="$2"
            shift 2
            ;;
        --squared_complexity| -sc)
            squared_complexity="$2"
            shift 2
            ;;
        --output_dir|-o)
            output_directory="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done


# Display the configuration
echo "Configuration:"
echo "Version in: $VERSION_IN"
echo "Version out: $VERSION_OUT"
echo "Instruction file: $instruction_file"
echo "Sub sample size: $sub_sample"
echo "Minimum number of barcode: $min_nBC"
echo "Outlier Z threshold: $Outlier_Z_threshold"
echo "Target time in seconds: $target_time_s"
echo "Constant time: $cst_time"
echo "Linear complexity: $linear_complexity"
echo "Squared complexity: $squared_complexity"
echo "output directory: $output_directory"

module purge
module load gcc/9.2.0
module load lapack/3.10.0 
module load R/4.1.0

 echo "running instruction: 
      Rscript $SCRIPT_DIR/02a_prerun_prepare_chunkfiles_core.R \\
        --version_in $VERSION_IN \\
        --version_out $VERSION_OUT \\
        --instruction_file $instruction_file \\
        --sub_sample $sub_sample \\
        --min_nBC $min_nBC \\
        --outlier_threshold $Outlier_Z_threshold \\
        --target_time_s $target_time_s \\
        --cst_time $cst_time \\
        --linear_complexity $linear_complexity \\
        --squared_complexity $squared_complexity \\
        --output_dir $output_directory \\        
        --instruction_line $SLURM_ARRAY_TASK_ID"

    Rscript $SCRIPT_DIR/02a_prerun_prepare_chunkfiles_core.R \
        --version_in $VERSION_IN \
        --version_out $VERSION_OUT \
        --instruction_file $instruction_file \
        --sub_sample $sub_sample \
        --min_nBC $min_nBC \
        --outlier_threshold $Outlier_Z_threshold \
        --target_time_s $target_time_s \
        --cst_time $cst_time \
        --linear_complexity $linear_complexity \
        --squared_complexity $squared_complexity \
        --output_dir $output_directory \
        --instruction_line $SLURM_ARRAY_TASK_ID

echo "All done!"