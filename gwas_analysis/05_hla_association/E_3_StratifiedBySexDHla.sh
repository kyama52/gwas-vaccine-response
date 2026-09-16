#! /bin/bash
# ==============================================================================
#  B_3_StratifiedBySexDHla.sh   Written by K.Yamazaki
#   CD: Aug 28 2023
#   - Support for apptainer usage environment at va-server
#------------------------------------------------------------------------------
#  Perform stratified analysis by sex to files exported by DEEP*HLA
#
#  Usage:
#   sbatch -J B_3_StratifiedBySexDHla ./B_3_StratifiedBySexDHla.sh
#------------------------------------------------------------------------------
#  Memo
#
#  ref)
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=4gb
#SBATCH -o logs/%x.%j
# #SBATCH -e logs/%x.%j
# #SBATCH -qos cpu-normal
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
STEP=B_3_StratifiedBySexDHla
SRC_DIR=$HOME/analysis/COVID-19
WK_DIR=${SRC_DIR}/2308_impHLA

APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${WK_DIR}/script/sub_${STEP}.R