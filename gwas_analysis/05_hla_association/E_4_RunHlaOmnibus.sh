#! /bin/bash
# ==============================================================================
#  B_4_RunHlaOmnibus.sh     Written by K.Yamazaki
#   CD: Aug 28 2023
#   - Support for apptainer usage environment at va-server
#------------------------------------------------------------------------------
#  Perform stratified analysis by sex to files exported by DEEP*HLA
#
#  Usage:
#   sbatch -J B_2_ReformatAssoc ./B_2_ReformatAssoc.sh
#------------------------------------------------------------------------------
#  Memo
#
#  ref)
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=16gb
#SBATCH -o logs/%x.%j
# #SBATCH -e logs/%x.%j
# #SBATCH -qos cpu-normal
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
STEP=B_4_RunHlaOmnibus
SRC_DIR=$HOME/analysis/COVID-19
WK_DIR=${SRC_DIR}/2308_impHLA

APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${WK_DIR}/script/sub_${STEP}.R