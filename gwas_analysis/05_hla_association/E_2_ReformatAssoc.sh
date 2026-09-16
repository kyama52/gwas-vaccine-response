#! /bin/bash
# ==============================================================================
#  B_2_ReformatAssoc.sh     Written by K.Yamazaki
#   CD: Aug 28 2023
#   - Support for apptainer usage environment at va-server
#------------------------------------------------------------------------------
#  Reformat files exported by DEEP*HLA and examine linear logistic regression
#
#  Usage:
#   sbatch -J B_2_ReformatAssoc ./B_2_ReformatAssoc.sh
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
STEP=B_2_ReformatAssoc
SRC_DIR=$HOME/analysis/COVID-19
WK_DIR=${SRC_DIR}/2308_impHLA

APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${WK_DIR}/script/sub_${STEP}.R