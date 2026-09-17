#! /bin/bash
#------------------------------------------------------------------------------
#  Perform sex-stratified omnibus analysis of HLA amino acid positions
#  [adverse reaction]
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=16gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

set -eu

# [Arguments] (1) Environmental args
SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
STEP=F_2_HlaOmnibus

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=$HOME/tools/container

apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${TOOL_DIR}/05_hla_association/sub_${STEP}.R