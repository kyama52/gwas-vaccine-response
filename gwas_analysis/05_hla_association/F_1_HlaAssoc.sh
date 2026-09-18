#! /bin/bash
#------------------------------------------------------------------------------
#  Perform sex-stratified analysis to exported files by DEEP*HLA
#  [adverse reaction]
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH --exclusive
#SBATCH -p cpu

set -eu

# [Arguments] (1) Environmental args
SRC_DIR=/path/to/project
STEP=F_1_HlaAssoc

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=/path/to/container

apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${TOOL_DIR}/05_hla_association/sub_${STEP}.R