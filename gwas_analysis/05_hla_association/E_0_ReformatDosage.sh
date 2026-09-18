#! /bin/bash
#------------------------------------------------------------------------------
#  Reformat exported files by DEEP*HLA
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=8gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
SRC_DIR=/path/to/project
STEP=E_0_ReformatDosage

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=/path/to/container

#------------------------------------------------------------------------------
# [Main script ]
apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${TOOL_DIR}/05_hla_association/sub_${STEP}.R