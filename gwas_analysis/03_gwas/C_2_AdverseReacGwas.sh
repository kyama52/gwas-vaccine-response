#! /bin/bash
#------------------------------------------------------------------------------
#  Perform sex-stratified analysis [adverse reaction]
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=8gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
# #SBATCH -p lm
# #SBATCH -qos cpu-normal
#SBATCH -p cpu
#SBATCH --qos normal

set -eu

# [Arguments] (2) Environmental args
SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
SMPL_NAME=covid-vac

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed

STEP=C_2_AdverseReacGwas
WK_DIR=${SRC_DIR}/${STEP}
[[ ! -d ${WK_DIR}/logs ]] && mkdir -p ${WK_DIR}/logs

PHENO_LIST=("seff_fever" "seff_arthralgia" "seff_fatigue" "seff_cold" "seff_headache" "seff_muscles")
# PHENO_LIST=("seff_fever")

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
SUB_SH_FILE=${TOOL_DIR}/03_gwas/sub_${STEP}.sh

#------------------------------------------------------------------------------
# [Main script ]
# (0) Prepare condition for analysis
# a. Start logging
echo -e "###########################################################"
echo -e "### [Association analysis]"
echo -e "###\t- SOFTWARE: OrdinalGWAS"
echo -e "###"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t- GENO_DIR_IMPQC:"
echo -e "###\t\t>> ${GENO_DIR_IMPQC/$GENO_DIR/.}"
echo -e "###"
echo -e "###\t- WK_DIR: $WK_DIR"
echo -e "###----------------------------------------------------###\n"

# (1) Perform OrdinalGWAS using by sub-routine script
echo -e "- Perform OrdinalGWAS using by sub-routine script\n"
array_num=${#PHENO_LIST[@]}

if [[ $array_num -gt 0 ]]; then
    sbatch \
        -J "ordinalGWAS" \
        --mem 32gb \
        -n 16 \
        -a 0-$((array_num - 1)) \
        -o "${WK_DIR}/logs/%x.%A_%a" \
        -e "${WK_DIR}/logs/%x.%A_%a" \
        "$SUB_SH_FILE"
fi