#! /bin/bash
#------------------------------------------------------------------------------
#  Perform sex-stratified analysis [antibody titer]
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=16gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
SRC_DIR=/path/to/project
SMPL_NAME=covid-vac

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed

STEP=C_1_AntibodyTiterGwas
WK_DIR=${SRC_DIR}/${STEP}
[[ ! -d ${WK_DIR}/logs ]] && mkdir -p ${WK_DIR}/logs

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
SUB_SH_FILE=${TOOL_DIR}/03_gwas/sub_${STEP}.sh
PREP_R_FILE=${TOOL_DIR}/03_gwas/subprep_${STEP}.R
APPC_DIR=/path/to/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Association analysis]"
echo -e "###\t- SOFTWARE: plink"
echo -e "###"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t- GENO_DIR_IMPQC:"
echo -e "###\t\t>> ${GENO_DIR_IMPQC/$GENO_DIR/.}"
echo -e "###\t- WK_DIR: $WK_DIR"
echo -e "###----------------------------------------------------###\n"

# (0) Prepare parameters
sex_list=(total male female)
array_num=${#sex_list[@]}
apptainer exec ${APPC_DIR}/gwas.sif Rscript ${PREP_R_FILE}

# (1) Perform analysis
if [[ $array_num -gt 0 ]]; then
    sbatch \
        -J "gwas" \
        --mem 32gb \
        -n 4 \
        -a 0-$((array_num - 1)) \
        -o "${WK_DIR}/logs/%x.%A_%a" \
        -e "${WK_DIR}/logs/%x.%A_%a" \
        --wrap="bash ${SUB_SH_FILE}"
fi

echo -e "\n\n### All finished: $(date +%y/%m/%d\ %H:%M:%S)###"