#! /bin/bash
#------------------------------------------------------------------------------
#  Perform conditional analysis based on HLA imputation association results
#  [antibody titer]
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=100G
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH --nodes=1               # -N, --node
#SBATCH --ntasks=1              # -n, --ntasks
#SBATCH --cpus-per-task=28
#SBATCH --exclusive
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
SRC_DIR=/path/to/project
SMPL_NAME=covid-vac

CHR=6
HLA_STR=24
HLA_END=36

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed
GENO_DIR_HLA=${GENO_DIR}/04_imputation_hla
MHC_GNAME=${GENO_DIR_HLA}/${SMPL_NAME}_chr${CHR}_hla_imputed

STEP=E_3_ConditionedOmnibus

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=/path/to/container

#------------------------------------------------------------------------------
# [Main script ]
# (1) Prepare imputed genotype data
if [[ ! -f ${MHC_GNAME}.bim.gz ]]; then
    echo -e "\n- Prepare BIM file around the HLA region\n"
    plink1.9 \
        --bfile "${GENO_DIR_IMPQC}/${SMPL_NAME}_chr${CHR}_imputed_qc" \
        --chr ${CHR} \
        --from-mb $HLA_STR --to-mb $HLA_END \
        --make-just-bim \
        --out "${MHC_GNAME}"
    gzip -f ${MHC_GNAME}.bim
    echo -e "\n>> ${MHC_GNAME/${GENO_DIR}/.}.bim\n"
fi

if [[ ! -f ${MHC_GNAME}.raw.gz ]]; then
    echo -e "\n- Prepare imputed genotype file around the HLA region\n"
    plink1.9 \
        --bfile "${GENO_DIR_IMPQC}/${SMPL_NAME}_chr${CHR}_imputed_qc" \
        --chr ${CHR} \
        --from-mb $HLA_STR --to-mb $HLA_END \
        --recode A \
        --out "${MHC_GNAME}"
    gzip -f ${MHC_GNAME}.raw
    echo -e "\n>> ${MHC_GNAME/${GENO_DIR}/.}.raw.gz\n"
fi

# (2) Condition analysis
echo -e "\n- Condition analysis\n"
apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${TOOL_DIR}/05_hla_association/sub_${STEP}.R

echo -e "\nFinished $(date +"%Y/%m/%d %H:%M")\n"