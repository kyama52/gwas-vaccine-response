#! /bin/bash
# ==============================================================================
#  B_5_ConditionedOmnibus.sh
#   CD: Apr 03 2025     Written by K.Yamazaki
#   - For multiprocessiong on R, add option --cpus-per-task=28 --mem=100G
#   UD: Apr 07 2025     Bug fix
#------------------------------------------------------------------------------
#  Condition analysis based on association results of HLA imputation
#
#  Usage:
#   sbatch -J B_5_ConditionedOmnibus ./B_5_ConditionedOmnibus.sh
#------------------------------------------------------------------------------
#  Memo
#
#  ref)
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
# #SBATCH --mem=4gb
#SBATCH --mem=100G
#SBATCH -o logs/%x.%j
# #SBATCH -e logs/%x.%j
# #SBATCH -qos cpu-normal
#SBATCH --nodes=1               # -N, --node
#SBATCH --ntasks=1              # -n, --ntasks
#SBATCH --cpus-per-task=28
#SBATCH --exclusive
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
SRC_DIR=$HOME/analysis/COVID-19
SMPL_NAME=CUH-GWAS_230215
GENO_DNAME=PLINK_160223_1213

CHR=6
HLA_STR=24
HLA_END=36

GENO_DIR=${SRC_DIR}/genotype/${GENO_DNAME}
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed/qc_2
GENO_DIR_HLA=${GENO_DIR}/04_imputation_hla
MHC_GNAME=${GENO_DIR_HLA}/${SMPL_NAME}_chr${CHR}_mhc_imputed

WK_DIR=${SRC_DIR}/2308_impHLA

STEP=B_5_ConditionedOmnibus

APPC_DIR=$HOME/tools/container

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
apptainer exec ${APPC_DIR}/gwas.sif Rscript ${WK_DIR}/script/sub_${STEP}.R

echo -e "\nFinished $(date +"%Y/%m/%d %H:%M")\n"