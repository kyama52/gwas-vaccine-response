#! /bin/bash
#------------------------------------------------------------------------------
#  Run SNP-by-sex interaction analysis for lead variants
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
GENO_DIR_CAND=${GENO_DIR}/06_cand_vars
[[ ! -d $GENO_DIR_CAND ]] && mkdir -p "$GENO_DIR_CAND"

STEP=G_1_SnpSexInteraction

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=/path/to/container

#------------------------------------------------------------------------------
# [Main script]
atexit() {
    [[ -n $tmpfile ]] && rm -f "$tmpfile"
}
tmpfile=$(mktemp)
trap atexit EXIT
trap 'trap - EXIT; atexit; exit -i' SIGHUP SIGINT SIGTERM

# Export genotype file for samples selection
echo -e "\n- Export genotype file for samples selection"

top_list=("rs2596506" "rs117342275")

cp /dev/null "$tmpfile"
for var in "${top_list[@]}"; do
    echo "$var" >>"$tmpfile"
done

apptainer exec "${APPC_DIR}/gwas.sif" \
    plink2 \
    --pfile "${GENO_DIR_IMPQC}/${SMPL_NAME}_chr6_imputed_qc" \
    --extract "$tmpfile" \
    --make-just-bim \
    --out "${GENO_DIR_CAND}/${SMPL_NAME}_toptag"

apptainer exec "${APPC_DIR}/gwas.sif" \
    plink2 \
    --pfile "${GENO_DIR_IMPQC}/${SMPL_NAME}_chr6_imputed_qc" \
    --extract "$tmpfile" \
    --export A \
    --out "${GENO_DIR_CAND}/${SMPL_NAME}_toptag"

[[ -f "${GENO_DIR_CAND}/${SMPL_NAME}_toptag.raw" ]] &&
    apptainer exec "${APPC_DIR}/gwas.sif" \
        Rscript "${TOOL_DIR}/06_additional_analysis/sub_${STEP}.R"