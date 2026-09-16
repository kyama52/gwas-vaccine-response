#! /bin/bash
# ==============================================================================
#  Concentrate and Convert to plink binary format
#------------------------------------------------------------------------------
#   Procedures:
#   (1) Cleaning and concentration VCF files
#    a. Cleanup dosage files [R2 > '${R2_VAL}']
#    b. Prepare temporary file and directory
#    c. Concentrate VCF files after cleaning
#    d. Annotate rsID using reference panel
#    e. Concentrate info files
#    f. Concentrate log files
#    g. Remove original files [chunkXX.dose.vcf.gz & chunkXX.info]
#    h. Convert to plink binary files
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=4gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
SMPL_NAME=covid-vac
R2_VAL=0.6

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMP=${GENO_DIR}/02_imputation
[ ! -d $GENO_DIR_IMP ] && mkdir -p $GENO_DIR_IMP
FAM_FILE=${GENO_DIR_IMP}/${SMPL_NAME}.fam

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Cleaning and Concentration]"
echo -e "###\t- UTILITY: bcftools"
echo -e "###"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t- GENO_DIR_IMP:"
echo -e "###\t\t>> ${GENO_DIR_IMP/$GENO_DIR/.}"
echo -e "###\t\t- FAM_FILE:"
echo -e "###\t\t>> ${FAM_FILE/$GENO_DIR/.}"
echo -e "###"
echo -e "###\t- THRESHOLD:"
echo -e "###\t\t- CUTOFF: R2 <= ${R2_VAL}"
echo -e "###----------------------------------------------------###\n"

# (1) Cleaning and concentration VCF files
echo -e "(1) Cleaning and concentration VCF files"

for chr in {1..22} X; do
    imp_dir=${GENO_DIR_IMP}/chr${chr}
    each_gname=${SMPL_NAME}_chr${chr}_imputed

    echo -e "\n[chr${chr}]\n"

    if [ ! -f "${imp_dir}/${each_gname}.vcf.gz" ]; then
        sbatch \
            -J "conc_chunk_chr${chr}" \
            -o "${imp_dir}/logs/%x-%j" \
            -e "${imp_dir}/logs/%x-%j" \
            --wrap="apptainer exec ${APPC_DIR}/gwas.sif \
                ${TOOL_DIR}/00_utils/subbatch/200_conc_chunks.sh \
                -i $imp_dir \
                -c $R2_VAL \
                -f $FAM_FILE"
    else
        echo -e "EXISTS imputed file :"
        echo -e "  : ${each_gname}.vcf.gz"
        continue
    fi
done