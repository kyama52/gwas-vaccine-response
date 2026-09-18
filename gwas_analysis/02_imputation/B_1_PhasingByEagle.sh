#! /bin/bash
#------------------------------------------------------------------------------
#  Phasing by Eagle
#------------------------------------------------------------------------------
#   Procedures:
#   (1) Extract common variants to references
#   a. Extract common variants between BIM and reference
#   b. Update alleles by exported files
#   c. Update ID by exported files
#   (2) Creating the input files for phasing and imputation.
#   a. Extract genotypes of overlapped SNPs
#   b. Split genotypes by chromosome and phasing
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=4gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

# [Arguments] (2) Environmental args
set -eu

RES_DIR=$HOME/resource
BUILD=hg19
FASTA_FILE=${RES_DIR}/UCSC/${BUILD}/${BUILD}.fa
GMAP_FILE=${RES_DIR}/Eagle/genetic_map_${BUILD}_withX.txt.gz
RES_DIR_EGL=${RES_DIR}/Eagle/BBJ1K
RES_DIR_MMC=${RES_DIR}/minimac/BBJ1K
RES_FNAME=BBJ1K_1KGP_RefPanel_b155

SRC_DIR=/path/to/project
SMPL_NAME=covid-vac

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_QC=${GENO_DIR}/01_qc
GENO_DIR_IMP=${GENO_DIR}/02_imputation
[ ! -d $GENO_DIR_IMP ] && mkdir -p $GENO_DIR_IMP

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
SUB_ROUTINE=${TOOL_DIR}/00_utils/subbatch/103_extract_cmn_vars.R

X_PAR_FILE=${GENO_DIR_IMP}/chrX_par.txt
if [ ! -f $X_PAR_FILE ]; then
    echo -e "X\t60001\t2699520\tPAR1" >${X_PAR_FILE}
    echo -e "X\t154931044\t155260560\tPAR2" >>${X_PAR_FILE}
fi

UPD_ID_FILE=${GENO_DIR_IMP}/snps-upd-id.txt
UPD_AL_FILE=${GENO_DIR_IMP}/snps-upd-allele.txt
CMN_SNP_FILE=${GENO_DIR_IMP}/snps-reference-dup.txt

APPC_DIR=/path/to/container

atexit() {
    [[ -n $tmpfile ]] && rm -f "$tmpfile"
}
tmpfile=$(mktemp)
trap atexit EXIT
trap 'trap - EXIT; atexit; exit -i' SIGHUP SIGINT SIGTERM

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Phasing and imputation]"
echo -e "###\t- PHASING: EAGLE_DIR"
echo -e "###"
echo -e "###\t- RES_DIR: ${RES_DIR}"
echo -e "###\t\t- EAGLE DIR: ${RES_DIR_EGL/$RES_DIR/.}"
echo -e "###\t\t- MINIMAC DIR: ${RES_DIR_MMC/$RES_DIR/.}"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t- APPLIED QC FILES:"
echo -e "###\t\t>> ${GENO_DIR_QC/$GENO_DIR/.}/${SMPL_NAME}.(bed|bim|fam)"
echo -e "###\t\t- GENO_DIR_IMP:"
echo -e "###\t\t>> ${GENO_DIR_IMP/$GENO_DIR/.}"
echo -e "###----------------------------------------------------###\n"

# (1) Extract common variants to references
echo -e "(1) Extract common variants to references"

# a. "Convert X chromosome from [23] to [X]"
# ref)  https://orebibou.com/ja/home/201601/20160114_001/ -> awk '{}1'
TGT_FNAME=${GENO_DIR_QC}/${SMPL_NAME}
if [ ! -e ${GENO_DIR_QC}.bimbak ]; then
    echo -e "- Convert X chromosome from [23] to [X]"
    echo -e "\t- Backup original BIM file:"
    echo -e "\t>> ${TGT_FNAME/${GENO_DIR}/.}.bimbak"
    mv ${TGT_FNAME}.bim ${TGT_FNAME}.bimbak
    awk -F "\t" \
        ' BEGIN {OFS = "\t"}  $1 == "23" { $1 = "X" }1' ${TGT_FNAME}.bimbak >${TGT_FNAME}.bim
    echo -e "\t- Converted BIM file:"
    echo -e "\t>> ${TGT_FNAME/${GENO_DIR}/.}.bim"
fi

# b. Strand flips by snpflip and PLINK
echo -e "\n- Strand flips by snpflip and PLINK"
[ ! -d ${GENO_DIR_IMP}/etc ] && mkdir -p ${GENO_DIR_IMP}/etc
SFLIP_FNAME=${GENO_DIR_IMP}/etc/${SMPL_NAME}
snpflip --fasta-genome=${FASTA_FILE} \
    --bim-file=${TGT_FNAME}.bim \
    -o ${SFLIP_FNAME}
echo -e ">> ${SFLIP_FNAME/${GENO_DIR}/.}.ambiguous"
echo -e ">> ${SFLIP_FNAME/${GENO_DIR}/.}.reverse\n"

plink1.9 \
    --bfile ${TGT_FNAME} \
    --exclude ${SFLIP_FNAME}.ambiguous \
    --flip ${SFLIP_FNAME}.reverse \
    --allow-no-sex \
    --make-bed \
    --out ${GENO_DIR_IMP}/${SMPL_NAME}_sflip
echo -e "\n>> ${GENO_DIR_IMP/${GENO_DIR}/.}/${SMPL_NAME}_sflip.(bed|bim|fam)\n"

# c. Extract common variants between BIM and reference
echo -e "- Using by sub-script [sub_C_extract_cmn_vars.R]"
apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${SUB_ROUTINE} \
    ${GENO_DIR_IMP}/${SMPL_NAME}_sflip.bim ${GENO_DIR_IMP}

# d. Update alleles by exported files
TGT_FNAME=${GENO_DIR_IMP}/${SMPL_NAME}
if [ -f $UPD_AL_FILE ]; then
    echo -e "\n\n- Update alleles by exported files"
    echo -e "\t>> ${UPD_AL_FILE/$GENO_DIR_IMP/.}\n"

    plink1.9 \
        --allow-no-sex \
        --bfile ${TGT_FNAME}_sflip \
        --update-alleles ${UPD_AL_FILE} \
        --make-bed \
        --out ${TGT_FNAME}_upd_al
fi

# e. Update ID by exported files
if [ -f $UPD_ID_FILE ]; then
    echo -e "\n\n- Update ID by exported files"
    echo -e "\t>> ${UPD_ID_FILE/$GENO_DIR_IMP/.}\n"

    if [ -f ${TGT_FNAME}_upd_al.bed ]; then
        src_fname=${TGT_FNAME}_upd_al
    else
        src_fname=${TGT_FNAME}_sflip
    fi
    plink1.9 \
        --allow-no-sex \
        --bfile ${src_fname} \
        --update-name ${UPD_ID_FILE} \
        --make-bed \
        --out ${TGT_FNAME}_upd_id
fi

# (2) Creating the input files for phasing and imputation.
# a. Extract genotypes of overlapped SNPs
echo -e "(2) Creating the input files for phasing and imputation"
echo -e "- Extract genotypes of overlapped SNPs\n"

if [ ! -f "${TGT_FNAME}.bed" ]; then
    if [ -f ${TGT_FNAME}_upd_id.bed ]; then
        src_fname=${TGT_FNAME}_upd_id
    else
        src_fname=${TGT_FNAME}_sflip
    fi
    plink1.9 \
        --allow-no-sex \
        --bfile ${src_fname} \
        --extract ${CMN_SNP_FILE} \
        --exclude range ${X_PAR_FILE} \
        --make-bed \
        --out ${TGT_FNAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
else
    echo -e "<< ${TGT_FNAME/$GENO_DIR/.}.bed\n"
fi

# b. Split genotypes by chromosome and phasing
echo -e "- Split genotypes by chromosome and phasing"

for chr in {1..22} X; do
    echo -e "\n[chr${chr}]\n"
    imp_dir=${GENO_DIR_IMP}/chr${chr}
    [ ! -d ${imp_dir}/logs ] && mkdir -p ${imp_dir}/logs

    each_geno_file=${imp_dir}/${SMPL_NAME}_chr${chr}
    each_ref_file=${RES_DIR_EGL}/${RES_FNAME}_chr${chr}.bcf

    # [UD: Apr 24 2023: Add --snps-only just-acgt]
    if [ ${chr} != "X" ]; then
        jobid=$(sbatch \
            --parsable \
            -J "split_chr${chr}" \
            -o "${imp_dir}/logs/%x-%j" \
            -e "${imp_dir}/logs/%x-%j" \
            --wrap="apptainer exec ${APPC_DIR}/gwas.sif bash -c ' \
                plink1.9 \
                    --allow-no-sex \
                    --snps-only just-acgt \
                    --bfile ${TGT_FNAME} \
                    --chr $chr \
                    --recode vcf-iid \
                    --out ${each_geno_file} && \
                bgzip -f ${each_geno_file}.vcf && \
                bcftools index -f ${each_geno_file}.vcf.gz'")
    else
        echo -e "23\tX" >${GENO_DIR_IMP}/nonauto_chr_conv.txt
        jobid=$(sbatch \
            --parsable \
            -J "split_chr${chr}" \
            -o "${imp_dir}/logs/%x-%j" \
            -e "${imp_dir}/logs/%x-%j" \
            --wrap="apptainer exec ${APPC_DIR}/gwas.sif bash -c ' \
                plink1.9 \
                    --allow-no-sex \
                    --snps-only just-acgt \
                    --bfile ${TGT_FNAME} \
                    --chr $chr \
                    --recode vcf-iid \
                    --out ${each_geno_file} && \
                bgzip -f ${each_geno_file}.vcf && \
                bcftools annotate \
                    --rename-chrs ${GENO_DIR_IMP}/nonauto_chr_conv.txt \
                    ${each_geno_file}.vcf.gz \
                    -O z -o $tmpfile && \
                mv $tmpfile ${each_geno_file}.vcf.gz && \
                bcftools index -f ${each_geno_file}.vcf.gz'")
    fi

    sbatch \
        -d afterany:$jobid \
        --parsable \
        -J "eagle_chr${chr}" \
        -o "${imp_dir}/logs/%x-%j" \
        -e "${imp_dir}/logs/%x-%j" \
        -c 4 \
        --mem 8G \
        --qos cpu-normal \
        --wrap="apptainer exec ${APPC_DIR}/gwas.sif eagle \
            --vcfTarget ${each_geno_file}.vcf.gz \
            --numThreads 4 \
            --geneticMapFile=${GMAP_FILE} \
            --vcfRef ${each_ref_file} \
            --vcfOutFormat z \
            --allowRefAltSwap \
            --outPrefix ${each_geno_file}_phased"
done