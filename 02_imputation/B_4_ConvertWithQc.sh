#! /bin/bash
#------------------------------------------------------------------------------
#  Convert imputed genotypes with QC
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
SRC_DIR=/path/to/project
SMPL_NAME=covid-vac
R2_VAL=0.6
MAF_VAL=0.05

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMP=${GENO_DIR}/02_imputation
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed
[ ! -d ${GENO_DIR_IMPQC}/logs ] && mkdir -p ${GENO_DIR_IMPQC}/logs
FAM_FILE=${GENO_DIR_IMP}/${SMPL_NAME}.fam

ETC_DIR=${SRC_DIR}/etc
PHENO_FILE=${ETC_DIR}/${SMPL_NAME}_pheno.txt
THRES_VAL_1=bef_titer
THRES_VAL_2=cov19_phx
RMV_SFILE=${ETC_DIR}/rmv_smpl_expvirus.txt

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
APPC_DIR=/path/to/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Cleaning and Concentration]"
echo -e "###\t- SRC_DIR: ${SRC_DIR}"
echo -e "###"
echo -e "###\t\t- PHENO_FILE:"
echo -e "###\t\t>> ${PHENO_FILE/$SRC_DIR/.}"
echo -e "###\t\t\t- ${THRES_VAL_1} > 0 && ${THRES_VAL_2} == 0"
echo -e "###\t\t\t>> ${RMV_SFILE/$SRC_DIR/.}"
echo -e "###"
echo -e "###\t- GENO_DIR: ${GENO_DIR/$SRC_DIR/.}"
echo -e "###\t\t- GENO_DIR_IMP:"
echo -e "###\t\t>> ${GENO_DIR_IMP/$GENO_DIR/.}"
echo -e "###\t\t- FAM_FILE:"
echo -e "###\t\t>> ${FAM_FILE/$GENO_DIR/.}"
echo -e "###----------------------------------------------------###\n"

# (0) Identify samples with prior SARS-CoV-2 exposure
if [ ! -f "$RMV_SFILE" ]; then
    echo "- Check samples exposured SARS-CoV-2 virus"
    if [ -f "$PHENO_FILE" ]; then
        echo -e "<< ${PHENO_FILE/$SRC_DIR/.}"

        awk -F "\t" -v cname_bef=$THRES_VAL_1 -v cname_cov19=$THRES_VAL_2 \
            ' BEGIN { OFS = "\t"; idx_bef = 0; idx_cov19 = 0 }
            {
                if(NR == 1) {
                    for (i=1; i<=NF; i++) {
                        if($i == cname_bef)  idx_bef = i
                        else if($i == cname_cov19)  idx_cov19 = i
                    }
                    print $1, $2
                }
                else {
                    if($idx_bef > 0 || $idx_cov19 == 1) print $1, $2
                }
            }' "$PHENO_FILE" | uniq >"$RMV_SFILE"
        echo -e ">> ${RMV_SFILE/$SRC_DIR/.}\n"

    else
        echo ">> There is NO PHENO_FILE: !"
        echo -e "${PHENO_FILE/$SRC_DIR/.}"
        exit 1
    fi
fi

# (1) Apply QC and convert to genotype data
echo -e "\n###----------------------------------------------------###"
echo -e "###\t- THRESHOLD:"
echo -e "###\t\t- R2 > ${R2_VAL}"
echo -e "###\t\t- MAF > ${MAF_VAL}"
echo -e "###----------------------------------------------------###"

for chr in {1..22} X; do
    imp_dir=${GENO_DIR_IMP}/chr${chr}
    each_gname=${SMPL_NAME}_chr${chr}_imputed

    if [ -f "${imp_dir}/${each_gname}.vcf.gz" ]; then
        # Support for apptainer usage environment
        sbatch \
            -J "cnv_chr${chr}" \
            -o "${GENO_DIR_IMPQC}/logs/%x.%A_%a" \
            -e "${GENO_DIR_IMPQC}/logs/%x.%A_%a" \
            --wrap="apptainer exec ${APPC_DIR}/gwas.sif \
                    bash ${TOOL_DIR}/00_utils/subbatch/001_conv_genofiles.sh \
                    -i ${imp_dir}/${each_gname}.vcf.gz \
                    -o $GENO_DIR_IMPQC \
                    -s $RMV_SFILE \
                    -f $FAM_FILE \
                    -m $MAF_VAL \
                    -c $R2_VAL"

    else
        echo -e "NOT EXISTS imputed file :"
        echo -e "  : ${each_gname}.vcf.gz"
        continue
    fi
done