#! /bin/bash
# ==============================================================================
#  Run HLA-imputation by DEEP*HLA
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=8gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -c 4
#SBATCH -p cpu

set -eu

#------------------------------------------------------------------------------
# [Arguments] (2) Environmental args
SRC_DIR=$HOME/analysis/COVID-19
SMPL_NAME=CUH-GWAS_230215
GENO_DNAME=PLINK_160223_1213

DHLA_DIR=$HOME/local/DEEP-HLA

GENO_DIR=${SRC_DIR}/genotype/${GENO_DNAME}
GENO_DIR_HLA=${GENO_DIR}/04_imputation_hla
PHASED_FNAME=${GENO_DIR_HLA}/${SMPL_NAME}_mhc_phased
IMPUTED_FNAME=${PHASED_FNAME/_phased/_imputed}

today=$(date +%y%m%d)

# (3) Prepare venv for DEEP*HLA
VENV_DIR=${HOME}/.venvs
VNAME=DEEP-HLA
if [ ! -d ${VENV_DIR} ]; then
    [ ! -d $VENV_DIR ] && mkdir $VENV_DIR
    if [ ! -d ${VENV_DIR}/${VNAME} ]; then
        bash ${SRC_DIR}/2308_impHLA/script/Z_2_ContructVenvDeepHLA.sh
    fi
fi

# atexit() {
#     [[ -n $tmpfile ]] && rm -f "$tmpfile"
# }
# tmpfile=$(mktemp)
# trap atexit EXIT
# trap 'trap - EXIT; atexit; exit -i' SIGHUP SIGINT SIGTERM

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [HLA imputation analysis]"
echo -e "###\t- SOFTWARE: DEEP*HLA"
echo -e "###\t\t>> ${DHLA_DIR}"
echo -e "###"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t>> ${GENO_DIR_HLA/$GENO_DIR/.}"
echo -e "###----------------------------------------------------###\n"

# (1) HLA imputation using by DEEP-HLA
echo -e "\n----------------------------------------"
echo -e "- HLA imputation using by DEEP*HLA"
echo -e "----------------------------------------\n"

# a. Check files for DEEP-HLA
echo -e "- Check files for DEEP-HLA"
if [ ! -f ${PHASED_FNAME}.haps ] || [ ! -f ${PHASED_FNAME}.sample ]; then
    [ ! -f ${PHASED_FNAME}.haps ] &&
        echo -e "\t- NOT Exist: ${PHASED_FNAME/$GENO_DIR/.}.haps"
    [ ! -f ${PHASED_FNAME}.sample ] &&
        echo -e "\t- NOT Exist: ${PHASED_FNAME/$GENO_DIR/.}.sample"
    return 1
fi

if [ ! -f ${PHASED_FNAME}.bim ]; then
    echo -e "\t- NOT exist BIM_FILE !"

    if [ -f ${PHASED_FNAME/_phased/}.bim ]; then
        cp ${PHASED_FNAME/_phased/}.bim ${PHASED_FNAME}.bim
        echo -e "\t>> COPY ${PHASED_FNAME/GENO_DIR/.}.bim\n"
    else
        awk 'BEGIN { OFS= "\t"} { print $1, $2, "0", $3, $4, $5}' \
            ${PHASED_FNAME}.haps >${PHASED_FNAME}.bim
        echo -e "\t>> ${PHASED_FNAME/$GENO_DIR/.}.bim"
    fi
fi
if [ ! -f ${PHASED_FNAME}.fam ]; then
    echo -e "\t- NOT exist FAM_FILE !"

    if [ -f ${PHASED_FNAME/_phased/}.fam ]; then
        cp ${PHASED_FNAME/_phased/}.fam ${PHASED_FNAME}.fam
        echo -e "\t>> COPY ${PHASED_FNAME/GENO_DIR/.}.fam\n"
    else
        awk 'BEGIN { OFS= "\t"} { if (NR>2) { print $1, $2, $4, $5, $6, $7}}' \
            ${PHASED_FNAME}.sample >${PHASED_FNAME}.fam
        echo -e "\t>> ${PHASED_FNAME/$GENO_DIR/.}.fam"
    fi
fi

# b. HLA imputation using by DEEP-HLA
echo -e "- HLA imputation using by DEEP-HLA\n"
source ${VENV_DIR}/${VNAME}/bin/activate

cd ${DHLA_DIR}

echo -e "\t- HLA Imputatetion\n"
python impute.py \
    --sample ${PHASED_FNAME} \
    --model ${GENO_DIR_HLA}/jp/jp \
    --hla ${GENO_DIR_HLA}/jp/jp \
    --model-dir ${GENO_DIR_HLA}/model \
    --phased-type haps \
    --out $IMPUTED_FNAME
mv ${DHLA_DIR}/imputation.${today}*.log ${GENO_DIR_HLA}

# c. Imputation of amino acid polymorphisms by DEEP-HLA
echo -e "\n\t- Imputation of amino acid polymorphisms"
python impute_aa.py \
    --dosage $IMPUTED_FNAME \
    --aa-table ${GENO_DIR_HLA}/jp/jp \
    --out $IMPUTED_FNAME
mv ${DHLA_DIR}/imputation_aa.${today}*.log ${GENO_DIR_HLA}

deactivate