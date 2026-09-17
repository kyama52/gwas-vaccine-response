#! /bin/bash
#------------------------------------------------------------------------------
#  Clean data by result of QC
#------------------------------------------------------------------------------
#   Procedures for QC:
#   (8) Identification of individuals of divergent ancestry [EAS Again]
#     a. Remove outliner by PC1 nad PC2
#     b. Remove samples using "${SMPL_NAME}_FailPcaSmpl.txt"
#         PC1 > 0.02
#         PC2 > 0.1
#     c. Performing PCA analysis Using EIGENSTRAT & CONVERTF
#   (9) Remove individuals failing QC (2)
#     a. List up all samples not satisfied with criteria
#     b. Apply sample QC
#     c. Performing PCA analysis Using EIGENSTRAT [ONLY samples]
#   (10) Check Markers
#     a. Identify all markers with an excessive missing data rate
#     b. Test markers for different genotype call rates between cases and contols
#     c. Remove all markers failing QC
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=4gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

# [Arguments] (2) Environmental args
SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
SMPL_NAME=covid-vac

STEP=A_2_CleanDataAfterQc
WK_DIR_QC=${SRC_DIR}/${STEP}
WK_DIR_PREQC=${SRC_DIR}/A_1_ApplyQc/03_pca

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_WK=${GENO_DIR}/00_working
GENO_DIR_QC=${GENO_DIR}/01_qc
[ ! -d $GENO_DIR_WK ] && mkdir -p $GENO_DIR_WK
[ ! -d $GENO_DIR_QC ] && mkdir -p $GENO_DIR_QC

COND_PC1="> 0.02"
COND_PC2="> 0.1"
SMPL_PTN="CUH-\\d{4}"

EIGEN_DIR=$HOME/local/EIG-6.1.4/bin
export PATH=$EIGEN_DIR:$PATH

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
. ${TOOL_DIR}/00_utils/gwas_tools.sh
APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Clean data after QC]"
echo -e "###\t- STEP: ${STEP}"
echo -e "###"
echo -e "###\t- WK_DIR: $WK_DIR"
echo -e "###\t\t- WK_DIR_PREQC:"
echo -e "###\t\t>> ${WK_DIR_PREQC/$WK_DIR/.}"
echo -e "###\t\t- WK_DIR_QC:"
echo -e "###\t\t>> ${WK_DIR_QC/$WK_DIR/.}"
echo -e "###----------------------------------------------------###\n"

# (8) Identification of individuals of divergent ancestry [EAS Again]
echo -e "(8) Identification of individuals of divergent ancestry [EAS Again]"
echo -e "- Remove outliner by PC1 nad PC2"
echo -e ">> ${WK_DIR_PREQC/${WK_DIR}/.}"

# a. Remove outliner by PC1 nad PC2
TGT_FNAME=${WK_DIR_PREQC}/${SMPL_NAME}_easpop
RMV_SLIST_FILE=${TGT_FNAME}_fail_pca.txt

[ -f $RMV_SLIST_FILE ] && rm $RMV_SLIST_FILE

apptainer exec ${APPC_DIR}/gwas.sif \
    Rscript ${TOOL_DIR}/00_utils/subbatch/101_remove_outlier.R \
    ${TGT_FNAME}.eigenvec \
    "${SMPL_PTN}" \
    "${COND_PC1}" \
    "${COND_PC2}"

# b. Remove samples using "$RMV_SLIST_FILE"
CHK_DIR=${WK_DIR_QC}/01_repca
[ ! -d $CHK_DIR ] && mkdir -p $CHK_DIR

echo -e "\n- Remove samples using ${RMV_SLIST_FILE##*/}"
echo -e ">> ${CHK_DIR/${WK_DIR_QC}/.}\n"

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_eashap_qc
plink1.9 \
    --allow-no-sex \
    --bfile ${WK_DIR_PREQC}/${SMPL_NAME}_eashap \
    --remove ${RMV_SLIST_FILE} \
    --make-bed \
    --out ${TGT_FNAME}
echo -e "  finished: $(date +%y/%m/%d\ %H:%M:%S)\n"

# c. Performing PCA analysis Using EIGENSTRAT & CONVERTF
if [ -f ${TGT_FNAME}.bed ]; then
    echo -e "- Start PCA analysis [EAS Again]\n"
    smartpca.perl \
        -i ${TGT_FNAME}.bed \
        -a ${TGT_FNAME}.bim \
        -b ${TGT_FNAME}.fam \
        -o ${TGT_FNAME}.pca \
        -p ${TGT_FNAME}.plot \
        -e ${TGT_FNAME}.eval \
        -l ${TGT_FNAME}.log \
        -k 10 -t 2 -m 0 -s 6.0

    convert_eigen2plink ${TGT_FNAME}.pca.evec
else
    echo -e "NOT Exists Genotype file :"
    echo -e "\t :${TGT_FNAME}.bed \n"
fi

# (9) Remove individuals failing QC (2)
echo -e "\n(9) Remove individuals failing QC (2)"

# a. List up all samples not satisfied with criteria
echo -e "- List up all samples not satisfied with criteria"
RMV_SLIST_FILE=${WK_DIR_QC}/${SMPL_NAME}_fail_all_smpl.txt
find ${WK_DIR_PREQC} ${WK_DIR_QC} -type f -name "*_fail_*.txt" |
    xargs cat |
    sed -e 's/\s/\t/g' | sort -k1n |
    uniq >$RMV_SLIST_FILE
echo -e ">> ${RMV_SLIST_FILE/$WK_DIR/.}\n"

# b. Apply sample QC
echo -e "- Apply sample QC"
GENO_FNAME=${GENO_DIR_QC}/${SMPL_NAME}
plink1.9 \
    --allow-no-sex \
    --bfile ${GENO_DIR_WK}/${SMPL_NAME} \
    --remove $RMV_SLIST_FILE \
    --chr 1-23 \
    --make-bed \
    --out $GENO_FNAME
echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"

# c. Performing PCA analysis Using EIGENSTRAT [ONLY samples]
PCA_FNAME=${CHK_DIR}/${SMPL_NAME}_eashap_qc
TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_csct

if [ ! -f ${TGT_FNAME}.bed ]; then
    echo -e "- Extract genotype data for PCA analysis [ONLY samples]\n"
    plink1.9 \
        --silent \
        --allow-no-sex \
        --bfile ${PCA_FNAME} \
        --keep-fam ${GENO_FNAME}.fam \
        --make-bed \
        --out ${TGT_FNAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

if [ -f ${TGT_FNAME}.bed ]; then
    echo -e "- Start PCA analysis [EAS]\n"
    smartpca.perl \
        -i ${TGT_FNAME}.bed \
        -a ${TGT_FNAME}.bim \
        -b ${TGT_FNAME}.fam \
        -o ${TGT_FNAME}.pca \
        -p ${TGT_FNAME}.plot \
        -e ${TGT_FNAME}.eval \
        -l ${TGT_FNAME}.log \
        -k 10 -t 2 -m 0 -s 6.0

    convert_eigen2plink ${TGT_FNAME}.pca.evec
else
    echo -e "Not Exists Genotype file :"
    echo -e "\t :${TGT_FNAME}.bed \n"
fi

# (10) Prepare each dataset and apply QC to markers
echo -e "\n(10) Prepare each dataset and apply QC to markers"
CHK_DIR=${WK_DIR_QC}/02_snp_check
[ ! -d $CHK_DIR ] && mkdir -p $CHK_DIR
echo -e ">> ${CHK_DIR/${WK_DIR_QC}/.}\n"

# a. Identify all markers with an excessive missing data rate
echo -e "-e Identify all markers with an excessive missing data rate\n"
plink1.9 \
    --bfile ${GENO_FNAME} \
    --missing \
    --out ${CHK_DIR}/${SMPL_NAME}
echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"