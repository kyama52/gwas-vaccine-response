#! /bin/bash
#------------------------------------------------------------------------------
#  Apply Qc to samples
#------------------------------------------------------------------------------
#  Procedures for QC:
#   (0) Prepare genotype files for appling minimal QC
#   (1) Identification of individuals with discordant sex information
#   (2) Identification of individuals with elevated missing data rates or
#       outlying heterozygosity rate
#   (3) Identification of duplicated or related individuals
#   (4) Remove individuals failing QC (1)
#   (5) Prepare genotype files for PCA
#   (6) Identification of individuals of divergent ancestry [ALL]
#   (7) Identification of individuals of divergent ancestry [EAS]
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=4gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
BUILD_VER=hg19-GRCh37
RES_DIR=$HOME/resource

RES_DIR_HAP=${RES_DIR}/Others/QcTools/PCA/HapmapForPCA
RESHAP=hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps
HIGH_LD=high-LD-regions.txt
NON_EAS=hapmap3r2_nonEAS.txt

SRC_DIR=/path/to/project
SMPL_NAME=covid-vac

STEP=A_1_ApplyQc
WK_DIR_QC=${SRC_DIR}/${STEP}
PRUNLD_FILE=${WK_DIR_QC}/${SMPL_NAME}_prunLD.prune.in

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_WK=${GENO_DIR}/00_working
[ ! -d $GENO_DIR_WK ] && mkdir -p $GENO_DIR_WK
GENO_FILE=${GENO_DIR_WK}/${SMPL_NAME}

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
. ${TOOL_DIR}/00_utils/gwas_tools.sh

EIGEN_DIR=$HOME/local/EIG-6.1.4/bin
export PATH=$EIGEN_DIR:$PATH

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Apply QC for genotypes]"
echo -e "###\t- STEP: ${STEP}"
echo -e "###\t- BUILD_VER: ${BUILD_VER}"
echo -e "###"
echo -e "###\t- HAP_DIR: ${RES_DIR_HAP}"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t- ORIGINAL FILES:"
echo -e "###\t\t>> ./${SMPL_NAME}.(ped|map)"
echo -e "###\t\t- GENO_DIR_WK:"
echo -e "###\t\t>> ${GENO_DIR_WK/$GENO_DIR/.}"
echo -e "###"
echo -e "###\t\ WK_DIR_QC:"
echo -e "###----------------------------------------------------###\n"

# (0) Prepare genotype files for appling minimal QC
plink1.9 \
    --bfile ${GENO_DIR}/${SMPL_NAME} \
    --chr 1-23 \
    --maf 0.01 \
    --geno 0.05 \
    --hwe 0.00001 \
    --allow-no-sex \
    --make-bed \
    --out $GENO_FILE
echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"

echo -e "\n----------------------------------------------------------\n"

# (1) Identification of individuals with discordant sex information
CHK_DIR=${WK_DIR_QC}/01_sex_check
[ ! -d $CHK_DIR ] && mkdir -p $CHK_DIR

echo -e "(1) Identification of individuals with discordant sex information"
echo -e ">> ${CHK_DIR/${WK_DIR_QC}/.}\n"

plink1.9 \
    --bfile $GENO_FILE \
    --check-sex \
    --out ${CHK_DIR}/${SMPL_NAME}
echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"

CHK_FILE=${CHK_DIR}/${SMPL_NAME}_fail_sex.txt
grep PROBLEM ${CHK_DIR}/${SMPL_NAME}.sexcheck |
    sed -e 's/^\s\+//g' | sed -e 's/\s\+/\t/g' |
    cut -f1-2 >$CHK_FILE

echo -e "\n----------------------------------------------------------\n"

# (2) Identification of individuals with elevated missing data rates or outlying heterozygosity rate
CHK_DIR=${WK_DIR_QC}/02_het_imiss_ibs
[ ! -d $CHK_DIR ] && mkdir -p $CHK_DIR

echo -e "(2) Identification of individuals with elevated missing data rates or outlying heterozygosity rate"
echo -e ">> ${CHK_DIR/${WK_DIR_QC}/.}\n"

plink1.9 \
    --bfile ${GENO_DIR_WK}/${SMPL_NAME} \
    --het \
    --out ${CHK_DIR}/${SMPL_NAME}
echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"

plink1.9 \
    --bfile ${GENO_DIR_WK}/${SMPL_NAME} \
    --missing \
    --out ${CHK_DIR}/${SMPL_NAME}
echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"

echo -e "\n----------------------------------------------------------\n"

# (3) Identification of duplicated or related individuals
CHK_DIR=${WK_DIR_QC}/02_het_imiss_ibs
[ ! -d $CHK_DIR ] && mkdir -p $CHK_DIR
echo -e "(3) Identification of duplicated or related individuals (1)"
echo -e "Starting ..."
echo -e ">> ${CHK_DIR/${WK_DIR_QC}/.}\n"

if [ ! -f $PRUNLD_FILE ]; then
    echo -e "- Variant pruning"
    plink1.9 \
        --bfile ${GENO_DIR_WK}/${SMPL_NAME} \
        --exclude ${RES_DIR_HAP}/${HIGH_LD} \
        --range --indep-pairwise 50 6 0.2 \
        --out ${PRUNLD_FILE/.prune.in/}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}.genome
if [ ! -f $TGT_FNAME ]; then
    echo -e "- Calculation IBS\n"
    plink1.9 \
        --bfile ${GENO_DIR_WK}/${SMPL_NAME} \
        --extract $PRUNLD_FILE \
        --genome --out ${CHK_DIR}/${SMPL_NAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

echo -e "- Extract duplicated individuals\n"
perl ${TOOL_DIR}/00_utils/plinkQC_IBD.pl ${CHK_DIR}/${SMPL_NAME}

echo -e ">> ./${SMPL_NAME}_fail_ibs.txt\n"

echo -e "\n----------------------------------------------------------\n"

# (4) Remove individuals failing QC (1)
echo -e "(4) Remove individuals failing QC (1)"
echo -e "- Exclude samples not satisfied minimal QC"
echo -e ">> ./${SMPL_NAME}_fail_min_smpl.txt"
find ${WK_DIR_QC} -type f -name "*_fail_*.txt" | xargs cat |
    sed -e 's/\s/\t/g' | sort -k1 |
    uniq >${WK_DIR_QC}/${SMPL_NAME}_fail_min_smpl.txt

TGT_FNAME=${GENO_DIR_WK}/${SMPL_NAME}_min_qc
if [ ! -f ${TGT_FNAME}.bed ]; then
    echo -e "- Apply minimal QC\n"
    plink1.9 \
        --allow-no-sex \
        --bfile ${GENO_DIR_WK}/${SMPL_NAME} \
        --remove ${WK_DIR_QC}/${SMPL_NAME}_fail_min_smpl.txt \
        --make-bed \
        --out $TGT_FNAME
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

echo -e "\n----------------------------------------------------------\n"

# (5) Prepare genotype files for PCA
CHK_DIR=${WK_DIR_QC}/03_pca
[ ! -d $CHK_DIR ] && mkdir -p $CHK_DIR
echo -e "(5) Prepare genotype files for PCA"
echo -e ">> ${CHK_DIR/${WK_DIR_QC}/.}\n"

# a. Merging data: For getting ${SMPL_NAME}.missnp
echo -e "- Merging data: For getting ${SMPL_NAME}.missnp"
TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_no_atcg
if [ ! -f ${TGT_FNAME}.bed ]; then
    plink1.9 \
        --allow-no-sex \
        --bfile ${GENO_DIR_WK}/${SMPL_NAME}_min_qc \
        --extract ${RES_DIR_HAP}/${RESHAP}.txt \
        --make-bed \
        --out $TGT_FNAME
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}
if [ ! -f ${TGT_FNAME}.bed ]; then
    plink1.9 \
        --allow-no-sex \
        --bfile ${CHK_DIR}/${SMPL_NAME}_no_atcg \
        --extract ${PRUNLD_FILE} \
        --make-bed \
        --out ${TGT_FNAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_hap
MISSNP_LIST=($(find ${CHK_DIR} -type f -name "${SMPL_NAME}*.missnp"))
if [ ${#MISSNP_LIST[@]} -eq 0 ]; then
    set +e
    plink1.9 \
        --allow-no-sex \
        --bfile ${CHK_DIR}/${SMPL_NAME} \
        --bmerge ${RES_DIR_HAP}/${RESHAP} \
        --flip-scan 'verbose' \
        --extract ${CHK_DIR}/${SMPL_NAME}.bim \
        --make-bed \
        --out ${TGT_FNAME}
    set -e
    echo -e "\n\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

# b. Merging data: Flipping by ${SMPL_NAME}.missnp
echo -e "- Merging data: Flipping by ${SMPL_NAME}.missnp"
MISSNP_LIST=($(find ${CHK_DIR} -type f -name "${SMPL_NAME}_hap*.missnp"))
if [ ${#MISSNP_LIST[@]} -eq 1 ]; then
    MISSNP_FILE=${MISSNP_LIST[0]}
else
    echo "Not ONLY ONE MISSNP FILE!"
    echo "  : ${CHK_DIR}/${SMPL_NAME}_hap*.missnp"
    exit 1
fi

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_flip
if [ ! -f ${TGT_FNAME}.bed ]; then
    plink1.9 \
        --allow-no-sex \
        --bfile ${CHK_DIR}/${SMPL_NAME} \
        --flip $MISSNP_FILE \
        --make-bed \
        --out $TGT_FNAME
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

if [ ! -f ${CHK_DIR}/${RESHAP}.bed ]; then
    echo -e ">> ./${SMPL_NAME}_flip.bed\n"
    plink1.9 \
        --allow-no-sex \
        --bfile ${RES_DIR_HAP}/${RESHAP} \
        --extract ${CHK_DIR}/${SMPL_NAME}_flip.bim \
        --make-bed \
        --out ${CHK_DIR}/${RESHAP}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

echo -e "\n----------------------------------------------------------\n"

# (6) Identification of individuals of divergent ancestry [ALL]
# a. Merge genotypes with HAPMAP using flipped SNPs
echo -e "(6) Identification of individuals of divergent ancestry [ALL]"
echo "- Prepare genotypes with HAPMAP using flipped SNPs\n"

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_allhap
if [ ! -f ${TGT_FNAME}.bed ]; then
    plink1.9 \
        --allow-no-sex \
        --bfile ${CHK_DIR}/${SMPL_NAME}_flip \
        --bmerge ${CHK_DIR}/${RESHAP} \
        --extract ${PRUNLD_FILE} \
        --chr 1-22 \
        --make-bed \
        --out ${TGT_FNAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_allhap_chk
if [ ! -f ${TGT_FNAME}.imiss ]; then
    plink1.9 \
        --bfile ${CHK_DIR}/${SMPL_NAME}_allhap \
        --missing \
        --out ${TGT_FNAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

# c. Performing PCA analysis Using EIGENSTRAT & CONVERTF
TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_allhap

echo -e "- Fill zero to cM column in BIM files"
echo -e ">> ${TGT_FNAME/${WK_DIR_QC}/.}.bim\n"
overwrite_bim_cmcol ${TGT_FNAME}.bim

echo -e "- Replace phenotype column [unknown: 0/-9 -> 99]"
echo -e ">> ${TGT_FNAME/${WK_DIR_QC}/.}.fam\n"
mv ${TGT_FNAME}.fam ${TGT_FNAME}.fambak
sed -e 's/-9$/99/' -e 's/0$/99/' ${TGT_FNAME}.fambak >${TGT_FNAME}.fam

if [ -f ${TGT_FNAME}.bed ]; then
    echo -e "- Start PCA analysis [ALL]\n"
    smartpca.perl \
        -i ${TGT_FNAME}.bed \
        -a ${TGT_FNAME}.bim \
        -b ${TGT_FNAME}.fam \
        -o ${TGT_FNAME/allhap/allpop}.pca \
        -p ${TGT_FNAME/allhap/allpop}.plot \
        -e ${TGT_FNAME/allhap/allpop}.eval \
        -l ${TGT_FNAME/allhap/allpop}.log \
        -k 10 -t 2 -m 0 -s 6.0

    convert_eigen2plink ${TGT_FNAME/allhap/allpop}.pca.evec

else
    echo -e "Not Exists Genotype file :"
    echo -e "\t :${TGT_FNAME}.bed \n"
fi

echo -e "\n----------------------------------------------------------\n"

# (7) Identification of individuals of divergent ancestry [EAS]
echo -e "(7) Identification of individuals of divergent ancestry [EAS]"

#   a. Extract genotype data for PCA analysis [EAS]
TGT_FNAME=${CHK_DIR}/${SMPL_NAME}_eashap
if [ ! -f ${TGT_FNAME}.bed ]; then
    echo "- Extract genotype data for PCA analysis [EAS]\n"
    plink1.9 \
        --allow-no-sex \
        --bfile ${TGT_FNAME/eas/all} \
        --remove ${RES_DIR_HAP}/${NON_EAS} \
        --make-bed \
        --out ${TGT_FNAME}
    echo -e "\tFinished: $(date +%y/%m/%d\ %H:%M:%S)\n"
fi

#   b. Performing PCA analysis Using EIGENSTRAT & CONVERTF
if [ -f ${TGT_FNAME}.bed ]; then
    echo -e "- Start PCA analysis [EAS]\n"
    smartpca.perl \
        -i ${TGT_FNAME}.bed \
        -a ${TGT_FNAME}.bim \
        -b ${TGT_FNAME}.fam \
        -o ${TGT_FNAME/eashap/easpop}.pca \
        -p ${TGT_FNAME/eashap/easpop}.plot \
        -e ${TGT_FNAME/eashap/easpop}.eval \
        -l ${TGT_FNAME/eashap/easpop}.log \
        -k 10 -t 2 -m 0 -s 6.0

    convert_eigen2plink ${TGT_FNAME/eashap/easpop}.pca.evec

else
    echo -e "Not Exists Genotype file :"
    echo -e "\t :${TGT_FNAME}.bed \n"
fi