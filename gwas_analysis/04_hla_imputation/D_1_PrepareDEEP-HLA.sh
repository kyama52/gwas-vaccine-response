#! /bin/bash
#------------------------------------------------------------------------------
#  Prepare for running DEEP-HLA
#------------------------------------------------------------------------------
#  Memo
# (1) Prepare genotype data
# (2) HLA imputation using by DEEP-HLA
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=8gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -c 12
#SBATCH -p cpu

set -eu

#------------------------------------------------------------------------------
# [Arguments] (2) Environmental args
SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
SMPL_NAME=covid-vac

CHR=6
HLA_STR=24
HLA_END=36

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_QC=${GENO_DIR}/01_qc
GENO_DIR_HLA=${GENO_DIR}/04_imputation_hla
GENO_DIR_ETC=${GENO_DIR_HLA}/etc

BUILD=b37
GMAP_FILE=$HOME/resource/genetic_maps/chr${CHR}.${BUILD}.gmap.gz

APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [HLA imputation analysis]"
echo -e "###\t- SOFTWARE(PRE-PHASING): shapeit2:"
echo -e "###\t\t- BUILD_VER: ${BUILD}"
echo -e "###\t\t- GENETIC_MAP_FILE:"
echo -e "###\t\t>> ${GMAP_FILE}"
echo -e "###"
echo -e "###\t- SOFTWARE(IMPUTATION): DEEP*HLA:"
echo -e "###"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t>> ${GENO_DIR_QC/$GENO_DIR/.}"
echo -e "###\t\t>> ${GENO_DIR_HLA/$GENO_DIR/.}"
echo -e "###----------------------------------------------------###\n"

# (1) Prepare genotype data
hla_fname=${SMPL_NAME}_hla
hla_fpath=${GENO_DIR_HLA}/${hla_fname}

# a. Remove duplicated variants and flip
if [ ! -f ${hla_fpath}.bed ]; then
    echo -e "\n- Check duplicated variants\n"
    plink1.9 \
        --bfile ${GENO_DIR_QC}/${SMPL_NAME} \
        --chr ${CHR} \
        --flip ${GENO_DIR_ETC}/${hla_fname}.flip \
        --exclude ${GENO_DIR_ETC}/${hla_fname}.exclude \
        --from-mb $HLA_STR --to-mb $HLA_END \
        --list-duplicate-vars ids-only suppress-first \
        --out ${GENO_DIR_HLA}/${SMPL_NAME}

    exclude_vfile=${hla_fpath}.dupvar
    cat ${GENO_DIR_ETC}/${hla_fname}.exclude ${GENO_DIR_HLA}/${SMPL_NAME}.dupvar >$exclude_vfile
    echo -e "\n>> ${exclude_vfile/${GENO_DIR_HLA}/.}"

    echo -e "\n- Prepare genotype data"
    plink1.9 \
        --bfile ${GENO_DIR_QC}/${SMPL_NAME} \
        --chr ${CHR} \
        --flip ${GENO_DIR_ETC}/${hla_fname}.flip \
        --exclude ${exclude_vfile} \
        --from-mb $HLA_STR --to-mb $HLA_END \
        --make-bed \
        --out ${hla_fpath}

    echo -e "\n\t>> ${hla_fpath/${GENO_DIR}/.}.(bim|bed|fam)"
fi

# b. Reformat BIM file for DEEP-HLA
if [ ! -f ${hla_fpath}.bimbak ]; then
    echo -e "\n- Reformat BIM file for DEEP-HLA"
    mv ${hla_fpath}.bim ${hla_fpath}.bimbak
    awk '{print $1,"chr6:"$4,$3,$4,$5,$6}' ${hla_fpath}.bimbak >${hla_fpath}.bim
    echo -e "\t ${hla_fpath/${GENO_DIR_HLA}/.}.bim"
    echo -e "\t -> ${hla_fpath/${GENO_DIR_HLA}/.}.bimbak"
fi

# (2) Prephasing by shapeit
# a. Run shapeit
echo -e "\n- Prephasing by shapeit"
apptainer exec ${APPC_DIR}/gwas.sif shapeit \
    --input-bed ${hla_fpath}.bed ${hla_fpath}.bim ${hla_fpath}.fam \
    --input-map ${GMAP_FILE} \
    --output-max ${hla_fpath}_phased \
    --thread 12 \
    --output-log ${hla_fpath}_phased