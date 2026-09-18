#! /bin/bash
#------------------------------------------------------------------------------
#  Sub batch for performing sex-stratified analysis [antibody titer]
#------------------------------------------------------------------------------

set -eu

# [Arguments] (1) Environmental args
SRC_DIR=/path/to/project
SMPL_NAME=covid-vac

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed

STEP=C_1_AntibodyTiterGwas
WK_DIR=${SRC_DIR}/${STEP}

ETC_DIR=${SRC_DIR}/etc
COV_FILE=${ETC_DIR}/${SMPL_NAME}_plink.covar

PHENOTYPE=post_titer_norm
SEX_LIST=(total male female)

#------------------------------------------------------------------------------
# [Main script ]
# (0) Prepare parameters
sex_num=${#SEX_LIST[@]}
sex_idx=$((SLURM_ARRAY_TASK_ID % sex_num))
sex=${SEX_LIST[$sex_idx]}
sex_name=$(printf "%02d" $sex_idx)_${sex}

echo -e "\n-----------------------------------------------------------"
echo -e "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID"
echo -e "[$sex_name]"
echo -e "-----------------------------------------------------------\n"

# (1) Perform analysis
imp_dir=$GENO_DIR_IMPQC
stat_dir=${WK_DIR}/${sex_name}
[[ ! -d $stat_dir ]] && mkdir -p $stat_dir
pheno_file=${ETC_DIR}/${SMPL_NAME}_${sex_name}.txt

if [[ $sex_idx -eq 0 ]]; then
    glm_cond="--glm sex"
elif [[ $sex_idx -eq 1 ]] || [[ $sex_idx -eq 2 ]]; then
    glm_cond="--glm --filter-${sex}s"
fi

for chr in {1..22} X; do
    imp_gname=${SMPL_NAME}_chr${chr}_imputed_qc

    if [ ! -f ${imp_dir}/${imp_gname}.pgen ]; then
        echo -e "\t!! NOT EXIST plink binary files !!"
        echo -e "\t>> ${imp_dir/${GENO_DIR}/.}/${imp_gname}.(pgen|pvar|psam)"
        continue
    fi

    # a. Association study [plink]
    echo -e "\n- Association study [chr${chr}]"
    echo -e "\t<< ${imp_dir/$GENO_DIR/.}/${imp_gname}.(pgen|pvar|psam)"
    stat_fname=${SMPL_NAME}_${PHENOTYPE}_chr${chr}_${sex_name}
    plink2 \
        --silent \
        --pfile ${imp_dir}/${imp_gname} \
        --pheno $pheno_file \
        --pheno-name $PHENOTYPE \
        --covar-variance-standardize \
        --require-pheno $PHENOTYPE \
        --covar $COV_FILE \
        $glm_cond \
        --ci 0.95 \
        --out ${stat_dir}/${stat_fname}
    echo -e "\t>> ${stat_dir/$WK_DIR/.}/${stat_fname}.glm.linear"
done

# (2) Summarize results [glm]
echo -e "\n- Summarize results [glm]"
stat_flist=($(find ${stat_dir} -maxdepth 1 -type f \
    -name "${SMPL_NAME}_${PHENOTYPE}_chr*_${sex_name}.${PHENOTYPE}.glm.linear" | sort -V))
if [ ${#stat_flist[@]} -eq 0 ]; then
    echo -e "\t!! NOT EXIST association result files !!"
    echo -e "\t>> ${stat_dir/${WK_DIR}/.}"
else
    [ ! -d ${stat_dir}/archives ] && mkdir -p ${stat_dir}/archives
fi

assoc_file=${stat_dir}/${SMPL_NAME}_${PHENOTYPE}_${sex_name}.glm
cp /dev/null $assoc_file
for stat_file in "${stat_flist[@]}"; do
    if [ ! -s $assoc_file ]; then
        head -n 1 $stat_file | sed -e "s/^#//" >$assoc_file
    fi
    grep -w "ADD" $stat_file >>$assoc_file
    gzip -f $stat_file
    mv ${stat_file}.gz ${stat_dir}/archives
done
gzip -f $assoc_file
echo -e "\t>> ${assoc_file/${WK_DIR}/.}.gz"