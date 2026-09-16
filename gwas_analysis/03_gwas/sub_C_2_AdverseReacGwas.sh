#! /bin/bash
#------------------------------------------------------------------------------
#  Sub batch for performing sex-stratified analysis [adverse reaction]
#------------------------------------------------------------------------------

# [Arguments] (0) Environmental args
SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
SMPL_NAME=covid-vac

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMPQC=${GENO_DIR}/03_qc_imputed

STEP=C_2_AdverseReacGwas
WK_DIR=${SRC_DIR}/${STEP}
[[ ! -d "$WK_DIR" ]] && mkdir -p "$WK_DIR"

ETC_DIR=${SRC_DIR}/etc
PHENO_FILE=${ETC_DIR}/${SMPL_NAME}_seff.txt
ORD_PFILE=${ETC_DIR}/${SMPL_NAME}_seff.csv

PHENO_LIST=("seff_fever" "seff_arthralgia" "seff_fatigue" "seff_cold" "seff_headache" "seff_muscles")

SEX_LIST=(total male female)

TOOL_DIR=${SRC_DIR}/scripts/gwas_analysis
SUB_JL_FILE=${TOOL_DIR}/03_gwas/sub_${STEP}.jl
APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Prepare condition for analysis
if [ ! -f $ORD_PFILE ]; then
    awk -F "\t" -v OFS=',' \
        'NR == 1 {
            for (i = 1; i <= NF; i++) {
                $i = tolower($i)
                if ($i ~ /^seff_/) {
                    indexList[i] = 1
                }
            }
        }
        NR > 1 {
            for (i in indexList) {
                $i = $i + 1
            }
        }
        { print $0 }' $PHENO_FILE >$ORD_PFILE
    echo -e ">> ${ORD_PFILE/$SRC_DIR/.}\n"
fi
sex_idx=$(head -n 1 $ORD_PFILE | awk -F',' '{for (i=1; i<=NF; i++) if ($i == "sex") print i}')

# (1) Perform analysis
pheno_num=${#PHENO_LIST[@]}
phenotype=${PHENO_LIST[$((SLURM_ARRAY_TASK_ID % pheno_num))]}

imp_dir=$GENO_DIR_IMPQC
stat_dir=${WK_DIR}/${phenotype}
[[ ! -d "$stat_dir" ]] && mkdir -p "$stat_dir"

echo -e "\n-----------------------------------------------------------"
echo -e "[$phenotype]"
echo -e "-----------------------------------------------------------\n"

for ((i = 0; i < ${#SEX_LIST[@]}; i++)); do
    sex=${SEX_LIST[$i]}
    sex_name=$(printf "%02d" $i)_${sex}

    if [ $i -eq 0 ]; then
        ord_file=$ORD_PFILE
    elif [ $i -eq 1 ] || [ $i -eq 2 ]; then
        ord_file=$(echo $ORD_PFILE | sed -e "s/.csv/_${sex_name}.csv/")

        if [[ ! -f $ord_file && -n $sex_idx ]]; then
            awk -F ',' -v OFS="," -v sex_idx=$sex_idx -v sex_val=$i \
                'NR == 1 {
                        print $0
                    }
                    NR > 1 {
                        if($sex_idx == sex_val) { print $0 }
                } ' $ORD_PFILE | cut -d, -f1-$((sex_idx - 1)),$((sex_idx + 1))- >$ord_file
        fi
    fi
    cov_list=$(awk -F ',' 'NR==1 {
            for(i=1; i<=NF; i++) {
                if ($i ~ /sex/ || $i ~ /age/ || $i ~ /pc[0-9]+/) {
                    cnames = (cnames ? cnames "+" : "") $i
                }
            }
            print cnames
        }' "${ord_file}")

    echo -e "\n- ${sex}"

    edir_stat=${stat_dir}/${sex_name}
    [ ! -d $edir_stat ] && mkdir -p $edir_stat
    assoc_fname="${edir_stat}/${SMPL_NAME}_${phenotype}_${sex_name}"

    if [[ -f "${assoc_fname}_pval.txt.gz" ]]; then
        echo -e "\tAlready performed this step"
        echo -e "\t>> ${assoc_fname/${WK_DIR}/.}_pval.txt.gz\n"
        continue
    fi

    for chr in {1..22} X; do
        ord_fname=${SMPL_NAME}_${phenotype}_chr${chr}_${sex_name}
        csv_file=${edir_stat}/${ord_fname}_pval.csv
        if [[ -f $csv_file ]]; then
            echo -e "\t<< ${csv_file/${WK_DIR}/.}"
            continue
        fi

        # a. Prepare and check files
        imp_gname=${imp_dir}/${SMPL_NAME}_chr${chr}_imputed_qc
        if [ ! -f ${imp_gname}.bim ]; then
            echo -e "\t!! NOT EXIST plink binary files !!"
            echo -e "\t>> ${imp_gname/${GENO_DIR}/.}.(fam|bim|bed)"
            continue
        fi
        echo -e "\t<< ${imp_gname/$GENO_DIR/.}.(bim|bed)"

        # b. Association study [OrdinalGWAS]
        # ${APPC_DIR}/OrdinalGWAS.sif \
        apptainer exec \
            ${APPC_DIR}/OrdinalGWAS \
            julia -e 'include("'"$SUB_JL_FILE"'"); \
                ordinal_lrt(
                    phenotype="'"$phenotype"'",
                    cov_list="'"$cov_list"'",
                    ord_file="'"$ord_file"'",
                    bim_name="'"$imp_gname"'",
                    out_prefix="'"${edir_stat}/${ord_fname}"'")'
    done

    # (2) Summarize results
    # a. Check whether error logs are exported or not
    elogs_fnum=$(find $edir_stat -maxdepth 1 -type f -name "*.error" | wc -l)

    if [[ $elogs_fnum -gt 0 ]]; then
        echo -e "\n\t- Exported ERROR_LOG !!"
        echo -e "\t>> ${edir_stat/${WK_DIR}/.}"
        continue
    fi

    # a. Merge result files to ONE file
    [ ! -d ${edir_stat}/archives ] && mkdir -p ${edir_stat}/archives
    echo -e "\n\t- Summarize results"
    stat_flist=($(find ${edir_stat} -maxdepth 1 -type f \
        -name "${SMPL_NAME}_${phenotype}_chr*_${sex_name}_pval.csv" | sort -V))
    if [ ${#stat_flist[@]} -eq 0 ]; then
        echo -e "\t!! NOT EXIST association result files !!"
        echo -e "\t>> ${edir_stat/${WK_DIR}/.}"
        continue
    fi

    assoc_file="${assoc_fname}_pval.txt"
    cp /dev/null $assoc_file
    for stat_file in "${stat_flist[@]}"; do
        if [ ! -s $assoc_file ]; then
            sed -e "s/^#//" -e "s/,/\\t/g" $stat_file >$assoc_file
        else
            sed "1d" $stat_file | sed -e "s/,/\\t/g" >>$assoc_file
        fi
        gzip -f $stat_file
        mv ${stat_file}.gz ${edir_stat}/archives
    done
    gzip -f $assoc_file
    echo -e "\t>> ${assoc_file/${WK_DIR}/.}.gz"
done