#! /bin/bash
# ==============================================================================
#  Sub routine script for concentration chunk files
#------------------------------------------------------------------------------
#  [Causion!]
#   Gloval variables in this script
#   - $HOME/resource/Eagle/BBJ1K_1KGP_RefPanel_b155
#------------------------------------------------------------------------------
#   Procedures:
#    a. Prepare temporary file and directory
#    b. Cleanup dosage files [R2 > '${R2_VAL}']
#    c. Concentrate VCF files after cleaning
#    d. Annotate rsID using reference panel
#    e. Concentrate info files
#    f. Concentrate log files
#    g. Remove original files [chunkXX.dose.vcf.gz & chunkXX.info]
#    [DISCARD] h. Convert to plink binary files
#===============================================================================
# [Arguments]
# (1) Environmental args
IN_DIR=
FAM_FILE=
R2_VAL=0.3
OUT_DIR=
REF_DIR=
ANNOT_DIR=$HOME/resource/minimac/BBJ1K
ANNOT_FNAME=BBJ1K_1KGP_RefPanel_b155

# (2) Get options
usage() {
    echo "Usage: ${0} -i IN_DIR -c R2_VAL -f FAM_FILE [-o OUT_DIR -r REF_DIR]"
    echo "  This script is submodule to concentrate and clean chunk files"
    echo ""
    echo "Options: "
    echo "  -i: INPUT directory with chunk files"
    echo "  -c: CUTOFF value of R2 [DEFAULT=0.3]"
    echo "  -f: FAM file"
    echo "  -o: OUTPUT path [optional]"
    echo "  -r: REFERENCE path [optional]"
    exit 1
}

while getopts "i:c:f:o:r:h" OPT; do
    case "$OPT" in
    i) IN_DIR=$OPTARG ;;
    c) R2_VAL=${OPTARG:=0.3} ;;
    f) FAM_FILE=$OPTARG ;;
    o) OUT_DIR=$OPTARG ;;
    r) REF_DIR=$OPTARG ;;
    h | *) usage ;;
    esac
done

if [ ! -d $IN_DIR ]; then
    echo -e "NOT exists indicated directory !"
    echo -e "\t >> $IN_DIR"
    exit 1
fi

if [ ! -f $FAM_FILE ]; then
    echo -e "NOT Exist FAM_FILE"
    echo -e "\t>> ${FAM_FILE}"
fi

[ -z $OUT_DIR ] && OUT_DIR=$IN_DIR
[ ! -d "$OUT_DIR" ] && mkdir -p $OUT_DIR
if [ -n "$REF_DIR" ] && [ ! -d "$REF_DIR" ]; then
    REF_DIR=""
fi

# [Main script ] -------------------------------------------------------
# (1) Get information
chunk_flist=($(find $IN_DIR -maxdepth 1 -type f -name "*chunk*" |
    grep -E ".*_chr([0-9]+|X)_imputed_chunk[0-9]+.dose.vcf.gz" | sort -V))

if [ ${#chunk_flist[@]} -eq 0 ]; then
    echo -e "The indicated directory NOT including chunk files!!"
    echo -e ">> $IN_DIR"
    exit 1
fi

sname_list=($(find $IN_DIR -maxdepth 1 -type f -name "*chunk*" |
    xargs -n1 basename |
    grep -v "chunk.tar.gz" |
    sed -E "s/(.*)_imputed_chunk[0-9]+(.dose.vcf.gz|.dose.log|.info)/\1_imputed/" | sort -V | uniq))

echo -e "----------------------------------------------------\n"
echo -e " Cleaning and concentration to chunk files"
echo -e "- INPUT_DIR:"
echo -e "\t>> $IN_DIR"
echo -e "- FAM_FILE:"
echo -e "\t>> $FAM_FILE"
echo -e "- THRESHOLD: R2_VAL < ${R2_VAL}"
echo -e ""
echo -e - ANNOTATION_FILES:
echo -e "\t>> ${ANNOT_DIR}/${ANNOT_FNAME}_chr[XX]"
echo -e "----------------------------------------------------\n"

# (2) Main procedures
for sname in "${sname_list[@]}"; do
    # a. Prepare temporary file and directory
    tmpfile=$(mktemp)
    trap 'rm -f $tmpfile' EXIT
    tmp_dir=$(mktemp -d -t concat-XXXXXXXX)
    trap 'rm -rf -- "$tmp_dir"' EXIT

    chr=$(echo $sname | sed -E "s/(.*)_chr([0-9]+|X)_imputed/\2/")

    echo -e "[$sname]\n"

    # b. Cleanup dosage files [R2 > '${R2_VAL}']
    echo -e "\n- Cleanup dosage files [R2 > '${R2_VAL}']"

    dose_flist=($(find $IN_DIR -maxdepth 1 -type f -name "${sname}_chunk*.dose.vcf.gz" | sort -V))

    if [ ${#dose_flist[@]} -eq 0 ]; then
        echo -e "\tNOT exist dosage file in this directory!"
        echo -e "\t>> ${IN_DIR/${REF_DIR}/.}"
        continue
    fi

    for dose_file in "${dose_flist[@]}"; do
        dose_fname=$(echo ${dose_file##*/} | sed -e 's/.vcf.gz//')
        plink2 \
            --silent \
            --memory 2000 require \
            --vcf ${dose_file} dosage=HDS \
            --exclude-if-info "R2<=${R2_VAL}" \
            --fam $FAM_FILE \
            --export vcf id-paste=iid vcf-dosage=HDS \
            --out ${tmp_dir}/${dose_fname}
        mv ${tmp_dir}/${dose_fname}.log $IN_DIR
        bgzip ${tmp_dir}/${dose_fname}.vcf
        echo -e "\t>> ${tmp_dir}/${dose_fname}.vcf.gz"
    done

    # c. Concentrate VCF files after cleaning
    echo -e "\n- Concentrate VCF files after cleaning"
    qc_flist=($(find ${tmp_dir} -maxdepth 1 -type f -name "${sname}_chunk*.dose.vcf.gz" | sort -V))
    if [ ${#qc_flist[@]} -eq 0 ]; then
        echo -e "\tNOT exist cleaned file in this directory!"
        echo -e "\t>> ${tmp_dir}"
        continue
    fi
    bcftools concat --threads 4 "${qc_flist[@]}" --output-type z --output $tmpfile

    # d. Annotate rsID using reference panel
    res_file=${ANNOT_DIR}/${ANNOT_FNAME}_chr${chr}.vcf.gz
    if [ -f $res_file ]; then
        echo -e "\n- Annotate rsID using reference panel"
        [ ! -f ${res_file}.csi ] && bcftools index $res_file
        echo -e "\t<< ${res_file/${RES_DIR}/.}"

        bcftools index $tmpfile
        bcftools annotate \
            --annotations ${res_file} \
            --columns ID \
            --output-type z \
            --output ${IN_DIR}/${sname}.vcf.gz \
            $tmpfile
    else
        mv $tmpfile ${IN_DIR}/${sname}.vcf.gz
    fi
    bcftools index -f ${IN_DIR}/${sname}.vcf.gz
    echo -e "\n\t>> ${IN_DIR/${REF_DIR}/.}/${sname}.vcf.gz"

    # e. Concentrate info files
    echo -e "\n- Concentrate info files"
    info_flist=($(find ${IN_DIR} -maxdepth 1 -type f -name "${sname}_chunk*.info" | sort -V))
    info_cfile=${IN_DIR}/${sname}.info
    cp /dev/null $info_cfile
    for info_file in "${info_flist[@]}"; do
        [ ! -s $info_cfile ] && head -n1 $info_file >$info_cfile
        sed -n '2,$p' $info_file >>$info_cfile
    done
    if [ -s $info_cfile ]; then
        gzip $info_cfile
        echo -e "\t>> ${info_cfile/${REF_DIR}/.}.gz"
    else
        rm $info_cfile
    fi

    # f. Concentrate log files
    echo -e "\n- Concentrate log files"
    log_flist=($(find ${IN_DIR} -maxdepth 1 -type f -name "${sname}_chunk*.info" | sort -V))
    log_cfile=${IN_DIR}/${sname}.log
    cat "${log_flist[@]}" >$log_cfile
    echo -e "\t>> ${log_cfile/${REF_DIR}/.}"

    # g. Archive and remove chunk files [chunkXX.dose.vcf.gz & chunkXX.info]
    echo -e "\n- Archive and remove chunk files [chunkXX.dose.vcf.gz & chunkXX.info]"

    cd $IN_DIR || exit 1
    find . -maxdepth 1 -type f -name "*chunk*" |
        grep -E "${sname}_chunk[0-9]+.(dose.vcf.gz|dose.log|info)" |
        sort -V |
        xargs tar cvzf ${sname}_chunk.tar.gz -C ${OUT_DIR}
    echo -e "\t>> ${OUT_DIR/$REF_DIR/.}/${sname}_chunk.tar.gz"
    rm ${IN_DIR}/${sname}_chunk*.dose.vcf.gz
    rm ${IN_DIR}/${sname}_chunk*.dose.log
    rm ${IN_DIR}/${sname}_chunk*.info

    # [2023.4.5]
    # # h. Convert to plink binary files
    # echo -e "\n- Convert to plink binary files"
    # plink2 \
    #     --silent \
    #     --vcf ${IN_DIR}/${sname}.vcf.gz dosage=HDS \
    #     --fam $FAM_FILE \
    #     --maf 0.01 \
    #     --hwe 0.00001 \
    #     --geno 0.05 \
    #     --make-pgen \
    #     --out ${IN_DIR}/${sname}

    # echo -e "\t>> ${IN_DIR/$REF_DIR/.}/${sname}.(pgen|psam|pvar)"

    echo -e "\n ...Finished $(date +"%Y/%m/%d %H:%M")\n"
done