#! /bin/bash
# ==============================================================================
#  sub_F_conc_snptest.sh
#
#  CD: Apr 12 2023      K.Yamazaki
#  UD: Aug 24 2024      Add option to archived files "tar -f"
#-------------------------------------------------------------------------------
#  Sub routine script for concentrate files exportd by snptest
#===============================================================================
# [Arguments]
# (1) Environmental args
IN_DIR=
REF_DIR=

# (2) Get options
usage() {
    echo "Usage: ${0} -i IN_DIR [REF_DIR]"
    echo "  This script is submodule to concentrate files exported by snptest"
    echo ""
    echo "Options: "
    echo "  -i: INPUT directory with chunk files"
    echo -e "  -r: REFERENCE path [OPTIONAL]"
    exit 1
}

while getopts "i:r:h" OPT; do
    case "$OPT" in
    i) IN_DIR=$OPTARG ;;
    r) REF_DIR=$OPTARG ;;
    h | *) usage ;;
    esac
done

if [ ! -d $IN_DIR ]; then
    echo -e "NOT exists indicated directory !"
    echo -e "\t >> $IN_DIR"
    exit 1
fi

if [ -n "$REF_DIR" ] && [ ! -d "$REF_DIR" ]; then
    REF_DIR=""
fi

# [Main script ] -------------------------------------------------------
# (1) Get information
stat_flist=($(find $IN_DIR -maxdepth 1 -type f -name "*_imputed" |
    grep -E ".*_chr([0-9]+|X)_imputed" | sort -V))

if [ ${#stat_flist[@]} -eq 0 ]; then
    echo -e "The indicated directory NOT including chunk files!!"
    echo -e ">> $IN_DIR"
    exit 1
fi

echo -e "----------------------------------------------------\n"
echo -e "[Concentrate files exportd by snptest]"
echo -e "- INPUT_DIR:"
echo -e "\t>> $IN_DIR"
echo -e "----------------------------------------------------\n"

# (2) Main procedures
assoc_file=$(echo ${stat_flist[0]} | sed -E 's/(.*)_chr(.*)_imputed/\1.assoc/')

cp /dev/null $assoc_file
for stat_file in "${stat_flist[@]}"; do
    if [[ ! $stat_file =~ "chrX" ]]; then
        if [ ! -s $assoc_file ]; then
            grep -v "^#" $stat_file | head -n 1 |
                cut -f 1-9,14-31,39- >$assoc_file
        fi
        grep -v "^#" $stat_file | sed -n '2,$p' |
            cut -f 1-9,14-31,39- >>$assoc_file
    fi
done
gzip -f $assoc_file
echo -e "\t>> ${assoc_file/$REF_DIR/.}.gz"

[ ! -d ${IN_DIR}/archives ] && mkdir ${IN_DIR}/archives
mv "${stat_flist[@]}" ${IN_DIR}/archives

cd ${IN_DIR}/archives
echo "${stat_flist[@]}" | xargs -n1 basename | xargs -n1 gzip
echo -e "\t>> ${IN_DIR/$REF_DIR/.}/archives"