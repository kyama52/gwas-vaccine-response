#! /bin/bash
# ==============================================================================
#  tools_gwas.sh
#===============================================================================
# Gloval Argument

# TODAY=`date +%y%m%d`

atexit() {
    [[ -n $tmpfile ]] && rm -f "$tmpfile"
}
tmpfile=$(mktemp)
trap atexit EXIT
trap 'trap - EXIT; atexit; exit -i' SIGHUP SIGINT SIGTERM

#===============================================================================
# overwrite_bim_cmcol
#
# ARGV[1] : BIM file
#-------------------------------------------------------------------------------
function overwrite_bim_cmcol() {
    local input_file
    local output_file

    if [ $# -ne 1 ]; then
        echo "ERROR: Input Files (BIM)"
        return 1
    fi

    input_file=$1
    if [ ! -f $input_file ]; then
        echo "ERROR: Not exist indicated files"
        echo "  : $input_file"
        return 1
    fi
    output_file=$input_file
    input_file=${output_file//bim/bimbak}
    mv $output_file $input_file

    sed -e 's/\s\+/\t/g' $input_file |
        awk -F "\t" \
            ' BEGIN {OFS="\t"} \
            { \
                if(NF == 6) print $1, $2, 0, $4, $5, $6; \
                else $0; \
            }' >$output_file

    return 0
}
#===============================================================================
# convert_eigen2plink
#
# ARGV[1] : evec file exported by EIGEN
# ARGV[2](option) : out_dir
#-------------------------------------------------------------------------------
function convert_eigen2plink() {
    local input_file
    local out_dir
    local title

    if [ $# -ne 1 ] && [ $# -ne 2 ]; then
        echo "ERROR: Input File and out_dir [optional] "
        return 1
    fi

    if [ $# -ne 1 ]; then
        echo "ERROR: Input filename exported by EIGEN"
        return 1
    fi

    input_file=$1
    if [ ! -f $input_file ]; then
        echo "ERROR: Not exist indicated file: $input_file"
        return 1
    fi

    out_dir=""
    if [ $# -eq 2 ]; then
        out_dir=$2
        if [ ! -d $out_dir ]; then
            mkdir -p $out_dir
        fi
    else
        out_dir="${input_file%/*}"
    fi

    output_file=${input_file##*/}
    output_file=${out_dir}/${output_file//.pca.evec/.eigenvec}

    sed -e 's/^\s\+//g' $input_file | sed -e 's/\s\+/\t/g' | sed -e 's/:/\t/g' |
        awk -F "\t" \
            ' BEGIN {OFS="\t"} \
            { \
                if (NR == 1) {
                    for (i=1; i<=(NF-2); i++){ pc_list=pc_list "\tPC"i } \
                    print "#FID", "IID" pc_list; \
                } else { \
                    print $0 \
                } \
            }' |
        cut -f 1-12 >$output_file

    input_file=${input_file//.pca.evec/.pca}
    output_file=${input_file##*/}
    output_file=${out_dir}/${output_file//.pca/.eigenval}
    [ -f $input_file ] && sed -n 2,11p $input_file >$output_file

    return 0
}