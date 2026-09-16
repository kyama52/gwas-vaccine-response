#! /bin/bash
# ==============================================================================
#  001_conv_genofiles.sh
#-------------------------------------------------------------------------------
#  Sub routine script for convert genotype data from VCF files
#------------------------------------------------------------------------------
#   Procedures:
#   a. Apply QC to VCF file
#   b. Convert to plink binary files
#   c. Convert to Oxford files (gen + gzip)
#------------------------------------------------------------------------------
#   [Attension]
#   - Convert vcf to Oxford genotype format
#       [qctool] NOT supported "VCFv4.3"
#           - CANNOT convert VCF from 4.3 to 4.2
#       [plink2]
#           - convert VCF to bgen (--export bgen-1.3) -> NOT working
#           - convert VCF to gen (--export oxford-v2) -> OK
#           - convert VCF to gen (--export oxford-v2) + gzip -> OK
#       [plink2 + qctool]
#           - VCF 4.3 to oxford-v2 by plink2
#           - Convert gen to bgen (v1.2) -> OK
#
#===============================================================================
# [Arguments]

set -eu

# (1) Environmental args
IN_FILE=
FAM_FILE=
R2_VAL=0.3
MAF_VAL=0.01
REF_DIR=""

OUT_DIR=

# (2) Get options
usage() {
    echo -e "Usage: ${0}"
    echo -e "\t-i IN_FILE -f FAM_FILE -m MAF_VAL -c R2_VAL"
    echo -e "\t [-o OUT_DIR -r REF_DIR]"
    echo -e "  This script is submodule to convert VCF files to various formats"
    echo -e ""
    echo -e "Options: "
    echo -e "  -i: INPUT VCF file exported by PLINK2"
    echo -e "  -f: FAM file"
    echo -e "  -m: CUTOFF value of MAF [DEFAULT=0.01]"
    echo -e "  -c: CUTOFF value of R2 [DEFAULT=0.3]"
    echo -e "  -o: OUTPUT path [OPTIONAL]"
    echo -e "  -s: REMOVE sample file [OPTIONAL]"
    echo -e "  -r: REFERENCE path [OPTIONAL]"
    exit 1
}

while getopts "i:f:m:c:o:s:r:h" OPT; do
    case "$OPT" in
    i) IN_FILE=$OPTARG ;;
    f) FAM_FILE=$OPTARG ;;
    m) MAF_VAL=$OPTARG ;;
    c) R2_VAL=$OPTARG ;;
    o) OUT_DIR=$OPTARG ;;
    s) RMV_SFILE=$OPTARG ;;
    r) REF_DIR=$OPTARG ;;
    h | *) usage ;;
    esac
done

if [ ! -f $IN_FILE ]; then
    echo -e "NOT exists indicated file !"
    echo -e "\t >> $IN_FILE"
    exit 1
fi
IN_DIR=${IN_FILE%/*}
IN_FNAME=$(echo ${IN_FILE##*/} | sed -E 's/.(vcf.gz|vcf)$//')
OUT_FNAME="${IN_FNAME}_qc"

if [ ! -f $FAM_FILE ]; then
    echo -e "NOT Exist FAM_FILE"
    echo -e "\t>> ${FAM_FILE}"
fi

[ -z $OUT_DIR ] && OUT_DIR=$IN_DIR
RMV_SCOND=""
if [ -n "$RMV_SFILE" ] && [ -f "$RMV_SFILE" ]; then
    RMV_SCOND="--remove $RMV_SFILE"
fi
[ -z $OUT_DIR ] && OUT_DIR=$IN_DIR
[ ! -d "$OUT_DIR" ] && mkdir -p $OUT_DIR
if [ -n "$REF_DIR" ] && [ ! -d "$REF_DIR" ]; then
    REF_DIR=""
fi

# [Main script ] -------------------------------------------------------
# (1) Logging
echo -e "----------------------------------------------------"
echo -e " Convert genotype data from VCF files"
echo -e "- IN_FILE:"
echo -e "\t>> $IN_FILE"
echo -e "- FAM_FILE:"
echo -e "\t>> $FAM_FILE"
echo -e "- THRESHOLD:"
echo -e "\t- R2_VAL < ${R2_VAL}"
echo -e "\t- MAF_VAL < ${MAF_VAL}"
echo -e "----------------------------------------------------"

# (2) Main procedures
# a. Apply QC to VCF file
echo -e "\n- Apply QC to VCF file and exported to plink2 format\n"

if [ -n "$RMV_SCOND" ]; then
    # --vcf $IN_FILE dosage=HDS \
    plink2 \
        --vcf $IN_FILE dosage=DS \
        --exclude-if-info "R2<=${R2_VAL}" \
        --fam $FAM_FILE \
        $RMV_SCOND \
        --maf $MAF_VAL \
        --make-pgen \
        --out ${OUT_DIR}/${OUT_FNAME}
else
    plink2 \
        --vcf $IN_FILE dosage=DS \
        --exclude-if-info "R2<=${R2_VAL}" \
        --fam $FAM_FILE \
        --maf $MAF_VAL \
        --make-pgen \
        --out ${OUT_DIR}/${OUT_FNAME}
fi
echo -e "\n\t>> ${OUT_DIR/$REF_DIR/.}/${OUT_FNAME}.(bed|bim|fam)\n"

# b. Convert plink2.0 binary files to various formats
echo -e "\n- Convert plink2.0 binary files to various formats\n"
echo -e "[VCF format]\n"
plink2 \
    --pfile ${OUT_DIR}/${OUT_FNAME} \
    --export vcf id-paste=iid vcf-dosage=DS \
    --out ${OUT_DIR}/${OUT_FNAME}
bgzip -f ${OUT_DIR}/${OUT_FNAME}.vcf
bcftools index -f -t ${OUT_DIR}/${OUT_FNAME}.vcf.gz
echo -e "\t>> ${OUT_DIR/$REF_DIR/.}/${OUT_FNAME}.vcf.gz"

echo -e "[plink1.9 format]\n"
plink2 \
    --pfile ${OUT_DIR}/${OUT_FNAME} \
    --hard-call-threshold 0.1 \
    --make-bed \
    --out ${OUT_DIR}/${OUT_FNAME}
echo -e "\n\t>> ${OUT_DIR/$REF_DIR/.}/${OUT_FNAME}.(bed|bim|fam)\n"

# echo -e "[Oxford binary files]\n"
# plink2 \
#     --pfile ${OUT_DIR}/${OUT_FNAME} \
#     --export oxford-v2 \
#     --out ${OUT_DIR}/${OUT_FNAME}
# gzip ${OUT_DIR}/${OUT_FNAME}.gen

# echo -e "\t>> ${OUT_DIR/$REF_DIR/.}/${OUT_FNAME}.(gen.gz|sample)"

echo -e "\n ...Finished $(date +"%Y/%m/%d %H:%M")\n"