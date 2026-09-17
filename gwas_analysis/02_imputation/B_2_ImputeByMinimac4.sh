#! /bin/bash
#------------------------------------------------------------------------------
#  Imputing by minimac4
#------------------------------------------------------------------------------
#   Procedures:
#   (1) Create a file with all the SNP names that are in the reference set
#   (2) Creating the input files for phasing and imputation.
#     a. List up SNP list for imputation
#     b. Extract Genotypes of overlapped SNPs
#     c. By chromosome
#------------------------------------------------------------------------------
# [Arguments] (1) SLURM arguments
#SBATCH --mem=4gb
#SBATCH -o logs/%x.%j
#SBATCH -e logs/%x.%j
#SBATCH -p cpu

set -eu

# [Arguments] (2) Environmental args
RES_DIR=$HOME/resource
GMAP_FILE=${RES_DIR}/Eagle/genetic_map_hg19_withX.txt.gz
RES_DIR_EGL=${RES_DIR}/Eagle/BBJ1K
RES_DIR_MMC=${RES_DIR}/minimac/BBJ1K
RES_FNAME=BBJ1K_1KGP_RefPanel_b155

SRC_DIR=$HOME/analysis/COVID-19/2609_PrepGitRepo
SMPL_NAME=covid-vac

GENO_DIR=${SRC_DIR}/genotype
GENO_DIR_IMP=${GENO_DIR}/02_imputation
[ ! -d $GENO_DIR_IMP ] && mkdir -p $GENO_DIR_IMP

APPC_DIR=$HOME/tools/container

#------------------------------------------------------------------------------
# [Main script ]
# (0) Start logging
echo -e "###########################################################"
echo -e "### [Phasing and imputation]"
echo -e "###\t- IMPUTATION: minimac4"
echo -e "###"
echo -e "###\t- RES_DIR: ${RES_DIR}"
echo -e "###\t\t- EAGLE DIR: ${RES_DIR_EGL/$RES_DIR/.}"
echo -e "###\t\t- MINIMAC DIR: ${RES_DIR_MMC/$RES_DIR/.}"
echo -e "###\t- GENO_DIR: ${GENO_DIR}"
echo -e "###\t\t- GENO_DIR_IMP:"
echo -e "###\t\t>> ${GENO_DIR_IMP/$GENO_DIR/.}"
echo -e "###----------------------------------------------------###\n"

# (1) Imputation with minimac4. For every chromosome, perform imputations in chunks of 5 Mb
echo -e "(1) Imputation with minimac4."
echo -e "\tFor every chromosome, perform imputations in chunks of 5 Mb"

chunk_len=25000000
chunk_om=5000000
for chr in {1..22} X; do
    echo -e "\n[chr${chr}]\n"
    imp_dir=${GENO_DIR_IMP}/chr${chr}
    log_dir=${imp_dir}/logs

    each_gname=${imp_dir}/${SMPL_NAME}_chr${chr}_phased
    each_ref_file=${RES_DIR_MMC}/${RES_FNAME}_chr${chr}.m3vcf.gz

    if [ -f "${each_gname}.vcf.gz" ]; then
        [ ! -d ${log_dir} ] && mkdir -p ${log_dir}

        chr_len=$(zcat ${GMAP_FILE} | sed -e 's/^23/X/g' |
            awk -v chr=$chr '$1==chr {print $2}' |
            sort -n | tail -n 1)
        chunk_nr=$((chr_len / chunk_len))

        for chunk in $(seq 0 $chunk_nr); do
            str_pos=$(((chunk * chunk_len) + 1))
            end_pos=$(((chunk + 1) * chunk_len))
            [ ${end_pos} -gt ${chr_len} ] && end_pos=${chr_len}

            echo -e "- chunk(${chunk}): ${str_pos}-${end_pos}"
            que_type="cpu-normal"
            # [ $chr -le 8 ] && que_type="cpu-middle"

            sbatch \
                -J "mmc4_chr${chr}_${chunk}" \
                -o "${imp_dir}/logs/%x-%j" \
                -e "${imp_dir}/logs/%x-%j" \
                -c 4 \
                --mem 8G \
                --qos ${que_type} \
                --wrap="apptainer exec ${APPC_DIR}/gwas.sif minimac4 \
                    --cpus 4 \
                    --refHaps ${each_ref_file} \
                    --haps ${each_gname}.vcf.gz \
                    --start ${str_pos} \
                    --end ${end_pos} \
                    --window ${chunk_om} \
                    --format GT,DS,GP,HDS  \
                    --prefix ${each_gname/phased/imputed}_chunk${chunk} \
                    --chr ${chr} \
                    --map ${GMAP_FILE} \
                    --allTypedSites \
                    --minRatio 0.00001"
        done
    else
        echo -e "NOT Exists Phased file :"
        echo "  : ${each_gname}.vcf.gz"
        break
    fi
done