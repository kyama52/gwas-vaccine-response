# ==============================================================================
# Sub routine script for preparing phenotype files for stratified GWAS by sex
# ==============================================================================
# [Argumnets]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "data.table", "R.utils",
    update = F
)

# (1) Environmental arguments
src_dir <- "~/analysis/COVID-19/2609_PrepGitRepo"
smpl_name <- "covid-vac"
etc_dir <- file.path(src_dir, "etc")

geno_dir <- file.path(src_dir, "genotype", "02_imputation")
fam_file <- file.path(geno_dir, str_c(smpl_name, ".fam"))

pca_file <- file.path(
    src_dir, "A_2_CleanDataAfterQc", "01_repca", str_c(smpl_name, "_csct.eigenvec")
)

pheno_org_file <- file.path(etc_dir, str_c(smpl_name, ".tsv"))
pheno_file <- file.path(etc_dir, str_c(smpl_name, ".txt"))
cov_file <- file.path(etc_dir, str_c(smpl_name, "_plink.covar"))

sex_tlist <- c("00_total", "01_male", "02_female")

# (2) Load packages and in-house script
source("~/tools/script/GWAS/plink.tools.R")

# [Main] ----------------------------------------------------------------------
# (1) Check distribution phenotype and covariate files
message("- Check distribution phenotype and covariate files\n")
if (!file.exists(fam_file)) {
    message("\tNOT exist FAM file!:\n")
    fam_qflist <- list.files(
        geno_dir, str_c(smpl_name, "_chr(\\d+|X|Y)_imputed_qc.fam"),
        recursive = T, full.names = T
    ) %>% sort()
    if (length(fam_qflist)) {
        system(paste("cp", fam_qflist[1], fam_file))
        fam_file <- fam_qflist[1] %>%
            str_replace("(.*03_qc_imputed)", geno_dir)
        message("\t<< ", str_replace(fam_file, src_dir, "."), "\n")
    } else {
        stop("\t>> ", str_replace(fam_file, src_dir, "."), "\n")
    }
}
fam_data <- fread(fam_file, header = F) %>%
    select(-ncol(.)) %>%
    setNames(c("FID", "IID", "FATID", "MATID", "sex"))
message("\tFAM_FILE:")
message("\t<< ", str_replace(fam_file, src_dir, "."), "\n")

if (!file.exists(pca_file)) {
    message("\tNOT exist PCA file!:\n")
    stop("\t>> ", str_replace(pca_file, src_dir, "."), "\n")
}
message("\tPCA_FILE:")
pca_data <- fread(pca_file, header = T, sep = "\t")
message("\t<< ", str_replace(pca_file, src_dir, "."), "\n")

if (!file.exists(pheno_org_file)) {
    message("\tNOT exist ORIGINAL PHENOTYPE file!:")
    stop("\t>> ", str_replace(pheno_org_file, src_dir, "."), "\n")
}
message("\tORG_PHENO_FILE:")
pheno_odata <- fread(pheno_org_file, header = T, sep = "\t")
message("\t<< ", str_replace(pheno_org_file, src_dir, "."), "\n")

# (2) Prepare phenotype and covariate files stratified by sex
message("- Prepare phenotype and covariate files stratified by sex\n")
base_data <- fam_data %>%
    left_join(
        pheno_odata %>% select(IID, age, post_titer_log),
        by = "IID"
    ) %>%
    left_join(
        pca_data %>% select(IID, starts_with("PC")),
        by = "IID"
    )

if (!file.exists(cov_file)) {
    cov_data <- base_data %>%
        select(FID, IID, age, starts_with("PC")) %>%
        rename(`#FID` = FID)
    fwrite(cov_data, cov_file, row.names = FALSE, col.names = TRUE, sep = " ")
    message("\t>> ", str_replace(cov_file, src_dir, "."))
} else {
    message("\t[EXIST] ", str_replace(cov_file, src_dir, "."))
}

for (i in seq_along(sex_tlist)) {
    sex_type <- sex_tlist[i]
    epheno_file <- str_replace(pheno_file, "\\.txt", str_c("_", sex_type, ".txt"))

    if (!file.exists(epheno_file)) {
        epheno_data <- base_data
        if (sex_type == "01_male") {
            epheno_data <- epheno_data %>% filter(sex == 1)
        } else if (sex_type == "02_female") {
            epheno_data <- epheno_data %>% filter(sex == 2)
        }

        epheno_data <- epheno_data %>%
            mutate(post_titer_norm = inv_nor_trans(post_titer_log)) %>%
            select(-post_titer_log)
        fwrite(epheno_data, epheno_file, row.names = F, col.names = T, sep = "\t")
        message("\t>> ", str_replace(epheno_file, src_dir, "."))
    } else {
        message("\t[EXIST] ", str_replace(epheno_file, src_dir, "."))
    }
}