#------------------------------------------------------------------------------
# Subroutine for reformatting exported files by DEEP*HLA
#------------------------------------------------------------------------------
# [Argumnets]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "data.table", "R.utils", "broom",
    update = F
)

# (1) Environmental arguments
src_dir <- "~/analysis/COVID-19/2609_PrepGitRepo"
smpl_name <- "covid-vac"

geno_dir <- file.path(src_dir, "genotype", "04_imputation_hla")
geno_dir_etc <- file.path(geno_dir, "etc")
phased_fname <- file.path(geno_dir, str_c(smpl_name, "_hla_phased"))
imputed_fname <- file.path(geno_dir, str_c(smpl_name, "_hla_imputed"))

fam_file <- str_c(phased_fname, ".fam")
dhla_ext <- "deephla.dosage"
dhla_cv_file <- str_c(imputed_fname, ".", dhla_ext)
dhla_aa_file <- str_c(imputed_fname, ".aa.", dhla_ext)
dhla_gflist <- c(dhla_cv_file, dhla_aa_file)

sex_tlist <- c("total", "male", "female")

hla_pfile <- file.path(geno_dir_etc, "jp.hla.txt")
acc_ext <- "accuracy.txt"
acc_cv_file <- file.path(geno_dir_etc, str_c(smpl_name, "_hla.cv.", acc_ext))
acc_aa_file <- file.path(geno_dir_etc, str_c(smpl_name, "_hla.cv.aa.", acc_ext))

etc_dir <- file.path(src_dir, "etc")
dhla_pfile <- file.path(etc_dir, str_c(smpl_name, "_dhla_pheno.txt"))
pheno_file <- file.path(etc_dir, str_c(smpl_name, ".txt"))
cov_file <- file.path(etc_dir, str_c(smpl_name, "_plink.covar"))
rmv_sfile <- file.path(etc_dir, "rmv_smpl_expvirus.txt")

# (2) Load packages and in-house script
tool_dir <- file.path(src_dir, "scripts", "gwas_analysis", "00_utils")
source(file.path(tool_dir, "gwas_tools.R"))

# [Main] ----------------------------------------------------------------------
# (0) Logging
cat("[Environmental argmunets]\n\n")
cat("- Root directory:\n\t>> ", src_dir, "\n", sep = "")
cat("- Genotype directory:\n\t>> ")
cat(str_replace(geno_dir, src_dir, "."), "\n", sep = "")
cat("\t- Dosage files:\n")
cat("\t\t- ", str_replace(dhla_cv_file, geno_dir, "."), "\n")
cat("\t\t- ", str_replace(dhla_aa_file, geno_dir, "."), "\n")
cat("\t- Fam file:", str_replace(fam_file, geno_dir, "."), "\n")
cat("\t- Accuracy calculated by DEEP*HLA: \n")
cat("\t\t- ", str_replace(acc_cv_file, geno_dir, "."), "\n")
cat("\t\t- ", str_replace(acc_aa_file, geno_dir, "."), "\n")
cat("\t- Remove sample file:\n")
cat("\t\t>> ", str_replace(rmv_sfile, src_dir, "."), "\n")
cat("\t- Phenotype file:\n")
cat("\t\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (0) Check file whether to exist or not.
cat("- Check file whether to exist or not\n")
if (!file.exists(fam_file)) {
    cat("\t-NOT EXIST FAM FILE !\n")
    cat("\t>>", str_replace(fam_file, geno_dir, "."), "\n", sep = "")
    stop("Required input file(s) not found.")
} else {
    fam_data <- fread(fam_file, header = F, showProgress = F) %>%
        magrittr::set_colnames(c("fid", "iid", "", "", "sex", "phenotype")) %>%
        select(1, 2, 5)
    cat("<< ", str_replace(fam_file, src_dir, "."), "\n", sep = "")
}

if (!file.exists(hla_pfile)) {
    cat("\t-NOT EXIST HLA_POS FILE !\n")
    cat("\t>>", str_replace(hla_pfile, geno_dir, "."), "\n", sep = "")
    stop("Required input file(s) not found.")
} else {
    hla_pdata <- fread(hla_pfile, header = F, showProgress = F) %>%
        rename("name" = 1, "pos" = 2)
    cat("<< ", str_replace(hla_pfile, src_dir, "."), "\n", sep = "")
}

if (!file.exists(dhla_pfile)) {
    cat("\t-NOT EXIST PHENOTYPE FILE FOR DEEP*HLA !\n")
    cat("\t>>", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    if (file.exists(pheno_file) && file.exists(cov_file)) {
        pheno_data <- fread(pheno_file, header = T, showProgress = F) %>%
            select(-1) %>%
            rename("iid" = "IID")
        cov_data <- fread(cov_file, header = T, showProgress = F) %>%
            select(-1) %>%
            rename("iid" = "IID") %>%
            select(iid, starts_with(("PC")))
        dhla_pdata <- pheno_data %>%
            left_join(cov_data, by = "iid")
        fwrite(dhla_pdata, dhla_pfile, row.names = F, col.names = T, sep = "\t")
        cat("\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    } else {
        if (!file.exists(pheno_file)) {
            cat("\t- NOT EXIST PHENO_FILE !\n")
            cat("\t >> ", str_replace(pheno_file, src_dir, "."), "\n", sep = "")
        }
        if (!file.exists(cov_file)) {
            cat("\t- NOT EXIST COV_FILE !\n")
            cat("\t >> ", str_replace(cov_file, src_dir, "."), "\n", sep = "")
        }
        stop("Required input file(s) not found.")
    }
} else {
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F)
    cat("<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
}
cov_nlist <- colnames(dhla_pdata) %>% str_subset("^PC\\d")

rmv_slist <- NULL
if (file.exists(rmv_sfile)) {
    rmv_slist <- fread(rmv_sfile, header = T, sep = "\t") %>% pull(2)
    cat("<< ", str_replace(rmv_sfile, src_dir, "."), "\n", sep = "")
}

# (1) File check and allele information
for (dhla_gfile in dhla_gflist) {
    cat("\n- File check and allele information\n")
    if (!file.exists(dhla_gfile)) {
        cat("- NOT EXIST DOSAGE FILE !\n")
        cat(">>", str_replace(dhla_gfile, geno_dir, "."), "\n", sep = "")
        next
    }
    cat("<< ", str_replace(dhla_gfile, geno_dir, "."), "\n", sep = "")
    dhla_data <- fread(dhla_gfile, header = F, showProgress = F) %>%
        rename("hla_name" = 1, "al_1" = 2, "al_2" = 3)
    hinfo_data <- dhla_data %>% select(seq_len(3))

    # a. Extract information from dosaga.file
    cat("\n- Extract information from dosaga.file\n")
    hinfo_file <- file.path(
        etc_dir, str_replace(dhla_gfile, dhla_ext, "info") %>% basename()
    )

    if (dhla_gfile == dhla_cv_file) {
        hinfo_data <- hinfo_data %>%
            mutate(
                name =
                    str_replace(hla_name, "(HLA_.*)_(.*)", "\\1") %>%
                        str_replace("(MICA|MICB|TAP\\d)_(.*)", "\\1")
            ) %>%
            left_join(hla_pdata, by = "name") %>%
            rename("gene" = "name")

        acc_file <- acc_cv_file
    } else if (dhla_gfile == dhla_aa_file) {
        aa_ptn <- "^(AA|INS)_(.*)_(-\\d+|\\d+|\\d+\\w\\d+)_(\\d{8}|\\d{8}_\\W+)"
        hinfo_data <- hinfo_data %>%
            mutate(type = str_replace(hla_name, aa_ptn, "\\1")) %>%
            mutate(gene = str_replace(hla_name, aa_ptn, "\\2")) %>%
            mutate(aa_pos = str_replace(hla_name, aa_ptn, "\\3")) %>%
            mutate(new_aa = if_else(
                str_detect(aa_pos, "_"), str_replace(aa_pos, "(.*)_(.*)", "\\2"),
                ""
            )) %>%
            mutate(pos = str_replace(hla_name, aa_ptn, "\\4")) %>%
            mutate_at(
                vars(type:aa_pos, pos),
                ~ ifelse(str_detect(.x, "_"), str_replace(.x, "(.*)_(.*)", "\\1"), .x)
            )

        acc_file <- acc_aa_file
    }

    # (2) Load accuracy information
    cat("\n- Load accuracy information\n")
    acc_data <- data.frame()
    if (file.exists(acc_file)) {
        acc_data <- fread(acc_file, header = T, showProgress = F)
        cat("\t<< ", str_replace(acc_file, src_dir, "."), "\n", sep = "")
    } else {
        cat("\t-NOT EXIST ACCURACY FILE !\n")
        cat("\t>>", str_replace(acc_file, src_dir, "."), "\n", sep = "")
    }
    if (nrow(acc_data) > 0) {
        hinfo_data <- hinfo_data %>%
            left_join(acc_data, by = c("hla_name" = "id"))
    }
    if (nrow(hinfo_data) > 0) {
        hinfo_data %>% fwrite(hinfo_file, row.names = F, col.names = T, sep = "\t")
        cat("\t>> ", str_replace(hinfo_file, src_dir, "."), "\n", sep = "")
    }

    # (3) Convert dosage file exported by DEEP*HLA
    cat("- Convert dosage file exported by DEEP*HLA\n")

    geno_data <- dhla_data %>%
        select(-c(2, 3)) %>%
        column_to_rownames("hla_name") %>%
        as.matrix() %>%
        t() %>%
        as.data.table()
    rownames(geno_data) <- fam_data %>% pull(iid)
    geno_data <- dhla_pdata %>%
        right_join(geno_data %>% rownames_to_column("iid"), by = "iid")
    if (length(rmv_slist)) {
        cat("\n- Apply sample QC which remove samples with detectable bef_titer\n")
        cat("\t- Samples removed: ", length(rmv_slist), "\n", sep = "")
        geno_data <- geno_data %>% filter(!(iid %in% rmv_slist))
    }

    for (i in seq_along(sex_tlist)) {
        sex_val <- i - 1
        sex_name <- str_c(sprintf("%02d", sex_val), "_", sex_tlist[i])
        geno_dfile <- str_replace(dhla_gfile, dhla_ext, str_c("dhla_", sex_name, ".txt"))

        if (!file.exists(geno_dfile)) {
            if (i == 1) {
                wk_gdata <- geno_data %>% mutate(age = scale(age))
            } else {
                wk_gdata <- geno_data %>%
                    filter(sex == sex_val) %>%
                    mutate(age_cname = scale(age)) %>%
                    mutate(post_titer_norm = inv_nor_trans(post_titer_log))
            }

            fwrite(wk_gdata, geno_dfile, row.names = F, col.names = T, sep = "\t")
            cat("\t>> ", str_replace(geno_dfile, src_dir, "."), "\n", sep = "")
        } else {
            cat("\t[EXIST] ", str_replace(geno_dfile, src_dir, "."), "\n", sep = "")
        }
    }
}