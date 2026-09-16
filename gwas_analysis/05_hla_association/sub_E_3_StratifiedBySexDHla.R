# ==============================================================================
# sub_B_3_StratifiedBySexDHla.R
#
#  CD: Apr 02 2024   K.Yamazaki
#  UD: Apr 09 2026
#  - Add QC condition (maf > 0.01)
#  UD: Jan 23 2025
#  - Change from sex_vlist(0:2) to sex_tlist("total", "male", "female")
#  - Add sex as covariate in all samples
#  - Add procedure to calculate frequency and MAF.
#  - Add normalize for two columns (sex as.factor, scale(age)) in each group
#  - Comment out the description to calculate ORs,
#       because ORs CANNOT be calculated essentially from the results of "lm".
#  UD: Dec 09 2024
#  - Add procedure to apply normalization to each sex-stratified group
#  UD: Aug 28 2024
#  - Change the environment from windows to linux on va-server
#  - Change procedure from to run this script to load this script in shell script with apptainer
#       -> Rename filename from "B_3_StratifiedBySexDHla.R" to "sub_B_3_StratifiedBySexDHla.R"
#  - Comment out description of logging
#  - Bug fix
#  UD: May 09 2024
#  - Change PHENOTYPE from "post_titer_log" to "post_titer_norm"
#  - Remove samples to be detectable antibody titer before anti SARS-CoV-2 vaccination
#  UD: Apr 03 2024   Unify files' name
#------------------------------------------------------------------------------
# Memo
#  Perform stratified analysis by sex to files exported by DEEP*HLA
#------------------------------------------------------------------------------
# ref)
# ==============================================================================
# [Argumnets]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "data.table", "R.utils", "broom",
    update = F
)

# (1) Environmental arguments
src_dir <- "~/analysis/COVID-19"
smpl_name <- "CUH-GWAS_230215"
geno_dname <- "PLINK_160223_1213"
step <- "B_3_StratifiedBySexDHla"

id_cname <- "iid"
sex_cname <- "sex"
age_cname <- "age"
# pheno_cname <- "post_titer_log"
pheno_cname <- "post_titer_norm"
r2_cfvalue <- 0.7
today <- as.character(format(Sys.time(), "%Y%m%d"))
sex_tlist <- c("total", "male", "female")
sig_p <- 5e-8

geno_dir <- file.path(src_dir, "genotype", geno_dname)
geno_dir_hla <- file.path(geno_dir, "04_imputation_hla")
imputed_fname <- file.path(geno_dir_hla, str_c(smpl_name, "_mhc_imputed"))
hla_pfile <- file.path(geno_dir, "etc", "jp.hla.txt")

dgeno_ext <- "dhla.txt"
dgeno_cv_file <- str_c(imputed_fname, ".", dgeno_ext)
dgeno_aa_file <- str_c(imputed_fname, ".aa.", dgeno_ext)
dgeno_gflist <- c(dgeno_cv_file, dgeno_aa_file)

wk_dir <- file.path(src_dir, "2308_impHLA")
wk_dir_hla <- file.path(wk_dir, "B_2_ReformatAssoc")
wk_dir_sex <- file.path(wk_dir, step)
if (!file.exists(wk_dir_sex)) dir.create(wk_dir_sex, recursive = T)

# pheno_file <- file.path(src_dir, "etc", str_c(smpl_name, ".txt"))
# dhla_pfile <- file.path(wk_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))
dhla_pfile <- file.path(src_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))
dinfo_cv_file <- file.path(
    wk_dir_hla, basename(dgeno_cv_file) %>% str_replace(dgeno_ext, "info")
)
dinfo_aa_file <- str_replace(dinfo_cv_file, "info$", "aa.info")

# log_file <- file.path(wk_dir, "script", "logs", str_c(step, "_", today, ".log"))

# (2) Load packages and in-house script
source("~/tools/script/GWAS/plink.tools.R")
# source("D:/Dropbox/Tools/script/GWAS/plink.tools.R")

# [Main] ----------------------------------------------------------------------
# (0) Logging
# sink(file = log_filed/)
cat("[Environmental argmunets]\n\n")
cat("- Root directory:\n\t>> ", src_dir, "\n", sep = "")
cat("- Genotype directory:\n\t>> ")
cat(str_replace(geno_dir, src_dir, "."), "\n", sep = "")
cat("\t- Reformatted dosage files:\n")
cat("\t\t- ", str_replace(dgeno_cv_file, geno_dir, "."), "\n")
cat("\t\t- ", str_replace(dgeno_aa_file, geno_dir, "."), "\n")
cat("- Working directory:\n\t>> ")
cat(str_replace(wk_dir, src_dir, "."), "\n", sep = "")
cat("\t- Phenotype:", pheno_cname, "\n")
cat("\t- Phenotype file:\n")
cat("\t\t>> ", str_replace(dhla_pfile, wk_dir, "."), "\n", sep = "")
cat("\t- Draft assoc directory:\n")
cat("\t\t>> ", str_replace(wk_dir_hla, wk_dir, "."), "\n", sep = "")
cat("\t- Information files related DEEP*HLA: \n")
cat("\t\t- ", str_replace(dinfo_cv_file, wk_dir, "."), "\n")
cat("\t\t- ", str_replace(dinfo_aa_file, wk_dir, "."), "\n")
cat("\t- Output directory:\n")
cat("\t\t>> ", str_replace(wk_dir_sex, wk_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (0) Check file whether to exist or not.
if (file.exists(dhla_pfile)) {
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F)
    cat("<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t-NOT EXIST PHENOTYPE FILE FOR DEEP*HLA !\n")
    cat("\t>>", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    stop()
}
cov_pclist <- str_subset(colnames(dhla_pdata), "PC\\d+")

for (dgeno_gfile in dgeno_gflist) {
    cat("\n[", basename(dgeno_gfile), "]\n", sep = "")
    # (1) File check and allele information
    cat("\n- File check and allele information\n")

    # a. Check and Load reformatted dosage file
    if (!file.exists(dgeno_gfile)) {
        cat("-NOT EXIST RE_FORMATTED DOSAGE FILE !\n")
        cat(">>", str_replace(dgeno_gfile, geno_dir, "."), "\n", sep = "")
        next
    }
    cat("<< ", str_replace(dgeno_gfile, geno_dir, "."), "\n", sep = "")
    dgeno_data <- fread(dgeno_gfile, header = T, showProgress = F)

    # b. Check and load HLA allele information file
    if (dgeno_gfile == dgeno_cv_file) {
        dinfo_file <- dinfo_cv_file
    } else if (dgeno_gfile == dgeno_aa_file) {
        dinfo_file <- dinfo_aa_file
    }
    if (file.exists(dinfo_file)) {
        dinfo_data <- fread(dinfo_file, header = T, showProgress = F)
        cat("<< ", str_replace(dinfo_file, src_dir, "."), "\n", sep = "")
    } else {
        cat("\t- NOT EXIST HLA-ALLELE INFORMATION FILE !\n")
        cat("\t>>", str_replace(dinfo_file, src_dir, "."), "\n", sep = "")
    }

    all_covlist <- NULL
    qc_flist <- NULL
    sig_vlist <- NULL

    for (i in seq_along(sex_tlist)) {
        sex_val <- i - 1
        sex_name <- str_c(sprintf("%02d", sex_val), "_", sex_tlist[i])

        all_covlist <- NULL
        wk_gdata <- dgeno_data %>%
            select(-c("bef_titer", "post_titer")) %>%
            mutate(across(all_of(sex_cname), as.factor))

        if (i == 1) {
            wk_gdata <- wk_gdata %>%
                mutate(across(all_of(age_cname), scale))
            all_covlist <- c(sex_cname, age_cname, cov_pclist)
        } else {
            ehla_pfile <- str_replace(dgeno_gfile, "\\.txt", str_c("_", sex_name, ".txt"))
            if (!file.exists(ehla_pfile)) {
                wk_gdata <- wk_gdata %>%
                    filter(sex == sex_val) %>%
                    mutate(across(all_of(age_cname), scale)) %>%
                    mutate(post_titer_norm = inv_nor_trans(post_titer_log))
                fwrite(wk_gdata, ehla_pfile, row.names = F, col.names = T, sep = "\t")
                cat("\t>> ", str_replace(ehla_pfile, src_dir, "."), "\n", sep = "")
            } else {
                wk_gdata <- fread(ehla_pfile, header = T, showProgress = F)
                cat("\t<< ", str_replace(ehla_pfile, src_dir, "."), "\n", sep = "")
            }
            wk_gdata <- wk_gdata %>% select(-sex)
            all_covlist <- c(age_cname, cov_pclist)
        }
        wk_gdata <- wk_gdata %>% select(-matches("_titer"), all_of(pheno_cname))

        # (2) Linear logistic regression by single marker
        cat("\n\t[", sex_name, "]\n", sep = "")
        cat("\n\t- Linear logistic regression by single marker\n")
        stat_data <- wk_gdata %>%
            pivot_longer(
                -all_of(c(id_cname, pheno_cname, all_covlist)),
                names_to = "genotype",
            ) %>%
            group_by(genotype) %>%
            nest() %>%
            mutate(
                freq = map_dbl(
                    data,
                    ~ sum(.x$value) / (nrow(.x) * 2)
                ),
                maf = if_else(freq > 0.5, 1 - freq, freq)
            ) %>%
            mutate(
                model = map(
                    data,
                    ~ lm(
                        str_c(
                            pheno_cname, "~ ", "value + ", str_flatten(all_covlist, collapse = "+")
                        ),
                        data = .,
                        na.action = na.omit
                    )
                )
            ) %>%
            select(-data) %>%
            # mutate(summary = list(tidy(model[[1]], conf.int = TRUE))) %>%
            mutate(
                results = map(
                    model,
                    ~ tidy(.x, conf.int = TRUE)
                )
            ) %>%
            unnest(results) %>%
            select(-model) %>%
            # mutate_at(vars(estimate, conf.low, conf.high), exp) %>%
            # rename(
            #     OR = estimate,
            #     L95 = conf.low,
            #     U95 = conf.high
            # ) %>%
            filter(term != "(Intercept)") %>%
            # select(-std.error) %>%
            # select(genotype:term, statistic, everything()) %>%
            relocate(p.value, .after = last_col())

        # (3) Export statistic results and apply QC
        cat("\n\t- Export statistic results and apply QC\n")
        stat_file <- basename(dgeno_gfile) %>%
            str_remove("_imputed") %>%
            str_replace(str_c(".", dgeno_ext), str_c(".", sex_name, "_rawlm.txt")) %>%
            file.path(wk_dir_sex, .)
        stat_data <- dinfo_data %>%
            select(-freq) %>%
            right_join(stat_data, by = c("hla_name" = "genotype"))
        fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
        system(paste("gzip -f", stat_file))
        cat("\t>> ", str_replace(stat_file, src_dir, "."), ".gz\n", sep = "")

        qc_file <- str_replace(stat_file, "_rawlm.txt", "_qclm.txt")
        qc_data <- stat_data %>%
            filter(term == "value") %>%
            select(-term) %>%
            filter(
                r2 >= 0.7 & !is.na(r2) & maf > 0.01
                # r2 >= 0.7 & !is.na(r2)
                # if (dhla_gfile == dhla_aa_file) {
                #     r2 >= 0.7 & !is.na(r2)
                # } else {
                #     !is.na(r2)
                # }
            )
        fwrite(qc_data, qc_file, row.names = F, col.names = T, sep = "\t")
        cat("\t>> ", str_replace(qc_file, src_dir, "."), "\n", sep = "")
        qc_flist <- c(qc_flist, qc_file)

        buf <- qc_data %>%
            filter(p.value < sig_p) %>%
            pull(hla_name)
        sig_vlist <- c(sig_vlist, buf)
    }

    # (4) Summarize statistics with significant association
    cat("\n- Summarize statistics with significant association\n")

    # a. Extract data from each file applied QC
    cat("\n\t- Extract data from each file applied QC\n")
    sig_vlist <- sig_vlist %>%
        sort() %>%
        unique()
    cat("\t- No. of significant variants: ", length(sig_vlist), "\n", sep = "")

    geno_fname <- basename(dgeno_gfile) %>%
        str_remove("_imputed") %>%
        str_remove(str_c(".", dgeno_ext))

    sig_data <- data.frame()
    for (qc_file in qc_flist) {
        sex_name <- qc_file %>%
            basename() %>%
            str_replace(str_c(geno_fname, "\\.(\\d{2})_(.*)_qclm.txt"), "\\1_\\2")
        buf <- fread(qc_file, header = T) %>%
            mutate(study = sex_name) %>%
            filter(hla_name %in% sig_vlist) %>%
            select(-c(maf, sensitivity:NPV))
        sig_data <- bind_rows(sig_data, buf)

        cat("\t<< ", str_replace(qc_file, src_dir, "."), "\n", sep = "")
    }

    # b. Reformat and add minimal p-value
    cat("\n\t- Reformat and add minimal p-value\n")

    buf <- sig_data %>%
        pivot_wider(
            names_from = study, values_from = c(freq:p.value),
            names_sep = ":", names_vary = "slowest"
        ) %>%
        mutate(min.p = pmap_dbl(dplyr::select(., starts_with("p.value:")), pmin, na.rm = T)) %>%
        dplyr::select(hla_name:r2, min.p, everything())
    colnames(buf) <- colnames(buf) %>%
        str_replace("(.*):(\\d{2}_.*)", "\\2:\\1")

    sig_file <- file.path(wk_dir_sex, str_c(geno_fname, "_siglm.txt"))
    fwrite(buf, sig_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(sig_file, wk_dir, "."), "\n\n", sep = "")
}

# sink()