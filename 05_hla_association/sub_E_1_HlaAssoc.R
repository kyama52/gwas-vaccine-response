#------------------------------------------------------------------------------
#  Subroutine for performing sex-stratified analysis to exported files by DEEP*HLA
#  [antibody titer]
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
step <- "E_1_HlaAssoc"

pheno_cname <- "post_titer_norm"
freq_val <- 0.01
r2_val <- 0.7
today <- as.character(format(Sys.time(), "%Y%m%d"))
sex_tlist <- c("total", "male", "female")
sig_p <- 5e-8

geno_dir <- file.path(src_dir, "genotype", "04_imputation_hla")

imputed_fname <- str_c(smpl_name, "_hla_imputed")
dgeno_cv_fname <- file.path(geno_dir, str_c(imputed_fname, ".dhla"))
dgeno_aa_fname <- file.path(geno_dir, str_c(imputed_fname, ".aa.dhla"))
dgeno_flist <- c(dgeno_cv_fname, dgeno_aa_fname)

wk_dir <- file.path(src_dir, step)
if (!file.exists(wk_dir)) dir.create(wk_dir, recursive = T)

etc_dir <- file.path(src_dir, "etc")
dhla_pfile <- file.path(etc_dir, str_c(smpl_name, "_dhla_pheno.txt"))

# [Main] ----------------------------------------------------------------------
# (0) Logging
cat("[Environmental argmunets]\n\n")
cat("- Root directory:\n\t>> ", src_dir, "\n", sep = "")
cat("- Genotype directory:\n\t>> ")
cat(str_replace(geno_dir, src_dir, "."), "\n", sep = "")
cat("- Phenotype:", pheno_cname, "\n")
cat("- Phenotype file:\n")
cat("\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
cat("- Output directory:\n")
cat("\t>> ", str_replace(wk_dir, src_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (0) Check file whether to exist or not.
if (file.exists(dhla_pfile)) {
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F)
    cat("<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t-NOT EXIST PHENOTYPE FILE FOR DEEP*HLA !\n")
    cat("\t>>", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    stop("Required input file(s) not found.")
}
cov_pclist <- str_subset(colnames(dhla_pdata), "PC\\d+")

for (dgeno_fname in dgeno_flist) {
    # (1) File check and allele information
    cat("\n[", basename(dgeno_fname), "]\n", sep = "")
    cat("\n- File check and allele information\n")

    # a. Check and load HLA allele information file
    if (dgeno_fname == dgeno_cv_fname) {
        dinfo_file <- file.path(etc_dir, str_c(imputed_fname, ".info"))
    } else if (dgeno_fname == dgeno_aa_fname) {
        dinfo_file <- file.path(etc_dir, str_c(imputed_fname, ".aa.info"))
    }
    if (file.exists(dinfo_file)) {
        dinfo_data <- fread(dinfo_file, header = T, showProgress = F)
        cat("<< ", str_replace(dinfo_file, src_dir, "."), "\n", sep = "")
    } else {
        cat("\t- NOT EXIST HLA-ALLELE INFORMATION FILE !\n")
        cat("\t>>", str_replace(dinfo_file, src_dir, "."), "\n", sep = "")
    }

    qc_flist <- NULL
    sig_vlist <- NULL

    for (i in seq_along(sex_tlist)) {
        sex_name <- str_c(sprintf("%02d", i - 1), "_", sex_tlist[i])
        cat("\n\t[", sex_name, "]\n", sep = "")

        # a. Check and Load reformatted dosage file
        cat("\n\t- Check and Load reformatted dosage file\n")
        dgeno_file <- str_c(dgeno_fname, "_", sex_name, ".txt")
        if (!file.exists(dgeno_file)) {
            cat("-NOT EXIST RE_FORMATTED DOSAGE FILE !\n")
            cat(">>", str_replace(dgeno_file, geno_dir, "."), "\n", sep = "")
            next
        }
        geno_gdata <- fread(dgeno_file, header = T, showProgress = F) %>%
            select(-c("bef_titer", "post_titer")) %>%
            mutate(sex = as.factor(sex))
        cat("\t<< ", str_replace(dgeno_file, geno_dir, "."), "\n", sep = "")

        if (i == 1) {
            all_covlist <- c("sex", "age", cov_pclist)
            wk_gdata <- geno_gdata
        } else {
            all_covlist <- c("age", cov_pclist)
            wk_gdata <- geno_gdata %>% select(-sex)
        }
        wk_gdata <- wk_gdata %>% select(-matches("_titer"), all_of(pheno_cname))

        # (2) Linear logistic regression by single marker
        cat("\n\t- Linear logistic regression by single marker\n")
        stat_data <- wk_gdata %>%
            pivot_longer(
                -all_of(c("iid", pheno_cname, all_covlist)),
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
            mutate(
                results = map(
                    model,
                    ~ tidy(.x, conf.int = TRUE)
                )
            ) %>%
            unnest(results) %>%
            select(-model) %>%
            filter(term != "(Intercept)") %>%
            relocate(p.value, .after = last_col())

        # (3) Export statistic results and apply QC
        cat("\n\t- Export statistic results and apply QC\n")
        stat_file <- file.path(wk_dir, str_c(smpl_name, "_hla.", sex_name, "_rawlm.txt"))
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
            filter(r2 >= r2_val & !is.na(r2) & maf > freq_val)
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

    sig_fname <- basename(dgeno_fname) %>% str_remove_all("_imputed|\\.dhla")
    sig_data <- data.frame()
    for (qc_file in qc_flist) {
        sex_name <- qc_file %>%
            basename() %>%
            str_replace(str_c(sig_fname, "\\.(\\d{2})_(.*)_qclm.txt"), "\\1_\\2")
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

    sig_file <- file.path(wk_dir, str_c(sig_fname, "_siglm.txt"))
    fwrite(buf, sig_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(sig_file, wk_dir, "."), "\n\n", sep = "")
}