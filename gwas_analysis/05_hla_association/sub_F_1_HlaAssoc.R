#------------------------------------------------------------------------------
#  Subroutine for performing sex-stratified analysis to exported files by DEEP*HLA
#  [adverse reaction]
#------------------------------------------------------------------------------
# [Argumnets]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "data.table", "R.utils", "broom", "conflicted",
    update = F
)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# (1) Environmental arguments
src_dir <- "~/analysis/COVID-19/2609_PrepGitRepo"
smpl_name <- "covid-vac"
step <- "F_1_HlaAssoc"

freq_val <- 0.01
r2_val <- 0.7
today <- as.character(format(Sys.time(), "%Y%m%d"))
# sig_p <- 5e-8
sig_p <- 1e-5
sex_tlist <- c("total", "male", "female")

geno_dir <- file.path(src_dir, "genotype", "04_imputation_hla")

imputed_fname <- str_c(smpl_name, "_hla_imputed")
dgeno_cv_fname <- file.path(geno_dir, str_c(imputed_fname, ".dhla"))
dgeno_aa_fname <- file.path(geno_dir, str_c(imputed_fname, ".aa.dhla"))
dgeno_flist <- c(dgeno_cv_fname, dgeno_aa_fname)

etc_dir <- file.path(src_dir, "etc")
dhla_pfile <- file.path(etc_dir, str_c(smpl_name, "_seff.txt"))

wk_dir <- file.path(src_dir, step)
if (!file.exists(wk_dir)) dir.create(wk_dir, recursive = T)

# [Main] ----------------------------------------------------------------------
# (0) Logging
cat("[Environmental argmunets]\n\n")
cat("- Root directory:\n\t>> ", src_dir, "\n", sep = "")
cat("- Genotype directory:\n\t>> ")
cat(str_replace(geno_dir, src_dir, "."), "\n", sep = "")
cat("- Phenotype file:\n")
cat("\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
cat("- Output directory:\n")
cat("\t>> ", str_replace(src_dir, wk_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (0) Check file whether to exist or not.
cat("- Check file whether to exist or not\n")

if (file.exists(dhla_pfile)) {
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F) %>%
        mutate(across(c(starts_with("seff_")), as.factor))
    cat("<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t-NOT EXIST PHENOTYPE FILE !\n")
    cat("\t>>", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    stop("Required input file(s) not found.")
}

cov_pclist <- str_subset(colnames(dhla_pdata), "PC\\d+")

for (dgeno_fname in dgeno_flist) {
    # (1) File check and allele information
    cat("\n-----------------------------------------------------------\n")
    cat("\n[", basename(dgeno_fname), "]\n", sep = "")
    cat("-----------------------------------------------------------\n\n")

    cat("\n- File check and allele information\n")

    #  Check and load HLA allele information file
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

    stat_fname <- basename(dgeno_fname) %>% str_remove_all("_imputed|\\.dhla$")

    for (i in seq_along(sex_tlist)) {
        sex_name <- str_c(sprintf("%02d", i - 1), "_", sex_tlist[i])
        cat("\n\t[", sex_name, "]\n", sep = "")

        stat_file <- file.path(wk_dir, str_c(stat_fname, ".", sex_name, "_rawlm.txt"))
        qc_file <- str_replace(stat_file, "_rawlm.txt", "_qclm.txt")

        if (file.exists(qc_file)) {
            cat("- ALREADY FINISHED THIS STEP !\n")
            cat(">>", str_replace(qc_file, src_dir, "."), "\n", sep = "")
            next
        }

        # (2) Load reformatted dosage file
        cat("\n- Load reformatted dosage file\n")
        dgeno_file <- str_c(dgeno_fname, "_", sex_name, ".txt")
        if (!file.exists(dgeno_file)) {
            cat("- NOT EXIST RE_FORMATTED DOSAGE FILE !\n")
            cat(">>", str_replace(dgeno_file, src_dir, "."), "\n", sep = "")
            next
        }

        dgeno_data <- fread(dgeno_file, header = T, showProgress = F) %>%
            select(-contains("titer")) %>%
            left_join(
                dhla_pdata %>% select(IID, starts_with("seff_")),
                by = setNames("IID", "iid")
            ) %>%
            mutate(sex = as.factor(sex))
        cat("<< ", str_replace(dgeno_file, geno_dir, "."), "\n", sep = "")

        if (i == 1) {
            all_covlist <- c("sex", "age", cov_pclist)
            dgeno_data <- dgeno_data
        } else {
            all_covlist <- c("age", cov_pclist)
            dgeno_data <- dgeno_data %>% select(-sex)
        }
        dgeno_data <- dgeno_data %>% mutate(age = scale(age))

        # (2) Linear logistic regression by single marker
        cat("\n- Linear logistic regression by single marker\n")
        pheno_clist <- str_subset(colnames(dgeno_data), "^seff_")
        geno_clist <- colnames(dgeno_data) %>%
            setdiff(c("iid", all_covlist, pheno_clist))

        wk_data <- dgeno_data %>%
            pivot_longer(
                cols = all_of(geno_clist),
                names_to = "genotype",
                values_to = "value"
            ) %>%
            pivot_longer(
                cols = all_of(pheno_clist),
                names_to = "phenotype",
                values_to = "pheno_value"
            ) %>%
            group_by(phenotype, genotype) %>%
            nest() %>%
            mutate(
                freq = map_dbl(
                    data,
                    ~ sum(.x$value) / (nrow(.x) * 2)
                ),
                maf = if_else(freq > 0.5, 1 - freq, freq)
            )

        stat_data <- wk_data %>%
            mutate(
                model = map2(
                    data,
                    maf,
                    ~ {
                        if (.y <= freq_val) {
                            return(NULL)
                        }
                        tryCatch(
                            {
                                formula <- as.formula(
                                    str_c("pheno_value ~ value + ", paste(all_covlist, collapse = " + "))
                                )
                                MASS::polr(formula, data = .x, Hess = T)
                            },
                            error = function(e) {
                                message("Error in model fitting: ", conditionMessage(e))
                                NULL
                            }
                        )
                    }
                )
            ) %>%
            select(-data)

        buf <- stat_data %>%
            mutate(
                results = map(
                    model,
                    ~ {
                        if (is.null(.x)) {
                            return(NULL)
                        }
                        tryCatch(
                            {
                                tbl <- tidy(.x, conf.int = TRUE, exponentiate = TRUE)
                                if (!nrow(tbl)) {
                                    return(NULL)
                                }
                                coef_info <- coef(summary(.x))
                                tbl <- tbl %>%
                                    mutate(
                                        p.value = pnorm(abs(coef_info[term, "t value"]), lower.tail = FALSE) * 2
                                    ) %>%
                                    select(-coef.type)
                                return(tbl)
                            },
                            error = function(e) {
                                message("Error in model fitting: ", conditionMessage(e))
                                NULL
                            }
                        )
                    }
                )
            )
        na_tbl <- pull(buf, results) %>%
            discard(is.null) %>%
            .[[1]] %>%
            slice(1) %>%
            mutate(across(everything(), ~NA))
        res_data <- buf %>%
            select(-model) %>%
            mutate(
                results = map(
                    results,
                    ~ {
                        if (is.null(.x)) {
                            return(na_tbl)
                        } else {
                            .x
                        }
                    }
                )
            ) %>%
            unnest(results)

        # (3) Export statistic results and apply QC
        cat("\n- Export statistic results and apply QC\n")
        stat_data <- dinfo_data %>%
            select(-freq) %>%
            right_join(res_data, by = c("hla_name" = "genotype")) %>%
            select(phenotype, hla_name:pos, freq, maf, everything())
        if ("p.value" %in% colnames(stat_data)) {
            stat_data <- stat_data %>% relocate(p.value, .after = last_col())
        }
        fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
        system(paste("gzip -f", stat_file))
        cat(">> ", str_replace(stat_file, src_dir, "."), ".gz\n", sep = "")

        qc_data <- stat_data %>%
            filter(term == "value") %>%
            select(-term) %>%
            filter(r2 >= r2_val & !is.na(r2))
        # filter(r2 >= r2_val & !is.na(r2) & maf > freq_val)
        fwrite(qc_data, qc_file, row.names = F, col.names = T, sep = "\t")
        cat(">> ", str_replace(qc_file, src_dir, "."), "\n\n", sep = "")
        qc_flist <- c(qc_flist, qc_file)

        buf <- qc_data %>%
            filter(p.value < sig_p) %>%
            select(phenotype, hla_name)
        sig_vlist <- bind_rows(sig_vlist, buf)
    }

    # (4) Summarize statistics with significant association
    cat("\n- Summarize statistics with significant association\n")

    # a. Extract data from each file applied QC
    cat("\n- Extract data from each file applied QC\n")

    if (!length(qc_flist)) {
        qc_flist <- list.files(
            wk_dir, str_c(stat_fname, "\\.\\d{2}_.*_qclm.txt"),
            recursive = T, full.names = T
        )
        for (qc_file in qc_flist) {
            buf <- fread(qc_file, header = T) %>%
                filter(p.value < sig_p) %>%
                select(phenotype, hla_name)
            sig_vlist <- bind_rows(sig_vlist, buf)
        }
    }

    if (!nrow(sig_vlist)) {
        cat("\t\t>> NOT exist variants satisfied p.value <= ", sig_p, "\n", sep = "")
        next
    }
    sig_vlist <- sig_vlist %>% distinct()
    cat("\t- No. of significant variants: ", nrow(sig_vlist), "\n", sep = "")

    sig_data <- data.frame()
    for (qc_file in qc_flist) {
        sex_name <- qc_file %>%
            basename() %>%
            str_replace(str_c(stat_fname, "\\.(\\d{2})_(.*)_qclm.txt"), "\\1_\\2")
        buf <- fread(qc_file, header = T) %>%
            mutate(study = sex_name) %>%
            right_join(sig_vlist, join_by(phenotype, hla_name)) %>%
            select(-c(maf, sensitivity:NPV)) %>%
            select(phenotype, hla_name:pos, r2, study, everything())

        if ("p.value" %in% colnames(buf)) {
            buf <- buf %>% relocate(all_of(c("p.value", "study")), .after = last_col())
        }
        sig_data <- bind_rows(sig_data, buf)

        cat("\t<< ", str_replace(qc_file, src_dir, "."), "\n", sep = "")
    }

    # b. Reformat and add minimal p-value
    cat("\n\t- Reformat and add minimal p-value\n")

    buf <- sig_data %>%
        pivot_wider(
            names_from = study,
            values_from = c(freq:p.value),
            names_sep = ":",
            names_vary = "slowest",
            names_glue = "{study}:{.value}"
        ) %>%
        mutate(min.p = pmap_dbl(dplyr::select(., ends_with("p.value")), pmin, na.rm = T)) %>%
        dplyr::select(phenotype, hla_name:r2, min.p, everything())

    sig_file <- file.path(wk_dir, str_c(stat_fname, "_siglm.txt"))
    fwrite(buf, sig_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(sig_file, wk_dir, "."), "\n\n", sep = "")
}