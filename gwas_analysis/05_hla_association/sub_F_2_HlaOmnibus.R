#------------------------------------------------------------------------------
# Subroutine for rerforming sex-stratified omnibus analysis of HLA amino acid positions
#  [adverse reaction]
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
step <- "F_2_HlaOmnibus"

freq_val <- 0.01
r2_val <- 0.7
today <- as.character(format(Sys.time(), "%Y%m%d"))
# sig_p <- 5e-8
sig_p <- 1e-5
sex_tlist <- c("total", "male", "female")

geno_dir <- file.path(src_dir, "genotype", "04_imputation_hla")

imputed_fname <- str_c(smpl_name, "_hla_imputed")
dgeno_aa_fname <- file.path(geno_dir, str_c(imputed_fname, ".aa.dhla"))

etc_dir <- file.path(src_dir, "etc")
dhla_pfile <- file.path(etc_dir, str_c(smpl_name, "_seff.txt"))
dinfo_aa_file <- file.path(etc_dir, str_c(imputed_fname, ".aa.info"))

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
cat("- Information files related DEEP*HLA: \n")
cat("\t- ", str_replace(dinfo_aa_file, src_dir, "."), "\n")
cat("- Output directory:\n")
cat("\t>> ", str_replace(wk_dir, src_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (1) Check file whether to exist or not.
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
if (file.exists(dinfo_aa_file)) {
    dinfo_data <- fread(dinfo_aa_file, header = T, showProgress = F)
    cat("<< ", str_replace(dinfo_aa_file, src_dir, "."), "\n", sep = "")
} else {
    cat("\t- NOT EXIST HLA-ALLELE INFORMATION FILE !\n")
    cat("\t>>", str_replace(dinfo_aa_file, src_dir, "."), "\n", sep = "")
}

sig_vlist <- NULL
stat_flist <- NULL
cov_pclist <- str_subset(colnames(dhla_pdata), "PC\\d+")

for (i in seq_along(sex_tlist)) {
    sex_name <- str_c(sprintf("%02d", i - 1), "_", sex_tlist[i])
    cat("\n\t[", sex_name, "]\n", sep = "")

    stat_file <- file.path(wk_dir, str_c(smpl_name, "_dhla.", sex_name, "_omni.txt"))
    if (file.exists(stat_file)) {
        cat("- ALREADY FINISHED THIS STEP !\n")
        cat(">>", str_replace(stat_file, src_dir, "."), "\n", sep = "")
        next
    }

    # (1) Load reformatted dosage file
    cat("\n- Load reformatted dosage file\n")
    dgeno_aa_file <- str_c(dgeno_aa_fname, "_", sex_name, ".txt")
    if (!file.exists(dgeno_aa_file)) {
        cat("\t- NOT EXIST RE_FORMATTED DOSAGE FILE !\n")
        cat("\t>>", str_replace(dgeno_aa_file, geno_dir, "."), "\n", sep = "")
        next
    }

    dgeno_data <- fread(dgeno_aa_file, header = T, showProgress = F) %>%
        select(-contains("titer")) %>%
        left_join(
            dhla_pdata %>% select(IID, starts_with("seff_")),
            by = setNames("IID", "iid")
        ) %>%
        mutate(sex = as.factor(sex))
    cat("\t<< ", str_replace(dgeno_aa_file, src_dir, "."), "\n", sep = "")

    if (i == 1) {
        all_covlist <- c("sex", "age", cov_pclist)
    } else {
        all_covlist <- c("age", cov_pclist)
        dgeno_data <- dgeno_data %>% select(-sex)
    }
    dgeno_data <- dgeno_data %>% mutate(age = scale(age))

    # (2) Prepare for OMNIBUS test
    # a. Prepare variants' information for OMNIBUS test
    omni_idata <- dinfo_data %>%
        mutate(aa_id = str_c(type, gene, aa_pos, sep = "_")) %>%
        filter(r2 > r2_val) %>%
        filter(freq >= freq_val & freq <= (1 - freq_val)) %>%
        select(hla_name, aa_id)

    # b. Reformat genotype data for OMNIBUS test
    pheno_clist <- str_subset(colnames(dgeno_data), "^seff_")
    geno_clist <- colnames(dgeno_data) %>%
        setdiff(c("iid", all_covlist, pheno_clist))

    omni_gdata <- dgeno_data %>%
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
        right_join(
            omni_idata,
            by = c("genotype" = "hla_name")
        ) %>%
        group_by(phenotype, aa_id) %>%
        nest() %>%
        mutate(
            data = map(
                data,
                ~ .x %>% pivot_wider(names_from = genotype, values_from = value)
            )
        )

    # (2) Linear logistic regression by single marker
    # a. Extract and reformat data for OMNIBUS test
    wk_gdata <- omni_gdata %>%
        mutate(
            vars_name = map2(
                aa_id, data,
                ~ colnames(.y) %>% str_subset(.x)
            )
        ) %>%
        mutate(
            var_num = length(unlist(vars_name))
        ) %>%
        mutate(
            idx_name = map(var_num, ~ str_c("genotype_", seq_len(.x)))
        ) %>%
        mutate(
            data = map(
                data,
                ~ .x %>%
                    rename_with(~ unlist(idx_name), .cols = all_of(unlist(vars_name)))
            )
        ) %>%
        mutate(vars_name = str_flatten_comma(unlist(vars_name))) %>%
        filter(var_num > 1) %>%
        select(phenotype, aa_id, var_num, vars_name, data)

    # b. OMNIBUS test to AA with over two variants (remove ONLY single variant)
    cat("\t\n- Linear logistic regression by  ominibus test\n")
    buf <- wk_gdata %>%
        mutate(
            anova = map(
                data,
                ~ {
                    tryCatch(
                        {
                            bformula <- as.formula(
                                str_c(
                                    "pheno_value ~ ",
                                    str_c(str_flatten(all_covlist, collapse = "+"), sep = "+")
                                )
                            )
                            formula <- as.formula(
                                str_c(
                                    "pheno_value ~ ",
                                    str_c(
                                        str_flatten(
                                            str_subset(colnames(.x), "genotype_\\d+"),
                                            collapse = "+"
                                        ),
                                        str_flatten(all_covlist, collapse = "+"),
                                        sep = "+"
                                    )
                                )
                            )
                            bmodel <- MASS::polr(bformula, data = .x, Hess = T)
                            fmodel <- MASS::polr(formula, data = .x, Hess = T)
                            anova(bmodel, fmodel) %>% as.data.frame()
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

    stat_data <- buf %>%
        mutate(results = map(anova, ~ tail(.x, n = 1))) %>%
        unnest(results, keep_empty = T) %>%
        select(-any_of(c("anova", "Test", "Model"))) %>%
        rename_with(
            ~ ifelse(.x %in% c("Pr(>Chi)", "Pr(Chi)"), "p.value", .x)
        ) %>%
        rename_with(~ gsub(" ", "", tolower(.x))) %>%
        rename_with(~ str_remove(., "\\.$")) %>%
        select(phenotype:vars_name, df, everything())

    stat_data <- omni_idata %>%
        select(aa_id, hla_name) %>%
        mutate(gene = str_replace(aa_id, "^(AA|INS)_(.*)_(-\\d+|\\d+)", "\\2")) %>%
        mutate(aa_pos = str_replace(aa_id, "^(AA|INS)_(.*)_(-\\d+|\\d+)", "\\3")) %>%
        mutate(pos = str_replace(hla_name, str_c(aa_id, "_(-\\d+|\\d+)_(.*)"), "\\1")) %>%
        select(-hla_name) %>%
        unique() %>%
        right_join(stat_data, by = "aa_id") %>%
        select(phenotype, everything())

    # (3) Export statistic results
    cat("\t\n- Export statistic results and apply QC\n")

    fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(stat_file, src_dir, "."), "\n", sep = "")

    buf <- stat_data %>%
        filter(p.value < sig_p) %>%
        select(phenotype, aa_id)
    sig_vlist <- bind_rows(sig_vlist, buf)
    stat_flist <- c(stat_flist, stat_file)
}

# (4) Summarize statistics with significant association
cat("\n- Summarize statistics with significant association\n")

# a. Extract data from each file applied QC
cat("\n\t- Extract data from each file applied QC\n")

if (!length(stat_flist)) {
    stat_flist <- list.files(
        wk_dir, str_c(smpl_name, "_dhla.\\d{2}_.*_omni.txt"),
        recursive = T, full.names = T
    )
    for (stat_file in stat_flist) {
        buf <- fread(stat_file, header = T) %>%
            filter(p.value < sig_p) %>%
            select(phenotype, aa_id)
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
for (stat_file in stat_flist) {
    sex_name <- stat_file %>%
        basename() %>%
        str_replace(str_c(smpl_name, "_dhla\\.(\\d{2})_(.*)_omni.txt"), "\\1_\\2")
    buf <- fread(stat_file, header = T) %>%
        mutate(study = sex_name) %>%
        right_join(sig_vlist, join_by(phenotype, aa_id)) %>%
        select(phenotype, aa_id:pos, study, df:p.value)
    sig_data <- bind_rows(sig_data, buf)

    cat("\t<< ", str_replace(stat_file, src_dir, "."), "\n", sep = "")
}

# b. Reformat and add minimal p-value
cat("\n\t- Reformat and add minimal p-value\n")

buf <- sig_data %>%
    pivot_wider(
        names_from = study,
        values_from = c(df:p.value),
        names_sep = ":",
        names_vary = "slowest",
        names_glue = "{study}:{.value}"
    ) %>%
    mutate(min.p = pmap_dbl(dplyr::select(., ends_with("p.value")), pmin, na.rm = T)) %>%
    dplyr::select(phenotype, aa_id:pos, min.p, everything()) %>%
    arrange(min.p)

sig_file <- file.path(wk_dir, str_c(smpl_name, "_dhla_sigomni.txt"))
fwrite(buf, sig_file, row.names = F, col.names = T, sep = "\t")
cat("\t>> ", str_replace(sig_file, wk_dir, "."), "\n\n", sep = "")