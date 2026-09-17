#------------------------------------------------------------------------------
# Subroutine for rerforming sex-stratified omnibus analysis of HLA amino acid positions
#  [antibody titer]
#------------------------------------------------------------------------------
# [Functions]
# Use global args(pheno_cname, all_covlist) as default value
fit_model_value <- function(
  df, geno_vars = NULL, pheno_cname = NULL, all_covlist = NULL
) {
    if (is.null(pheno_cname)) pheno_cname <- get("pheno_cname", envir = .GlobalEnv)
    if (is.null(all_covlist)) all_covlist <- get("all_covlist", envir = .GlobalEnv)

    if (length(geno_vars)) {
        lm(
            formula = as.formula(
                str_c(
                    pheno_cname, " ~ cond_matrix + ", str_flatten(all_covlist, collapse = " + ")
                )
            ),
            data = mutate(df, cond_matrix = I(as.matrix(df[, geno_vars]))),
            na.action = na.omit
        )
    } else {
        lm(
            formula = as.formula(
                str_c(
                    pheno_cname, " ~ ", str_flatten(all_covlist, collapse = " + ")
                )
            ),
            data = df,
            na.action = na.omit
        )
    }
}

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
step <- "E_2_HlaOmnibus"

# id_cname <- "iid"
# sex_cname <- "sex"
# age_cname <- "age"
pheno_cname <- "post_titer_norm"
freq_val <- 0.01
r2_val <- 0.7
today <- as.character(format(Sys.time(), "%Y%m%d"))
sex_tlist <- c("total", "male", "female")
sig_p <- 5e-8

geno_dir <- file.path(src_dir, "genotype", "04_imputation_hla")

imputed_fname <- str_c(smpl_name, "_hla_imputed")
# dgeno_ext <- "dhla.txt"
# dgeno_aa_file <- str_c(imputed_fname, ".aa.", dgeno_ext)
dgeno_aa_fname <- file.path(geno_dir, str_c(imputed_fname, ".aa.dhla"))

wk_dir <- file.path(src_dir, step)
if (!file.exists(wk_dir)) dir.create(wk_dir, recursive = T)

# dhla_pfile <- file.path(wk_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))
etc_dir <- file.path(src_dir, "etc")
dhla_pfile <- file.path(etc_dir, str_c(smpl_name, "_dhla_pheno.txt"))
dinfo_aa_file <- file.path(etc_dir, str_c(imputed_fname, ".aa.info"))

# [Main] ----------------------------------------------------------------------
# (0) Logging
cat("[Environmental argmunets]\n\n")
cat("- Root directory:\n\t>> ", src_dir, "\n", sep = "")
cat("- Genotype directory:\n\t>> ")
cat(str_replace(geno_dir, src_dir, "."), "\n", sep = "")
cat("- Phenotype:", pheno_cname, "\n")
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
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F)
    cat("<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t-NOT EXIST PHENOTYPE FILE FOR DEEP*HLA !\n")
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
    cat("\n[", sex_name, "]\n", sep = "")

    dgeno_aa_file <- str_c(dgeno_aa_fname, "_", sex_name, ".txt")
    if (i == 1) {
        all_covlist <- c("sex", "age", cov_pclist)
    } else {
        all_covlist <- c("age", cov_pclist)
    }

    # (1) Load reformatted dosage file
    cat("\n- Load reformatted dosage file\n")

    if (!file.exists(dgeno_aa_file)) {
        cat("\t- NOT EXIST RE_FORMATTED DOSAGE FILE !\n")
        cat("\t>>", str_replace(dgeno_aa_file, geno_dir, "."), "\n", sep = "")
        next
    }
    dgeno_data <- fread(dgeno_aa_file, header = T, showProgress = F)
    cat("\t<< ", str_replace(dgeno_aa_file, geno_dir, "."), "\n", sep = "")

    # (2) Prepare for OMNIBUS test
    # a. Prepare variants' information for OMNIBUS test
    cat("\n- Prepare variants' information for OMNIBUS test\n")
    omni_idata <- dinfo_data %>%
        mutate(aa_id = str_c(type, gene, aa_pos, sep = "_")) %>%
        filter(r2 > r2_val) %>%
        filter(freq >= freq_val & freq <= (1 - freq_val)) %>%
        select(hla_name, aa_id)
    cat("\t- (Overall): ", nrow(dinfo_data), " variants\n", sep = "")
    cat("\t- (after QC [MAF > ", freq_val, ", r2 > ", r2_val, "]): ", sep = "")
    cat(nrow(omni_idata), " variants\n", sep = "")

    # b. Reformat genotype data for OMNIBUS test
    cat("\n- Reformat genotype data for OMNIBUS test\n")
    omni_gdata <- dgeno_data %>%
        pivot_longer(
            -all_of(c("iid", all_covlist, pheno_cname)),
            names_to = "genotype",
        ) %>%
        right_join(
            omni_idata,
            by = c("genotype" = "hla_name")
        ) %>%
        group_by(aa_id) %>%
        nest() %>%
        mutate(
            data = map(
                data,
                ~ .x %>% pivot_wider(names_from = genotype, values_from = value)
            )
        )

    # c. Extract and reformat data for OMNIBUS test
    #   OMNIBUS test to AA with over two variants (remove ONLY single variant)
    wk_gdata <- omni_gdata %>%
        mutate(
            vars_name = map2(
                aa_id, data,
                ~ colnames(.y) %>% str_subset(pattern = str_c(.x, "_"))
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
        select(aa_id, var_num, vars_name, idx_name, data)
    cat("\t- (Overall): ", nrow(omni_gdata), " AA\n", sep = "")
    cat("\t- (var_num > 1): ", nrow(wk_gdata), " AA\n", sep = "")

    # (2) Linear logistic regression by single marker
    cat("\t\n- Linear logistic regression by  ominibus test\n")
    stat_data <- wk_gdata %>%
        mutate(
            bmodel = map(
                data, ~ fit_model_value(df = .x)
            )
        ) %>%
        mutate(
            fmodel = map2(
                data, idx_name, ~ fit_model_value(df = .x, geno_vars = .y)
            )
        ) %>%
        mutate(
            anova = map2(bmodel, fmodel, ~ anova(.x, .y))
            # anova = map2(bmodel, fmodel, ~ anova(.x, .y, test = "Chisq"))
        ) %>%
        mutate(results = map(
            anova,
            ~ tidy(.x) %>% slice(2)
        )) %>%
        select(-c(data, idx_name, bmodel, fmodel, anova)) %>%
        unnest(cols = c(results)) %>%
        select(-c(term:rss))
    stat_data <- omni_idata %>%
        select(aa_id, hla_name) %>%
        mutate(gene = str_replace(aa_id, "^(AA|INS)_(.*)_(-\\d+|\\d+)", "\\2")) %>%
        mutate(aa_pos = str_replace(aa_id, "^(AA|INS)_(.*)_(-\\d+|\\d+)", "\\3")) %>%
        mutate(pos = str_replace(hla_name, str_c(aa_id, "_(-\\d+|\\d+)_(.*)"), "\\1")) %>%
        select(-hla_name) %>%
        unique() %>%
        right_join(stat_data, by = "aa_id")

    # (3) Export statistic results
    cat("\t\n- Export statistic results and apply QC\n")

    stat_file <- file.path(
        wk_dir,
        str_c(smpl_name, "_dhla.", sex_name, "_omni.txt")
    )
    fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(stat_file, src_dir, "."), "\n", sep = "")

    buf <- stat_data %>%
        filter(p.value < sig_p) %>%
        pull(aa_id)
    sig_vlist <- c(sig_vlist, buf)
    stat_flist <- c(stat_flist, stat_file)
}

# (4) Summarize statistics with significant association
cat("\n-----------------------------------------------------------\n")
cat("\n- Summarize statistics with significant association\n")

# a. Extract data from each file applied QC
cat("\n\t- Extract data from each file applied QC\n")
sig_vlist <- sig_vlist %>%
    sort() %>%
    unique()
cat("\t- No. of significant variants: ", length(sig_vlist), "\n", sep = "")

sig_data <- data.frame()
for (stat_file in stat_flist) {
    sex_name <- stat_file %>%
        basename() %>%
        str_replace(str_c(smpl_name, "_dhla\\.(\\d{2})_(.*)_omni.txt"), "\\1_\\2")
    buf <- fread(stat_file, header = T) %>%
        mutate(study = sex_name) %>%
        filter(aa_id %in% sig_vlist) %>%
        select(aa_id:pos, study, df:p.value)
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
    dplyr::select(aa_id:pos, min.p, everything())

sig_file <- file.path(wk_dir, str_c(smpl_name, "_dhla_sigomni.txt"))
fwrite(buf, sig_file, row.names = F, col.names = T, sep = "\t")
cat("\t>> ", str_replace(sig_file, wk_dir, "."), "\n\n", sep = "")