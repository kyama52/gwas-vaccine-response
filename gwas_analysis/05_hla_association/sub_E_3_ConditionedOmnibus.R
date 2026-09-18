#------------------------------------------------------------------------------
#  Subroutine for performing conditional analysis based on HLA imputation association results
#  [antibody titer]
#------------------------------------------------------------------------------
# [Functions]
# Use global args(pheno_cname, all_covlist) as default value
fit_model_value <- function(
  df, cond_vars = NULL, pheno_cname = NULL, all_covlist = NULL
) {
    if (is.null(pheno_cname)) pheno_cname <- get("pheno_cname", envir = .GlobalEnv)
    if (is.null(all_covlist)) all_covlist <- get("all_covlist", envir = .GlobalEnv)

    if (length(cond_vars)) {
        lm(
            formula = as.formula(
                str_c(
                    pheno_cname, " ~ value + cond_matrix + ",
                    str_flatten(all_covlist, collapse = " + ")
                )
            ),
            data = mutate(df, cond_matrix = I(as.matrix(df[, cond_vars]))),
            na.action = na.omit
        )
    } else {
        lm(
            formula = as.formula(
                str_c(
                    pheno_cname, " ~ value + ",
                    str_flatten(all_covlist, collapse = " + ")
                )
            ),
            data = df,
            na.action = na.omit
        )
    }
}

# Use global args(pheno_cname, all_covlist) as default value
fit_model_varlist <- function(df, hla_vars = NULL, cond_vars, pheno_cname = NULL, all_covlist = NULL) {
    if (is.null(pheno_cname)) pheno_cname <- get("pheno_cname", envir = .GlobalEnv)
    if (is.null(all_covlist)) all_covlist <- get("all_covlist", envir = .GlobalEnv)

    lm(
        formula = as.formula(
            str_c(
                pheno_cname, " ~ cond_matrix +", str_flatten(all_covlist, collapse = " + ")
            )
        ),
        data = mutate(df, cond_matrix = I(as.matrix(df[, c(hla_vars, cond_vars)]))),
        na.action = na.omit
    )
}


fit_model_lm <- function(df, cond_vars = NULL, pheno_cname = NULL, all_covlist = NULL) {
    if (is.null(pheno_cname)) pheno_cname <- get("pheno_cname", envir = .GlobalEnv)
    if (is.null(all_covlist)) all_covlist <- get("all_covlist", envir = .GlobalEnv)

    if (length(cond_vars)) {
        lm(
            formula = as.formula(
                str_c(
                    pheno_cname, " ~ value + cond_matrix + ",
                    str_flatten(all_covlist, collapse = " + ")
                )
            ),
            data = mutate(df, cond_matrix = I(as.matrix(df[, ..cond_vars]))),
            na.action = na.omit
        )
    } else {
        lm(
            formula = as.formula(
                str_c(
                    pheno_cname, " ~ value + ",
                    str_flatten(all_covlist, collapse = " + ")
                )
            ),
            na.action = na.omit
        )
    }
}

fit_model_matrix <- function(df, cond_vars = NULL, pheno_cname = NULL, all_covlist = NULL) {
    if (is.null(pheno_cname)) pheno_cname <- get("pheno_cname", envir = .GlobalEnv)
    if (is.null(all_covlist)) all_covlist <- get("all_covlist", envir = .GlobalEnv)

    df_filtered <- df[complete.cases(df[, c("value", cond_vars, all_covlist, pheno_cname), with = FALSE]), ]

    # Create design-matrix using `model.matrix()` (split by gene, condition variables and covariate)
    mat_geno <- model.matrix(~ value - 1, data = df_filtered)
    mat_cond <- model.matrix(~ . - 1, data = df_filtered[, ..cond_vars, drop = FALSE])
    mat_covs <- model.matrix(~., data = df_filtered[, ..all_covlist, drop = FALSE]) # including intercept

    # Combine all matrices
    X <- cbind(mat_geno, mat_cond, mat_covs)

    # Perform linear regression
    lm(df_filtered[[pheno_cname]] ~ X - 1, na.action = na.omit)
}

# [Argumnets]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "data.table", "R.utils", "broom", "progressr",
    "future", "furrr",
    update = F
)

# (1) Environmental arguments
src_dir <- "~/analysis/COVID-19/2609_PrepGitRepo"
smpl_name <- "covid-vac"
step <- "E_3_ConditionedOmnibus"

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
raw_pfile <- file.path(geno_dir, str_c(smpl_name, "_chr6_hla_imputed.raw.gz"))
bim_pfile <- file.path(geno_dir, str_c(smpl_name, "_chr6_hla_imputed.bim.gz"))

imputed_fname <- str_c(smpl_name, "_hla_imputed")
dgeno_cv_fname <- file.path(geno_dir, str_c(imputed_fname, ".dhla"))
dgeno_aa_fname <- file.path(geno_dir, str_c(imputed_fname, ".aa.dhla"))

wk_dir_omn <- file.path(src_dir, "E_2_HlaOmnibus")
wk_dir <- file.path(src_dir, step)
if (!file.exists(wk_dir)) dir.create(wk_dir, recursive = T)

etc_dir <- file.path(src_dir, "etc")
dhla_pfile <- file.path(etc_dir, str_c(smpl_name, "_dhla_pheno.txt"))
dinfo_cv_file <- file.path(etc_dir, str_c(imputed_fname, ".info"))

# Get the number of CPU cores allocated by Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores) || n_cores <= 0) {
    n_cores <- parallel::detectCores() - 1
    # n_cores <- 1
}
plan(multisession, workers = n_cores)

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
cat("\t- ", str_replace(dinfo_cv_file, src_dir, "."), "\n")
cat("\t- ", str_replace(dinfo_aa_file, src_dir, "."), "\n")
cat("- Output directory:\n")
cat("\t>> ", str_replace(wk_dir, src_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (1) Check file whether to exist or not.
cat("- Check file whether to exist or not\n")
if (file.exists(dhla_pfile)) {
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F)
    cat("\t<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t-NOT EXIST PHENOTYPE FILE FOR DEEP*HLA !\n")
    cat("\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    stop("Required input file(s) not found.")
}
if (file.exists(dinfo_cv_file)) {
    dinfo_cdata <- fread(dinfo_cv_file, header = T, showProgress = F)
    cat("\t<< ", str_replace(dinfo_cv_file, src_dir, "."), "\n", sep = "")
} else {
    cat("\t- NOT EXIST HLA-ALLELE INFORMATION FILE [DIGIT]!\n")
    cat("\t>> ", str_replace(dinfo_cv_file, src_dir, "."), "\n", sep = "")
}
if (file.exists(dinfo_aa_file)) {
    dinfo_adata <- fread(dinfo_aa_file, header = T, showProgress = F)
    cat("\t<< ", str_replace(dinfo_aa_file, src_dir, "."), "\n", sep = "")
} else {
    cat("\t- NOT EXIST HLA-ALLELE INFORMATION FILE [AA]!\n")
    cat("\t>> ", str_replace(dinfo_aa_file, src_dir, "."), "\n", sep = "")
}
dinfo_data <- dinfo_adata %>% bind_rows(dinfo_cdata)

raw_data <- data.frame()
if (file.exists(raw_pfile)) {
    raw_data <- fread(raw_pfile, header = T, showProgress = F) %>%
        rename("iid" = "IID") %>%
        select(-c(1, 3:6))
    colnames(raw_data) <- colnames(raw_data) %>% str_remove("_[^_]+$")
    cat("\t<< ", str_replace(raw_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t- NOT EXIST PLINK-DOSAGE FILE!\n")
    cat("\t>> ", str_replace(raw_pfile, src_dir, "."), "\n", sep = "")
}
bim_data <- data.frame()
if (file.exists(bim_pfile)) {
    bim_data <- fread(bim_pfile, header = F, showProgress = F)
    colnames(bim_data) <- c("chr", "var_id", "dis", "pos", "ref", "alt")
    bim_data <- bim_data %>%
        select(-dis) %>%
        select(var_id, everything())
    cat("\t<< ", str_replace(bim_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t- NOT EXIST PLINK-BIM FILE!\n")
    cat("\t>> ", str_replace(bim_pfile, src_dir, "."), "\n", sep = "")
}

# (2) Prepare and exaimine conditioned analysis
cov_pclist <- str_subset(colnames(dhla_pdata), "PC\\d+")

for (i in seq_along(sex_tlist)) {
    # (A) Prepare variables and dataset for analysis
    # [i] Prepare Variables
    sex_name <- str_c(sprintf("%02d", i - 1), "_", sex_tlist[i])

    cat("\n------------------------------\n")
    cat("[", sex_name, "]\n", sep = "")
    cat("------------------------------\n\n")

    stat_dir <- file.path(wk_dir, sex_name)
    if (!file.exists(stat_dir)) dir.create(stat_dir, recursive = T)

    dgeno_cv_file <- str_c(dgeno_cv_fname, "_", sex_name, ".txt")
    dgeno_aa_file <- str_c(dgeno_aa_fname, "_", sex_name, ".txt")

    if (i == 1) {
        all_covlist <- c("sex", "age", cov_pclist)
    } else {
        all_covlist <- c("age", cov_pclist)
    }

    # [ii] Prepare data for analysis
    # a. Load imputed allele files
    cat("- Load imputed allele files\n")
    if (file.exists(dgeno_cv_file)) {
        rmv_vlist <- dinfo_cdata %>%
            filter(freq < freq_val | freq > (1 - freq_val) | r2 < r2_val) %>%
            pull(hla_name)
        dgeno_cdata <- fread(dgeno_cv_file, header = T, showProgress = F) %>%
            select(setdiff(colnames(.), rmv_vlist))
        cat("\t<<", str_replace(dgeno_cv_file, geno_dir, "."), "\n", sep = "")
    } else {
        cat("\t- NOT EXIST RE_FORMATTED DOSAGE FILE [DIGIT] !\n")
        cat("\t>>", str_replace(dgeno_cv_file, geno_dir, "."), "\n", sep = "")
        next
    }
    if (file.exists(dgeno_aa_file)) {
        rmv_vlist <- dinfo_adata %>%
            filter(freq < freq_val | freq > (1 - freq_val) | r2 < r2_val) %>%
            pull(hla_name)
        dgeno_adata <- fread(dgeno_aa_file, header = T, showProgress = F) %>%
            select(setdiff(colnames(.), rmv_vlist))
        cat("\t<<", str_replace(dgeno_aa_file, geno_dir, "."), "\n", sep = "")
    } else {
        cat("\t- NOT EXIST RE_FORMATTED DOSAGE FILE [AA] !\n")
        cat("\t>>", str_replace(dgeno_aa_file, geno_dir, "."), "\n", sep = "")
        next
    }

    # b. Prepare data for analysis
    dgeno_data <- dgeno_cdata %>%
        left_join(.,
            dgeno_adata %>% select("iid", starts_with("AA_"), starts_with("INS_")),
            by = "iid"
        )
    nouse_cnames <- dgeno_data %>%
        colnames() %>%
        str_subset("titer") %>%
        str_subset(pheno_cname, negate = T)
    wk_gdata <- dgeno_data %>% select(-any_of(nouse_cnames))

    if (i == 1) {
        wk_gdata <- wk_gdata %>% mutate(sex = as.factor(sex))
    }

    # (B) Perform conditioning analysis
    j <- 1
    flag <- 0

    while (flag == 0) {
        rtime <- sprintf("%02d", j)
        bef_rtime <- sprintf("%02d", j - 1)
        cat("\n<", rtime, ">\n", sep = "")
        cat("\tStarting ", format(Sys.time(), "%Y.%m.%d %H:%M"), "\n\n", sep = "")

        #----------------------------------------------------------------------
        # [i] Check whether there is significant association < 5e-8
        cat("- Check whether there is significant association < 5e-8\n")
        if (j == 1) {
            stat_file <- file.path(wk_dir_omn, str_c(smpl_name, "_dhla.", sex_name, "_omni.txt"))
        } else {
            stat_file <- file.path(
                stat_dir, str_c(smpl_name, "_dhla.", sex_name, "_", bef_rtime, "_omni.txt")
            )
        }
        if (!file.exists(stat_file)) {
            cat("\tNOT EXIST OMNIBUS TEST FILE!\n")
            cat("\t>>", str_replace(stat_file, src_dir, "."), "\n", sep = "")
            flag <- 1
            break
        }

        cat("\t<< ", str_replace(stat_file, wk_dir, "."), "\n", sep = "")
        buf <- fread(stat_file, header = T, showProgress = F) %>%
            arrange(p.value) %>%
            slice(1) %>%
            select(aa_id, vars_name, p.value)
        pval <- pull(buf, p.value)[1]

        cat("\t- minimal P: ", sprintf("%.2e", pval), sep = "")
        cat(" (", pull(buf, aa_id)[1], ")\n", sep = "")
        if (pval > sig_p) {
            flag <- 1
            j <- j + 1
            next
        }
        cond_vlist <- unlist(str_split(pull(buf, vars_name), ", "))

        #----------------------------------------------------------------------
        # [ii] Conditioned analysis [single marker (1)]
        cat("\n- Conditioned analysis [single marker(1)]\n")

        # a. Perform conditioning analysis [single marker]
        stat_data <- wk_gdata %>%
            pivot_longer(
                -all_of(c("iid", pheno_cname, all_covlist, cond_vlist)),
                names_to = "genotype",
            ) %>%
            group_by(genotype) %>%
            nest() %>%
            mutate(
                model = map(
                    data,
                    ~ fit_model_value(df = .x, cond_vars = cond_vlist)
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

        # b. Export result to file
        res_data <- dinfo_data %>%
            right_join(stat_data %>% rename(hla_name = genotype), by = "hla_name")
        res_file <- file.path(stat_dir, str_c(smpl_name, "_hla_", sex_name, "_", rtime, "_rawlm.txt"))
        fwrite(res_data, res_file, row.names = F, col.names = T, sep = "\t")
        system(paste("gzip -f", res_file))
        cat("\t>> ", str_replace(res_file, src_dir, "."), ".gz\n", sep = "")

        qc_file <- str_replace(res_file, "_rawlm\\.txt", "_qclm.txt")
        res_data %>%
            filter(term == "value") %>%
            select(-term) %>%
            fwrite(qc_file, row.names = F, col.names = T, sep = "\t")
        cat("\t>> ", str_replace(qc_file, src_dir, "."), ".gz\n", sep = "")

        #----------------------------------------------------------------------
        # [iii] Conditioned analysis [single marker (2)]
        cat("\n- Conditioned analysis [single marker(2)]\n")
        if (nrow(raw_data)) {
            wk_pdata <- wk_gdata %>%
                select(all_of(c("iid", pheno_cname, all_covlist, cond_vlist))) %>%
                left_join(raw_data, by = "iid")

            # a. Perform conditioning analysis [single marker]
            geno_vlist <- setdiff(
                names(wk_pdata), c("iid", pheno_cname, all_covlist, cond_vlist)
            )

            stat_data <- future_map_dfr(geno_vlist, function(snp_name) {
                # Create temporary column 'value'
                df <- copy(wk_pdata)
                df[, value := get(snp_name)]

                model <- fit_model_matrix(
                    df = df,
                    cond_vars = cond_vlist,
                    pheno_cname = pheno_cname,
                    all_covlist = all_covlist
                )

                broom::tidy(model, conf.int = TRUE) %>%
                    mutate(var_id = snp_name) %>%
                    select(var_id, everything()) %>%
                    relocate(p.value, .after = last_col())
            })
            # }, .progress = T)

            # b. Export result to file
            res_data <- bim_data %>% right_join(stat_data, by = "var_id")
            res_file <- file.path(stat_dir, str_c(smpl_name, "_gwas_", sex_name, "_", rtime, "_rawlm.txt"))
            fwrite(stat_data, res_file, row.names = F, col.names = T, sep = "\t")
            system(paste("gzip -f", res_file))
            cat("\t>> ", str_replace(res_file, src_dir, "."), ".gz\n", sep = "")

            qc_file <- str_replace(res_file, "_rawlm\\.txt", "_qclm.txt")
            res_data %>%
                filter(term == "Xvalue") %>%
                select(-term) %>%
                fwrite(qc_file, row.names = F, col.names = T, sep = "\t")
            cat("\t>> ", str_replace(qc_file, src_dir, "."), ".gz\n", sep = "")
        } else {
            cat("\tSKIP THIS PROCEDURE")
        }

        #----------------------------------------------------------------------
        # [iv] Prepare for OMNIBUS test
        # a. Prepare aa_id list for OMNIBUS test
        aa_data <- dinfo_adata %>%
            mutate(aa_id = str_c(type, gene, aa_pos, sep = "_")) %>%
            select(hla_name, aa_id)

        # b. Reformat genotype data for OMNIBUS test
        wk_adata <- wk_gdata %>%
            pivot_longer(
                -all_of(c("iid", pheno_cname, all_covlist, cond_vlist)),
                names_to = "genotype",
            ) %>%
            inner_join(
                aa_data,
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
        wk_adata <- wk_adata %>%
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

        # c. OMNIBUS test to AA with over two variants (remove ONLY single variant)
        cat("\t\n- Linear logistic regression by  ominibus test\n")
        stat_data <- wk_adata %>%
            mutate(
                bmodel = map(
                    data,
                    ~ fit_model_varlist(df = .x, hla_vars = NULL, cond_vars = cond_vlist)
                )
            ) %>%
            mutate(
                fmodel = map2(
                    data, idx_name,
                    ~ fit_model_varlist(df = .x, hla_vars = .y, cond_vars = cond_vlist)
                )
            ) %>%
            mutate(
                anova = map2(bmodel, fmodel, ~ anova(.x, .y))
            ) %>%
            mutate(
                results = map(anova, ~ tidy(.x) %>% slice(2))
            ) %>%
            select(-c(data, bmodel, fmodel, anova)) %>%
            unnest(cols = c(results)) %>%
            select(-c(idx_name, term:rss))

        stat_data <- aa_data %>%
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
            stat_dir, str_c(smpl_name, "_dhla.", sex_name, "_", rtime, "_omni.txt")
        )
        fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
        cat("\t>> ", str_replace(stat_file, src_dir, "."), "\n", sep = "")

        j <- j + 1
        cat("\tFinished ", format(Sys.time(), "%Y.%m.%d %H:%M"), "\n", sep = "")
    }
}