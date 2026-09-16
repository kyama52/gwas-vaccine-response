# ==============================================================================
# sub_B_5_ConditionedOmnibus.R
#
# CD: Apr 04 2025   K.Yamazaki
# UD: Apr 08 2025   Chage rtime when start `1`
# UD: Apr 07 2025   Bugfix
# - Add procedures as below;
#   1) scaling age in total
#   2) remove sex column in sex-specific groups
#   3) correct file name exported results by omnibus test
#   4) print start and end times.
#------------------------------------------------------------------------------
# Memo
#  Perform stratified analysis by sex to files exported by DEEP*HLA
#------------------------------------------------------------------------------
# ref)
# - It was made as reference supplied from Dr Naito HLA omnibus test as below;
#   "~/analysis/COVID-19/DEEP-HLA/240424_Omnibus/HLAassoc.share.R"
# - https://github.com/immunogenomics/HLA_analyses_tutorial/blob/main/tutorial_association.md
#
# [ATTENTION]
# - If the ANOVA test is performed with the VARIANT name,
#   an ERROR is generated with a variable name with a MINUS amino acid position.
#------------------------------------------------------------------------------
#  Result of "microbenchmark"   record == 1000, times = 10
#   1: fit_model_matrix()[model.matrix] + future_map_dfr
#   2: fit_model_matrix()[model.matrix] + map_dfr
#   3: fit_model_lm() + map_dfr
#       min       lq     mean   median       uq     max neval
#  20.86666 21.02499 21.25963 21.26920 21.47104 21.6798    10
#  36.11426 36.22089 36.35553 36.29796 36.54038 36.6264    10
#  42.72598 42.87361 43.12754 43.03333 43.25617 44.1262    10
# ==============================================================================
# [Functions]
# Use global args(pheno_cname, all_covlist) as default value
fit_model_value <- function(
    df, cond_vars = NULL, pheno_cname = NULL, all_covlist = NULL) {
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
fit_model_varlist <- function(
    df, hla_vars = NULL, cond_vars, pheno_cname = NULL, all_covlist = NULL) {
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
src_dir <- "~/analysis/COVID-19"
smpl_name <- "CUH-GWAS_230215"
geno_dname <- "PLINK_160223_1213"
step <- "B_5_ConditionedOmnibus"

id_cname <- "iid"
sex_cname <- "sex"
age_cname <- "age"
pheno_cname <- "post_titer_norm"

# chr <- 6
# str_pos <- 30
# end_pos <- 34
freq_val <- 0.01
r2_val <- 0.7
today <- as.character(format(Sys.time(), "%Y%m%d"))
sex_tlist <- c("total", "male", "female")
sig_p <- 5e-8

geno_dir <- file.path(src_dir, "genotype", geno_dname)
geno_dir_hla <- file.path(geno_dir, "04_imputation_hla")
imputed_fname <- file.path(geno_dir_hla, str_c(smpl_name, "_mhc_imputed"))
raw_pfile <- file.path(geno_dir_hla, str_c(smpl_name, "_chr6_mhc_imputed.raw.gz"))
bim_pfile <- file.path(geno_dir_hla, str_c(smpl_name, "_chr6_mhc_imputed.bim.gz"))
hla_pfile <- file.path(geno_dir, "etc", "jp.hla.txt")

dgeno_ext <- "dhla.txt"
dgeno_cv_file <- str_c(imputed_fname, ".", dgeno_ext)
dgeno_aa_file <- str_c(imputed_fname, ".aa.", dgeno_ext)
dgeno_gflist <- c(dgeno_cv_file, dgeno_aa_file)

wk_dir <- file.path(src_dir, "2308_impHLA")
wk_dir_hla <- file.path(wk_dir, "B_2_ReformatAssoc")
wk_dir_omn <- file.path(wk_dir, "B_4_RunHlaOmnibus")
wk_dir_cnd <- file.path(wk_dir, step)
if (!file.exists(wk_dir_cnd)) dir.create(wk_dir_cnd, recursive = T)

# dhla_pfile <- file.path(wk_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))
dhla_pfile <- file.path(src_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))
dinfo_cv_file <- file.path(
    wk_dir_hla, basename(dgeno_cv_file) %>% str_replace(dgeno_ext, "info")
)
dinfo_aa_file <- str_replace(dinfo_cv_file, "info$", "aa.info")

# Get the number of CPU cores allocated by Slurm
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(n_cores) || n_cores <= 0) {
    n_cores <- parallel::detectCores() - 1
    # n_cores <- 1
}
plan(multisession, workers = n_cores)

# log_file <- file.path(wk_dir, "script", "logs", str_c(step, "_", today, ".log"))

# (2) Load packages and in-house script
# source("~/tools/script/GWAS/plink.tools.R")
# source("D:/Dropbox/Tools/script/GWAS/plink.tools.R")

# [Main] ----------------------------------------------------------------------
# (0) Logging
# sink(file = log_file)
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
cat("\t- Information files related DEEP*HLA: \n")
cat("\t\t- ", str_replace(dinfo_cv_file, wk_dir, "."), "\n")
cat("\t\t- ", str_replace(dinfo_aa_file, wk_dir, "."), "\n")
cat("\t- Output directory:\n")
cat("\t\t>> ", str_replace(wk_dir_cnd, wk_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (1) Check file whether to exist or not.
cat("- Check file whether to exist or not\n")
if (file.exists(dhla_pfile)) {
    dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F)
    cat("\t<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
} else {
    cat("\t-NOT EXIST PHENOTYPE FILE FOR DEEP*HLA !\n")
    cat("\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    stop()
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

dinfo_data <- dinfo_adata %>% bind_rows(dinfo_cdata)

# (2) Prepare and exaimine conditioned analysis
sig_vlist <- NULL
stat_flist <- NULL
cov_pclist <- str_subset(colnames(dhla_pdata), "PC\\d+")

for (i in seq_along(sex_tlist)) {
    # (A) Prepare variables and dataset for analysis
    # [i] Prepare Variables
    # i <- 1
    sex_val <- i - 1
    sex_name <- str_c(sprintf("%02d", sex_val), "_", sex_tlist[i])
    edir_stat <- file.path(wk_dir_cnd, sex_name)
    if (!file.exists(edir_stat)) dir.create(edir_stat, recursive = T)

    if (i == 1) {
        egeno_cv_file <- dgeno_cv_file
        egeno_aa_file <- dgeno_aa_file
        all_covlist <- c(sex_cname, age_cname, cov_pclist)
    } else {
        egeno_cv_file <- str_replace(dgeno_cv_file, "\\.txt", str_c("_", sex_name, ".txt"))
        egeno_aa_file <- str_replace(dgeno_aa_file, "\\.txt", str_c("_", sex_name, ".txt"))
        all_covlist <- c(age_cname, cov_pclist)
    }

    cat("\n------------------------------\n")
    cat("[", sex_name, "]\n", sep = "")
    cat("------------------------------\n\n")

    sig_file <- file.path(wk_dir_omn, str_c("sig_vlist_", sex_name, ".txt"))
    sig_vdata <- data.table()

    # [ii] Prepare data for analysis
    # a. Load imputed allele files
    cat("- Load imputed allele files\n")
    if (file.exists(egeno_cv_file)) {
        rmv_vlist <- dinfo_cdata %>%
            filter(freq < freq_val | freq > (1 - freq_val) | r2 < r2_val) %>%
            pull(hla_name)
        egeno_cdata <- fread(egeno_cv_file, header = T, showProgress = F) %>%
            select(setdiff(colnames(.), rmv_vlist))
        cat("\t<<", str_replace(egeno_cv_file, geno_dir, "."), "\n", sep = "")
    } else {
        cat("\t- NOT EXIST RE_FORMATTED DOSAGE FILE [DIGIT] !\n")
        cat("\t>>", str_replace(egeno_cv_file, geno_dir, "."), "\n", sep = "")
        next
    }
    if (file.exists(egeno_aa_file)) {
        rmv_vlist <- dinfo_adata %>%
            filter(freq < freq_val | freq > (1 - freq_val) | r2 < r2_val) %>%
            pull(hla_name)
        egeno_adata <- fread(egeno_aa_file, header = T, showProgress = F) %>%
            select(setdiff(colnames(.), rmv_vlist))
        cat("\t<<", str_replace(egeno_aa_file, geno_dir, "."), "\n", sep = "")
    } else {
        cat("\t- NOT EXIST RE_FORMATTED DOSAGE FILE [AA] !\n")
        cat("\t>>", str_replace(egeno_aa_file, geno_dir, "."), "\n", sep = "")
        next
    }

    # b. Prepare data for analysis
    dgeno_data <- egeno_cdata %>%
        left_join(.,
            egeno_adata %>% select(all_of(id_cname), starts_with("AA_"), starts_with("INS_")),
            by = "iid"
        )
    nouse_cnames <- dgeno_data %>%
        colnames() %>%
        str_subset("titer") %>%
        str_subset(pheno_cname, negate = T)
    wk_gdata <- dgeno_data %>% select(-any_of(nouse_cnames))

    if (i == 1) {
        wk_gdata <- wk_gdata %>%
            mutate(across(all_of(sex_cname), as.factor)) %>%
            mutate(across(all_of(age_cname), scale))
    } else {
        wk_gdata <- wk_gdata %>% select(-any_of(nouse_cnames))
    }

    # (B) Perform conditioning analysis
    j <- 1
    flag <- 0

    while (flag == 0) {
        rtime <- sprintf("%02d", j)
        bef_rtime <- sprintf("%02d", j-1)
        cat("\n<", rtime, ">\n", sep = "")
        cat("\tStarting ", format(Sys.time(), "%Y.%m.%d %H:%M"), "\n\n", sep = "")

        #----------------------------------------------------------------------
        # [i] Check whether there is significant association < 5e-8
        cat("- Check whether there is significant association < 5e-8\n")
        if (j == 1) {
            stat_file <- file.path(wk_dir_omn, str_c(smpl_name, "_dhla.", sex_name, "_omni.txt"))
        } else {
            stat_file <- file.path(
                edir_stat, str_c(smpl_name, "_dhla.", sex_name, "_", bef_rtime, "_omni.txt")
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
                -all_of(c(id_cname, pheno_cname, all_covlist, cond_vlist)),
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
            right_join(stat_data %>% rename(hla_name = genotype), by = "hla_name") # <<< Check
        res_file <- file.path(edir_stat, str_c(smpl_name, "_hla_", sex_name, "_", rtime, "_rawlm.txt"))
        fwrite(stat_data, res_file, row.names = F, col.names = T, sep = "\t")
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
                select(all_of(c(id_cname, pheno_cname, all_covlist, cond_vlist))) %>%
                left_join(raw_data, by = id_cname)

            # a. Perform conditioning analysis [single marker]
            geno_vlist <- setdiff(
                names(wk_pdata), c(id_cname, pheno_cname, all_covlist, cond_vlist)
            )

            if (FALSE) {
                microbenchmark::microbenchmark(
                    purrr_data <- future_map_dfr(geno_vlist, function(snp_name) {
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
                    }),
                    matric_data <- map_dfr(geno_vlist, function(snp_name) {
                        df <- copy(wk_pdata)
                        df[, value := get(snp_name)] # `value` に遺伝子データをセット

                        model <- fit_model_matrix(
                            df = df,
                            cond_vars = cond_vlist,
                            pheno_cname = pheno_cname,
                            all_covlist = all_covlist
                        )

                        tidy(model, conf.int = TRUE) %>%
                            mutate(var_id = snp_name) %>%
                            select(var_id, everything()) %>%
                            relocate(p.value, .after = last_col())
                    }),
                    lm_data <- map_dfr(geno_vlist, function(snp_name) {
                        df <- copy(wk_pdata)
                        df[, value := get(snp_name)] # `value` に遺伝子データをセット

                        model <- fit_model_lm(
                            df = df,
                            cond_vars = cond_vlist,
                            pheno_cname = pheno_cname,
                            all_covlist = all_covlist
                        )

                        tidy(model, conf.int = TRUE) %>%
                            mutate(var_id = snp_name) %>%
                            select(var_id, everything()) %>%
                            relocate(p.value, .after = last_col())
                    }),
                    times = 10
                )
            }

            if (FALSE) {
                stat_data <- with_progress({
                    p <- progressor(along = geno_vlist)

                    map_dfr(geno_vlist, function(snp_name) {
                        p(message = snp_name)

                        # Create temporary column 'value'
                        df <- copy(wk_pdata)
                        df[, value := get(snp_name)]

                        model <- fit_model_matrix(
                            df = df,
                            cond_vars = cond_vlist,
                            pheno_cname = pheno_cname,
                            all_covlist = all_covlist
                        )

                        tidy(model, conf.int = TRUE) %>%
                            mutate(var_id = snp_name) %>%
                            select(var_id, everything()) %>%
                            relocate(p.value, .after = last_col())
                    })
                })
            }
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
            res_file <- file.path(edir_stat, str_c(smpl_name, "_gwas_", sex_name, "_", rtime, "_rawlm.txt"))
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
                -all_of(c(id_cname, pheno_cname, all_covlist, cond_vlist)),
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
                # anova = map2(bmodel, fmodel, ~ anova(.x, .y, test = "Chisq"))
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
            edir_stat, str_c(smpl_name, "_dhla.", sex_name, "_", rtime, "_omni.txt")
        )
        fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
        cat("\t>> ", str_replace(stat_file, src_dir, "."), "\n", sep = "")

        j <- j + 1
        cat("\tFinished ", format(Sys.time(), "%Y.%m.%d %H:%M"), "\n", sep = "")
    }

    # buf <- stat_data %>%
    #     filter(p.value < sig_p) %>%
    #     pull(aa_id)
    # sig_vlist <- c(sig_vlist, buf)
    # stat_flist <- c(stat_flist, stat_file)
}

if (FALSE) {
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
    # dplyr::select(aa_id:pos, min.p, everything()) %>%
    # arrange(min.p)
    # colnames(buf) <- colnames(buf) %>%
    #     str_replace("(.*):(\\d{2}_.*)", "\\2:\\1")

    sig_file <- file.path(wk_dir_cnd, str_c(smpl_name, "_dhla_sigomni.txt"))
    fwrite(buf, sig_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(sig_file, wk_dir, "."), "\n\n", sep = "")
}

# sink()