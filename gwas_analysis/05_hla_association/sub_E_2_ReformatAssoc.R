# ==============================================================================
# sub_B_2_ReformatAssoc.R
#
# CD: Mar 28 2024   K.Yamazaki
# UD: Aug 28 2024
#   - Change the environment from windows to linux on va-server
#   - Change procedure from to run this script to load this script in shell script with apptainer
#     -> Rename filename from "B_2_ReformatAssoc.R" to "sub_B_2_ReformatAssoc.R"
#   - Comment out description of logging
#   - Bug fix
# UD: May 08 2024
#   - Bug fix [combine genotype with phenotype using by "right_join (<- left_join)"]
#   - Change PHENO_FILE from "${SMPL_NAME}.tsv" to "${SMPL_NAME}.txt"
#   - Change PHENOTYPE from "post_titer_log" to "post_titer_norm"
#   - Remove samples to be detectable antibody titer before anti SARS-CoV-2 vaccination
# UD: Apr 01 2024   Change exporting filename
#------------------------------------------------------------------------------
# Memo
# Reformat files exported by DEEP*HLA and examine linear logistic regression
#------------------------------------------------------------------------------
# ref)
# - https://stackoverflow.com/questions/67694428/using-mutateacross-with-purrrmap
# - https://stackoverflow.com/questions/33668912/lm-across-many-columns-in-a-dataframe-in-r
# - https://stackoverflow.com/questions/65091970/multiple-linear-models-on-each-column-of-dataframe
# - https://sudori.info/stat/stat_tidyverse_06.html
# - https://bookdown.org/palmjulia/r_intro_script/inference-ii.html#regression
#
# [glm]
# - https://www2.kobe-u.ac.jp/~bunji/files/lecture/MVA/mva-05-regression.pdf
# - https://statmath.wu.ac.at/courses/heather_turner/glmCourse_001.pdf
# - https://cran.r-project.org/web/packages/jtools/vignettes/summ.html
#
# [Odds]
# - https://stats.stackexchange.com/questions/8661/logistic-regression-in-r-odds-ratio
# - https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS1.pdf
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
step <- "B_2_ReformatAssoc"

sex_cname <- "sex"
age_cname <- "age"
# pheno_cname <- "post_titer_log"
pheno_cname <- "post_titer_norm"
today <- as.character(format(Sys.time(), "%Y%m%d"))

geno_dir <- file.path(src_dir, "genotype", geno_dname)
geno_dir_hla <- file.path(geno_dir, "04_imputation_hla")
geno_dir_etc <- file.path(geno_dir_hla, "etc")
phased_fname <- file.path(geno_dir_hla, str_c(smpl_name, "_mhc_phased"))
imputed_fname <- file.path(geno_dir_hla, str_c(smpl_name, "_mhc_imputed"))

fam_file <- str_c(phased_fname, ".fam")
dhla_ext <- "deephla.dosage"
dhla_cv_file <- str_c(imputed_fname, ".", dhla_ext)
dhla_aa_file <- str_c(imputed_fname, ".aa.", dhla_ext)
dhla_gflist <- c(dhla_cv_file, dhla_aa_file)

hla_pfile <- file.path(geno_dir_etc, "jp.hla.txt")
acc_ext <- "accuracy.txt"
acc_cv_file <- file.path(geno_dir_etc, str_c(smpl_name, "_mhc.cv.", acc_ext))
acc_aa_file <- file.path(geno_dir_etc, str_c(smpl_name, "_mhc.cv.aa.", acc_ext))

wk_dir <- file.path(src_dir, "2308_impHLA")
wk_dir_hla <- file.path(wk_dir, step)
# wk_dir_etc <- file.path(wk_dir, "etc")
if (!file.exists(wk_dir_hla)) dir.create(wk_dir_hla, recursive = T)
# if (!file.exists(wk_dir_etc)) dir.create(wk_dir_etc, recursive = T)
# dhla_pfile <- file.path(wk_dir_etc, str_c(smpl_name, "_dhla_pheno.txt"))
dhla_pfile <- file.path(src_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))
# log_file <- file.path(wk_dir, "script", "logs", str_c(step, "_", today, ".log"))

# wk_dir_gwas <- file.path(src_dir, "2304_applyQC")
# pheno_file <- file.path(wk_dir_gwas, "etc", str_c(smpl_name, ".tsv"))
pheno_file <- file.path(src_dir, "etc", str_c(smpl_name, ".txt"))
cov_file <- file.path(src_dir, "etc", str_c(smpl_name, "_plink.covar"))
rmv_sfile <- file.path(src_dir, "etc", "rmv_smpl_expvirus.txt")

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
cat("\t- Dosage files:\n")
cat("\t\t- ", str_replace(dhla_cv_file, geno_dir, "."), "\n")
cat("\t\t- ", str_replace(dhla_aa_file, geno_dir, "."), "\n")
cat("\t- Fam file:", str_replace(fam_file, geno_dir, "."), "\n")
cat("\t- Accuracy calculated by DEEP*HLA: \n")
cat("\t\t- ", str_replace(acc_cv_file, geno_dir, "."), "\n")
cat("\t\t- ", str_replace(acc_aa_file, geno_dir, "."), "\n")
cat("- Working directory:\n\t>> ")
cat(str_replace(wk_dir, src_dir, "."), "\n", sep = "")
cat("\t- Remove sample file:\n")
cat("\t\t>> ", str_replace(rmv_sfile, src_dir, "."), "\n")
cat("\t- Phenotype:", pheno_cname, "\n")
cat("\t- Phenotype file:\n")
cat("\t\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
cat("\t- Output directory:\n")
cat("\t\t>> ", str_replace(wk_dir_hla, src_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n\n")

# (0) Check file whether to exist or not.
cat("- Check file whether to exist or not\n")
if (!file.exists(fam_file)) {
    cat("\t-NOT EXIST FAM FILE !\n")
    cat("\t>>", str_replace(fam_file, geno_dir, "."), "\n", sep = "")
    stop()
} else {
    fam_data <- fread(fam_file, header = F, showProgress = F) %>%
        magrittr::set_colnames(c("fid", "iid", "", "", "sex", "phenotype")) %>%
        select(1, 2, 5)
    cat("<< ", str_replace(fam_file, src_dir, "."), "\n", sep = "")
}

if (!file.exists(hla_pfile)) {
    cat("\t-NOT EXIST HLA_POS FILE !\n")
    cat("\t>>", str_replace(hla_pfile, geno_dir, "."), "\n", sep = "")
    stop()
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
        stop()
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
        wk_dir_hla,
        str_replace(dhla_gfile, dhla_ext, "info") %>% basename()
    )

    if (dhla_gfile == dhla_cv_file) {
        hinfo_data <- hinfo_data %>%
            mutate(
                name =
                    str_replace(hla_name, "(HLA_.*)_(.*)", "\\1") %>%
                        str_replace("(MICA|MICB|TAP\\d)_(.*)", "\\1")
            ) %>%
            left_join(hla_pdata, by = "name") %>%
            # group_by(name) %>%
            # mutate(mgn = seq_len(n())) %>%
            # mutate(assm_pos = pos + mgn - 1) %>%
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

    geno_dfile <- str_replace(dhla_gfile, dhla_ext, "dhla.txt")
    cat("\t<< ", str_replace(dhla_gfile, src_dir, "."), "\n", sep = "")
    geno_data <- dhla_data %>%
        select(-c(2, 3)) %>%
        column_to_rownames("hla_name") %>%
        as.matrix() %>%
        t() %>%
        as.data.table()
    rownames(geno_data) <- fam_data %>% pull(iid)
    geno_data <- dhla_pdata %>%
        right_join(geno_data %>% rownames_to_column("iid"), by = "iid")
    fwrite(geno_data, geno_dfile, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(geno_dfile, src_dir, "."), "\n", sep = "")

    # (4) Apply sample QC which remove samples with detectable bef_titer
    if (length(rmv_slist)) {
        cat("\n- Apply sample QC which remove samples with detectable bef_titer\n")
        cat("\t- Remove sample No: ", length(rmv_slist), "\n", sep = "")
        geno_data <- geno_data %>% filter(!(iid %in% rmv_slist))
    }

    # (5) Linear logistic regression by single marker
    cat("\n- Linear logistic regression by single marker\n")
    cat("\t- Sample No: ", nrow(geno_data), "\n", sep = "")
    all_covlist <- c(sex_cname, age_cname, cov_nlist)
    stat_data <- geno_data %>%
        select(-c("bef_titer", "post_titer")) %>%
        pivot_longer(
            -all_of(c("iid", pheno_cname, all_covlist)),
            names_to = "genotype",
        ) %>%
        group_by(genotype) %>%
        nest() %>%
        mutate(
            model = map(
                data,
                ~ lm(
                    str_c(
                        pheno_cname, "~ ",
                        str_c(
                            "value", sex_cname, age_cname, str_flatten(cov_nlist, collapse = "+"),
                            sep = "+"
                        )
                    ),
                    data = .,
                    na.action = na.omit
                )
            )
        ) %>%
        select(-data) %>%
        mutate(summary = list(tidy(model[[1]]))) %>%
        unnest(summary) %>%
        select(-model) %>%
        mutate(or = exp(estimate)) %>%
        mutate(low95_or = exp(estimate - 1.96 * std.error)) %>%
        mutate(up95_or = exp(estimate + 1.96 * std.error)) %>%
        filter(term != "(Intercept)") %>%
        rename(test = term)

    # (5) Export statistic results and apply QC
    cat("\n- Export statistic results and apply QC\n")
    stat_file <- basename(dhla_gfile) %>%
        str_remove("_imputed") %>%
        str_replace(dhla_ext, "rawlm.txt") %>%
        file.path(wk_dir_hla, .)
    stat_data <- hinfo_data %>%
        right_join(stat_data, by = c("hla_name" = "genotype"))
    fwrite(stat_data, stat_file, row.names = F, col.names = T, sep = "\t")
    system(paste("gzip -f", stat_file))
    cat(">> ", str_replace(stat_file, src_dir, "."), ".gz\n", sep = "")

    qc_file <- str_replace(stat_file, "rawlm.txt", "qclm.txt")
    qc_data <- stat_data %>%
        filter(test == "value") %>%
        select(-test) %>%
        filter(
            r2 >= 0.7 & !is.na(r2)
            # if (dhla_gfile == dhla_aa_file) {
            #     r2 >= 0.7 & !is.na(r2)
            # } else {
            #     !is.na(r2)
            # }
        )
    fwrite(qc_data, qc_file, row.names = F, col.names = T, sep = "\t")
    cat(">> ", str_replace(qc_file, src_dir, "."), "\n", sep = "")
}

# sink()
