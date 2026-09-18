#------------------------------------------------------------------------------
# Subroutine for running SNP-by-sex interaction analysis for lead variants
#------------------------------------------------------------------------------
# [Arguments]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "broom", "data.table", "R.utils",
    update = F
)

# (1) Environmental arguments
src_dir <- "/path/to/project"
smpl_name <- "covid-vac"
top_list <- c("rs2596506", "rs117342275")

today <- as.character(format(Sys.time(), "%y%m%d"))
phenotype <- "post_titer_norm"

geno_dir <- file.path(src_dir, "genotype", "06_cand_vars")
raw_file <- file.path(geno_dir, str_c(smpl_name, "_toptag.raw"))
cmb_file <- str_replace(raw_file, ".raw$", ".txt")

step <- "G_1_SnpSexInteraction"
wk_dir <- file.path(src_dir, step)
if (!file.exists(wk_dir)) dir.create(wk_dir, recursive = T)

dhla_pfile <- file.path(src_dir, "etc", str_c(smpl_name, "_dhla_pheno.txt"))

# [Main] ----------------------------------------------------------------------
# (0) Logging
cat("-----------------------------------------------------------\n")
cat("[Environmental argmunets]\n\n")
cat("- Root directory:\n\t>> ", src_dir, "\n", sep = "")
cat("- Working directory:\n\t>> ")
cat(str_replace(geno_dir, src_dir, "."), "\n", sep = "")
cat("\t\t- ", str_replace(raw_file, geno_dir, "."), "\n")
cat("- Phenotype file:\n")
cat("\t>> ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
cat("\t- Output directory:\n")
cat("\t>> ", str_replace(wk_dir, src_dir, "."), "\n", sep = "")
cat("\n-----------------------------------------------------------\n")

# (1) Load datasets
cat("\n- Load datasets\n")

if (!file.exists(cmb_file)) {
    cat("\t- Prepare dataset for ")
    # a. Check file whether to exist or not.
    cat("- Check file whether to exist or not\n")

    if (!file.exists(raw_file)) {
        cat("\t[NOT EXIST RAW FILE !]\n")
        cat("\t>>", str_replace(raw_file, src_dir, "."), "\n", sep = "")
        stop()
    } else {
        raw_gdata <- fread(raw_file, header = T, showProgress = F) %>%
            rename_with(
                .fn = ~ str_remove(., "_[^_]+$"), .cols = 7:ncol(.)
            ) %>%
            {
                select(., all_of(c(2, 7:ncol(.))))
            }
        cat("<< ", str_replace(raw_file, src_dir, "."), "\n", sep = "")
    }

    if (!file.exists(dhla_pfile)) {
        cat("\t[NOT EXIST DHLA_PHENO FILE !]\n")
        cat("\t>>", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
        stop()
    } else {
        dhla_pdata <- fread(dhla_pfile, header = T, showProgress = F) %>%
            select(1:3, post_titer_norm, starts_with("PC"))
        cat("<< ", str_replace(dhla_pfile, src_dir, "."), "\n", sep = "")
    }

    cat("\n\t- Export combined dataset\n")
    cmb_data <- dhla_pdata %>% right_join(raw_gdata, by = c("iid" = "IID"))
    fwrite(cmb_data, cmb_file, row.names = F, col.names = T, sep = "\t")
    cat("\t>> ", str_replace(cmb_file, src_dir, "."), "\n", sep = "")
} else {
    cmb_data <- fread(cmb_file, header = T, showProgress = F)
    cat("<< ", str_replace(cmb_file, src_dir, "."), "\n", sep = "")
}
cmb_data <- cmb_data %>% mutate(sex = as.factor(sex))
cov_list <- str_subset(colnames(cmb_data), "PC\\d+")

# (2) Run SNP-by-sex interaction analysis
cat("\n- Test SNP-by-sex interactions\n")

int_dt <- data.table()

for (vname in top_list) {
    fml <- as.formula(
        paste0(
            phenotype, " ~ sex + age + ", paste(cov_list, collapse = " + "),
            " + sex * ", vname
        )
    )
    model <- lm(fml, data = cmb_data)

    int_dt <- bind_rows(
        int_dt,
        tidy(model) %>%
            mutate(variant = vname, .before = term)
        # filter(str_detect(term, "sex.*:|:.*sex")) %>%
    )
}

# (3) Export results
cat("\n- Export results\n")
int_file <- file.path(wk_dir, str_c(smpl_name, "_sex_interaction_", today, ".txt"))
fwrite(int_dt, int_file, sep = "\t")
cat(">> ", str_replace(int_file, src_dir, "."), "\n", sep = "")