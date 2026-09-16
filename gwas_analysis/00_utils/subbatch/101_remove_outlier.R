# ==============================================================================
#  101_remove_outlier.R
#   Sub routine script for removing outlier
#
#   Arguments:
#       args[1](eigen_file) Eigenvec file exported or converted by EIGEN | plink
#       args[2](cond_pc1)   Condition for removing by PC1
#       args[3](cond_pc2)   Condition for removing by PC2
#   Usage
#       Rscript sub_B_remove_outlier.R eigen_file fam_file
#
# ==============================================================================
# [Argumnets]
# (1) Get arguments
args <- commandArgs(trailingOnly = T)
if (length(args) == 4) {
    eigen_file <- args[1]
    smpl_ptn <- args[2]
    cond_pc1 <- args[3]
    cond_pc2 <- args[4]

    if (!file.exists(eigen_file)) {
        message("NOT exists EIGEN_FILE:")
        stop("\t>> ", eigen_file, sep = "")
    }
} else {
    message("Please supply FOUR argumnets: EIGEN_FILE, SMPL_PTN, COND_PC1, COND_PC2")
    message("Usage:")
    message("\tRscript 101_remove_outlier.R")
    stop("\t>> [1]EIGEN_FILE [2]SMPL_PTN [3]CONDITION(PC1) [4]CONDITION(PC2)")
}
message("\neigen_file:\n\t", eigen_file)
message("SMPL_PTN: ", smpl_ptn)
message("COND_1: PC1 ", cond_pc1)
message("COND_2: PC2 ", cond_pc2, "\n")

# (2) Load packages and in-house script
tool_dir <- "~/analysis/COVID-19/2609_PrepGitRepo/scripts/gwas_analysis/00_utils"
source(file.path(tool_dir, "gwas_tools.R"))

# [Main]
rmv_file <- remove_smpl_by_pc(eigen_file, smpl_ptn, cond_pc1, cond_pc2)

if (file.exists(rmv_file)) {
    message(">> ", rmv_file, sep = "")
} else {
    message("NO OUTLIER !")
}