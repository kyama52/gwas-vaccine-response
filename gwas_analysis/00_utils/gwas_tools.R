#------------------------------------------------------------------------------
# Utility functions for sample QC and phenotype transformation
#------------------------------------------------------------------------------
# [Package management]
library(dplyr)
library(stringr)

# [Functions]
# remove_smpl_by_pc--------------------------------------------------------------
remove_smpl_by_pc <- function(eigen_file,
                              smpl_ptn,
                              cond_pc1 = "> 0.01",
                              cond_pc2 = "> 0.015") {
    out_file <- str_replace(eigen_file, "(.*)\\.(.*)", "\\1_fail_pca.txt")

    if (!file.exists(eigen_file)) {
        stop(paste("Not Exist this file:", eigen_file, sep = ""))
    }
    eigen_data <- read.table(eigen_file, header = T, comment.char = "")
    if (!is.null(smpl_ptn) && nzchar(smpl_ptn)) {
        eigen_data <- eigen_data %>%
            dplyr::filter(str_detect(IID, smpl_ptn))
    }
    rmv_slist <- eigen_data %>%
        dplyr::filter(
            !!rlang::parse_expr(str_c("PC1 ", cond_pc1)) |
                !!rlang::parse_expr(str_c("PC2 ", cond_pc2))
        ) %>%
        dplyr::select(1, 2)
    if (nrow(rmv_slist)) {
        write.table(rmv_slist,
            file = out_file, sep = "\t",
            col.names = F, row.names = F, quote = F
        )
    } else {
        out_file <- ""
    }

    return(out_file)
}
# inv_nor_trans  --------------------------------------------------------------
#   IVT: inverse-normal transformation, rank-based inverse normal transformation
#   ref)
#   - https://stackoverflow.com/questions/63969719/tidyr-separate-a-column-into-a-variable-number-of-columns
inv_nor_trans <- function(x) {
    return(qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}