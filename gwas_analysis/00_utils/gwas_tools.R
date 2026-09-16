# ==============================================================================
# plink tools (plink.tools.r)
# ==============================================================================
# [Package management]
# 0) Environmental arguments
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "geneplotter",
    "tidyverse", "data.table", "RColorBrewer", "R.utils", "scales",
    "cowplot", "khroma", "ggsci", "qqman", "plotly", "showtext",
    update = F
)

# loadfonts(quiet = TRUE)
# loadfonts(device = "win", quiet = TRUE)
showtext_auto()

# source("D:/Dropbox/Tools/script/Utils/tools.R")
# source("D:/Dropbox/Tools/script/Genome/genome.tools.R")
# source("D:/Dropbox/Tools/script/GWAS/manhattan.R")
source("~/tools/script/GWAS/manhattan.R")
# source("~/tools/script/Utils/tools.R")
# source("~/tools/script/Genome/genome.tools.R")
# source("~/tools/script/GWAS/manhattan.R")

# [Functions]
# batch_plink_plot  -----------------------------------------------------------
batch_plink_plot <- function(input_file_1,
                             input_file_2,
                             kind # 1==imiss_het, 2==plot_pca, 3==plot_lmiss_hist
) {
    outdir <- dirname(input_file_1)
    filebody <- basename(input_file_1) %>% stringr::str_replace("(.*).(.*)", "\\1")
    outfile <- paste(outdir, "/", filebody, ".pdf", sep = "")

    if (!file.exists(input_file_1)) stop(paste("NOT Exist this file:", input_file_1, sep = ""))
    if (!file.exists(input_file_2)) stop(paste("NOT Exist this file:", input_file_2, sep = ""))

    pdf(outfile, width = 8.5, height = 8)
    switch(kind,
        1 <- plot_imiss_het(input_file_1, input_file_2),
        2 <- plot_pca(input_file_1, input_file_2),
        3 <- plot_lmiss_hist(input_file_1)
    )
    dev.off()

    return(outfile)
}
# plot_imiss_het-----------------------------------------------------------------
plot_imiss_het <- function(imiss_file,
                           het_file) {
    imiss <- read.table(imiss_file, h = TRUE)
    imiss$logF_MISS <- log10(imiss[, 6])

    het <- read.table(het_file, h = TRUE)
    het$meanHet <- (het$N.NM. - het$O.HOM.) / het$N.NM.
    colors <- densCols(imiss$logF_MISS, het$meanHet)

    failid <- numeric(0)
    cutoff_lw <- mean(het$meanHet) - (2 * sd(het$meanHet))
    cutoff_up <- mean(het$meanHet) + (2 * sd(het$meanHet))
    failid <- which(het$meanHet < cutoff_lw | het$meanHet > cutoff_up)
    if (length(failid) > 0) {
        failsmpl <- het[failid, ]
        write.table(failsmpl,
            # file = stringr::str_c(removeFileExt(imiss_file), "_failsmpl.txt"),
            file = stringr::str_replace(imiss_file, ".imiss$", "_failsmpl.txt"),
            sep = "\t", col.names = T, row.names = F, quote = F
        )
    }

    plot(imiss$logF_MISS, het$meanHet,
        col = colors, xlim = c(-3, 0), ylim = c(0, 0.5),
        pch = 20, xlab = "Proportion of missing genotypes", ylab = "Heterozygosity rate", axes = F
    )
    axis(2, at = c(0, 0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5), tick = T)
    axis(1, at = c(-3, -2, -1, 0), labels = c(0.001, 0.01, 0.1, 1))
    abline(h = mean(het$meanHet) - (2 * sd(het$meanHet)), col = "GRAY32", lty = 2)
    abline(h = mean(het$meanHet) + (2 * sd(het$meanHet)), col = "GRAY32", lty = 2)
    abline(v = -1.522879, col = "GRAY32", lty = 2)
}
# plot_pca---------------------------------------------------------------------
# ref)    $HOME/UserData/Resource/ExampleData/NatureProtocol/raw-GWA-data.tgz
# -> plot-pca-results.Rscript
plot_pca <- function(eigen_file,
                     fam_file) {
    round_for_pca <- function(value) {
        buf <- trunc(value * 10)
        if (buf == 0) {
            buf <- trunc(value * 100)
            if (buf == 0) {
                if (value > 0) {
                    buf <- 0.01
                } else {
                    buf <- -0.01
                }
            } else {
                if (buf > 0) {
                    buf <- buf + 1
                } else {
                    buf <- buf - 1
                }
                buf <- buf / 100
            }
        } else {
            if (buf > 0) {
                buf <- buf + 1
            } else {
                buf <- buf - 1
            }
            buf <- buf / 10
        }

        return(buf)
    }
    unk_color <- "#D6D6D6"
    if (!file.exists(eigen_file)) stop(paste("Not Exist this file:", eigen_file, sep = ""))
    eigen_data <- read.table(eigen_file, h = F, skip = 1)

    if (!file.exists(fam_file)) stop(paste("Not Exist this file:", fam_file, sep = ""))
    fam_data <- read.table(fam_file, h = F)

    pop_list <- c("case", "control", "CEU", "CHB", "JPT", "YRI", "unknown")
    data <- data.frame(cbind(eigen_data, fam_data[, ncol(fam_data)]))
    colnames(data) <- c("FID", "IID", stringr::str_c("PC", seq(1, 10)), "POP")
    cpalette <- data.frame(
        POP = factor(pop_list, levels = pop_list),
        color = khroma::color("bright")(length(pop_list)) %>% as.character()
    ) %>%
        dplyr::mutate(color = if_else(POP == "unknown", unk_color, color))

    plot_data <- data %>%
        dplyr::select(dplyr::ends_with("ID"), "PC1", "PC2", "POP") %>%
        dplyr::mutate(POP = if_else(POP %in% c(-9, 0, 99), length(pop_list), POP)) %>%
        dplyr::mutate(POP = factor(pop_list[POP], levels = pop_list)) %>%
        dplyr::inner_join(cpalette, by = "POP") %>%
        dplyr::arrange(match(.[["POP"]], levels(.[["POP"]])))
    # dplyr::inner_join(cpalette, by = "POP") %>%
    # dplyr::mutate(alpha = as.factor(if_else(POP %in% c("case", "control", "unknown"), 1, 0.6)))

    xmin <- round_for_pca(plot_data %>% dplyr::pull("PC1") %>% min())
    xmax <- round_for_pca(plot_data %>% dplyr::pull("PC1") %>% max())
    ymin <- round_for_pca(plot_data %>% dplyr::pull("PC2") %>% min())
    ymax <- round_for_pca(plot_data %>% dplyr::pull("PC2") %>% max())

    pca_plot <- ggplot(
        plot_data,
        aes(x = PC1, y = PC2, label = IID, color = POP)
    ) +
        geom_point(size = 2, alpha = 0.8) +
        scale_color_manual(
            breaks = plot_data$POP,
            values = plot_data$color, name = ""
        ) +
        theme_cowplot(
            font_family = "Open Sans",
            line_size = 1
        ) +
        theme(
            # plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12),
            legend.position = c(0.85, 0.25)
        ) +
        xlim(c(xmin, xmax)) +
        ylim(c(ymin, ymax)) +
        xlab("Principal Component 1") +
        ylab("Principal Component 2") +
        guides(alpha = "none")

    return(pca_plot)
}

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
# plot_lmiss_hist---------------------------------------------------------------
# ref)    $HOME/UserData/Resource/ExampleData/NatureProtocol/raw-GWA-data.tgz
# -> lmiss-hist.Rscript
plot_lmiss_hist <- function(lmiss_file) {
    if (!file.exists(lmiss_file)) stop(paste("Not Exist this file:", lmiss_file, sep = ""))
    x <- read.table(lmiss_file, header = T)

    ymax <- 3 # ylim (0, 150000)
    scale <- 50000
    xlabels <- as.character(10^(rep(-4:0)))
    ylabels <- stringr::str_c(rep(0:ymax * (scale / 1000)), "K")
    par(mfrow = c(1, 1))

    hist(log10(x$F_MISS),
        axes = F, xlim = c(-4, 0), col = "RED",
        ylab = "Number of SNPs", xlab = "Fraction of missing data", main = "All SNPs",
        ylim = c(0, ymax * scale)
    )
    axis(side = 2, labels = F)
    mtext(ylabels, side = 2, las = 2, at = rep(0:ymax * scale), line = 1)
    axis(side = 1, labels = F)
    mtext(xlabels, side = 1, at = c(-4, -3, -2, -1, 0), line = 1)
    abline(v = log10(0.05), lty = 2)
}
# calc_lambda_gc -----------------------------------------------------------
calc_lambda_gc <- function(p_values, # vector
                           col_flag = FALSE # p_values[1]: 1 == cname, 0 == data
) {
    if (col_flag) {
        label <- p_values[0]
    } else {
        label <- ""
    }

    lambda <- 0

    lambda <- qchisq(median(p_values), df = 1, lower.tail = FALSE) / qchisq(0.5, df = 1, lower.tail = FALSE)
    sprintf("lamda GC [%s] is %.4f", label, lambda)

    chi <- numeric(length(p_values))
    chi <- qchisq(p_values, df = 1, lower.tail = FALSE)
    chi <- ifelse(is.nan(chi), 0, chi)
    lambda <- quantile(chi, 0.5) / qchisq(0.5, df = 1, lower.tail = FALSE) # Dr.Okada
    sprintf("lamda GC [%s] is %.4f", label, lambda)

    return(invisible(lambda))
}
# plot_manhattan_qq-------------------------------------------------------------
plot_manhattan_qq <- function(res_file,
                              auto_chr = FALSE) {
    outdir <- dirname(res_file)
    filebody <- stringr::str_remove(basename(res_file), ".(txt|txt.gz)$")
    if (!file.exists(res_file)) stop(paste("Not Exist this file:", res_file, sep = ""))
    x <- data.table::fread(res_file, header = T)

    if (auto_chr) {
        buf <- x[which(x$CHR <= 22), ]
        ychrlabs <- as.character(c(1:22))
    } else {
        buf <- x
        ychrlabs <- c(1:22, "X")
    }

    ymax <- ceiling(max(-log10(buf$P)))
    if (ymax <= 10) ymax <- 10

    outfile <- paste(outdir, "/", filebody, "_1.png", sep = "")
    png(outfile, height = 960, width = 1372)
    par(family = "Open Sans", mar = c(6, 6, 5, 2), mgp = c(3, 1, 0))

    buf <- manhattan(buf,
        ylim = c(0, ymax),
        cex = 0.6, cex.axis = 1.2, cex.lab = 1.6, cex.main = 2.0,
        suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
        chrlabs = ychrlabs
    )
    dev.off()

    outfile <- paste(outdir, "/", filebody, "_2.png", sep = "")
    png(outfile, height = 960, width = 960)
    par(family = "Open Sans", mar = c(6, 6, 5, 2), mgp = c(3, 1, 0))
    qq(buf$P,
        main = "Q-Q plot of GWAS p-values",
        cex = 0.6, cex.axis = 1.2, cex.lab = 1.6, cex.main = 2.0
    )
    dev.off()

    outfile <- paste(outdir, "/", filebody, "_lambda.txt", sep = "")
    lambda <- calc_lambda_gc(buf$P, col_flag = TRUE)
    write.table(lambda, file = outfile)

    return(lambda)
}

# ------------------------------------------------------------------------------
# plot_manhattan:
#   Arguments:
#       filename, # files of GWAS result (including SNP,CHR,BP,P)
#       chr = "CHR", # Column's name: CHR
#       bp = "BP", # Column's name: BP
#       p = "P", # Column's name: P
#       snp = "SNP", # Column's name: SNP
#       col = c("gray30", "gray60"),
#       auto_chr = TRUE,
#       sug_line = 1e-5,
#       sig_line = 5e-8,
#       mark_1 = NULL,
#       markcol_1 = "orangered2",
#       mark_2 = NULL,
#       markcol_2 = "blue3",
#       mark_3 = NULL,
#       markcol_3 = "green3",
#       title = NULL, ...) {
#       df          data.frame
#       add_df      data.frame with additional data
#       key_value   key value for joining
#   Attension
#       MUST initialize data.frame
#
plot_manhattan <- function(filename, # files of GWAS result (including SNP,CHR,BP,P)
                           chr = "CHR", # Column's name: CHR
                           bp = "BP", # Column's name: BP
                           p = "P", # Column's name: P
                           snp = "SNP", # Column's name: SNP
                           col = c("gray30", "gray60"),
                           auto_chr = TRUE,
                           ymax = NULL,
                           ylab = NULL,
                           sug_line = 1e-5,
                           sig_line = 5e-8,
                           mark_1 = NULL,
                           markcol_1 = "orangered2",
                           mark_2 = NULL,
                           markcol_2 = "blue3",
                           mark_3 = NULL,
                           markcol_3 = "green3",
                           title = NULL, ...) {
    # 0) Check Arguments
    if (!file.exists(filename)) {
        stop(paste("Not Exist GWAS result file: ",
            basename(filename),
            sep = ""
        ))
    }
    # stat_data <- fread(filename, header = TRUE, showProgress = FALSE)
    stat_data <- data.table::fread(
        cmd = paste("zcat", shQuote(path.expand(filename))),
        header = TRUE, na.strings = c("NA", "."), showProgress = FALSE
    )

    # Check for sensible dataset
    if (!(chr %in% colnames(stat_data))) stop("NOT found column:", chr)
    if (!(bp %in% colnames(stat_data))) stop("NOT found column:", bp)
    if (!(p %in% colnames(stat_data))) stop("NOT found column:", p)
    if (!(snp %in% colnames(stat_data))) stop("NOT found column:", snp)

    # 1) Prepare dataset [pos_data]
    pos_cols <- c(chr, bp, p, snp)
    pos_data <- stat_data %>%
        dplyr::select(all_of(pos_cols)) %>%
        dplyr::filter(!is.nan(!!sym(p)) & !is.na(!!sym(p)))
    colnames(pos_data) <- pos_cols
    if (!is.numeric(pull(pos_data, !!chr))) {
        pos_data <- pos_data %>%
            mutate(!!chr := case_when(
                !!sym(chr) == "X" ~ "23",
                !!sym(chr) == "Y" ~ "24",
                !!sym(chr) == "XY" ~ "25",
                !!sym(chr) == "MT" ~ "26",
                TRUE ~ !!sym(chr)
            ))
        buf <- pull(pos_data, !!chr) %>%
            as.numeric() %>%
            is.na() %>%
            sum()
        if (!buf) {
            pos_data <- pos_data %>% mutate(!!chr := as.integer(!!sym(chr)))
        } else {
            stop(paste(chr, "column should be [1-22,X,Y,XY,MT]"))
        }
    }
    if (auto_chr) pos_data <- pos_data %>% filter(!!sym(chr) <= 22)

    # 2) Preparing data for Manhattan plot [plot_data]
    plot_data <- pos_data %>%
        group_by(!!sym(chr)) %>%
        summarise(chr_len = max(!!sym(bp)), .groups = "keep") %>% # Compute chromosome size
        ungroup() %>%
        mutate(total_mgn = cumsum(as.numeric(chr_len)) - chr_len) %>% # Calculate cumulative position of each chromosome
        left_join(pos_data, ., by = chr) %>%
        mutate(pos = !!sym(bp) + total_mgn) %>% # Add a cumulative position of each SNP
        dplyr::select(-c("chr_len", "total_mgn"))

    # Setting "x-axis label"
    #   get chromosome center positions for x-axis
    axis_data <- plot_data %>%
        group_by(!!sym(chr)) %>%
        summarize(center = (max(pos) + min(pos)) / 2) %>%
        dplyr::mutate(labchr = case_when(
            !!sym(chr) == 23 ~ "X",
            !!sym(chr) == 24 ~ "Y",
            !!sym(chr) == 25 ~ "XY",
            !!sym(chr) == 26 ~ "MT",
            TRUE ~ as.character(!!sym(chr))
        ))

    # Setting "max of y-axis"
    data_ymax <- ceiling(max(-log10(pull(pos_data, !!p))))
    if (is.null(ymax)) {
        ymax <- data_ymax
    } else {
        ymax <- if_else(data_ymax > ymax, data_ymax, ymax)
    }
    if (ymax <= 10) ymax <- 10
    ymax <- ceiling(ymax / 2) * 2

    if (is.null(ylab)) {
        ylab <- "-log10(P)"
    }

    # 3) Plotting using by ggplot2
    man_plot <- ggplot(plot_data, aes(x = pos, y = -log10(!!sym(p)))) +
        theme_classic(base_family = "Open Sans") +
        geom_point(aes(color = as.factor(!!sym(chr))), size = 0.5) +
        scale_color_manual(values = rep(col, nrow(axis_data))) +
        scale_x_continuous(
            label = axis_data$labchr,
            breaks = axis_data$center
        ) +
        # expand=c(0,0) removes space between plot area and x axis
        # scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
        scale_y_continuous(
            expand = c(0, 0), limits = c(0, ymax),
            breaks = seq(0, ceiling(ymax / 2) * 2, by = 2)
        ) +
        labs(x = "Chromosome", y = ylab) +
        theme(
            legend.position = "none",
            plot.title = element_text(face = "bold", hjust = 0.5),
            axis.text.x = element_text(face = "bold")
        )
    if (!is.null(title)) man_plot <- man_plot + ggtitle(title)
    if (!is.null(mark_1)) {
        man_plot <- man_plot +
            geom_point(
                data = subset(plot_data, SNP %in% mark_1),
                col = markcol_1, size = 1
            )
    }
    if (!is.null(mark_2)) {
        man_plot <- man_plot +
            geom_point(
                data = subset(plot_data, SNP %in% mark_2),
                col = markcol_2, size = 1
            )
    }
    if (!is.null(mark_3)) {
        man_plot <- man_plot +
            geom_point(
                data = subset(plot_data, SNP %in% mark_3),
                col = markcol_3, size = 1
            )
    }
    man_plot <- man_plot +
        geom_hline(yintercept = -log10(sig_line), color = "gray40") +
        geom_hline(
            yintercept = -log10(sug_line), color = "gray40",
            linetype = "dotted"
        )

    return(man_plot)
}

# Plot_manhattan_with_mark-------------------------------------------------------
plot_manhattan_with_mark <- function(res_file,
                                     sig_snps = NULL,
                                     sugg_snps = NULL,
                                     auto_chr = FALSE) {
    outdir <- dirname(res_file)
    filebody <- stringr::str_remove(basename(res_file), ".(txt|txt.gz)$")
    if (!file.exists(res_file)) stop(paste("Not Exist this file:", res_file, sep = ""))
    x <- data.table::fread(res_file, header = T) %>%
        dplyr::filter(!is.nan(P) & !is.na(P))

    if (auto_chr) {
        buf <- x[which(x$CHR <= 22), ]
        ychrlabs <- as.character(c(1:22))
    } else {
        buf <- x
        ychrlabs <- c(1:22, "X")
    }

    ymax <- ceiling(max(-log10(buf$P)))
    if (ymax <= 10) ymax <- 10

    outfile <- paste(outdir, "/", filebody, "_3.png", sep = "")
    png(outfile, height = 960, width = 1372)
    par(family = "Open Sans", mar = c(6, 6, 5, 2), mgp = c(3, 1, 0))

    buf <- manhattan(buf,
        ylim = c(0, ymax),
        cex = 0.6, cex.axis = 1.2, cex.lab = 1.6, cex.main = 2.0,
        suggestiveline = F, genomewideline = F,
        chrlabs = ychrlabs
    )

    if (!is.null(sig_snps)) {
        marker_info <- buf[which(buf$SNP %in% sig_snps), ]
        with(marker_info, points(pos, logp, col = "red3", pch = 20, cex = 1))
    }
    if (!is.null(sugg_snps)) {
        marker_info <- buf[which(buf$SNP %in% sugg_snps), ]
        with(marker_info, points(pos, logp, col = "green3", pch = 20, cex = 1))
    }

    abline(h = -log10(1e-05), lty = 3, col = "gray40")
    abline(h = -log10(5e-08), col = "gray40")

    dev.off()
}
# extract_sig_signals------------------------------------------------------------
extract_sig_signals <- function(res_file,
                                threshold = 1e-5,
                                CHR = "CHR",
                                BP = "BP",
                                P = "P") {
    outdir <- dirname(res_file)
    filebody <- stringr::str_remove(basename(res_file), ".(txt|txt.gz)$")
    if (!file.exists(res_file)) stop(stringr::str_c("Not Exist this file:", res_file))
    x <- data.table::fread(res_file, header = T)

    idx_p <- match(P, colnames(x))
    idx_c <- match(CHR, colnames(x))
    idx_b <- match(BP, colnames(x))
    if (!length(idx_p)) {
        stop(stringr::str_c("Not Including Column [", P, "]:", res_file))
    }

    sig <- x[which(x[[idx_p]] < threshold), ]
    if (nrow(sig) > 0 && length(idx_c) && length(idx_b)) {
        # sig[, NEAREST_GENE:=mapply(search_gene_all_chr, CHR, BP, margin=5000)]
        sig[, NEAREST_GENE := mapply(search_gene_all_chr, sig[[idx_c]], sig[[idx_b]], margin = 5000)]
        outfile <- stringr::str_c(outdir, "/", filebody, "_", threshold, ".txt")
        write.table(sig, file = outfile, sep = "\t", col.names = T, row.names = F, quote = F)
    } else {
        message(stringr::str_c("No significant signals THRESHOLD=", threshold))
        message(stringr::str_c("\t-> Target file:", res_file))
    }

    return(sig)
}
# extract_region_by_snp  ----------------------------------------------------------
extract_region_by_snp <- function(stat_file,
                                  margin = 5,
                                  marker_list,
                                  ID = "ID",
                                  CHR = "CHR",
                                  BP = "BP",
                                  out_dir) {
    get_cname <- function(df, cnames) {
        idx <- which(colnames(df) %in% cnames) %>% max()
        if (idx) {
            cname <- colnames(df)[idx]
        } else {
            cname <- ""
        }
        return(cname)
    }
    extract_region <- function(chr,
                               bp) {
        reg_data <- stat_data %>%
            filter(!!sym(cname_c) == chr) %>%
            filter(!!sym(cname_b) >= (bp - margin * 1000000) &
                !!sym(cname_b) <= (bp + margin * 1000000))

        return(reg_data)
    }

    if (!file.exists(stat_file)) stop(stringr::str_c("Not Exist this file:", stat_file))
    if (out_dir == "") out_dir <- dirname(stat_file)
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = T)
    filebody <- stringr::str_remove(basename(stat_file), ".(txt|txt.gz|gz)$")
    # filebody <- stringr::str_remove(basename(stat_file), ".(txt|txt.gz|gz)$") %>%
    #     stringr::str_replace("(.*)\\.(.*)", "\\1")

    stat_data <- data.table::fread(stat_file, header = T, showProgress = F)

    cname_i <- get_cname(stat_data, c(ID, "ID", "SNP", "alternate_ids"))
    cname_c <- get_cname(stat_data, c(CHR, "CHR", "CHROM", "chromosome"))
    cname_b <- get_cname(stat_data, c(BP, "BP", "POS", "position"))
    if (cname_i == "") {
        stop(stringr::str_c("Not Including Column [", ID, "|ID|SNP]:\n\t", stat_file))
    }
    if (cname_c == "") {
        stop(stringr::str_c("Not Including Column [", CHR, "|CHR|CHROM]:\n\t", stat_file))
    }
    if (cname_b == "") {
        stop(stringr::str_c("Not Including Column [", BP, "|BP|POS]:\n\t", stat_file))
    }

    marker_info <- data.table("no" = seq_len(length(marker_list)), marker_list) %>%
        rename(!!cname_i := marker_list) %>%
        left_join(stat_data %>% filter(!!sym(cname_i) %in% marker_list), by = cname_i) %>%
        tibble() %>%
        mutate(reg = map2(!!sym(cname_c), !!sym(cname_b), extract_region)) %>%
        mutate(out_file = file.path(
            out_dir,
            str_c(
                filebody, "_", sprintf("%02d", no), "_",
                str_replace_all(!!sym(cname_i), ":", "-"), ".txt"
            )
        ))
    purrr::walk2(
        marker_info[["reg"]], marker_info[["out_file"]],
        ~ write.table(.x, file = .y, sep = "\t", quote = F, col.names = T, row.names = F)
    )

    return(pull(marker_info, out_file))
}

# extract_region  ----------------------------------------------------------
extract_region <- function(stat_file,
                           margin = 5,
                           marker_list,
                           ID = "ID",
                           CHR = "CHR",
                           BP = "BP",
                           out_dir) {
    get_cname <- function(df, cnames) {
        idx <- which(colnames(df) %in% cnames) %>% max()
        if (idx) {
            cname <- colnames(df)[idx]
        } else {
            cname <- ""
        }
        return(cname)
    }
    extract_region <- function(chr,
                               bp) {
        reg_data <- stat_data %>%
            filter(!!sym(cname_c) == chr) %>%
            filter(!!sym(cname_b) >= (bp - margin * 1000000) &
                !!sym(cname_b) <= (bp + margin * 1000000))

        return(reg_data)
    }

    if (!file.exists(stat_file)) stop(stringr::str_c("Not Exist this file:", stat_file))
    if (out_dir == "") out_dir <- dirname(stat_file)
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = T)
    filebody <- stringr::str_remove(basename(stat_file), ".(txt|txt.gz|gz)$")
    # filebody <- stringr::str_remove(basename(stat_file), ".(txt|txt.gz|gz)$") %>%
    #     stringr::str_replace("(.*)\\.(.*)", "\\1")

    stat_data <- data.table::fread(stat_file, header = T, showProgress = F)

    cname_i <- get_cname(stat_data, c(ID, "ID", "SNP", "alternate_ids"))
    cname_c <- get_cname(stat_data, c(CHR, "CHR", "CHROM", "chromosome"))
    cname_b <- get_cname(stat_data, c(BP, "BP", "POS", "position"))
    if (cname_i == "") {
        stop(stringr::str_c("Not Including Column [", ID, "|ID|SNP]:\n\t", stat_file))
    }
    if (cname_c == "") {
        stop(stringr::str_c("Not Including Column [", CHR, "|CHR|CHROM]:\n\t", stat_file))
    }
    if (cname_b == "") {
        stop(stringr::str_c("Not Including Column [", BP, "|BP|POS]:\n\t", stat_file))
    }

    marker_info <- data.table("no" = seq_len(length(marker_list)), marker_list) %>%
        rename(!!cname_i := marker_list) %>%
        left_join(stat_data %>% filter(!!sym(cname_i) %in% marker_list), by = cname_i) %>%
        tibble() %>%
        mutate(reg = map2(!!sym(cname_c), !!sym(cname_b), extract_region)) %>%
        mutate(out_file = file.path(
            out_dir,
            str_c(
                filebody, "_", sprintf("%02d", no), "_",
                str_replace_all(!!sym(cname_i), ":", "-"), ".txt"
            )
        ))
    purrr::walk2(
        marker_info[["reg"]], marker_info[["out_file"]],
        ~ write.table(.x, file = .y, sep = "\t", quote = F, col.names = T, row.names = F)
    )

    return(pull(marker_info, out_file))
}
# reformat_pvar  --------------------------------------------------------------
#   ref)
#   - https://stackoverflow.com/questions/63969719/tidyr-separate-a-column-into-a-variable-number-of-columns
reformat_pvar <- function(pvar_file,
                          out_dir = "") {
    add_genotyped <- function(str) {
        buf <- str %>%
            str_split(";") %>%
            unlist()
        if (str_detect(last(buf), "TYPED|IMPUTED")) {
            buf[length(buf)] <- str_c("Genotyped=", last(buf))
        }
        return(stringr::str_flatten(buf, ";"))
    }

    if (!file.exists(pvar_file)) stop(stringr::str_c("Not Exist this file:", pvar_file))
    if (out_dir == "") out_dir <- dirname(pvar_file)
    if (!file.exists(out_dir)) dir.create(out_dir, recursive = T)
    rfmt_file <- file.path(out_dir, str_replace(basename(pvar_file), "(\\.pvar|\\.pvar.gz)$", ".pvar.txt"))

    pvar_data <- fread(pvar_file, sep = "\t", header = T, showProgress = F)

    pvar_data <- pvar_data %>%
        mutate(INFO = map_chr(INFO, add_genotyped)) %>%
        separate_rows(INFO, sep = ";") %>%
        separate(INFO, c("field", "value"), sep = "=") %>%
        pivot_wider(names_from = field, values_from = value)

    fwrite(pvar_data, rfmt_file, row.names = F, col.names = T, sep = "\t")
    system(paste("gzip -f", rfmt_file))
    message(">> ", str_c(rfmt_file, ".gz\n"))
}
# inv_nor_trans  --------------------------------------------------------------
#   IVT: inverse-normal transformation, rank-based inverse normal transformation
#   ref)
#   - https://stackoverflow.com/questions/63969719/tidyr-separate-a-column-into-a-variable-number-of-columns
inv_nor_trans <- function(x) {
    return(qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x))))
}