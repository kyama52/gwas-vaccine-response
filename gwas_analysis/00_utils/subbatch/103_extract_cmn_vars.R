# ==============================================================================
#  sub_C_extract_cmn_vars.R
#
#  CD : 14 Feb, 2023    K.Yamazaki
#  - Add procedure for [chr X <-> chr 23]
#------------------------------------------------------------------------------
#  Memo
#   Sub routine script for extract common variants
#       between BIM(target) and VCF(reference)
#
#   [Causion!]
#   Gloval variables in this script
#   - $HOME/resource/Eagle/BBJ1K_1KGP_RefPanel_b155
#   - $HOME/resource/illumina/infinium-asian-screening-array-24v1-0/ASA-24v1-0_A1.csv
#
#   Arguments:
#       args[1](bim_file)   Bim file of targets
#       args[2](out_dir)    Output directory [OPTIONAL]
#   Usage
#       Rscript sub_C_extract_cmn_vars.R  bim_file [out_dir]
#
# ==============================================================================
# [Functions]
# base_revcmp ------------------------------------------------------------------
#   ref)
#   - https://stackoverflow.com/questions/28934185/reverse-complementary-base
base_revcmp <- function(seq) {
    char <- stringi::stri_reverse(seq) %>%
        chartr(old = "AaTtGgCc", new = "TTAACCGG") %>%
        return(char)
}

#------------------------------------------------------------------------------
# [Argumnets]
# (0) Check and load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    "tidyverse", "data.table", "R.utils", "Biostrings",
    update = F
)

# (1) Get and check arguments
args <- commandArgs(trailingOnly = T)
if (between(length(args), 1, 2)) {
    bim_file <- args[1] %>% path.expand()
    if (length(args) == 2) {
        out_dir <- args[2]
    } else {
        out_dir <- dirname(bim_file)
    }

    if (!file.exists(bim_file)) {
        message("NOT exists BIM_FILE:")
        stop("\t>> ", bim_file, sep = "")
    }

    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = T)
} else {
    message("Please supply TWO argumnets: BIM_FILE, [OUT_DIR](optional)")
    message("Usage:")
    message("\tRscript sub_C_extract_cmn_vars.R")
    stop("\t>> [1]BIM_FILE [2]OUT_DIR(optional)")
}
message("\nBIM_FILE:\n\t>> ", bim_file)
message("\nOUT_DIR:\n\t>> ", out_dir, "\n")

# (2) Global arguments
ref_dir <- "~/resource/Eagle/BBJ1K"
ref_fname <- "BBJ1K_1KGP_RefPanel_b155"
bcf_flist <- list.files(ref_dir, str_c(ref_fname, "_chr(\\d{1,2}|X)\\.bcf$"), full.names = T)
# asa_file <- "~/resource/illumina/infinium-asian-screening-array-24v1-0/ASA-24v1-0_E1.csv"
asa_file <- "~/resource/illumina/infinium-asian-screening-array-24v1-0/ASA-24v1-0_A1.csv"

clist_file <- file.path(out_dir, "snps-reference-cmn.txt")
cinfo_file <- file.path(out_dir, "snps-reference-cmn-info.txt")
dup_file <- file.path(out_dir, "snps-reference-dup.txt")
upd_id_file <- file.path(out_dir, "snps-upd-id.txt")
upd_al_file <- file.path(out_dir, "snps-upd-allele.txt")

tmpfile_1 <- tempfile()
tmpfile_2 <- tempfile()

# [Main script ] -------------------------------------------------------
# (1) Logging argmunets
cat("[Environmental argmunets]\n\n")
cat("- Reference directory:\n\t>> ", ref_dir, "\n", sep = "")
cat("\t- Reference name: ", ref_fname, "\n\n", sep = "")
cat("- Bim file:\n\t>> ", bim_file, "\n", sep = "")
cat("- Output directory:\n\t>> ", out_dir, "\n", sep = "")

# (1) Prepare datasets for combine
cat("\n[Processing data]\n\n")
cat("- Load BIM file:\n")
cat("\t<< ", bim_file, "\n", sep = "")
bim_data <- fread(bim_file, header = F, showProgress = F) %>%
    dplyr::select(1, 4, 2, 5, 6)
colnames(bim_data) <- c("chr", "pos", "id", "alt", "ref")
bim_data <- bim_data %>%
    dplyr::mutate(chr = if_else(chr == "23", "X", as.character(chr)))
cat("\t- ", nrow(bim_data), " records\n", sep = "")

cat("- Load BCF files:\n")
if (!file.exists(clist_file)) {
    bim_data %>%
        dplyr::select(chr, pos) %>%
        fwrite(tmpfile_1, row.names = F, col.names = F, sep = "\t")

    for (bcf_file in bcf_flist) {
        cat("\t<< ", str_replace(bcf_file, path.expand(ref_dir), "."), "\n", sep = "")
        system(
            paste(
                "bcftools view -R ", tmpfile_1, bcf_file,
                "| bcftools query -f '%CHROM %POS %ID %REF %ALT\n' >>", tmpfile_2
            )
        )
    }
    bcf_data <- fread(tmpfile_2, header = F, showProgress = F)
    colnames(bcf_data) <- c("chr", "pos", "id", "a0", "a1")
    fwrite(bcf_data, clist_file, row.names = F, col.names = T, sep = "\t")
    cat("\t- ", nrow(bcf_data), " records\n", sep = "")
    cat("\t\t>> ", str_replace(clist_file, out_dir, "."), "\n", sep = "")
} else {
    bcf_data <- fread(clist_file, header = T, showProgress = F)
    cat("\t- ", nrow(bcf_data), " records\n", sep = "")
    cat("\t\t<< ", str_replace(clist_file, out_dir, "."), "\n", sep = "")
}

# (2) Extract common variants
cat("- Extract common variants\n")
cmn_data <- inner_join(bim_data, bcf_data, by = c("chr", "pos")) %>%
    dplyr::rename(id_bim = 3, id_bcf = 6) %>%
    mutate(pos_name = str_c(chr, "_", pos)) %>%
    group_by(pos_name) %>%
    mutate(num = n()) %>%
    mutate_at(c("alt", "ref", "a0", "a1"), str_to_upper) %>%
    as.data.table()
fwrite(cmn_data, cinfo_file, row.names = F, col.names = T, sep = "\t")
cat("\t- ", nrow(cmn_data), " records\n", sep = "")
cat("\t>> ", cinfo_file, "\n", sep = "")

# (3) Check variant whether these variants are same or not (1) [ID, SNP]
cat("\n- Check variant whether these variants are same or not [ID, SNP]\n")

# a. SNP ID
cat("\t- Check ID: \n")
cmn_vlist <- NULL
cmn_vlist <- cmn_data %>%
    dplyr::filter(id_bim == id_bcf) %>%
    dplyr::pull(id_bcf)
cat("\t\t- Common ID: ", length(cmn_vlist), " variants\n", sep = "")

# b. Allele combination(same strand)
cat("\t- Check allele combination: \n")
tmp_data <- cmn_data %>% dplyr::filter(!(id_bcf %in% cmn_vlist))

tgt_data <- tmp_data %>%
    dplyr::filter(
        (str_c(pmin(alt, ref), "-", pmax(alt, ref)))
        == str_c(pmin(a0, a1), "-", pmax(a0, a1))
    ) %>%
    dplyr::group_by(pos_name) %>%
    dplyr::slice_min(id_bim) %>%
    dplyr::mutate(num = n())
cat("\t\t- Same strand: ", nrow(tgt_data), " variants\n")
if (nrow(tgt_data)) {
    upd_id_data <- tgt_data %>% dplyr::select(id_bim, id_bcf, pos_name)
    cmn_vlist <- c(cmn_vlist, pull(tgt_data, id_bcf))
}

# c. Allele combination(reverse complement strand)
tmp_data <- tmp_data %>%
    dplyr::filter(!(pos_name) %in% pull(upd_id_data, pos_name))
tgt_data <- tmp_data %>%
    dplyr::filter(
        (str_c(pmin(alt, ref), "-", pmax(alt, ref)))
        == str_c(
                pmin(base_revcmp(a0), base_revcmp(a1)),
                "-",
                pmax(base_revcmp(a0), base_revcmp(a1))
            )
    ) %>%
    dplyr::group_by(pos_name) %>%
    dplyr::slice_min(id_bim) %>%
    dplyr::mutate(num = n())
cat("\t\t- Reverse complement strand: ", nrow(tgt_data), " variants\n")
if (nrow(tgt_data)) {
    upd_id_data <- upd_id_data %>%
        bind_rows(tgt_data %>% dplyr::select(id_bim, id_bcf, pos_name))
    cmn_vlist <- c(cmn_vlist, pull(tgt_data, id_bcf))
}

# (4) Check variant whether these variants are same or not (2) [INDEL]
cat("\n- Check variant whether these variants are same or not [INDEL]\n")

# a. Convert allele pairs from manifest file
tmp_data <- tmp_data %>%
    dplyr::filter(!(pos_name) %in% pull(upd_id_data, pos_name)) %>%
    dplyr::filter(str_c(pmin(alt, ref), "-", pmax(alt, ref)) == "D-I") %>%
    dplyr::filter(str_length(str_c(a0, "-", a1)) > 3) %>%
    dplyr::group_by(pos_name) %>%
    dplyr::slice_min(id_bim) %>%
    dplyr::mutate(num = n())
cat("\t- INDEL: ", nrow(tmp_data), " recoreds\n\n", sep = "")

cat("\t- Load manifest file:\n\t\t<< ", asa_file, "\n", sep = "")
system(
    paste(
        "sed -e '/Controls/,$d' -e '1,7d'", asa_file, ">", tmpfile_1
    )
)
asa_data <- fread(tmpfile_1, header = T, showProgress = F) %>%
    dplyr::select(IlmnID, Name, Chr, MapInfo, SourceSeq)

buf <- asa_data %>%
    dplyr::filter(Name %in% pull(tmp_data, id_bim)) %>%
    dplyr::mutate(fseq = str_replace(SourceSeq, "(.*)\\[(.*)\\](.*)", "\\1")) %>%
    dplyr::mutate(revseq = stringi::stri_reverse(fseq)) %>%
    dplyr::mutate(inseq = str_replace(SourceSeq, "(.*)\\[-/(.*)\\](.*)", "\\2")) %>%
    dplyr::mutate(org = str_sub(str_remove_all(revseq, inseq), 1, 1)) %>%
    dplyr::mutate(ins = str_c(org, inseq)) %>%
    dplyr::mutate_at(c("org", "ins"), str_to_upper) %>%
    dplyr::select(Name, org, ins)
buf <- tmp_data %>%
    dplyr::inner_join(buf, by = c("id_bim" = "Name")) %>%
    dplyr::select(chr, pos, id_bim, alt, ref, org, ins, everything()) %>%
    dplyr::mutate(ins0 = if_else(alt == "I", org, ins)) %>%
    dplyr::mutate(ins1 = if_else(alt == "I", ins, org))
upd_al_data <- buf %>%
    ungroup() %>%
    dplyr::select(id_bim, alt, ref, ins0, ins1)
buf <- buf %>%
    dplyr::select(-c("alt", "ref", "org", "ins")) %>%
    dplyr::rename("alt" = "ins0", "ref" = "ins1") %>%
    dplyr::select(chr:id_bim, alt, ref, everything())

# b. Allele combination(same strand)
cat("\t- Check allele combination: \n")
tgt_data <- buf %>%
    dplyr::filter(
        (str_c(pmin(alt, ref), "-", pmax(alt, ref)))
        == str_c(pmin(a0, a1), "-", pmax(a0, a1))
    ) %>%
    dplyr::group_by(pos_name) %>%
    dplyr::slice_min(id_bim) %>%
    dplyr::mutate(num = n())

cat("\t\t- Same strand: ", nrow(tgt_data), " variants\n")
if (nrow(tgt_data)) {
    upd_id_data <- upd_id_data %>%
        bind_rows(tgt_data %>% dplyr::select(id_bim, id_bcf, pos_name))
    cmn_vlist <- c(cmn_vlist, pull(tgt_data, id_bcf))
}

# c. Allele combination(reverse complement strand)
tgt_data <- buf %>%
    dplyr::filter(!(pos_name) %in% pull(upd_id_data, pos_name)) %>%
    dplyr::filter(
        (str_c(pmin(alt, ref), "-", pmax(alt, ref)))
        == str_c(
                pmin(base_revcmp(a0), base_revcmp(a1)),
                "-",
                pmax(base_revcmp(a0), base_revcmp(a1))
            )
    ) %>%
    dplyr::group_by(pos_name) %>%
    dplyr::slice_min(id_bim) %>%
    dplyr::mutate(num = n())

cat("\t\t- Reverse complement strand: ", nrow(tgt_data), " variants\n")
if (nrow(tgt_data)) {
    upd_id_data <- upd_id_data %>%
        bind_rows(tgt_data %>% dplyr::select(id_bim, id_bcf, pos_name))
    cmn_vlist <- c(cmn_vlist, pull(tgt_data, id_bcf))
}

upd_al_data <- upd_al_data %>%
    dplyr::filter(id_bim %in% pull(upd_id_data, id_bim))

# d. NOT equol id_bim and Name in asa_file BUT id_bim including id_bcf
tgt_data <- tmp_data %>%
    dplyr::filter(!(pos_name) %in% pull(upd_id_data, pos_name)) %>%
    dplyr::filter(str_detect(id_bim, id_bcf)) %>%
    dplyr::mutate(a1_ins = if_else(
        str_length(a1) >= str_length(a0), 1, 0
    )) %>%
    dplyr::mutate(ins0 = if_else(
        alt == "I", if_else(a1_ins == 1, a1, a0), if_else(a1_ins == 1, a0, a1)
    )) %>%
    dplyr::mutate(ins1 = if_else(
        ref == "I", if_else(a1_ins == 1, a1, a0), if_else(a1_ins == 1, a0, a1)
    ))

cat("\t\t- Including id_bcf (eq. ilmnseq_rsXXXX): ", nrow(tgt_data), " variants\n")
if (nrow(tgt_data)) {
    upd_id_data <- upd_id_data %>%
        bind_rows(tgt_data %>% dplyr::select(id_bim, id_bcf, pos_name))
    upd_al_data <- upd_al_data %>%
        bind_rows(tgt_data %>% dplyr::select(id_bim, alt, ref, ins0, ins1))
    cmn_vlist <- c(cmn_vlist, tgt_data %>% as.data.table() %>% pull(id_bcf))
}

# (5) Export files for update ID and alleles
cat("\n- Export files for update ID and alleles\n")
if (nrow(upd_id_data)) {
    cat("\t- Update ID: ", nrow(upd_id_data), " variants\n", sep = "")
    fwrite(upd_id_data, upd_id_file, row.names = F, col.names = F, sep = "\t")
    cat("\t\t>> ", str_replace(upd_id_file, out_dir, "."), "\n", sep = "")
}

if (nrow(upd_al_data)) {
    cat("\t- Update alleles: ", nrow(upd_al_data), " variants\n", sep = "")
    fwrite(upd_al_data, upd_al_file, row.names = F, col.names = F, sep = "\t")
    cat("\t\t>> ", str_replace(upd_al_file, out_dir, "."), "\n", sep = "")
}

if (length(cmn_vlist)) {
    cat("\n- Extract variants commonly between BIM and reference\n")
    cat("\t", length(cmn_vlist), " variants\n", sep = "")
    write(cmn_vlist, dup_file)
    cat("\t\t>> ", str_replace(dup_file, out_dir, "."), "\n", sep = "")
}

cat("\n[Finished]\n")