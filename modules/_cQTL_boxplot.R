#!/usr/bin/env Rscript
###########################################
# _cQTL_boxplot.R
# Lightweight version adapted for circular RNA percentages
# Converts PLINK counted-allele dosage to ALT dosage for plotting
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript _cQTL_boxplot.R <pairs_file> <snp_file> <circ_file> <cov_file> <out_dir> <tissue>")
}

GSfile    <- args[1]
snp_file  <- args[2]
circ_file <- args[3]
cov_file  <- args[4]
out_dir   <- args[5]
tissue    <- args[6]

# TargetALS remains the default for legacy runs.
experiment <- Sys.getenv("EXPERIMENT", "target_ALS")
if (!experiment %in% c("target_ALS", "answer_ALS")) {
  stop("EXPERIMENT must be target_ALS or answer_ALS.")
}

# Existing files needed only to orient 0/1/2 as ALT dosage
snp_loc <- file.path(dirname(snp_file), "snp_location.txt")
RAW_FILE <- file.path(
  path.expand("~/donglab/data"),
  experiment,
  "QTL",
  "plink",
  "joint_all_chrs_matrixEQTL.raw"
)

out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_cQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_cQTL_boxplots_colored.pdf"))

# Read the pairs structure (expecting header containing circ_id and snpid)
pairs <- fread(GSfile, header=FALSE)
if (ncol(pairs) >= 2) {
  setnames(pairs, c("circ_id", "snpid"))
} else {
  stop("Error: Pairs file must contain at least two columns (circ_id and snpid).")
}

message("# Loading genomic matrices...")
snp_mat  <- fread(snp_file, header=TRUE)
circ_mat <- fread(circ_file, header=TRUE)
cov_mat  <- fread(cov_file, header=TRUE)

# Standardize internal matrix headers
setnames(circ_mat, names(circ_mat)[1], "circ_id")
setnames(snp_mat, names(snp_mat)[1], "snpid")

# Recover REF/ALT and PLINK-counted allele for all SNPs once
snp_map <- fread(snp_loc)
raw_cols <- scan(RAW_FILE, what=character(), nlines=1, quiet=TRUE)
raw_cols <- raw_cols[grepl("^[^:]+:[0-9]+:[ACGT]+:[ACGT]+_[ACGT]+$", raw_cols)]

alleles <- data.table(raw_col=raw_cols)
alleles[, c("variant", "counted") := tstrsplit(raw_col, "_", fixed=TRUE)]
alleles[, c("chr", "pos", "ref", "alt") := tstrsplit(variant, ":", fixed=TRUE)]
alleles[, chr := sub("^chr", "", chr)]
snp_map[, `:=`(chr_key=sub("^chr", "", as.character(chr)), pos_key=as.character(pos))]
alleles <- merge(
  snp_map[, .(snpid, chr_key, pos_key)],
  alleles[, .(chr_key=chr, pos_key=pos, ref, alt, counted)],
  by=c("chr_key", "pos_key"), all.x=TRUE
)

message(sprintf("# Allele mappings recovered: %d / %d SNPs",
                sum(!is.na(alleles$counted)), nrow(alleles)))

# Determine case/control status.
sample_ids <- names(circ_mat)[-1]

if (experiment == "answer_ALS") {
  status <- ifelse(
    grepl("^CASE-", sample_ids, ignore.case=TRUE), "CASE",
    ifelse(grepl("^CTRL-", sample_ids, ignore.case=TRUE), "CTRL", NA_character_)
  )

  if (anyNA(status)) {
    stop(
      "Could not assign AnswerALS CASE/CTRL status for: ",
      paste(sample_ids[is.na(status)], collapse=", ")
    )
  }

  als_vals <- as.integer(status == "CASE")
  status_labels <- c("CASE", "CTRL")
} else {
  # TargetALS keeps its existing is_als covariate behavior.
  cov_first_col <- names(cov_mat)[1]
  als_row <- cov_mat[get(cov_first_col) == "is_als", ]

  if (nrow(als_row) == 0) {
    als_vals <- rep(0, ncol(cov_mat) - 1)
  } else {
    als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
  }

  status_labels <- c("ALS", "Control")
}

if (length(als_vals) != length(sample_ids)) {
  stop("Case/control status count does not match the number of circRNA-matrix samples.")
}

p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)

run_plotting <- function(pdf_path, use_status_colors = FALSE) {
    pdf(pdf_path, width=12, height=5)
    plots_made <- 0
    pairs_skipped <- 0

    for (i in seq_len(nrow(pairs))) {
        c_id <- pairs$circ_id[i]
        S    <- pairs$snpid[i]

        circ_row <- circ_mat[circ_id == c_id]
        snp_row  <- snp_mat[snpid == S]

        if (nrow(circ_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        circ_vals <- as.numeric(unlist(circ_row[, -1, with=FALSE], use.names=FALSE))
        # Input circ matrix stores fractions (0–1); display percentages (0–100).
        if (any(!is.na(circ_vals) & (!is.finite(circ_vals) | circ_vals < 0 | circ_vals > 1))) {
            stop("Expected circ fractions in [0, 1]; do not supply a percent-scaled matrix.")
        }
        circ_vals <- 100 * circ_vals
        snp_vals  <- as.numeric(unlist(snp_row[, -1, with=FALSE], use.names=FALSE))

        if (length(circ_vals) == 0 || length(snp_vals) == 0 || length(circ_vals) != length(snp_vals)) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        # Convert PLINK counted-allele dosage to ALT dosage
        a <- alleles[snpid == S]
        if (nrow(a) != 1 || is.na(a$counted[1]) || !(a$counted[1] %in% c(a$ref[1], a$alt[1]))) {
            message("# WARNING: allele orientation unavailable for ", S, "; skipping.")
            pairs_skipped <- pairs_skipped + 1
            next
        }
        if (a$counted[1] == a$ref[1]) snp_vals <- 2 - snp_vals

        df <- data.frame(
            circ_percent = circ_vals,
            SNP = snp_vals,
            p_col = if(use_status_colors) p_colors else rep("darkred", length(circ_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(circ_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(circ_vals)),
            stringsAsFactors = FALSE
        )

        df <- df[!is.na(df$SNP) & !is.na(df$circ_percent), ]
        if (nrow(df) < 5) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        get_p <- function(formula, data) {
            tryCatch({
                fit <- lm(formula, data=data)
                formatC(summary(fit)$coefficients[2,4], format="e", digits=3)
            }, error = function(e) "NA")
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))

        # Additive Model
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)

        boxplot(circ_percent ~ SNP_f, data=df, col="lightgreen", outline=FALSE,
                ylab="circRNA Percentage (%)", xlab="Genotype",
                main=paste0("Additive Model\n(p = ", get_p(circ_percent ~ SNP, df), ")"))
        points(jitter(as.numeric(df$SNP_f), amount=0.15), df$circ_percent,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        if (use_status_colors) legend("topleft", legend=status_labels, col=c("red", "#9932CC"),
                                      pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

        # Dominant Model
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)

        boxplot(circ_percent ~ dom_f, data=df, col="lightblue", outline=FALSE,
                ylab="circRNA Percentage (%)", xlab="Genotype",
                main=paste0("Dominant Model\n(p = ", get_p(circ_percent ~ dom_val, df), ")"))
        points(jitter(as.numeric(df$dom_f), amount=0.15), df$circ_percent,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        if (use_status_colors) legend("topleft", legend=status_labels, col=c("red", "#9932CC"),
                                      pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

        # Recessive Model
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)

        boxplot(circ_percent ~ rec_f, data=df, col="lightpink", outline=FALSE,
                ylab="circRNA Percentage (%)", xlab="Genotype",
                main=paste0("Recessive Model\n(p = ", get_p(circ_percent ~ rec_val, df), ")"))
        points(jitter(as.numeric(df$rec_f), amount=0.15), df$circ_percent,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)
        if (use_status_colors) legend("topleft", legend=status_labels, col=c("red", "#9932CC"),
                                      pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

        mtext(paste("circRNA:", c_id, "| SNP:", S, "| REF:", a$ref[1], "| ALT:", a$alt[1]),
              outer=TRUE, cex=1.0, font=2, line=0.5)

        plots_made <- plots_made + 1
    }

    dev.off()
    message(sprintf("# Plots made: %d, pairs skipped: %d", plots_made, pairs_skipped))
}

message(paste("# Processing", nrow(pairs), "circ-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
