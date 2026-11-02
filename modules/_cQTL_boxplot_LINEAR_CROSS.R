#!/usr/bin/env Rscript
###########################################
# _cQTL_boxplot_LINEAR_CROSS.R
# Interaction Layout: Split Boxplots + Regression Slopes (cQTL Edition)
# Converts PLINK counted-allele dosage to ALT dosage for plotting
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript _cQTL_boxplot_LINEAR_CROSS.R <pairs_file> <snp_file> <expr_file> <cov_file> <snp_loc> <out_dir> <tissue>")
}

GSfile    <- args[1]
snp_file  <- args[2]
expr_file <- args[3]
cov_file  <- args[4]
snp_loc   <- args[5]
out_dir   <- args[6]
tissue    <- args[7]

# TargetALS remains the default for legacy runs.
experiment <- Sys.getenv("EXPERIMENT", "target_ALS")
if (!experiment %in% c("target_ALS", "answer_ALS")) {
  stop("EXPERIMENT must be target_ALS or answer_ALS.")
}

RAW_FILE <- file.path(
  path.expand("~/donglab/data"),
  experiment,
  "QTL",
  "plink",
  "joint_all_chrs_matrixEQTL.raw"
)

out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_Interaction_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_Interaction_boxplots_colored.pdf"))

# Read file structure (explicitly handle headerless top-pairs file)
pairs <- fread(GSfile, header=FALSE)
setnames(pairs, c("circ_id", "snpid"))

message("# Loading genomic matrices...")
snp_mat  <- fread(snp_file, header=TRUE)
expr_mat <- fread(expr_file, header=TRUE)
cov_mat  <- fread(cov_file, header=TRUE)
snp_map  <- fread(snp_loc, header=TRUE)

setnames(expr_mat, names(expr_mat)[1], "circ_id")
setnames(snp_mat, names(snp_mat)[1], "snpid")

# Recover REF/ALT + PLINK-counted allele from RAW header
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
sample_ids <- names(expr_mat)[-1]

if (experiment == "answer_ALS") {
  status <- ifelse(
    grepl("^CASE-", sample_ids, ignore.case=TRUE), "ALS",
    ifelse(grepl("^CTRL-", sample_ids, ignore.case=TRUE), "Control", NA_character_)
  )

  if (anyNA(status)) {
    stop(
      "Could not assign AnswerALS CASE/CTRL status for: ",
      paste(sample_ids[is.na(status)], collapse=", ")
    )
  }

  als_vals <- as.integer(status == "ALS")
  case_label <- "CASE"
  control_label <- "CTRL"
} else {
  # TargetALS keeps its existing is_als covariate behavior.
  cov_first_col <- names(cov_mat)[1]
  als_row <- cov_mat[cov_mat[[1]] == "is_als", ]

  if (nrow(als_row) == 0) {
    warning("Warning: 'is_als' row not found in covariate matrix. Defaulting all values to 0.")
    als_vals <- rep(0, ncol(cov_mat) - 1)
  } else {
    als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
  }

  case_label <- "ALS"
  control_label <- "Control"
}

if (length(als_vals) != length(sample_ids)) {
  stop("Case/control status count does not match the number of phenotype-matrix samples.")
}

# Jitter plotting definitions
p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)

run_plotting <- function(pdf_path, use_status_colors=FALSE) {
    pdf(pdf_path, width=14, height=5.5)
    plots_made <- 0
    pairs_skipped <- 0

    for (i in seq_len(nrow(pairs))) {
        c_id <- pairs$circ_id[i]
        S    <- pairs$snpid[i]

        expr_row <- expr_mat[circ_id == c_id]
        snp_row  <- snp_mat[snpid == S]

        if (nrow(expr_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        expr_vals <- as.numeric(unlist(expr_row[, -1, with=FALSE], use.names=FALSE))
        # Input circ matrix stores fractions (0–1); display percentages (0–100).
        if (any(!is.na(expr_vals) & (!is.finite(expr_vals) | expr_vals < 0 | expr_vals > 1))) {
            stop("Expected circ fractions in [0, 1]; do not supply a percent-scaled matrix.")
        }
        expr_vals <- 100 * expr_vals
        snp_vals  <- as.numeric(unlist(snp_row[, -1, with=FALSE], use.names=FALSE))

        if (length(expr_vals) == 0 || length(snp_vals) == 0 || length(expr_vals) != length(snp_vals)) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        # Convert PLINK counted-allele dosage -> ALT dosage
        a <- alleles[snpid == S]
        if (nrow(a) != 1 || is.na(a$counted[1]) || !(a$counted[1] %in% c(a$ref[1], a$alt[1]))) {
            message("# WARNING: allele orientation unavailable for ", S, "; skipping.")
            pairs_skipped <- pairs_skipped + 1
            next
        }
        if (a$counted[1] == a$ref[1]) snp_vals <- 2 - snp_vals

        df <- data.frame(
            expression = expr_vals,
            SNP = snp_vals,
            Status = ifelse(als_vals > 0.5, "ALS", "Control"),
            p_col = if(use_status_colors) p_colors else rep("darkred", length(expr_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(expr_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(expr_vals)),
            stringsAsFactors = FALSE
        )

        df <- df[!is.na(df$SNP) & !is.na(df$expression), ]

        if (nrow(df[df$Status == "ALS", ]) < 3 || nrow(df[df$Status == "Control", ]) < 3) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        # Pull p-value helper specifically for the Interaction term (Genotype:Status)
        get_interaction_p <- function(data) {
            tryCatch({
                fit <- lm(expression ~ SNP * factor(Status), data=data)
                p_val <- summary(fit)$coefficients["SNP:factor(Status)Control", 4]
                formatC(p_val, format="e", digits=3)
            }, error=function(e) "NA")
        }

        # 3-Panel Layout: Controls, ALS, and Overlay Line Trend
        par(mfrow=c(1,3), mar=c(5,4,4,1), oma=c(0,0,3,0))

        # 1. Controls Boxplot
        df_ctrl <- df[df$Status == "Control", ]
        ctrl_counts <- table(factor(df_ctrl$SNP, levels=0:2))
        ctrl_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", ctrl_counts, ")")
        df_ctrl$SNP_f <- factor(df_ctrl$SNP, levels=0:2, labels=ctrl_labels)

        boxplot(expression ~ SNP_f, data=df_ctrl, col="lightblue", outline=FALSE,
                ylab="circRNA Percentage (%)", xlab="Genotype", main=paste(control_label, "Samples Only"))
        points(jitter(as.numeric(df_ctrl$SNP_f), amount=0.15), df_ctrl$expression,
               pch=df_ctrl$p_pch, col=df_ctrl$p_col, cex=df_ctrl$p_cex)

        # 2. ALS Boxplot
        df_als <- df[df$Status == "ALS", ]
        als_counts <- table(factor(df_als$SNP, levels=0:2))
        als_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", als_counts, ")")
        df_als$SNP_f <- factor(df_als$SNP, levels=0:2, labels=als_labels)

        boxplot(expression ~ SNP_f, data=df_als, col="lightpink", outline=FALSE,
                ylab="circRNA Percentage (%)", xlab="Genotype", main=paste(case_label, "Samples Only"))
        points(jitter(as.numeric(df_als$SNP_f), amount=0.15), df_als$expression,
               pch=df_als$p_pch, col=df_als$p_col, cex=df_als$p_cex)

        # 3. Interaction Slope Visualization
        plot(df$SNP, df$expression, type="n", xaxt="n", xlab="Genotype",
             ylab="circRNA Percentage (%)", main="Interaction Slope Context",
             sub=paste("Interaction p =", get_interaction_p(df)), col.sub="darkblue", font.sub=2)

        axis(1, at=0:2, labels=c("Ref/Ref", "Het", "Hom Alt"))
        points(jitter(df$SNP, amount=0.08), df$expression,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        tryCatch({
            abline(lm(expression ~ SNP, data=df_ctrl), col="#9932CC", lwd=3)
            abline(lm(expression ~ SNP, data=df_als), col="red", lwd=3)
        }, error=function(e) NULL)

        legend("topright", legend=c(paste(case_label, "Slope"), paste(control_label, "Slope")),
               col=c("red", "#9932CC"), lty=1, lwd=3, bty="n", cex=0.8)

        mtext(
            paste("circRNA:", c_id, "| SNP:", S, "| REF:", a$ref[1], "| ALT:", a$alt[1]),
            outer=TRUE, cex=1.1, font=2, line=0.5
        )

        plots_made <- plots_made + 1
    }

    dev.off()
    message(sprintf("# Interaction plots generated: %d, skipped: %d", plots_made, pairs_skipped))
}

message(paste("# Processing", nrow(pairs), "interaction circRNA-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
