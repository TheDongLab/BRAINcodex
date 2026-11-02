#!/usr/bin/env Rscript
###########################################
# _sQTL_boxplot.R
# Visualizes Genotype Groups vs Splicing Intron PSI Values
# Converts PLINK counted-allele dosage to ALT dosage for plotting
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript _sQTL_boxplot.R <pairs_file> <snp_file> <psi_file> <cov_file> <snp_loc> <out_dir> <tissue>")
}

GSfile    <- args[1]
snp_file  <- args[2]
psi_file  <- args[3]
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

# Un-hash these rows to generate boxplots for ALL significant SNPs:
out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_sQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_sQTL_boxplots_colored.pdf"))
pairs <- fread(GSfile, header=FALSE, col.names=c("junction_id","snpid"))

# Un-hash these rows to generate boxplots for the top 1000 subset:
# out_file_std <- file.path(out_dir, paste0(tissue, "_top1000_sig_sQTL_boxplots.pdf"))
# out_file_col <- file.path(out_dir, paste0(tissue, "_top1000_sig_sQTL_boxplots_colored.pdf"))
# pairs <- fread(GSfile, header=FALSE, col.names=c("junction_id","snpid"), nrows=1000)

message("# Loading genomic matrices...")

snp_mat <- fread(snp_file, header=TRUE)
psi_mat <- fread(psi_file, header=TRUE)
cov_mat <- fread(cov_file, header=TRUE)
snp_map <- fread(snp_loc, header=TRUE)

if ("geneid" %in% names(psi_mat)) setnames(psi_mat, "geneid", "junction_id")

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

message(sprintf("# Pairs to plot: %d", nrow(pairs)))
message(sprintf("# SNP matrix: %d SNPs x %d samples", nrow(snp_mat), ncol(snp_mat)-1))
message(sprintf("# Splicing (PSI) matrix: %d junctions x %d samples", nrow(psi_mat), ncol(psi_mat)-1))

# Get the first column name (the row identifier)
cov_first_col <- names(cov_mat)[1]

# Find the 'is_als' row
als_row <- cov_mat[get(cov_first_col) == "is_als", ]

if (nrow(als_row) == 0) {
  message("# WARNING: 'is_als' row not found. Defaulting to standard purple points.")
  als_vals <- rep(0, ncol(cov_mat) - 1)
} else {
  als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
}

p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)

run_plotting <- function(pdf_path, use_status_colors = FALSE) {
    pdf(pdf_path, width=12, height=5)

    plots_made <- 0
    pairs_skipped <- 0

    for (i in seq_len(nrow(pairs))) {
        G <- pairs$junction_id[i]
        S <- pairs$snpid[i]

        psi_row <- psi_mat[junction_id == G]
        snp_row <- snp_mat[snpid == S]

        if (nrow(psi_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        psi_vals <- as.numeric(unlist(psi_row[, -1, with=FALSE], use.names=FALSE))
        snp_vals <- as.numeric(unlist(snp_row[, -1, with=FALSE], use.names=FALSE))

        if (length(psi_vals) == 0 || length(snp_vals) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        if (length(psi_vals) != length(snp_vals)) {
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
            splicing = psi_vals,
            SNP = snp_vals,
            p_col = if(use_status_colors) p_colors else rep("darkred", length(psi_vals)),
            p_pch = if(use_status_colors) p_shapes else rep(16, length(psi_vals)),
            p_cex = if(use_status_colors) p_sizes else rep(0.7, length(psi_vals)),
            stringsAsFactors = FALSE
        )

        df <- df[!is.na(df$SNP) & !is.na(df$splicing), ]

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

        # Model 1: Additive
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)

        boxplot(splicing ~ SNP_f, data=df, col="lightgreen", outline=FALSE,
                ylab="Splicing Level (PSI)", xlab="Genotype",
                main=paste0("Additive Model\n(p = ", get_p(splicing ~ SNP, df), ")"))

        points(jitter(as.numeric(df$SNP_f), amount=0.15), df$splicing,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        if (use_status_colors) {
          legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
                 pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)
        }

        # Model 2: Dominant
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)

        boxplot(splicing ~ dom_f, data=df, col="lightblue", outline=FALSE,
                ylab="Splicing Level (PSI)", xlab="Genotype",
                main=paste0("Dominant Model\n(p = ", get_p(splicing ~ dom_val, df), ")"))

        points(jitter(as.numeric(df$dom_f), amount=0.15), df$splicing,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        if (use_status_colors) {
          legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
                 pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)
        }

        # Model 3: Recessive
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)

        boxplot(splicing ~ rec_f, data=df, col="lightpink", outline=FALSE,
                ylab="Splicing Level (PSI)", xlab="Genotype",
                main=paste0("Recessive Model\n(p = ", get_p(splicing ~ rec_val, df), ")"))

        points(jitter(as.numeric(df$rec_f), amount=0.15), df$splicing,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        if (use_status_colors) {
          legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
                 pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)
        }

        mtext(
            paste("Junction:", G, "| SNP:", S, "| REF:", a$ref[1], "| ALT:", a$alt[1]),
            outer=TRUE, cex=1.1, font=2, line=0.5
        )

        plots_made <- plots_made + 1
    }

    dev.off()
    message(sprintf("# Plots made: %d, pairs skipped: %d", plots_made, pairs_skipped))
}

# --- EXECUTION BLOCK ---
message(paste("# Processing", nrow(pairs), "junction-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
