#!/usr/bin/env Rscript
circ_mode <- identical(Sys.getenv("EQTL_MODE", "gene"), "circ")
###########################################
# _eQTL_boxplot.R
# Lightweight allele-aware eQTL boxplots
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript _eQTL_boxplot.R <pairs_file> <snp_file> <expr_file> <cov_file> <snp_loc> <out_dir> <tissue>")
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
if (circ_mode) RAW_FILE <- Sys.getenv("EQTL_RAW_FILE", RAW_FILE)

out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_eQTL_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_eQTL_boxplots_colored.pdf"))

pairs <- fread(GSfile, header=TRUE)
if (circ_mode) {
    if (nrow(pairs) == 0) { message("No significant circRNA pairs to plot"); quit(status=0) }
    annotation <- fread(file.path(dirname(expr_file), "circ_annotation.tsv"))
    pairs <- merge(pairs, annotation, by="geneid", all.x=TRUE, sort=FALSE)
    pairs[, gene_symbol := ifelse(is.na(host_gene) | host_gene == "", geneid, paste0(host_gene, " | ", geneid))]
}
phenotype_label <- if (circ_mode) "circRNA RPM (inverse-normal transformed)" else "Expression (Z-score)"

message("# Loading genomic matrices...")
snp_mat  <- fread(snp_file, header=TRUE)
expr_mat <- fread(expr_file, header=TRUE)
cov_mat  <- fread(cov_file, header=TRUE)
snp_map  <- fread(snp_loc, header=TRUE)

setnames(expr_mat, names(expr_mat)[1], "geneid")
setnames(snp_mat, names(snp_mat)[1], "snpid")

# ── Recover REF, ALT and PLINK-counted allele from RAW header ─────────
raw_cols <- scan(RAW_FILE, what=character(), nlines=1, quiet=TRUE)
raw_cols <- raw_cols[grepl("^[^:]+:[0-9]+:[ACGT]+:[ACGT]+_[ACGT]+$", raw_cols)]

alleles <- data.table(raw_col=raw_cols)
alleles[, c("variant", "counted") := tstrsplit(raw_col, "_", fixed=TRUE)]
alleles[, c("chr", "pos", "ref", "alt") := tstrsplit(variant, ":", fixed=TRUE)]
alleles[, chr := sub("^chr", "", chr)]
snp_map[, `:=`(chr_key=sub("^chr", "", as.character(chr)), pos_key=as.character(pos))]
alleles <- merge(snp_map[, .(snpid, chr_key, pos_key)],
                 alleles[, .(chr_key=chr, pos_key=pos, ref, alt, counted)],
                 by=c("chr_key", "pos_key"), all.x=TRUE)

message(sprintf("# Allele mappings recovered: %d / %d SNPs",
                sum(!is.na(alleles$counted)), nrow(alleles)))

# ── Process covariates ────────────────────────────────────────────────
cov_first_col <- names(cov_mat)[1]
als_row <- cov_mat[get(cov_first_col) == "is_als", ]

if (nrow(als_row) == 0) {
  als_vals <- rep(0, ncol(cov_mat) - 1)
} else {
  als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
}

p_colors <- ifelse(als_vals > 0.5, "red", "#9932CC")
p_shapes <- ifelse(als_vals > 0.5, 16, 18)
p_sizes  <- ifelse(als_vals > 0.5, 0.7, 1.1)

run_plotting <- function(pdf_path, use_status_colors=FALSE) {
    pdf(pdf_path, width=12, height=5)
    plots_made <- 0
    pairs_skipped <- 0

    for (i in seq_len(nrow(pairs))) {
        ens_id <- pairs$geneid[i]
        sym_id <- pairs$gene_symbol[i]
        S      <- pairs$snpid[i]

        expr_row <- expr_mat[geneid == ens_id]
        snp_row  <- snp_mat[snpid == S]

        if (nrow(expr_row) == 0 || nrow(snp_row) == 0) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        expr_vals <- as.numeric(unlist(expr_row[, -1, with=FALSE], use.names=FALSE))
        snp_vals  <- as.numeric(unlist(snp_row[, -1, with=FALSE], use.names=FALSE))

        if (length(expr_vals) == 0 || length(snp_vals) == 0 || length(expr_vals) != length(snp_vals)) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        # Convert PLINK counted-allele dosage to ALT dosage.
        a <- alleles[snpid == S]
        if (nrow(a) != 1 || is.na(a$counted[1]) || !(a$counted[1] %in% c(a$ref[1], a$alt[1]))) {
            message("# WARNING: allele orientation unavailable for ", S, "; skipping.")
            pairs_skipped <- pairs_skipped + 1
            next
        }
        if (a$counted[1] == a$ref[1]) snp_vals <- 2 - snp_vals

        df <- data.frame(
            expression=expr_vals,
            SNP=snp_vals,
            p_col=if(use_status_colors) p_colors else rep("darkred", length(expr_vals)),
            p_pch=if(use_status_colors) p_shapes else rep(16, length(expr_vals)),
            p_cex=if(use_status_colors) p_sizes else rep(0.7, length(expr_vals)),
            stringsAsFactors=FALSE
        )

        df <- df[!is.na(df$SNP) & !is.na(df$expression), ]

        if (nrow(df) < 5) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        get_p <- function(formula, data) {
            tryCatch({
                fit <- lm(formula, data=data)
                formatC(summary(fit)$coefficients[2,4], format="e", digits=3)
            }, error=function(e) "NA")
        }

        par(mfrow=c(1,3), mar=c(5,4,4,2), oma=c(0,0,3,0))

        # Additive Model
        add_counts <- table(factor(df$SNP, levels=0:2))
        add_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", add_counts, ")")
        df$SNP_f <- factor(df$SNP, levels=0:2, labels=add_labels)

        boxplot(expression ~ SNP_f, data=df, col="lightgreen", outline=FALSE,
                ylab=phenotype_label, xlab="Genotype",
                main=paste0("Additive Model\n(p = ", get_p(expression ~ SNP, df), ")"))

        points(jitter(as.numeric(df$SNP_f), amount=0.15), df$expression,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        if (use_status_colors)
            legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
                   pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

        # Dominant ALT Model
        df$dom_val <- ifelse(df$SNP > 0, 1, 0)
        dom_counts <- table(factor(df$dom_val, levels=0:1))
        dom_labels <- paste0(c("Ref/Ref", "Any Alt"), "\n(N=", dom_counts, ")")
        df$dom_f <- factor(df$dom_val, levels=0:1, labels=dom_labels)

        boxplot(expression ~ dom_f, data=df, col="lightblue", outline=FALSE,
                ylab=phenotype_label, xlab="Genotype",
                main=paste0("Dominant Model\n(p = ", get_p(expression ~ dom_val, df), ")"))

        points(jitter(as.numeric(df$dom_f), amount=0.15), df$expression,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        if (use_status_colors)
            legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
                   pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

        # Recessive ALT Model
        df$rec_val <- ifelse(df$SNP == 2, 1, 0)
        rec_counts <- table(factor(df$rec_val, levels=0:1))
        rec_labels <- paste0(c("Non-Hom Alt", "Hom Alt"), "\n(N=", rec_counts, ")")
        df$rec_f <- factor(df$rec_val, levels=0:1, labels=rec_labels)

        boxplot(expression ~ rec_f, data=df, col="lightpink", outline=FALSE,
                ylab=phenotype_label, xlab="Genotype",
                main=paste0("Recessive Model\n(p = ", get_p(expression ~ rec_val, df), ")"))

        points(jitter(as.numeric(df$rec_f), amount=0.15), df$expression,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        if (use_status_colors)
            legend("topleft", legend=c("ALS", "Control"), col=c("red", "#9932CC"),
                   pch=c(16,18), pt.cex=c(0.7,1.1), bty="n", cex=0.8)

        mtext(
            paste("Gene:", sym_id, "(", ens_id, ") | SNP:", S,
                  "| REF:", a$ref[1], "| ALT:", a$alt[1]),
            outer=TRUE, cex=1.0, font=2, line=0.5
        )

        plots_made <- plots_made + 1
    }

    dev.off()
    message(sprintf("# Plots made: %d, pairs skipped: %d", plots_made, pairs_skipped))
}

message(paste("# Processing", nrow(pairs), "gene-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
