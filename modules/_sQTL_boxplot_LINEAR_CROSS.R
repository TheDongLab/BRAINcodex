#!/usr/bin/env Rscript
###########################################
# _sQTL_boxplot_LINEAR_CROSS.R
# Visualizes Genotype Groups vs Splicing Intron PSI Values Across ALS/Control Status
# Interaction Layout: Split Boxplots (Controls, ALS) + Interaction Slope Overlay
# Converts PLINK counted-allele dosage to ALT dosage for plotting
###########################################

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript _sQTL_boxplot_LINEAR_CROSS.R <pairs_file> <snp_file> <psi_file> <cov_file> <snp_loc> <out_dir> <tissue>")
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

out_file_std <- file.path(out_dir, paste0(tissue, "_all_sig_sQTL_Interaction_boxplots.pdf"))
out_file_col <- file.path(out_dir, paste0(tissue, "_all_sig_sQTL_Interaction_boxplots_colored.pdf"))

# Read 2-column sQTL pairs file (junction_id, snpid)
pairs <- fread(GSfile, header=FALSE, col.names=c("junction_id", "snpid"))

message("# Loading genomic matrices...")
snp_mat <- fread(snp_file, header=TRUE)
psi_mat <- fread(psi_file, header=TRUE)
cov_mat <- fread(cov_file, header=TRUE)
snp_map <- fread(snp_loc, header=TRUE)

# Standardize internal matrix headers
if ("geneid" %in% names(psi_mat)) setnames(psi_mat, "geneid", "junction_id")
setnames(psi_mat, names(psi_mat)[1], "junction_id")
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

message(sprintf("# Pairs to plot: %d", nrow(pairs)))
message(sprintf("# SNP matrix: %d SNPs x %d samples", nrow(snp_mat), ncol(snp_mat)-1))
message(sprintf("# Splicing (PSI) matrix: %d junctions x %d samples", nrow(psi_mat), ncol(psi_mat)-1))

# Process clinical interaction covariate status
cov_first_col <- names(cov_mat)[1]
als_row <- cov_mat[get(cov_first_col) == "is_als", ]

if (nrow(als_row) == 0) {
  message("# WARNING: 'is_als' row not found in covariates file! Defaulting to Control status.")
  als_vals <- rep(0, ncol(cov_mat) - 1)
} else {
  als_vals <- as.numeric(unlist(als_row[, -1, with=FALSE], use.names=FALSE))
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

        if (length(psi_vals) == 0 || length(snp_vals) == 0 || length(psi_vals) != length(snp_vals)) {
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
            splicing=psi_vals,
            SNP=snp_vals,
            Status=ifelse(als_vals > 0.5, "ALS", "Control"),
            p_col=if(use_status_colors) p_colors else rep("darkred", length(psi_vals)),
            p_pch=if(use_status_colors) p_shapes else rep(16, length(psi_vals)),
            p_cex=if(use_status_colors) p_sizes else rep(0.7, length(psi_vals)),
            stringsAsFactors=FALSE
        )

        df <- df[!is.na(df$SNP) & !is.na(df$splicing), ]

        if (nrow(df[df$Status == "ALS", ]) < 3 || nrow(df[df$Status == "Control", ]) < 3) {
            pairs_skipped <- pairs_skipped + 1
            next
        }

        # Interaction term p-value calculation (Genotype x Status)
        get_interaction_p <- function(data) {
            tryCatch({
                fit <- lm(splicing ~ SNP * factor(Status), data=data)
                p_val <- summary(fit)$coefficients["SNP:factor(Status)Control", 4]
                formatC(p_val, format="e", digits=3)
            }, error=function(e) "NA")
        }

        # 3-Panel Layout: Controls Boxplot, ALS Boxplot, Regression Slope Comparison
        par(mfrow=c(1,3), mar=c(5,4,4,1), oma=c(0,0,3,0))

        # 1. Controls Boxplot
        df_ctrl <- df[df$Status == "Control", ]
        ctrl_counts <- table(factor(df_ctrl$SNP, levels=0:2))
        ctrl_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", ctrl_counts, ")")
        df_ctrl$SNP_f <- factor(df_ctrl$SNP, levels=0:2, labels=ctrl_labels)

        boxplot(splicing ~ SNP_f, data=df_ctrl, col="lightblue", outline=FALSE,
                ylab="Splicing Level (PSI)", xlab="Genotype", main="Controls Only")
        points(jitter(as.numeric(df_ctrl$SNP_f), amount=0.15), df_ctrl$splicing,
               pch=df_ctrl$p_pch, col=df_ctrl$p_col, cex=df_ctrl$p_cex)

        # 2. ALS Boxplot
        df_als <- df[df$Status == "ALS", ]
        als_counts <- table(factor(df_als$SNP, levels=0:2))
        als_labels <- paste0(c("Ref/Ref", "Het", "Hom Alt"), "\n(N=", als_counts, ")")
        df_als$SNP_f <- factor(df_als$SNP, levels=0:2, labels=als_labels)

        boxplot(splicing ~ SNP_f, data=df_als, col="lightpink", outline=FALSE,
                ylab="Splicing Level (PSI)", xlab="Genotype", main="ALS Only")
        points(jitter(as.numeric(df_als$SNP_f), amount=0.15), df_als$splicing,
               pch=df_als$p_pch, col=df_als$p_col, cex=df_als$p_cex)

        # 3. Interaction Slope Visualization
        plot(df$SNP, df$splicing, type="n", xaxt="n", xlab="Genotype",
             ylab="Splicing Level (PSI)", main="Interaction Slope Context",
             sub=paste("Interaction p =", get_interaction_p(df)),
             col.sub="darkblue", font.sub=2)

        axis(1, at=0:2, labels=c("Ref/Ref", "Het", "Hom Alt"))

        # Jitter points
        points(jitter(df$SNP, amount=0.08), df$splicing,
               pch=df$p_pch, col=df$p_col, cex=df$p_cex)

        # Add regression fits per cohort
        tryCatch({
            abline(lm(splicing ~ SNP, data=df_ctrl), col="#9932CC", lwd=3)
            abline(lm(splicing ~ SNP, data=df_als), col="red", lwd=3)
        }, error=function(e) NULL)

        legend("topright", legend=c("ALS Slope", "Control Slope"),
               col=c("red", "#9932CC"), lty=1, lwd=3, bty="n", cex=0.8)

        # Page Header
        mtext(
            paste("Junction:", G, "| SNP:", S, "| REF:", a$ref[1], "| ALT:", a$alt[1]),
            outer=TRUE, cex=1.1, font=2, line=0.5
        )

        plots_made <- plots_made + 1
    }

    dev.off()
    message(sprintf("# Interaction sQTL plots generated: %d, skipped: %d", plots_made, pairs_skipped))
}

message(paste("# Processing", nrow(pairs), "interaction junction-SNP pairs..."))
run_plotting(out_file_std, FALSE)
run_plotting(out_file_col, TRUE)
message("Done.")
