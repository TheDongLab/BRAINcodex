#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggrepel)

# ── Parse Command Line Arguments ──────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No tissue name provided. Usage: Rscript _sQTL_regional_zoom.R <TISSUE> [RUN_TYPE]", call. = FALSE)
}
tissue <- args[1]

# Safely capture run_type if passed by the bash script; default to standard
if (length(args) >= 2) {
  run_type <- args[2]
} else {
  run_type <- "standard"
}

# TargetALS remains the default for legacy runs.
experiment <- Sys.getenv("EXPERIMENT", "target_ALS")
if (!experiment %in% c("target_ALS", "answer_ALS")) {
  stop("EXPERIMENT must be target_ALS or answer_ALS.", call. = FALSE)
}

message(paste("## Running dynamic regional zoom plots for tissue:", tissue, paste0("(", run_type, " mode; ", experiment, ")")))

# Target window setup for Plot 1 & Plot 2
target_chr <- "17"
start_pos  <- 45465000
end_pos    <- 47055000

# Derive clean chromosome and position labels dynamically
chr_num_clean   <- gsub("chr", "", as.character(target_chr), ignore.case = TRUE)
chr_label_full  <- paste0("chr", chr_num_clean)
x_label_dynamic <- sprintf("Chromosome %s Position (Mb)", chr_num_clean)

# Base directory path setup for sQTL
base_dir <- path.expand(file.path("~/donglab/data", experiment, tissue, "sQTL"))
base_dir <- paste0(base_dir, "/")

# Dynamically set results folder depending on the run type
if (run_type == "interaction") {
  results_dir <- paste0(base_dir, "interaction_results/")
} else {
  results_dir <- paste0(base_dir, "results/")
}

message("## Reading data matrix tables...")
snp_raw  <- fread(paste0(results_dir, tissue, "_sQTL.cis.txt"))
lead_raw <- fread(paste0(results_dir, tissue, "_sQTL.lead_snps.txt"))
snp_loc  <- fread(paste0(base_dir, "snp_location.txt"), select = c("snpid", "chr", "pos"))

# Standardize headers for sQTL (handling phenotype, junction, feature, or gene ID columns)
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)

feature_columns <- c(
  "gene", "geneid", "junction_id", "phenotype",
  "junction", "feature", "intron", "source_gene"
)

for (col in feature_columns) {
  if (col %in% names(snp_raw)) {
    setnames(snp_raw, col, "featureid")
    break
  }
}

for (col in feature_columns) {
  if (col %in% names(lead_raw)) {
    setnames(lead_raw, col, "featureid")
    break
  }
}

# Calculate log p-values
snp_raw[, log10p := -log10(`p-value`)]
lead_raw[, log10p := -log10(`p-value`)]

# Filter location reference map first to save memory overhead
snp_loc[, chr_clean := gsub("chr", "", as.character(chr), ignore.case = TRUE)]
snp_loc_sub <- snp_loc[chr_clean == chr_num_clean & pos >= start_pos & pos <= end_pos]

# Keyed join
setkey(snp_raw, snpid); setkey(lead_raw, snpid); setkey(snp_loc_sub, snpid)
snp_zoom  <- snp_raw[snp_loc_sub, nomatch = NULL]
lead_zoom <- lead_raw[snp_loc_sub, nomatch = NULL]

if (nrow(snp_zoom) == 0) stop("## Zero variants found in this window.")

# Find the absolute top signal splicing feature/cluster in this window
top_feature <- snp_zoom[order(-log10p)][1, featureid]
message(paste("## Top signal splicing feature identified in this window:", top_feature))


# ========================================================================
# PLOT 1: THE REGIONAL ZOOM PLOT (ALL VARIANTS)
# ========================================================================
message("## Building Plot 1: Full variant sQTL zoom...")
snp_zoom[, color_cat := ifelse(featureid == top_feature & `p-value` <= 0.05, "Top Splicing Event Locus", "Other/Background")]
extreme_hits <- lead_zoom[order(-log10p)][1:min(3, .N)]

p1 <- ggplot() +
    geom_point(data = snp_zoom[color_cat == "Other/Background"], aes(x = pos / 1e6, y = log10p),
               colour = "#B0B0B0", size = 1.0, alpha = 0.4) +
    geom_point(data = snp_zoom[color_cat == "Top Splicing Event Locus"], aes(x = pos / 1e6, y = log10p),
               colour = "#E41A1C", size = 1.4, alpha = 0.8) +
    geom_point(data = lead_zoom, aes(x = pos / 1e6, y = log10p),
               shape = 18, size = 3.5, colour = "black") +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text_repel(data = extreme_hits, aes(x = pos / 1e6, y = log10p, label = featureid),
                    size = 3.5, colour = "black", fontface = "bold", box.padding = 0.5) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.01, 0.01)) +
    labs(title = sprintf("%s Regional cis-sQTL Zoom (All Variants)", gsub("_", " ", tissue)),
         subtitle = sprintf("Region: %s:%.2f-%.2f Mb | Highlighted Splicing Feature: %s", chr_label_full, start_pos/1e6, end_pos/1e6, top_feature),
         x = x_label_dynamic, y = expression(-log[10](p-value))) +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

out_file1 <- paste0(results_dir, tissue, "_sQTL.simple_zoom.png")
png(out_file1, width = 10, height = 6, units = "in", res = 300)
print(p1)
dev.off()


# ========================================================================
# PLOT 2: LEAD SNPS ONLY - PRUNED VIEW
# ========================================================================
message("## Building Plot 2: Lead variants summary structure...")

lead_snps_per_feature <- snp_zoom[order(featureid, -log10p), head(.SD, 1), by = featureid]
lead_snps_per_feature[, status := ifelse(featureid == top_feature, "Primary Splicing Locus Driver", "Secondary Splicing Event")]

p2 <- ggplot(data = lead_snps_per_feature, aes(x = pos / 1e6, y = log10p)) +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = status), size = 4.0, alpha = 0.9) +
    scale_colour_manual(values = c("Primary Splicing Locus Driver" = "#E41A1C", "Secondary Splicing Event" = "#377EB8")) +
    geom_text_repel(aes(label = featureid), size = 3.2, fontface = "bold",
                    box.padding = 0.4, max.overlaps = 20, cluster_groups = FALSE) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.02, 0.02)) +
    labs(title = sprintf("%s Lead sQTL Variant Summary", gsub("_", " ", tissue)),
         subtitle = sprintf("Pruned view: Showing top peak variants per unique splicing feature | Region: %s:%.2f-%.2f Mb", chr_label_full, start_pos/1e6, end_pos/1e6),
         x = x_label_dynamic, y = expression(-log[10](p-value))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.title = element_blank())

out_file2 <- paste0(results_dir, tissue, "_sQTL.lead_snps_summary.png")
png(out_file2, width = 10, height = 6, units = "in", res = 300)
print(p2)
dev.off()


# ========================================================================
# PLOT 3: SPLICING FEATURE DIAGNOSTIC REGIONAL ZOOM (Centered on anchor SNP)
# ========================================================================
message("## Building Plot 3: Splicing Diagnostic Zoom centered on anchor SNP...")

anchor_snp <- "rs62056809"
anchor_lookup <- snp_loc[snpid == anchor_snp]

if (nrow(anchor_lookup) > 0) {
    anchor_pos <- anchor_lookup[1, pos]
    diag_start <- anchor_pos - 500000
    diag_end   <- anchor_pos + 500000

    snp_diag <- snp_raw[snpid %in% snp_loc[chr_clean == chr_num_clean & pos >= diag_start & pos <= diag_end, snpid]]
    setkey(snp_diag, snpid)
    snp_diag <- snp_diag[snp_loc, nomatch = NULL]

    sig_features <- unique(snp_diag[log10p >= 5, featureid])
    snp_diag[, color_group := ifelse(featureid %in% sig_features & log10p >= 5, featureid, "Background / Non-Sig")]

    diag_leads <- lead_raw[snpid %in% snp_diag$snpid & featureid %in% sig_features]
    diag_leads <- diag_leads[snp_loc, nomatch = NULL, on = "snpid"]

    unique_sig_features <- sort(setdiff(unique(snp_diag$color_group), "Background / Non-Sig"))

    if (length(unique_sig_features) == 0) {
        message("## No significant splicing features near anchor SNP. Skipping Plot 3.")
    } else {
    color_palette <- setNames(scales::hue_pal()(length(unique_sig_features)), unique_sig_features)
    color_palette["Background / Non-Sig"] <- "#B0B0B0"

    p3 <- ggplot() +
        geom_point(data = snp_diag[color_group == "Background / Non-Sig"], aes(x = pos / 1e6, y = log10p),
                   colour = "#B0B0B0", size = 1.0, alpha = 0.3) +
        geom_point(data = snp_diag[color_group != "Background / Non-Sig"], aes(x = pos / 1e6, y = log10p, colour = color_group),
                   size = 1.5, alpha = 0.8) +
        geom_point(data = diag_leads, aes(x = pos / 1e6, y = log10p, colour = featureid),
                   shape = 23, fill = "white", size = 3.0, stroke = 1.5) +
        geom_vline(xintercept = anchor_pos / 1e6, linetype = "dotted", colour = "#984EA3", linewidth = 1.0) +
        geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
        geom_text_repel(data = diag_leads, aes(x = pos / 1e6, y = log10p, label = featureid),
                        size = 3.2, fontface = "bold", box.padding = 0.5, max.overlaps = 15) +
        scale_colour_manual(values = color_palette) +
        scale_x_continuous(limits = c(diag_start / 1e6, diag_end / 1e6), expand = c(0.01, 0.01)) +
        labs(title = sprintf("Splicing Feature Diagnostic Regional Zoom - %s sQTL", gsub("_", " ", tissue)),
             subtitle = sprintf("Region: Locus window around %s (%s:%.2f-%.2f Mb) | Colored by unique significant splicing feature",
                                anchor_snp, chr_label_full, diag_start/1e6, diag_end/1e6),
             x = x_label_dynamic, y = expression(-log[10](p-value))) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 8)) +
        guides(colour = guide_legend(ncol = 4, byrow = TRUE))

    out_file3 <- paste0(results_dir, tissue, "_sQTL.gene_diagnostics.png")
    png(out_file3, width = 11, height = 7, units = "in", res = 300)
    print(p3)
    dev.off()
        message(paste("## Saved Layout 3:", out_file3))
    }
} else {
    message(sprintf("## Warning: Anchor SNP %s not found in location maps. Skipping Plot 3.", anchor_snp))
}

message("## Execution completed cleanly!")
