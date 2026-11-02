#!/usr/bin/env Rscript
library(data.table)
library(ggplot2)
library(ggrepel)

# ── Parse Command Line Arguments ──────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Error: No tissue name provided. Usage: Rscript _cQTL_regional_zoom.R <TISSUE> [RUN_TYPE]", call. = FALSE)
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
target_chr <- "chr5"
start_pos  <- 80050000
end_pos    <- 83500000

# Derive clean chromosome label and position labels dynamically
chr_num_clean <- gsub("chr", "", target_chr, ignore.case = TRUE)
x_label_dynamic <- sprintf("Chromosome %s Position (Mb)", chr_num_clean)

# Base directory path setup
base_dir <- path.expand(file.path("~/donglab/data", experiment, tissue, "cQTL"))
base_dir <- paste0(base_dir, "/")

# Dynamically set results folder depending on the run type
if (run_type == "interaction") {
  results_dir <- paste0(base_dir, "interaction_results/")
} else {
  results_dir <- paste0(base_dir, "results/")
}

message("## Reading data matrix tables...")
snp_raw  <- fread(paste0(results_dir, tissue, "_cQTL.cis.txt"))
lead_raw <- fread(paste0(results_dir, tissue, "_cQTL.lead_snps.txt"))
snp_loc  <- fread(paste0(base_dir, "snp_location.txt"), select = c("snpid", "chr", "pos"))

# Standardize headers for circular RNA structures
setnames(snp_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(lead_raw, "SNP", "snpid", skip_absent = TRUE)
setnames(snp_raw, "gene", "circ_id", skip_absent = TRUE)
setnames(lead_raw, "gene", "circ_id", skip_absent = TRUE)

# Calculate log p-values using base assignment
snp_raw$log10p <- -log10(snp_raw$`p-value`)
lead_raw$log10p <- -log10(lead_raw$`p-value`)

# Filter location reference map straight against the clean chr column
snp_loc_sub <- snp_loc[chr == target_chr & pos >= start_pos & pos <= end_pos]

# Keyed join
setkey(snp_raw, snpid); setkey(lead_raw, snpid); setkey(snp_loc_sub, snpid)
snp_zoom  <- snp_raw[snp_loc_sub, nomatch = NULL]
lead_zoom <- lead_raw[snp_loc_sub, nomatch = NULL]

if (nrow(snp_zoom) == 0) stop("## Zero variants found in this window.")

# Find the absolute top signal circRNA in this window
top_circ <- snp_zoom[order(-log10p)][1, circ_id]
message(paste("## Top signal circRNA identified in this window:", top_circ))


# ========================================================================
# PLOT 1: THE REGIONAL ZOOM PLOT (ALL VARIANTS)
# ========================================================================
message("## Building Plot 1: Full variant zoom...")
snp_zoom$color_cat <- ifelse(snp_zoom$circ_id == top_circ & snp_zoom$`p-value` <= 0.05, "Top circRNA Locus", "Other/Background")
extreme_hits <- lead_zoom[order(-log10p)][1:min(3, .N)]

p1 <- ggplot() +
    geom_point(data = snp_zoom[color_cat == "Other/Background"], aes(x = pos / 1e6, y = log10p),
               colour = "#B0B0B0", size = 1.0, alpha = 0.4) +
    geom_point(data = snp_zoom[color_cat == "Top circRNA Locus"], aes(x = pos / 1e6, y = log10p),
               colour = "#E41A1C", size = 1.4, alpha = 0.8) +
    geom_point(data = lead_zoom, aes(x = pos / 1e6, y = log10p),
               shape = 18, size = 3.5, colour = "black") +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text_repel(data = extreme_hits, aes(x = pos / 1e6, y = log10p, label = circ_id),
                    size = 3.5, colour = "black", fontface = "bold", box.padding = 0.5) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.01, 0.01)) +
    labs(title = sprintf("%s Regional cis-cQTL Zoom (All Variants)", gsub("_", " ", tissue)),
         subtitle = sprintf("Region: %s:%.2f-%.2f Mb | Highlighted: %s", target_chr, start_pos/1e6, end_pos/1e6, top_circ),
         x = x_label_dynamic, y = expression(-log[10](p-value))) +
    theme_bw() + theme(panel.grid.minor = element_blank(), legend.position = "none")

out_file1 <- paste0(results_dir, tissue, "_cQTL.simple_zoom.png")
png(out_file1, width = 10, height = 6, units = "in", res = 300)
print(p1)
dev.off()


# ========================================================================
# PLOT 2: LEAD SNPS ONLY - PRUNED VIEW
# ========================================================================
message("## Building Plot 2: Lead variants summary structure...")

lead_snps_per_circ <- snp_zoom[order(circ_id, -log10p), head(.SD, 1), by = circ_id]
lead_snps_per_circ$status <- ifelse(lead_snps_per_circ$circ_id == top_circ, "Primary Locus Driver", "Secondary Target circRNA")

p2 <- ggplot(data = lead_snps_per_circ, aes(x = pos / 1e6, y = log10p)) +
    geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_point(aes(colour = status), size = 4.0, alpha = 0.9) +
    scale_colour_manual(values = c("Primary Locus Driver" = "#E41A1C", "Secondary Target circRNA" = "#377EB8")) +
    geom_text_repel(aes(label = circ_id), size = 3.2, fontface = "bold",
                    box.padding = 0.4, max.overlaps = 50) +
    scale_x_continuous(limits = c(start_pos/1e6, end_pos/1e6), expand = c(0.02, 0.02)) +
    labs(title = sprintf("%s Lead cQTL Variant Summary", gsub("_", " ", tissue)),
         subtitle = sprintf("Pruned view: Showing top peak variants per unique circRNA structure | Region: %s:%.2f-%.2f Mb", target_chr, start_pos/1e6, end_pos/1e6),
         x = x_label_dynamic, y = expression(-log[10](p-value))) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.title = element_blank())

out_file2 <- paste0(results_dir, tissue, "_cQTL.lead_snps_summary.png")
png(out_file2, width = 10, height = 6, units = "in", res = 300)
print(p2)
dev.off()


# ========================================================================
# PLOT 3: circRNA DIAGNOSTIC REGIONAL ZOOM (CENTERED ON rs430943)
# ========================================================================
message("## Building Plot 3: circRNA Diagnostic Zoom centered on rs430943...")

anchor_snp <- "rs430943"
anchor_lookup <- snp_loc[snpid == anchor_snp]

if (nrow(anchor_lookup) > 0) {
    anchor_pos <- anchor_lookup[1, pos]
    diag_start <- anchor_pos - 500000
    diag_end   <- anchor_pos + 500000

    snp_diag <- snp_raw[snpid %in% snp_loc[chr == target_chr & pos >= diag_start & pos <= diag_end, snpid]
    setkey(snp_diag, snpid)
    snp_diag <- snp_diag[snp_loc, nomatch = NULL]

    sig_circs <- unique(snp_diag[log10p >= 5, circ_id])
    snp_diag$color_group <- ifelse(snp_diag$circ_id %in% sig_circs & snp_diag$log10p >= 5, snp_diag$circ_id, "Background / Non-Sig")

    diag_leads <- lead_raw[snpid %in% snp_diag$snpid & circ_id %in% sig_circs]
    diag_leads <- diag_leads[snp_loc, nomatch = NULL, on = "snpid"]

    unique_sig_circs <- sort(setdiff(unique(snp_diag$color_group), "Background / Non-Sig"))
    color_palette <- setNames(scales::hue_pal()(length(unique_sig_circs)), unique_sig_circs)
    color_palette["Background / Non-Sig"] <- "#B0B0B0"

    p3 <- ggplot() +
        geom_point(data = snp_diag[color_group == "Background / Non-Sig"], aes(x = pos / 1e6, y = log10p),
                   colour = "#B0B0B0", size = 1.0, alpha = 0.3) +
        geom_point(data = snp_diag[color_group != "Background / Non-Sig"], aes(x = pos / 1e6, y = log10p, colour = color_group),
                   size = 1.5, alpha = 0.8) +
        geom_point(data = diag_leads, aes(x = pos / 1e6, y = log10p, colour = circ_id),
                   shape = 23, fill = "white", size = 3.0, stroke = 1.5) +
        geom_vline(xintercept = anchor_pos / 1e6, linetype = "dotted", colour = "#984EA3", linewidth = 1.0) +
        geom_hline(yintercept = 5, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
        geom_text_repel(data = diag_leads, aes(x = pos / 1e6, y = log10p, label = circ_id),
                        size = 3.2, fontface = "bold", box.padding = 0.5, max.overlaps = 15) +
        scale_colour_manual(values = color_palette) +
        scale_x_continuous(limits = c(diag_start / 1e6, diag_end / 1e6), expand = c(0.01, 0.01)) +
        labs(title = sprintf("circRNA Diagnostic Regional Zoom - %s cQTL", gsub("_", " ", tissue)),
             subtitle = sprintf("Region: Locus window around %s (%s:%.2f-%.2f Mb) | Colored by unique significant circular RNA structural locus",
                                anchor_snp, target_chr, diag_start/1e6, diag_end/1e6),
             x = x_label_dynamic, y = expression(-log[10](p-value))) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank(),
              legend.text = element_text(size = 8)) +
        guides(colour = guide_legend(ncol = 4, byrow = TRUE))

    out_file3 <- paste0(results_dir, tissue, "_cQTL.gene_diagnostics.png")
    png(out_file3, width = 11, height = 7, units = "in", res = 300)
    print(p3)
    dev.off()
    message(paste("## Saved Layout 3:", out_file3))
} else {
    message(sprintf("## Warning: Anchor SNP %s not found in location maps. Skipping Plot 3.", anchor_snp))
}

message("## Execution completed cleanly!")
