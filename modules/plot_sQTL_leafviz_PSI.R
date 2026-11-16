#!/usr/bin/env Rscript

# Static LeafViz-style plot. Existing calls default to genotype mode unchanged.
# Opt in: --group-by tdp43 (no genotype or variant arguments required).
# Optional --tdp43 <aligned CSV> and --pathology <original pathology CSV>; defaults are aligned_rnaseq_tdp43_by_subject.csv and collections.postmortem_tissue_core.semiquantitative_tdp43_data.csv under --data-root.
# TDP43 mode matches subject + RNA sample + tissue, verifies pathology provenance, selects one highest-valid-RIN sample per subject before inspecting event coverage, excludes missing/conflicting/unsupported scores, and writes mapping/count audits.
# Arcs show unweighted mean original junction PSI, never normalized to displayed junctions.
# Optional --psi-matrix: override matching junctions using QTL matrix subject columns.
# The matrix must contain the anchor; other displayed junctions use source-file PSI.
# PSI means junction usage under its ORIGINAL denominator, not inferred exon inclusion.
# No installed leafcutter R package is required.

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(x) {
  out <- list()
  i <- 1
  while (i <= length(x)) {
    key <- x[i]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    if (i == length(x)) stop("Missing value for ", key)
    out[[sub("^--", "", key)]] <- x[i + 1]
    i <- i + 2
  }
  out
}

a <- parse_args(args)
a$prefix <- paste0(a$prefix, ".PSI")
group_by <- if (is.null(a[["group-by"]])) "genotype" else a[["group-by"]]
if (!group_by %in% c("genotype", "tdp43")) stop("--group-by must be genotype or tdp43")
pathology_mode <- group_by == "tdp43"

required <- c(
  "tissue", "tissue-regex", "anchor", "gene",
  "variant-chr", "variant-pos", "variant-id", "variant-ref", "variant-alt",
  "genotypes", "metadata", "data-root", "gtf", "outdir", "prefix"
)

if (pathology_mode) {
  required <- setdiff(required, c("genotypes", "variant-chr", "variant-pos",
                                 "variant-id", "variant-ref", "variant-alt"))
}
missing_args <- required[!required %in% names(a)]
if (length(missing_args)) {
  stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
}

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("R package 'gridExtra' is required.")
  }
})

library(ggplot2)
library(gridExtra)
library(grid)

dir.create(a$outdir, recursive = TRUE, showWarnings = FALSE)

cat("Loading metadata...\n")
meta <- read.csv(a$metadata, stringsAsFactors = FALSE, check.names = FALSE)

need_meta <- c("externalsampleid", "externalsubjectid", "tissue")
if (!all(need_meta %in% colnames(meta))) {
  stop(
    "Metadata must contain columns: ",
    paste(need_meta, collapse = ", ")
  )
}

meta <- meta[
  !is.na(meta$tissue) &
    grepl(a[["tissue-regex"]], meta$tissue, ignore.case = TRUE, perl = TRUE),
  ,
  drop = FALSE
]

if (!nrow(meta)) {
  stop("No metadata rows matched tissue pattern for ", a$tissue)
}

# File names/directories use underscores rather than dashes.
meta$sample_dir <- gsub("-", "_", meta$externalsampleid, fixed = TRUE)

# Only pathology mode enters this branch. Genotype joins/selection stay unchanged.
if (pathology_mode) {
  canonical_tissue <- function(x) {
    x <- tolower(gsub("[^a-z0-9]", "", tolower(trimws(x))))
    x[grepl("cerebell", x)] <- "cerebellum"
    x[grepl("frontal", x)] <- "frontalcortex"
    x[grepl("motor.*cortex|cortex.*motor|ba4", x)] <- "motorcortex"
    x[grepl("cervical", x)] <- "cervicalspinalcord"
    x[grepl("lumbar|lumbosacral", x)] <- "lumbarspinalcord"
    x[grepl("thoracic", x)] <- "thoracicspinalcord"
    x
  }
  score_levels <- c("Absent", "Sparse", "Moderate", "Frequent")
  score_value <- function(x) score_levels[match(tolower(trimws(x)), tolower(score_levels))]
  key <- function(...) do.call(paste, c(list(...), sep="\034"))
  need_columns <- function(x, cols, label) {
    if (!all(cols %in% names(x))) stop(label, " missing columns: ", paste(setdiff(cols,names(x)),collapse=", "))
  }
  if (is.null(a$tdp43)) a$tdp43 <- file.path(a[["data-root"]], "aligned_rnaseq_tdp43_by_subject.csv")
  if (is.null(a$pathology)) a$pathology <- file.path(a[["data-root"]], "collections.postmortem_tissue_core.semiquantitative_tdp43_data.csv")
  aligned <- read.csv(a$tdp43, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
  raw <- read.csv(a$pathology, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8-BOM")
  need_columns(aligned,c("Subject_ID","RNAseq_Sample_ID","RNAseq_Tissue_Name","Pathology_Sample_ID","Pathology_Tissue_Name","Neuronal_TDP43_Score"),"Aligned pathology table")
  need_columns(raw,c("Subject ID","Sample ID","Tissue Source","P Tdp 43 Inclusions Neuronal"),"Original pathology table")
  for (n in names(aligned)) aligned[[n]] <- trimws(as.character(aligned[[n]]))
  for (n in names(raw)) raw[[n]] <- trimws(as.character(raw[[n]]))
  meta$externalsubjectid <- trimws(meta$externalsubjectid)
  meta$externalsampleid <- trimws(meta$externalsampleid)
  meta$sample_dir <- gsub("-", "_", meta$externalsampleid, fixed=TRUE)
  wanted_tissue <- canonical_tissue(a$tissue)
  meta <- meta[canonical_tissue(meta$tissue)==wanted_tissue, , drop=FALSE]
  aligned <- aligned[!is.na(aligned$RNAseq_Tissue_Name) & canonical_tissue(aligned$RNAseq_Tissue_Name)==wanted_tissue, , drop=FALSE]
  aligned$normalized_score <- score_value(aligned$Neuronal_TDP43_Score)
  aligned$reason <- "verified"
  raw_key <- key(raw[["Subject ID"]], raw[["Sample ID"]], canonical_tissue(raw[["Tissue Source"]]))
  raw_groups <- split(seq_len(nrow(raw)),raw_key)
  meta_key <- key(meta$externalsubjectid,meta$externalsampleid,canonical_tissue(meta$tissue))
  aligned_key <- key(aligned$Subject_ID,aligned$RNAseq_Sample_ID,canonical_tissue(aligned$RNAseq_Tissue_Name))
  for (i in seq_len(nrow(aligned))) {
    r <- aligned[i, ]
    if (is.na(r$normalized_score)) { aligned$reason[i] <- "missing_or_unsupported_score"; next }
    if (any(is.na(r[c("Subject_ID","RNAseq_Sample_ID","Pathology_Sample_ID","Pathology_Tissue_Name")])) ||
        any(!nzchar(unlist(r[c("Subject_ID","RNAseq_Sample_ID","Pathology_Sample_ID","Pathology_Tissue_Name")])))) {
      aligned$reason[i] <- "missing_mapping_identifier"; next
    }
    if (canonical_tissue(r$Pathology_Tissue_Name)!=wanted_tissue) { aligned$reason[i] <- "pathology_tissue_mismatch"; next }
    if (!aligned_key[i] %in% meta_key) { aligned$reason[i] <- "no_exact_metadata_match"; next }
    ri <- raw_groups[[key(r$Subject_ID,r$Pathology_Sample_ID,canonical_tissue(r$Pathology_Tissue_Name))]]
    if (!length(ri)) { aligned$reason[i] <- "pathology_record_not_found"; next }
    raw_scores <- unique(score_value(raw[["P Tdp 43 Inclusions Neuronal"]][ri]))
    if (length(raw_scores)!=1L || is.na(raw_scores[1]) || raw_scores[1]!=r$normalized_score) {
      aligned$reason[i] <- "pathology_score_mismatch_or_conflict"
    }
  }
  # Different scores within a subject/tissue cannot become independent subjects.
  by_subject <- split(seq_len(nrow(aligned)),aligned$Subject_ID)
  for (ii in by_subject) {
    v <- unique(aligned$normalized_score[ii][!is.na(aligned$normalized_score[ii])])
    if (length(v)>1 || any(aligned$reason[ii]=="pathology_score_mismatch_or_conflict")) {
      aligned$reason[ii] <- "conflicting_subject_tissue_scores"
    }
  }
  write.table(aligned,file.path(a$outdir,paste0(a$prefix,".pathology_alignment.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  meta$group <- NA_character_
  meta$Pathology_Sample_ID <- NA_character_
  meta$Pathology_Tissue_Name <- NA_character_
  meta$mapping_status <- "no_verified_pathology_score"
  for (i in seq_len(nrow(meta))) {
    ii <- which(aligned_key==meta_key[i])
    if (!length(ii)) next
    good <- ii[aligned$reason[ii]=="verified"]
    if (!length(good)) { meta$mapping_status[i] <- paste(sort(unique(aligned$reason[ii])),collapse=";"); next }
    meta$group[i] <- aligned$normalized_score[good[1]]
    meta$Pathology_Sample_ID[i] <- paste(sort(unique(aligned$Pathology_Sample_ID[good])),collapse=";")
    meta$Pathology_Tissue_Name[i] <- paste(sort(unique(aligned$Pathology_Tissue_Name[good])),collapse=";")
    meta$mapping_status[i] <- "eligible"
  }
  meta$.mapping_row <- seq_len(nrow(meta))
  pathology_mapping <- meta
  meta <- meta[!is.na(meta$group), , drop=FALSE]
  meta$GT <- NA_character_ # Internal compatibility only; omitted from pathology output.
  if (!nrow(meta)) {
    write.table(pathology_mapping,file.path(a$outdir,paste0(a$prefix,".pathology_mapping.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    stop("No RNA-seq samples have a verified, unambiguous tissue-matched neuronal TDP43 score. See pathology audits.")
  }
} else {
cat("Loading VCF genotypes...\n")
gt <- read.delim(a$genotypes, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c("externalsubjectid", "GT") %in% colnames(gt))) {
  stop("Genotype table must contain externalsubjectid and GT.")
}

normalize_gt <- function(x) {
  x <- gsub("\\|", "/", x)
  out <- rep(NA_character_, length(x))
  out[x == "0/0"] <- "Ref/Ref"
  out[x %in% c("0/1", "1/0")] <- "Het"
  out[x == "1/1"] <- "Hom Alt"
  out
}

gt$group <- normalize_gt(gt$GT)

meta <- merge(
  meta,
  gt[, c("externalsubjectid", "GT", "group")],
  by = "externalsubjectid",
  all.x = TRUE
)

meta <- meta[!is.na(meta$group), , drop = FALSE]

if (!nrow(meta)) {
  stop("No tissue-matched RNA-seq samples had usable 0/0, 0/1, or 1/1 genotypes.")
}

}

# Select one RNA-seq sample per subject before checking PSI/event coverage.
# Highest RIN wins; ties (including missing RIN) use sample ID, then row order.
rin_headers <- tolower(gsub("[^A-Za-z0-9]", "", colnames(meta)))
rin_col <- match("rin", rin_headers)
if (is.na(rin_col)) rin_col <- match("rnaintegritynumber", rin_headers)
sample_rin <- rep(-Inf, nrow(meta))
if (!is.na(rin_col)) {
  sample_rin <- suppressWarnings(as.numeric(as.character(meta[[rin_col]])))
  sample_rin[!is.finite(sample_rin)] <- -Inf
  if (pathology_mode) sample_rin[sample_rin<0 | sample_rin>10] <- -Inf
}
selection_order <- order(
  meta$externalsubjectid, -sample_rin, meta$externalsampleid,
  seq_len(nrow(meta)), na.last = TRUE
)
selected_rows <- selection_order[
  !duplicated(meta$externalsubjectid[selection_order])
]
cat(
  "Subject-level sample selection: kept", length(selected_rows), "of",
  nrow(meta), "metadata rows (one sample per subject).\n"
)
if (pathology_mode) {
  pathology_mapping$mapping_status[pathology_mapping$mapping_status=="eligible"] <- "alternate_rna_sample"
  pathology_mapping$mapping_status[meta$.mapping_row[selected_rows]] <- "selected"
  write.table(pathology_mapping,file.path(a$outdir,paste0(a$prefix,".pathology_mapping.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
}
# Preserve the original ordering of the retained rows.
meta <- meta[sort(selected_rows), , drop = FALSE]

group_levels <- if (pathology_mode) score_levels else c("Ref/Ref", "Het", "Hom Alt")
meta$group <- factor(meta$group, levels = group_levels)

# Build an index of all per-sample PSI files across tissue folders.
cat("Indexing LeafCutter PSI files...\n")
psi_files <- Sys.glob(
  file.path(
    a[["data-root"]],
    "*",
    "RNAseq",
    "Processed",
    "*",
    "leafcutter",
    "psi",
    "*.leafcutter.PSI.tsv"
  )
)

if (!length(psi_files)) {
  stop("No LeafCutter PSI files found beneath ", a[["data-root"]])
}

# Prefer non-sorted PSI file if both exist; both have equivalent fields, but the
# exact *.leafcutter.PSI.tsv path is the canonical source used here.
psi_files <- psi_files[!grepl("\\.leafcutter\\.PSI\\.sorted\\.tsv$", psi_files)]

sample_from_path <- function(p) {
  # .../Processed/SAMPLE/leafcutter/psi/SAMPLE.leafcutter.PSI.tsv
  parts <- strsplit(normalizePath(p, mustWork = FALSE), .Platform$file.sep, fixed = TRUE)[[1]]
  k <- which(parts == "Processed")
  if (!length(k) || k[length(k)] == length(parts)) return(NA_character_)
  parts[k[length(k)] + 1]
}

psi_index <- data.frame(
  sample_dir = vapply(psi_files, sample_from_path, character(1)),
  psi_file = psi_files,
  stringsAsFactors = FALSE
)

psi_index <- psi_index[!is.na(psi_index$sample_dir), , drop = FALSE]
if (pathology_mode) {
  path_tissue <- vapply(strsplit(psi_index$psi_file, .Platform$file.sep, fixed=TRUE), function(parts) {
    k <- which(parts=="Processed"); if (!length(k) || k[1]<3) return(NA_character_)
    parts[k[1]-2]
  }, character(1))
  psi_index <- psi_index[!is.na(path_tissue) & canonical_tissue(path_tissue)==wanted_tissue, , drop=FALSE]
  pathology_selected <- meta
}

meta <- merge(meta, psi_index, by = "sample_dir", all.x = TRUE)

# A selected sample must not multiply into several rows during the PSI join.
if (anyDuplicated(meta$externalsubjectid)) {
  duplicate_subjects <- unique(meta$externalsubjectid[
    duplicated(meta$externalsubjectid) |
      duplicated(meta$externalsubjectid, fromLast = TRUE)
  ])
  stop(
    "Multiple PSI paths matched selected samples for subject(s): ",
    paste(duplicate_subjects, collapse = ", "),
    ". Resolve duplicate PSI paths before plotting; subjects must contribute once."
  )
}

missing_psi <- meta[is.na(meta$psi_file), c(
  "externalsampleid", "externalsubjectid", "tissue", "sample_dir", "GT", "group"
), drop = FALSE]

if (pathology_mode) missing_psi$GT <- NULL
if (nrow(missing_psi)) {
  write.table(
    missing_psi,
    file = file.path(a$outdir, paste0(a$prefix, ".missing_psi.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

if (is.null(a[["psi-matrix"]])) meta <- meta[!is.na(meta$psi_file), , drop = FALSE]

if (!nrow(meta)) {
  stop("None of the genotype-matched metadata samples had a LeafCutter PSI file.")
}

# Parse anchor: chr19:-:17641556-17642845
anchor_parts <- strsplit(a$anchor, ":", fixed = TRUE)[[1]]
if (length(anchor_parts) != 3) {
  stop("Anchor must look like chr19:-:17641556-17642845")
}

anchor_chr <- anchor_parts[1]
anchor_strand <- anchor_parts[2]
anchor_se <- strsplit(anchor_parts[3], "-", fixed = TRUE)[[1]]

if (length(anchor_se) != 2) {
  stop("Could not parse anchor start-end: ", anchor_parts[3])
}

anchor_start <- as.numeric(anchor_se[1])
anchor_end <- as.numeric(anchor_se[2])

if (is.na(anchor_start) || is.na(anchor_end)) {
  stop("Anchor coordinates are not numeric.")
}

read_psi <- function(path) {
  x <- tryCatch(
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) {
      warning("Could not read ", path, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(x)) return(NULL)

  need <- c("chrom", "strand", "start", "end", "reads", "cluster_reads", "PSI")
  if (!all(need %in% colnames(x))) {
    warning("Skipping malformed PSI file: ", path)
    return(NULL)
  }

  x$start <- as.numeric(x$start)
  x$end <- as.numeric(x$end)
  x$reads <- as.numeric(x$reads)
  x$PSI <- suppressWarnings(as.numeric(x$PSI))
  if (any(is.finite(x$PSI) & (x$PSI < 0 | x$PSI > 1))) stop("PSI must be a fraction in [0,1]: ",path)
  x
}

# -------------------------------------------------------------------------
# Event discovery
#
# Stable identity is chr:start:end:strand. We intentionally ignore sample-local
# clu_N identifiers. Starting from the sQTL anchor junction, recover the union
# across subjects of all refined LeafCutter junctions that share either the
# anchor donor/start site or anchor acceptor/end site.
#
# For UNC13A this recovers:
#   17641556-17642414
#   17641556-17642845  (anchor)
#   17642541-17642845
# -------------------------------------------------------------------------
cat("Discovering anchor-connected splice junctions...\n")

event_rows <- list()

for (i in seq_len(nrow(meta))) {
  x <- read_psi(meta$psi_file[i])
  if (is.null(x)) next

  hit <- x[
    x$chrom == anchor_chr &
      x$strand == anchor_strand &
      (x$start == anchor_start | x$end == anchor_end),
    c("chrom", "strand", "start", "end"),
    drop = FALSE
  ]

  if (nrow(hit)) event_rows[[length(event_rows) + 1]] <- hit
}

# Always retain the anchor itself.
event_rows[[length(event_rows) + 1]] <- data.frame(
  chrom = anchor_chr,
  strand = anchor_strand,
  start = anchor_start,
  end = anchor_end,
  stringsAsFactors = FALSE
)

event <- unique(do.call(rbind, event_rows))
event <- event[order(event$start, event$end), , drop = FALSE]
event$junction_id <- paste(event$chrom, event$strand,
                           paste0(event$start, "-", event$end), sep = ":")

if (!nrow(event)) {
  stop("No junctions discovered for anchor ", a$anchor)
}

cat("Discovered ", nrow(event), " anchor-connected junction(s):\n", sep = "")
cat(paste0("  ", event$junction_id, collapse = "\n"), "\n")

# -------------------------------------------------------------------------
# Build sample x junction read-count matrix.
# Missing PSI stays NA; observed zero PSI remains valid. Read counts are audit-only.
# -------------------------------------------------------------------------
psi_mat <- matrix(NA_real_, nrow=nrow(meta), ncol=nrow(event),
                  dimnames=list(meta$externalsampleid,event$junction_id))
psi_source <- rep("sample_file_PSI",nrow(event))
count_mat <- matrix(
  0,
  nrow = nrow(meta),
  ncol = nrow(event),
  dimnames = list(meta$externalsampleid, event$junction_id)
)

for (i in seq_len(nrow(meta))) {
  x <- read_psi(meta$psi_file[i])
  if (is.null(x)) next

  x <- x[x$chrom == anchor_chr & x$strand == anchor_strand, , drop = FALSE]
  if (!nrow(x)) next

  xid <- paste(x$chrom, x$strand, paste0(x$start, "-", x$end), sep = ":")
  m <- match(event$junction_id, xid)
  ok <- !is.na(m)
  count_mat[i, ok] <- x$reads[m[ok]]
  if (anyDuplicated(xid[xid %in% event$junction_id])) stop("Duplicate junctions in ",meta$psi_file[i])
  psi_mat[i, ok] <- x$PSI[m[ok]]
}

event_total <- rowSums(count_mat, na.rm = TRUE)
# Stream the optional matrix: only event rows are retained in memory.
if (!is.null(a[["psi-matrix"]])) {
  con <- file(a[["psi-matrix"]],open="r")
  header <- strsplit(readLines(con,n=1),"\t",fixed=TRUE)[[1]][-1]
  if (anyDuplicated(header)) stop("Duplicate subjects in PSI matrix")
  subject_idx <- match(meta$externalsubjectid,header)
  if (anyNA(subject_idx)) stop("PSI matrix lacks selected subjects; refusing to change the ranking cohort")
  found <- character()
  repeat {
    lines <- readLines(con,n=1000,warn=FALSE)
    if (!length(lines)) break
    for (line in lines) {
      id <- sub("\t.*$","",line)
      j <- match(id,event$junction_id)
      if (is.na(j)) next
      if (id %in% found) stop("Duplicate junction in PSI matrix: ",id)
      values <- strsplit(line,"\t",fixed=TRUE)[[1]][-1]
      if (length(values)!=length(header)) stop("Wrong PSI matrix row length: ",id)
      v <- suppressWarnings(as.numeric(values[subject_idx]))
      if (any(is.finite(v) & (v<0 | v>1))) stop("Matrix PSI must be in [0,1]")
      psi_mat[,j] <- v
      psi_source[j] <- "QTL_matrix_PSI"
      found <- c(found,id)
    }
  }
  close(con)
  if (!a$anchor %in% found) stop("Anchor missing from --psi-matrix")
}
psi_mat[!is.finite(psi_mat)] <- NA_real_
anchor_index <- match(a$anchor,event$junction_id)
usable <- is.finite(psi_mat[,anchor_index])


sample_qc <- data.frame(
  externalsampleid = meta$externalsampleid,
  externalsubjectid = meta$externalsubjectid,
  tissue = meta$tissue,
  sample_dir = meta$sample_dir,
  GT = meta$GT,
  group = as.character(meta$group),
  event_total_reads = event_total,
  anchor_PSI = psi_mat[,anchor_index],
  anchor_PSI_source = psi_source[anchor_index],
  usable_event = usable,
  psi_file = meta$psi_file,
  stringsAsFactors = FALSE
)

if (pathology_mode) {
  sample_qc$GT <- NULL
  sample_qc$Neuronal_TDP43_Score <- as.character(meta$group)
  sample_qc$Pathology_Sample_ID <- meta$Pathology_Sample_ID
  sample_qc$Pathology_Tissue_Name <- meta$Pathology_Tissue_Name
}
write.table(
  sample_qc,
  file = file.path(a$outdir, paste0(a$prefix, ".samples.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (pathology_mode) {
  counts <- do.call(rbind,lapply(group_levels,function(g) data.frame(
    tissue=a$tissue, Neuronal_TDP43_Score=g,
    mapped_unique_subjects=length(unique(pathology_mapping$externalsubjectid[!is.na(pathology_mapping$group) & pathology_mapping$group==g])),
    selected_subjects=sum(as.character(pathology_selected$group)==g),
    subjects_with_psi_file=sum(as.character(meta$group)==g),
    subjects_with_event_coverage=sum(as.character(meta$group)==g & usable))))
  counts$small_group <- counts$subjects_with_event_coverage<5
  write.table(counts,file.path(a$outdir,paste0(a$prefix,".group_counts.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
}
if (!any(usable)) {
  stop("No selected samples had finite original PSI for the anchor.")
}

write.table(data.frame(externalsubjectid=meta$externalsubjectid, group=meta$group, psi_mat,check.names=FALSE),
  file.path(a$outdir,paste0(a$prefix,".subject_PSI.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")

group_mean <- matrix(
  NA_real_,
  nrow = length(group_levels),
  ncol = ncol(psi_mat),
  dimnames = list(group_levels, colnames(psi_mat))
)

junction_n <- matrix(0L,nrow=length(group_levels),ncol=nrow(event),dimnames=list(group_levels,event$junction_id))
group_n <- setNames(integer(length(group_levels)), group_levels)

for (g in group_levels) {
  idx <- usable & as.character(meta$group) == g
  group_n[g] <- sum(idx)

  if (any(idx)) {
    group_mean[g, ] <- colMeans(psi_mat[idx, , drop = FALSE], na.rm = TRUE)
    junction_n[g, ] <- colSums(is.finite(psi_mat[idx, , drop = FALSE]))
    group_mean[g, !is.finite(group_mean[g, ])] <- NA_real_
  }
}

# Diagnostic long table.
usage_long <- do.call(rbind, lapply(group_levels, function(g) {
  data.frame(
    group = g,
    n = as.integer(junction_n[g, ]),
    anchor_group_n = group_n[g],
    PSI_source = psi_source,
    junction_id = event$junction_id,
    chrom = event$chrom,
    strand = event$strand,
    start = event$start,
    end = event$end,
    mean_PSI = as.numeric(group_mean[g, ]),
    stringsAsFactors = FALSE
  )
}))

# -------------------------------------------------------------------------
# GTF annotation and local exon model.
# -------------------------------------------------------------------------
cat("Loading local exon annotation...\n")

region_min <- max(1, min(event$start) - 2000)
region_max <- max(event$end) + 2000

gtf_cmd <- sprintf(
  "awk -F '\\t' 'BEGIN{OFS=\"\\t\"} $1==\"%s\" && $3==\"exon\" && $4<=%d && $5>=%d {print}' %s",
  anchor_chr,
  as.integer(region_max),
  as.integer(region_min),
  shQuote(a$gtf)
)

gtf <- tryCatch(
  read.delim(
    pipe(gtf_cmd),
    header = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  ),
  error = function(e) data.frame()
)

if (nrow(gtf)) {
  colnames(gtf)[1:9] <- c(
    "chr", "source", "feature", "start", "end", "score",
    "strand", "frame", "attributes"
  )

  extract_attr <- function(x, key) {
    pat <- paste0(key, ' "([^"]+)"')
    m <- regexec(pat, x, perl = TRUE)
    z <- regmatches(x, m)
    vapply(z, function(y) if (length(y) >= 2) y[2] else NA_character_, character(1))
  }

  gtf$gene_name <- extract_attr(gtf$attributes, "gene_name")
  gtf$gene_id <- extract_attr(gtf$attributes, "gene_id")
  gtf$transcript_id <- extract_attr(gtf$attributes, "transcript_id")

  gene_exons <- gtf[gtf$gene_name == a$gene, , drop = FALSE]
  if (!nrow(gene_exons)) {
    warning("No local GTF exons found for gene ", a$gene,
            "; using all local exons in the plotting region.")
    gene_exons <- gtf
  }
} else {
  gene_exons <- data.frame()
}

# Infer exact annotated introns from consecutive exons within transcripts.
annotated_keys <- character(0)

if (nrow(gene_exons) && any(!is.na(gene_exons$transcript_id))) {
  txs <- split(gene_exons, gene_exons$transcript_id)

  for (tx in txs) {
    tx <- tx[order(tx$start, tx$end), , drop = FALSE]
    if (nrow(tx) < 2) next

    left <- tx$end[-nrow(tx)]
    right <- tx$start[-1]
    annotated_keys <- c(
      annotated_keys,
      paste(left, right, sep = ":")
    )
  }

  annotated_keys <- unique(annotated_keys)
}

event$key <- paste(event$start, event$end, sep = ":")
event$verdict <- ifelse(event$key %in% annotated_keys, "annotated", "cryptic")
event$letter <- letters[seq_len(nrow(event))]

usage_long$verdict <- event$verdict[match(usage_long$junction_id, event$junction_id)]
usage_long$letter <- event$letter[match(usage_long$junction_id, event$junction_id)]

write.table(
  usage_long,
  file = file.path(a$outdir, paste0(a$prefix, ".junctions.tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -------------------------------------------------------------------------
# LeafViz-style coordinate compression.
# -------------------------------------------------------------------------
s <- sort(unique(c(event$start, event$end)))
if (length(s) < 2) stop("Not enough unique splice sites to plot.")

length_transform <- function(g) log(g + 1)
d <- diff(s)
coords <- c(0, cumsum(length_transform(d)))
names(coords) <- as.character(s)

invert_mapping <- function(pos) {
  if (pos <= min(s)) return(coords[1])
  if (pos >= max(s)) return(coords[length(coords)])
  if (as.character(pos) %in% names(coords)) return(coords[as.character(pos)])

  w <- which(pos > s[-length(s)] & pos < s[-1])
  if (length(w) != 1) return(NA_real_)

  coords[w] +
    (coords[w + 1] - coords[w]) *
    (pos - s[w]) / (s[w + 1] - s[w])
}

total_length <- max(coords)
xlim_vals <- c(-0.05 * total_length, 1.05 * total_length)

# Exons touching event splice sites, mimicking LeafViz's culling logic.
exon_df <- data.frame()

if (nrow(gene_exons)) {
  ex <- unique(gene_exons[, c("start", "end", "gene_name", "strand"), drop = FALSE])

  ex <- ex[
    (ex$end %in% event$start | ex$start %in% event$end) &
      (
        (ex$end - ex$start) <= 500 |
          ex$end == min(event$start) |
          ex$start == max(event$end)
      ),
    ,
    drop = FALSE
  ]

  if (nrow(ex)) {
    exon_df <- data.frame(
      x = vapply(ex$start, invert_mapping, numeric(1)),
      xend = vapply(ex$end, invert_mapping, numeric(1)),
      gene_name = ex$gene_name,
      strand = ex$strand,
      stringsAsFactors = FALSE
    )

    # Ensure tiny exons remain visible.
    min_exon_length <- 0.5
    too_short <- (exon_df$xend - exon_df$x) < min_exon_length
    exon_df$xend[too_short] <- exon_df$x[too_short] + min_exon_length
    exon_df <- exon_df[!duplicated(paste(exon_df$x, exon_df$xend)), , drop = FALSE]
  }
}

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------
make_panel <- function(g, show_title = FALSE, show_legend = FALSE) {
  if (group_n[g]==0) {
    p <- ggplot() + annotate("text",x=0,y=0,label=paste(g,"— no data"),size=5) +
      theme_void() + ylab(sprintf("%s (n=0)",g))
    if (show_title) p <- p + ggtitle(a$gene,subtitle=paste(a$anchor,"Neuronal TDP43 score",sep=" | "))
    return(p)
  }
  vals <- group_mean[g, ]
  # Missing means stay NA and are explicitly labelled rather than shown as zero.

  ed <- event
  ed$mean_PSI <- as.numeric(vals)
  ed$x <- coords[as.character(ed$start)]
  ed$xend <- coords[as.character(ed$end)]
  ed$xmid <- (ed$x + ed$xend) / 2
  ed$span <- pmax(ed$xend - ed$x, 0.01)

  # Alternate arcs above/below exactly as the original LeafViz plot does.
  ed$side <- ifelse(seq_len(nrow(ed)) %% 2 == 1, 1, -1)
  ed$ymid <- ed$side * (ed$span^0.65 / 2 + 0.25)

  # PSI view: linear width on a common 0–100% scale across all groups.
  ed$line_size <- ifelse(is.na(ed$mean_PSI),0.25,0.25 + 9.75 * ed$mean_PSI)
  ed$label <- ifelse(is.na(ed$mean_PSI),paste0("NA  ",ed$letter," (n=0)"),
                      sprintf("%.2f%%  %s (n=%d)",100*ed$mean_PSI,ed$letter,junction_n[g,]))

  ymax <- max(abs(ed$ymid), na.rm = TRUE) * 1.35
  if (!is.finite(ymax) || ymax <= 0) ymax <- 1

  p <- ggplot() +
    theme_bw(base_size = 13) +
    theme(
      panel.background = element_rect(fill = "white", colour = "white"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = if (pathology_mode) 10 else 12),
      panel.border = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      plot.margin = margin(5, 10, 5, 10)
    ) +
    geom_hline(yintercept = 0, colour = "white", size = 3) +
    geom_hline(yintercept = 0, colour = "black", size = 0.6)

  # Draw each arc as two curves meeting at the midpoint, following LeafViz.
  for (i in seq_len(nrow(ed))) {
    curvature <- if (ed$side[i] > 0) -0.1 else 0.1

    seg1 <- data.frame(
      x = ed$x[i], xend = ed$xmid[i],
      y = 0, yend = ed$ymid[i],
      verdict = ed$verdict[i]
    )
    seg2 <- data.frame(
      x = ed$xmid[i], xend = ed$xend[i],
      y = ed$ymid[i], yend = 0,
      verdict = ed$verdict[i]
    )

    p <- p +
      geom_curve(
        data = seg1,
        aes(x = x, xend = xend, y = y, yend = yend, colour = verdict),
        curvature = curvature,
        angle = 90,
        lineend = "round",
        size = ed$line_size[i],
        linetype = if (is.na(ed$mean_PSI[i])) "dotted" else "solid",
        inherit.aes = FALSE
      ) +
      geom_curve(
        data = seg2,
        aes(x = x, xend = xend, y = y, yend = yend, colour = verdict),
        curvature = curvature,
        angle = 90,
        lineend = "round",
        size = ed$line_size[i],
        linetype = if (is.na(ed$mean_PSI[i])) "dotted" else "solid",
        inherit.aes = FALSE
      ) +
      geom_label(
        data = ed[i, , drop = FALSE],
        aes(x = xmid, y = 0.92 * ymid, label = label),
        size = 3.2,
        label.size = NA,
        fill = "white",
        colour = "black",
        inherit.aes = FALSE
      )
  }

  if (nrow(exon_df)) {
    p <- p +
      geom_segment(
        data = exon_df,
        aes(x = x, xend = xend, y = 0, yend = 0),
        size = 6,
        colour = "black",
        inherit.aes = FALSE
      )

    gene_strand <- unique(exon_df$strand)
    gene_strand <- gene_strand[!is.na(gene_strand)]
    if (length(gene_strand) == 1) {
      # Add simple strand arrows in long gaps between visible exons.
      exs <- exon_df[order(exon_df$x), , drop = FALSE]
      if (nrow(exs) >= 2) {
        for (j in seq_len(nrow(exs) - 1)) {
          gap_start <- exs$xend[j]
          gap_end <- exs$x[j + 1]
          if ((gap_end - gap_start) > 0.025 * total_length) {
            mid <- (gap_start + gap_end) / 2
            if (gene_strand == "+") {
              xa <- gap_start
              xb <- mid
            } else {
              xa <- mid
              xb <- gap_end
            }

            p <- p + geom_segment(
              aes(x = xa, xend = xb, y = 0, yend = 0),
              arrow = arrow(
                ends = if (gene_strand == "+") "last" else "first",
                type = "open",
                angle = 30,
                length = unit(0.08, "inches")
              ),
              size = 0.7,
              colour = "black",
              inherit.aes = FALSE
            )
          }
        }
      }
    }
  }

  # LeafViz colors: annotated red, unannotated/cryptic blue.
  p <- p +
    scale_colour_manual(
      values = c(annotated = "red", cryptic = "dodgerblue"),
      breaks = c("annotated", "cryptic"),
      labels = c("annotated", "cryptic"),
      name = NULL,
      drop = FALSE
    ) +
    coord_cartesian(xlim = xlim_vals, ylim = c(-ymax, ymax), clip = "off") +
    ylab(sprintf("%s (n=%d)%s", g, group_n[g],
                 if (pathology_mode && group_n[g]<5) "\nsmall group" else ""))

  if (show_title) {
    if (pathology_mode) {
      variant_label <- "Neuronal TDP43 score"
    } else {
      variant_label <- paste0(
        a[["variant-chr"]], ":", a[["variant-pos"]],
        " ", a[["variant-ref"]], ">", a[["variant-alt"]]
      )
      if (!is.na(a[["variant-id"]]) && nzchar(a[["variant-id"]]) && a[["variant-id"]] != ".") {
        variant_label <- paste0(a[["variant-id"]], " | ", variant_label)
      }
    }
    p <- p + ggtitle(a$gene, subtitle=paste0(a$anchor, " | ", variant_label," | Original junction PSI (%)"))
  }

  p
}

plots <- list()
for (i in seq_along(group_levels)) {
  g <- group_levels[i]
  plots[[i]] <- make_panel(
    g,
    show_title = (i == 1),
    show_legend = (i == length(group_levels))
  )
}

# Junction legend/table grob.
table_df <- event[, c("letter", "chrom", "strand", "start", "end", "verdict"), drop = FALSE]
colnames(table_df) <- c("ID", "chr", "strand", "start", "end", "annotation")

table_grob <- gridExtra::tableGrob(
  table_df,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    base_size = 9,
    core = list(bg_params = list(fill = c("white", "grey95"), col = NA)),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

combined <- gridExtra::arrangeGrob(
  grobs = c(plots, list(table_grob)),
  ncol = 1,
  heights = c(rep(1, length(plots)), 0.8)
)

pdf_file <- file.path(a$outdir, paste0(a$prefix, ".pdf"))
png_file <- file.path(a$outdir, paste0(a$prefix, ".png"))
svg_file <- file.path(a$outdir, paste0(a$prefix, ".svg"))

ggsave(pdf_file, combined, width = 11, height = if (pathology_mode) 12 else 10, units = "in")
ggsave(png_file, combined, width = 11, height = if (pathology_mode) 12 else 10, units = "in", dpi = 300)
ggsave(svg_file, combined, width = 11, height = if (pathology_mode) 12 else 10, units = "in")

cat(if (pathology_mode) "\nNeuronal TDP43 sample counts with finite anchor PSI:\n" else "\nGenotype sample counts with finite anchor PSI:\n")
print(group_n)

cat("\nOutputs:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")
cat("  ", svg_file, "\n", sep = "")
cat("  ", file.path(a$outdir, paste0(a$prefix, ".samples.tsv")), "\n", sep = "")
cat("  ", file.path(a$outdir, paste0(a$prefix, ".junctions.tsv")), "\n", sep = "")
