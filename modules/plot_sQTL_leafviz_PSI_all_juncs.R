#!/usr/bin/env Rscript

# Static LeafViz-style plot of the full discovered LeafCutter cluster.
# Existing calls default to genotype mode unchanged.
# Opt in: --group-by tdp43 (no genotype or variant arguments required).
# Optional --tdp43 <aligned CSV> and --pathology <original pathology CSV>;
# defaults are aligned_rnaseq_tdp43_by_subject.csv and
# collections.postmortem_tissue_core.semiquantitative_tdp43_data.csv under --data-root.
# TDP43 mode matches subject + RNA sample + tissue, verifies pathology provenance,
# selects one highest-valid-RIN sample per subject before inspecting event coverage, excludes missing/conflicting/unsupported scores, and writes mapping/count audits.
# Arcs show mean per-subject PSI using ONE full, authoritative cluster denominator.
# Default: use each sample's leafcutter/clusters/leafcutter_clusters_refined.
# Display the union of full anchor-containing clusters, preserving each sample denominator.
# Outside-cluster junctions contribute no share of that sample's anchor cluster.
# Missing files/anchor exclude the whole sample; they are not zero measurements.
# Optional --cluster-counts: a common unfiltered LeafCutter count matrix instead.
# --psi-matrix is rejected: filtered, sample-local PSI cannot define a full cluster.
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

if (!is.null(a[["psi-matrix"]])) stop("--psi-matrix cannot supply full-cluster PSI. Use refined cluster files (default) or --cluster-counts instead.")

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

# Cluster counts define both membership and measurements; sample-local PSI is not used.
if (pathology_mode) pathology_selected <- meta
meta$psi_file <- NA_character_

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

# Read one authoritative cluster, including members absent from filtered QTL matrices.
# Standard LeafCutter IDs: chr:start:end:clu_N[_strand]. Counts or numerator/denominator.
read_cluster_counts <- function(path, meta, anchor_chr, anchor_strand, anchor_start, anchor_end) {
  open_counts <- function() if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  fields <- function(line) strsplit(trimws(line), "[[:space:]]+")[[1]]
  parse_id <- function(id) {
    z <- strsplit(id, ":", fixed = TRUE)[[1]]
    if (length(z) != 4L) stop("Expected LeafCutter chr:start:end:cluster ID; got ", id)
    se <- suppressWarnings(as.numeric(z[2:3]))
    if (any(!is.finite(se)) || se[1] >= se[2]) stop("Invalid count-file coordinates: ", id)
    strand <- if (grepl("_[+-]$", z[4])) sub(".*_([+-])$", "\\1", z[4]) else NA_character_
    list(chrom=z[1], start=se[1], end=se[2], cluster=z[4], strand=strand)
  }
  con <- open_counts()
  header <- fields(readLines(con, n=1L, warn=FALSE))
  first <- readLines(con, n=1L, warn=FALSE)
  close(con)
  if (!length(first)) stop("Empty cluster count file: ", path)
  ncols <- length(fields(first)) - 1L
  if (length(header) == ncols + 1L) header <- header[-1L]
  if (length(header) != ncols || anyDuplicated(header)) stop("Invalid/duplicate sample header in ", path)
  normalize_sample <- function(x) {
    x <- basename(x)
    x <- sub("\\.gz$", "", x)
    x <- sub("\\.junc$", "", x)
    x <- sub("\\.bam$", "", x)
    gsub("-", "_", x, fixed=TRUE)
  }
  normalized <- normalize_sample(header)
  columns <- vapply(seq_len(nrow(meta)), function(i) {
    aliases <- normalize_sample(c(meta$externalsampleid[i], meta$externalsubjectid[i]))
    hits <- which(normalized %in% aliases)
    if (length(hits) != 1L) stop("Expected exactly one count column for sample ", meta$externalsampleid[i],
                               "; found ", length(hits), ". Use matching sample or subject IDs in count headers.")
    hits
  }, integer(1))
  if (anyDuplicated(columns)) stop("Multiple selected samples map to the same count column")
  scan_rows <- function(callback) {
    con <- open_counts(); on.exit(close(con))
    readLines(con, n=1L, warn=FALSE)
    repeat {
      lines <- readLines(con, n=2000L, warn=FALSE)
      if (!length(lines)) break
      for (line in lines) {
        if (!nzchar(trimws(line))) next
        id <- sub("[[:space:]].*$", "", line)
        if (!startsWith(id, paste0(anchor_chr, ":"))) next
        callback(parse_id(id), line)
      }
    }
  }
  clusters <- character()
  scan_rows(function(z, line) {
    if (z$start == anchor_start && z$end == anchor_end &&
        (is.na(z$strand) || z$strand == anchor_strand)) clusters <<- c(clusters, z$cluster)
  })
  clusters <- unique(clusters)
  if (length(clusters) != 1L) stop("Anchor must identify exactly one cluster in --cluster-counts; found ", length(clusters))
  event_rows <- list(); count_rows <- list(); denominators <- list()
  scan_rows(function(z, line) {
    if (z$cluster != clusters[1]) return(invisible(NULL))
    values <- fields(line)[-1L]
    if (length(values) != length(header)) stop("Wrong number of count columns")
    values <- values[columns]
    ratio <- grepl("/", values, fixed=TRUE)
    numer <- suppressWarnings(as.numeric(sub("/.*$", "", values)))
    missing <- values %in% c("NA", "NaN", "nan", ".")
    if (any(!missing & (!is.finite(numer) | numer < 0 | numer != floor(numer)))) stop("Invalid junction counts")
    denom <- rep(NA_real_, length(values))
    if (any(ratio)) {
      if (any(!grepl("^[0-9]+/[0-9]+$", values[ratio]))) stop("Invalid numerator/denominator count")
      denom[ratio] <- as.numeric(sub("^.*/", "", values[ratio]))
      if (any(numer[ratio] > denom[ratio])) stop("Junction count exceeds cluster denominator")
    }
    event_rows[[length(event_rows)+1L]] <<- data.frame(chrom=z$chrom, strand=anchor_strand,
        start=z$start, end=z$end, cluster_id=z$cluster)
    count_rows[[length(count_rows)+1L]] <<- numer
    denominators[[length(denominators)+1L]] <<- denom
  })
  event <- do.call(rbind, event_rows)
  event$junction_id <- paste(event$chrom, event$strand, paste0(event$start, "-", event$end), sep=":")
  if (anyDuplicated(event$junction_id)) stop("Duplicate junction in authoritative cluster")
  ord <- order(event$start, event$end)
  event <- event[ord, , drop=FALSE]
  counts <- t(do.call(rbind, count_rows))[, ord, drop=FALSE]
  denom <- t(do.call(rbind, denominators))[, ord, drop=FALSE]
  total <- rowSums(counts)
  complete <- is.finite(total)
  for (i in which(complete)) {
    d <- denom[i, is.finite(denom[i, ])]
    if (length(d) && any(d != total[i])) stop("Count denominators disagree with full cluster sum for ", meta$externalsampleid[i],
                                            "; supply an unfiltered cluster count file")
  }
  dimnames(counts) <- list(meta$externalsampleid, event$junction_id)
  list(event=event, counts=counts, total=total, cluster=clusters[1], sample_columns=header[columns])
}

# A refined-file line is the authoritative sample-local cluster membership.
# Never infer membership from shared splice sites or equal cluster_reads values.
read_refined_cluster <- function(path, chrom, strand, start, end) {
  if (!file.exists(path)) return(list(reason="missing_refined_file"))
  con <- file(path, "rt"); on.exit(close(con))
  locus <- paste(chrom, strand, sep=":")
  anchor <- paste(start, end, sep=":")
  hit <- NULL; line_number <- 0L
  repeat {
    lines <- readLines(con, n=2000L, warn=FALSE)
    if (!length(lines)) break
    for (line in lines) {
      line_number <- line_number + 1L
      fields <- strsplit(trimws(line), "[[:space:]]+")[[1]]
      if (fields[1] != locus) next
      tokens <- fields[-1L]
      ids <- sub(":[^:]+$", "", tokens)
      if (!anchor %in% ids) next
      if (!is.null(hit)) stop("Anchor occurs in multiple refined clusters: ", path)
      z <- strsplit(tokens, ":", fixed=TRUE)
      if (any(lengths(z) != 3L)) stop("Malformed refined cluster: ", path)
      vals <- matrix(suppressWarnings(as.numeric(unlist(z))), ncol=3L, byrow=TRUE)
      if (any(!is.finite(vals)) || any(vals < 0) || any(vals != floor(vals)) ||
          any(vals[,1] >= vals[,2]) || anyDuplicated(ids)) stop("Invalid refined cluster counts/coordinates: ", path)
      hit <- list(reason="included", line=line_number,
                  event=data.frame(chrom=chrom, strand=strand, start=vals[,1], end=vals[,2]),
                  counts=vals[,3])
    }
  }
  if (is.null(hit)) list(reason="anchor_absent_from_refined_clusters") else hit
}

read_sample_refined_clusters <- function(meta, root, tissue, chrom, strand, start, end) {
  paths <- file.path(root, tissue, "RNAseq", "Processed", meta$sample_dir,
                     "leafcutter", "clusters", "leafcutter_clusters_refined")
  samples <- lapply(paths, read_refined_cluster, chrom=chrom, strand=strand, start=start, end=end)
  valid <- vapply(samples, function(x) x$reason == "included", logical(1))
  if (!any(valid)) stop("No selected sample has the anchor in leafcutter_clusters_refined. Checked ", paths[1])
  key <- function(x) paste(x$chrom, x$strand, paste0(x$start, "-", x$end), sep=":")
  event <- unique(do.call(rbind, lapply(samples[valid], `[[`, "event")))
  event <- event[order(event$start, event$end), , drop=FALSE]
  event$junction_id <- key(event)
  event$cluster_id <- "union_of_sample_anchor_clusters"
  counts <- matrix(NA_real_, nrow=nrow(meta), ncol=nrow(event),
                   dimnames=list(meta$externalsampleid, event$junction_id))
  membership <- matrix(FALSE, nrow=nrow(meta), ncol=nrow(event), dimnames=dimnames(counts))
  audit <- list()
  for (i in which(valid)) {
    # Zero outside this sample's anchor cluster means no share of THIS cluster;
    # it is not an imputed zero measurement or a claim of no reads at that locus.
    counts[i, ] <- 0
    j <- match(key(samples[[i]]$event), event$junction_id)
    counts[i, j] <- samples[[i]]$counts
    membership[i, j] <- TRUE
    audit[[length(audit)+1L]] <- data.frame(externalsampleid=meta$externalsampleid[i],
        junction_id=event$junction_id, member_of_sample_cluster=membership[i, ],
        cluster_assigned_reads=counts[i, ], refined_file=paths[i], cluster_line=samples[[i]]$line)
  }
  list(event=event, counts=counts, total=rowSums(counts),
       cluster="union_of_sample_anchor_clusters", sample_columns=meta$externalsampleid,
       paths=paths, reasons=vapply(samples, `[[`, character(1), "reason"),
       membership_audit=do.call(rbind, audit), memberships=membership)
}

if (is.null(a[["cluster-counts"]])) {
  cluster_data <- read_sample_refined_clusters(meta, a[["data-root"]], a$tissue,
                                              anchor_chr, anchor_strand, anchor_start, anchor_end)
  meta$psi_file <- cluster_data$paths
  write.table(cluster_data$membership_audit,
              file.path(a$outdir, paste0(a$prefix, ".sample_cluster_membership.tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
} else {
  cluster_data <- read_cluster_counts(a[["cluster-counts"]], meta, anchor_chr, anchor_strand, anchor_start, anchor_end)
  meta$psi_file <- a[["cluster-counts"]]
  cluster_data$reasons <- rep("included", nrow(meta))
}

event <- cluster_data$event
count_mat <- cluster_data$counts
event_total <- cluster_data$total
# A missing measurement excludes that subject from ALL arcs; measured zeros stay zero.
usable <- is.finite(event_total) & event_total > 0
psi_mat <- count_mat
psi_mat[,] <- NA_real_
psi_mat[usable, ] <- count_mat[usable, , drop=FALSE] / event_total[usable]
psi_source <- rep("full_cluster_counts", nrow(event))
anchor_index <- match(a$anchor, event$junction_id)
write.table(event, file.path(a$outdir, paste0(a$prefix, ".cluster_junctions.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("Full cluster ", cluster_data$cluster, ": ", nrow(event), " junctions spanning ",
    min(event$start), "-", max(event$end), "\n", sep="")

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
  count_file = meta$psi_file,
  count_sample_column = cluster_data$sample_columns,
  cluster_id = cluster_data$cluster,
  exclusion_reason = ifelse(usable, "included", ifelse(cluster_data$reasons != "included", cluster_data$reasons,
      ifelse(is.finite(event_total), "zero_cluster_coverage", "incomplete_cluster_counts"))),
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
    subjects_with_count_column=sum(as.character(meta$group)==g),
    subjects_with_event_coverage=sum(as.character(meta$group)==g & usable))))
  counts$small_group <- counts$subjects_with_event_coverage<5
  write.table(counts,file.path(a$outdir,paste0(a$prefix,".group_counts.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
}
if (!any(usable)) {
  stop("No selected samples had complete counts and positive coverage for the full cluster.")
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

# Enforce the same cohort and a 100% sum in every non-empty group.
for (g in group_levels[group_n > 0]) {
  if (abs(sum(group_mean[g, ]) - 1) > 1e-8) stop("Cluster PSI does not sum to 100% for ", g)
  if (any(junction_n[g, ] != group_n[g])) stop("Inconsistent sample cohort across cluster arcs")
}
write.table(data.frame(group=group_levels, n=as.integer(group_n),
                      sum_PSI_percent=100*rowSums(group_mean)),
            file.path(a$outdir, paste0(a$prefix, ".cluster_totals.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)

# Diagnostic long table.
usage_long <- do.call(rbind, lapply(group_levels, function(g) {
  data.frame(
    group = g,
    n = as.integer(junction_n[g, ]),
    cluster_group_n = group_n[g],
    cluster_id = cluster_data$cluster,
    subjects_with_junction_in_cluster = if (is.null(cluster_data$memberships)) as.integer(junction_n[g, ]) else
      colSums(cluster_data$memberships[usable & as.character(meta$group)==g, , drop=FALSE]),
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
# Stable IDs remain defined beyond the 26th junction.
event$letter <- vapply(seq_len(nrow(event)), function(i) {
  label <- ""
  while (i > 0L) {
    i <- i - 1L
    label <- paste0(letters[i %% 26L + 1L], label)
    i <- i %/% 26L
  }
  label
}, character(1))

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
  # Give nearby labels separate vertical lanes on each side of the baseline.
  for (side in c(-1, 1)) {
    ii <- which(ed$side == side)
    ii <- ii[order(abs(ed$ymid[ii]))]
    placed <- integer()
    for (j in ii) {
      height <- abs(ed$ymid[j])
      nearby <- placed[abs(ed$xmid[placed] - ed$xmid[j]) < 0.20 * total_length]
      if (length(nearby)) height <- max(height, max(abs(ed$ymid[nearby])) + 1.2)
      ed$ymid[j] <- side * height
      placed <- c(placed, j)
    }
  }

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
      axis.text.y = element_blank(),
      axis.text.x = element_text(size=8, angle=45, hjust=1),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size=9),
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
    scale_x_continuous(breaks=coords[unique(round(seq(1, length(s), length.out=min(8, length(s)))))],
                       labels=format(s[unique(round(seq(1, length(s), length.out=min(8, length(s)))))], scientific=FALSE, trim=TRUE),
                       name=paste0(anchor_chr, " genomic position (compressed spacing)")) +
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
    p <- p + ggtitle(a$gene, subtitle=paste0(a$anchor, " | ", variant_label," | Full cluster PSI (%)",
      if (is.null(a[["cluster-counts"]])) " | Union of sample anchor clusters" else ""))
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

panel_height <- max(2.5, 1.5 + 0.16 * nrow(event))
table_height <- max(1.5, 0.24 * (nrow(event) + 1))
figure_height <- length(plots) * panel_height + table_height

combined <- gridExtra::arrangeGrob(
  grobs = c(plots, list(table_grob)),
  ncol = 1,
  heights = c(rep(panel_height, length(plots)), table_height)
)

pdf_file <- file.path(a$outdir, paste0(a$prefix, ".pdf"))
png_file <- file.path(a$outdir, paste0(a$prefix, ".png"))
svg_file <- file.path(a$outdir, paste0(a$prefix, ".svg"))

ggsave(pdf_file, combined, width = 11, height = figure_height, units = "in", limitsize = FALSE)
ggsave(png_file, combined, width = 11, height = figure_height, units = "in", dpi = min(300, 30000 / figure_height), limitsize = FALSE)
ggsave(svg_file, combined, width = 11, height = figure_height, units = "in", limitsize = FALSE)

cat(if (pathology_mode) "\nNeuronal TDP43 sample counts with complete cluster coverage:\n" else "\nGenotype sample counts with complete cluster coverage:\n")
print(group_n)

cat("\nOutputs:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")
cat("  ", svg_file, "\n", sep = "")
cat("  ", file.path(a$outdir, paste0(a$prefix, ".samples.tsv")), "\n", sep = "")
cat("  ", file.path(a$outdir, paste0(a$prefix, ".junctions.tsv")), "\n", sep = "")
