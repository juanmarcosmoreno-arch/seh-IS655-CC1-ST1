#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else "scripts/25_run_mlst.R"
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = "") {
  hit <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", name, "="), "", hit[1])
}

results_dir <- file.path(project_dir, "results")
out_dir <- file.path(results_dir, "mlst")
medaka_dir <- file.path(results_dir, "medaka")
master_table <- file.path(results_dir, "master_table", "tabla_maestra_staphylococcus.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "")
}

write_tsv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

clean <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  trimws(x)
}

col_or_empty <- function(df, col) {
  if (nrow(df) == 0 || !(col %in% names(df))) return(rep("", nrow(df)))
  clean(df[[col]])
}

parse_sample_label <- function(sample) {
  hit <- regmatches(sample, regexec("^cepa_([^_]+)_(barcode[0-9]+)$", sample))[[1]]
  if (length(hit) == 3) {
    return(list(strain = hit[2], barcode = hit[3]))
  }
  list(strain = sample, barcode = "")
}

join <- function(x) {
  x <- unique(clean(x))
  x <- x[x != ""]
  paste(x, collapse = ";")
}

first_existing_executable <- function(candidates) {
  for (candidate in candidates) {
    candidate <- clean(candidate)
    if (candidate == "") next
    resolved <- Sys.which(candidate)
    if (resolved != "") return(unname(resolved))
    if (file.exists(candidate) && file.access(candidate, mode = 1) == 0) return(candidate)
  }
  ""
}

mlst_bin <- first_existing_executable(c(
  get_arg("--mlst-bin", ""),
  Sys.getenv("MLST_BIN", ""),
  file.path(project_dir, ".micromamba", "envs", "gridion_mlst", "bin", "mlst"),
  "mlst"
))

if (mlst_bin == "") {
  stop("No encontre el ejecutable mlst. Instalar o pasar --mlst-bin=/ruta/a/mlst")
}

mlst_bin_dir <- dirname(mlst_bin)
Sys.setenv(PATH = paste(mlst_bin_dir, Sys.getenv("PATH"), sep = .Platform$path.sep))

threads <- clean(get_arg("--threads", Sys.getenv("MLST_THREADS", "4")))
threads_num <- suppressWarnings(as.integer(threads))
use_threads <- !is.na(threads_num) && threads_num > 0

consensus_paths <- sort(Sys.glob(file.path(medaka_dir, "*", "consensus.fasta")))
if (length(consensus_paths) == 0) {
  stop("No encontre consensos en ", medaka_dir)
}

raw_path <- file.path(out_dir, "mlst_raw.tsv")
stderr_path <- file.path(out_dir, "mlst.stderr.log")

run_mlst <- function(extra_args) {
  if (file.exists(raw_path)) unlink(raw_path)
  if (file.exists(stderr_path)) unlink(stderr_path)
  status <- system2(
    mlst_bin,
    args = c(extra_args, consensus_paths),
    stdout = raw_path,
    stderr = stderr_path
  )
  if (is.null(status)) status <- 0L
  status
}

message("MLST bin: ", mlst_bin)
message("Consensos: ", length(consensus_paths))

mlst_args <- character(0)
if (use_threads) mlst_args <- c("--threads", as.character(threads_num))
status <- run_mlst(mlst_args)

if (!identical(as.integer(status), 0L) && use_threads) {
  err <- if (file.exists(stderr_path)) paste(readLines(stderr_path, warn = FALSE), collapse = "\n") else ""
  if (grepl("thread|Unknown option|invalid option", err, ignore.case = TRUE)) {
    message("Reintentando mlst sin --threads por compatibilidad.")
    status <- run_mlst(character(0))
  }
}

if (!identical(as.integer(status), 0L)) {
  err <- if (file.exists(stderr_path)) paste(readLines(stderr_path, warn = FALSE), collapse = "\n") else ""
  stop("mlst termino con codigo ", status, "\n", err)
}

if (!file.exists(raw_path) || file.info(raw_path)$size == 0) {
  stop("mlst no produjo salida en ", raw_path)
}

raw_lines <- readLines(raw_path, warn = FALSE)
raw_lines <- raw_lines[clean(raw_lines) != ""]
parts <- strsplit(raw_lines, "\t", fixed = TRUE)
max_cols <- max(lengths(parts))
raw_matrix <- do.call(rbind, lapply(parts, function(x) c(x, rep("", max_cols - length(x)))))
raw_df <- as.data.frame(raw_matrix, stringsAsFactors = FALSE, check.names = FALSE)
names(raw_df) <- c("fasta", "scheme", "st", paste0("allele_", seq_len(max(0, max_cols - 3))))

sample_from_fasta <- function(path) {
  clean_path <- clean(path)
  parent <- basename(dirname(clean_path))
  if (parent == "." || parent == "" || parent == basename(clean_path)) {
    parent <- sub("^.*(cepa_[^/]+_barcode[0-9]+).*$", "\\1", clean_path)
  }
  parent
}

master <- read_tsv(master_table)
if (nrow(master) > 0) {
  master_sample <- col_or_empty(master, "wgs_sample")
  fallback_sample <- ifelse(
    col_or_empty(master, "cepa_id") != "" & col_or_empty(master, "barcode_id") != "",
    paste0("cepa_", col_or_empty(master, "cepa_id"), "_", col_or_empty(master, "barcode_id")),
    ""
  )
  master$sample_key <- ifelse(master_sample != "", master_sample, fallback_sample)
} else {
  master$sample_key <- character(0)
}

master_lookup <- function(sample, col) {
  if (nrow(master) == 0 || !(col %in% names(master))) return("")
  hit <- which(master$sample_key == sample)
  if (length(hit) == 0) return("")
  clean(master[[col]][hit[1]])
}

allele_cols <- setdiff(names(raw_df), c("fasta", "scheme", "st"))

is_missing_allele <- function(x) {
  x <- clean(x)
  x == "" | x == "-" | x == "?" | grepl("\\((-|\\?|0)?\\)$", x)
}

valid_scheme <- function(x) {
  x <- tolower(clean(x))
  !(x %in% c("", "-", "na", "nan", "unknown", "none"))
}

valid_st <- function(x) {
  x <- tolower(clean(x))
  !(x %in% c("", "-", "0", "na", "nan", "unknown", "none", "?"))
}

summary_rows <- lapply(seq_len(nrow(raw_df)), function(i) {
  fasta <- clean(raw_df$fasta[i])
  sample <- sample_from_fasta(fasta)
  label <- parse_sample_label(sample)
  alleles <- if (length(allele_cols) > 0) clean(unlist(raw_df[i, allele_cols], use.names = FALSE)) else character(0)
  alleles <- alleles[alleles != ""]
  loci_detected <- sum(!is_missing_allele(alleles))
  loci_total <- length(alleles)
  scheme <- clean(raw_df$scheme[i])
  st <- clean(raw_df$st[i])
  quality <- if (!valid_scheme(scheme)) {
    "sin_esquema"
  } else if (valid_st(st)) {
    "ok"
  } else {
    "perfil_sin_ST"
  }

  species <- master_lookup(sample, "wgs_species_call")
  if (species == "") species <- master_lookup(sample, "ani_best_species")

  data.frame(
    sample = sample,
    strain = if (master_lookup(sample, "cepa_id") != "") master_lookup(sample, "cepa_id") else label$strain,
    barcode = if (master_lookup(sample, "barcode_id") != "") master_lookup(sample, "barcode_id") else label$barcode,
    species_wgs = species,
    scheme = scheme,
    st = st,
    mlst_quality = quality,
    loci_detected = loci_detected,
    loci_total = loci_total,
    allele_profile = join(alleles),
    fasta = fasta,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
summary_df <- summary_df[order(suppressWarnings(as.integer(sub("^barcode", "", summary_df$barcode))), summary_df$sample), ]
write_tsv(summary_df, file.path(out_dir, "mlst_summary.tsv"))

typed <- summary_df[valid_scheme(summary_df$scheme) & valid_st(summary_df$st), , drop = FALSE]
if (nrow(typed) > 0) {
  split_key <- paste(typed$species_wgs, typed$scheme, typed$st, sep = "\r")
  count_rows <- lapply(split(typed, split_key), function(x) {
    data.frame(
      species_wgs = x$species_wgs[1],
      scheme = x$scheme[1],
      st = x$st[1],
      n = nrow(x),
      strains = join(x$strain),
      samples = join(x$sample),
      stringsAsFactors = FALSE
    )
  })
  st_counts <- do.call(rbind, count_rows)
  st_counts <- st_counts[order(st_counts$species_wgs, st_counts$scheme, suppressWarnings(as.integer(st_counts$st)), st_counts$st), ]
} else {
  st_counts <- data.frame(
    species_wgs = character(0),
    scheme = character(0),
    st = character(0),
    n = integer(0),
    strains = character(0),
    samples = character(0),
    stringsAsFactors = FALSE
  )
}
write_tsv(st_counts, file.path(out_dir, "mlst_st_counts.tsv"))

quality_counts <- as.data.frame(table(summary_df$mlst_quality), stringsAsFactors = FALSE)
names(quality_counts) <- c("mlst_quality", "n")
write_tsv(quality_counts, file.path(out_dir, "mlst_quality_counts.tsv"))

species_lines <- character(0)
for (species in sort(unique(summary_df$species_wgs))) {
  x <- summary_df[summary_df$species_wgs == species, , drop = FALSE]
  x_typed <- x[valid_scheme(x$scheme) & valid_st(x$st), , drop = FALSE]
  st_label <- if (nrow(x_typed) > 0) {
    groups <- split(x_typed, paste(x_typed$scheme, x_typed$st, sep = "\r"))
    paste(vapply(groups, function(group) {
      paste0("ST", group$st[1], " (cepas ", join(group$strain), ")")
    }, character(1)), collapse = "; ")
  } else {
    "sin ST asignado"
  }
  species_lines <- c(species_lines, paste0("- ", species, ": ", nrow(x), " cepas; ", st_label))
}

md <- c(
  "# MLST",
  "",
  paste0("- Fecha: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- Herramienta: ", mlst_bin),
  paste0("- Cepas procesadas: ", nrow(summary_df)),
  paste0("- ST asignados: ", nrow(typed)),
  paste0("- Perfiles/especies sin ST exacto: ", nrow(summary_df) - nrow(typed)),
  "",
  "## Distribucion",
  species_lines,
  "",
  "## Nota",
  "El resultado depende de los esquemas incluidos en la base de `mlst`. Si una especie no tiene esquema disponible o el perfil de alelos no coincide con un ST definido, se conserva el perfil de alelos y se marca como `sin_esquema` o `perfil_sin_ST`."
)
writeLines(md, file.path(out_dir, "mlst_interpretation.md"))

message("Escribi: ", file.path(out_dir, "mlst_summary.tsv"))
message("Escribi: ", file.path(out_dir, "mlst_st_counts.tsv"))
message("Escribi: ", file.path(out_dir, "mlst_interpretation.md"))
