#!/usr/bin/env Rscript
# Figure 3: synteny comparison of the rlmH-cstA interval in 23EV612, cepa_18_barcode02, and GCA_031889545.1.
# CDS coordinates: Bakta GFF3 (results/seh_st1_bovine/phylogenomics_set/gff3/{23EV612,cepa_18_barcode02}.gff3)
#   and the curated ancestral-block CDS table for the ST2187 seh-negative genome
#   (results/seh_st1_bovine/st2187_ancestral_block/st2187_ancestral_block_cds_functional_table.tsv).
# Panel (a) 23EV612 is plotted in native genome coordinates (cstA already at the left edge).
# Panels (b) and (c) are mirrored (x -> length+1-x, strand flipped) so cstA is on the left / rlmH on the
# right in every panel, matching panel (a)'s reading orientation.
# Breakpoint uncertainty intervals: results/seh_st1_bovine/st2187_breakpoints/st2187_breakpoint_coordinates.tsv
# (measured on cepa_18_barcode02; the cstA-proximal offset is shared across panels because x = 0 is defined
# as the cstA-proximal outer edge of the extracted interval in all three panels).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(root, "results", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"

col_flank  <- "#999999"
col_minor  <- "#C6C6C6"
col_is655  <- "#0072B2"
col_seh    <- "#D55E00"
col_pyc    <- "#009E73"
col_rm     <- "#CC79A7"
col_break  <- "#F0E442"

panel_levels <- c(
  "a" = "(a) 23EV612 -- ST1, seh+ (8,571 bp)",
  "b" = "(b) cepa_18_barcode02 -- ST2187, seh+ (8,570 bp)",
  "c" = "(c) GCA_031889545.1 -- ST2187, seh- (12,453 bp)"
)

# ---- panel (a): 23EV612, contig_1:1880276-1888846, native orientation (cstA already left) ----
a_off <- 1880276
genes_a <- tibble::tribble(
  ~gene,                 ~start,  ~end,    ~strand, ~grp,   ~show_label,
  "cstA",                1880276, 1881343, "-",     "flank", TRUE,
  "cstR",                1881480, 1881740, "+",     "flank", FALSE,
  "tauE",                1881740, 1882495, "+",     "flank", FALSE,
  "Acetyltransferase",   1883046, 1883519, "-",     "minor", FALSE,
  "seh",                 1884415, 1885140, "-",     "seh",   TRUE,
  "IS655-like",          1885794, 1886423, "+",     "is655", TRUE,
  "NO VEIN, C-term.",    1886895, 1888151, "-",     "minor", FALSE,
  "rlmH",                1888367, 1888846, "-",     "flank", TRUE
) %>% dplyr::mutate(start = start - a_off, end = end - a_off, panel = "a")

# ---- panel (b): cepa_18_barcode02, contig_5:945979-954548, mirrored so cstA is on the left ----
b_len <- 954548 - 945979 + 1
mirror_b <- function(pos) (954548 + 945979) - pos
genes_b_raw <- tibble::tribble(
  ~gene,                 ~start,  ~end,    ~strand, ~grp,   ~show_label,
  "rlmH",                945979,  946458,  "+",     "flank", TRUE,
  "NO VEIN, C-term.",    946674,  947930,  "+",     "minor", FALSE,
  "IS655-like",          948402,  949031,  "-",     "is655", TRUE,
  "seh",                 949685,  950410,  "+",     "seh",   TRUE,
  "Acetyltransferase",   951306,  951779,  "+",     "minor", FALSE,
  "hypothetical",        951926,  952114,  "-",     "minor", FALSE,
  "tauE",                952329,  953084,  "-",     "flank", FALSE,
  "cstR",                953084,  953344,  "-",     "flank", FALSE,
  "cstA",                953481,  954548,  "+",     "flank", TRUE
)
genes_b <- genes_b_raw %>%
  dplyr::mutate(
    new_start = mirror_b(end), new_end = mirror_b(start),
    strand = ifelse(strand == "+", "-", "+"),
    start = new_start - 945979, end = new_end - 945979,
    panel = "b"
  ) %>%
  dplyr::select(gene, start, end, strand, grp, show_label, panel)

# ---- panel (c): GCA_031889545.1, ancestral-block coordinates (1..12453), mirrored so cstA is on the left ----
c_len <- 12453
mirror_c <- function(pos) (c_len + 1) - pos
genes_c_raw <- tibble::tribble(
  ~gene,          ~start, ~end,  ~strand, ~grp,
  "rlmH",         1,      480,   "+",     "flank",
  "pycC",         521,    1792,  "+",     "pyc",
  "pycTM",        1798,   2274,  "+",     "pyc",
  "RM_mtase",     2313,   3827,  "+",     "rm",
  "RM_spec",      3817,   5049,  "+",     "rm",
  "RM_endo",      5033,   8161,  "+",     "rm",
  "Membrane",     8273,   8980,  "+",     "minor",
  "hypothetical", 8973,   9476,  "+",     "minor",
  "dus",          9567,   10019, "-",     "flank",
  "tauE",         10234,  10989, "-",     "flank",
  "cstR",         10989,  11249, "-",     "flank",
  "cstA",         11386,  12453, "+",     "flank"
)
genes_c <- genes_c_raw %>%
  dplyr::mutate(
    new_start = mirror_c(end), new_end = mirror_c(start),
    strand = ifelse(strand == "+", "-", "+"),
    start = new_start, end = new_end,
    panel = "c",
    show_label = gene %in% c("cstA", "rlmH", "pycC", "pycTM", "dus")
  ) %>%
  dplyr::select(gene, start, end, strand, grp, show_label, panel)

genes <- dplyr::bind_rows(genes_a, genes_b, genes_c) %>%
  dplyr::mutate(panel = factor(panel, levels = names(panel_levels), labels = panel_levels))

pal <- c(flank = col_flank, minor = col_minor, is655 = col_is655, seh = col_seh, pyc = col_pyc, rm = col_rm)

gene_arrow_poly <- function(df) {
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    x0 <- df$start[i]; x1 <- df$end[i]
    head_len <- min(max((x1 - x0) * 0.10, 25), (x1 - x0) * 0.6)
    if (df$strand[i] == "+") {
      x <- c(x0, x1 - head_len, x1, x1 - head_len, x0)
    } else {
      x <- c(x1, x0 + head_len, x0, x0 + head_len, x1)
    }
    y <- c(-0.38, -0.38, 0, 0.38, 0.38)
    out[[i]] <- data.frame(x = x, y = y, id = i, gene = df$gene[i], grp = df$grp[i], panel = df$panel[i])
  }
  dplyr::bind_rows(out)
}
gene_poly <- gene_arrow_poly(genes)
gene_poly$uid <- interaction(gene_poly$panel, gene_poly$id, drop = TRUE)

label_df <- genes %>% dplyr::filter(show_label) %>%
  dplyr::mutate(xmid = (start + end) / 2)

# ---- breakpoint uncertainty zones ----
# cstA-proximal ("downstream/opposite") zone: same offset (x=0 is the cstA-proximal outer edge in all panels)
cstA_zone <- data.frame(xmin = -1045, xmax = 2505)
# rlmH-proximal zone: offset measured from each panel's own rlmH left edge (365-545 bp past that edge)
rlmH_left <- c(a = 8091, b = 8090, c = 11974)
rlmH_zone <- data.frame(
  panel = factor(names(panel_levels), levels = names(panel_levels), labels = panel_levels),
  xmin = rlmH_left + 365, xmax = rlmH_left + 545
)
cstA_zone_rep <- rlmH_zone %>% dplyr::mutate(xmin = cstA_zone$xmin, xmax = cstA_zone$xmax)

x_lo <- -1250
x_hi <- c_len + 750

p <- ggplot() +
  geom_rect(data = cstA_zone_rep, aes(xmin = xmin, xmax = xmax, ymin = -0.62, ymax = 0.62),
            fill = col_break, alpha = 0.55, color = NA) +
  geom_rect(data = rlmH_zone, aes(xmin = xmin, xmax = xmax, ymin = -0.62, ymax = 0.62),
            fill = col_break, alpha = 0.55, color = NA) +
  geom_polygon(data = gene_poly, aes(x = x, y = y, group = uid, fill = grp), color = "grey20", linewidth = 0.25) +
  geom_text(data = label_df, aes(x = xmid, y = 0.72, label = gene), size = 2.35, family = font_family,
            color = "grey15", fontface = "bold") +
  scale_fill_manual(values = pal, guide = "none") +
  facet_wrap(~panel, ncol = 1, scales = "fixed") +
  coord_cartesian(xlim = c(x_lo, x_hi), ylim = c(-1.15, 1.0), clip = "off") +
  theme_void(base_family = font_family) +
  theme(
    strip.text = element_text(hjust = 0, size = 8.2, face = "bold", family = font_family, color = "grey15",
                               margin = margin(b = 4, t = 8)),
    panel.spacing.y = unit(10, "mm")
  )

# panel-b identity annotation + panel-a / panel-c breakpoint labels (single instance each, avoids repetition)
p <- p +
  geom_text(data = data.frame(panel = factor(panel_levels["b"], levels = panel_levels), x = c_len / 2, y = 0.72),
            aes(x = x, y = y), label = "99.638% identity to ST1", size = 2.3, family = font_family,
            color = "grey35", fontface = "italic", vjust = 0, nudge_y = 0.16) +
  geom_text(data = data.frame(panel = factor(panel_levels["a"], levels = panel_levels),
                               x = mean(c(cstA_zone$xmin, cstA_zone$xmax)), y = -0.62),
            aes(x = x, y = y), label = "3,551 bp", size = 2.2, family = font_family, color = "#8a7600",
            vjust = 1, nudge_y = -0.05) +
  geom_text(data = data.frame(panel = factor(panel_levels["a"], levels = panel_levels),
                               x = rlmH_left[["a"]] + 455, y = -0.62),
            aes(x = x, y = y), label = "181 bp", size = 2.2, family = font_family, color = "#8a7600",
            vjust = 1, nudge_y = -0.05) +
  # scale bar (shared, drawn once under the bottom panel)
  geom_segment(data = data.frame(panel = factor(panel_levels["c"], levels = panel_levels)),
               aes(x = x_hi - 2000, xend = x_hi, y = -0.85, yend = -0.85), linewidth = 0.6, color = "black") +
  geom_text(data = data.frame(panel = factor(panel_levels["c"], levels = panel_levels)),
            aes(x = x_hi - 1000, y = -1.00), label = "2 kb", size = 2.5, family = font_family, color = "black") +
  # Type I R-M group label (panel c only)
  geom_segment(data = {
    rm_x <- genes_c %>% dplyr::filter(grp == "rm") %>% dplyr::summarise(x0 = min(start), x1 = max(end))
    data.frame(panel = factor(panel_levels["c"], levels = panel_levels), x = rm_x$x0, xend = rm_x$x1, y = 0.72, yend = 0.72)
  }, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4, color = "grey15") +
  geom_text(data = {
    rm_x <- genes_c %>% dplyr::filter(grp == "rm") %>% dplyr::summarise(x0 = min(start), x1 = max(end))
    data.frame(panel = factor(panel_levels["c"], levels = panel_levels), x = (rm_x$x0 + rm_x$x1) / 2, y = 0.72)
  }, aes(x = x, y = y), label = "Type I R-M", size = 2.35, family = font_family, color = "grey15",
  fontface = "bold", vjust = 0, nudge_y = 0.05)

svg_path <- file.path(out_dir, "manuscript_fig3_rlmH_cstA_synteny.svg")
png_path <- file.path(out_dir, "manuscript_fig3_rlmH_cstA_synteny.png")

ggsave(svg_path, p, device = svglite::svglite, width = 7.5, height = 6.2, units = "in")
ggsave(png_path, p, device = ragg::agg_png, width = 7.5, height = 6.2, units = "in", res = 300)

cat("Figure 3 SVG:", svg_path, "\n")
cat("Figure 3 PNG:", png_path, "\n")
