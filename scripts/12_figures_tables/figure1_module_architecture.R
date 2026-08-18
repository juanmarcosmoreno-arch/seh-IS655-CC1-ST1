#!/usr/bin/env Rscript
# Figure 1: proportional schematic of the seh-IS655-like module (23EV612 reference locus).
# Coordinates come from:
#   results/seh_st1_bovine/is3_coding_capacity/is655_like_coding_capacity.tsv (ORF1/ORF2 stops, overlap, longest intact ORF, DDN residues)
#   results/seh_st1_bovine/is3_boundaries/is3_boundary_best_candidate.tsv (element length, TIR)
#   results/seh_st1_bovine/seh_integrity/seh_upstream_regulatory_audit.tsv (SD/-35/-10 positions, spacer length)
#   results/seh_st1_bovine/st1_st2187_phylogenomics/locus_state_matrix.tsv (seh variant position)

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(root, "results", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"

# Okabe-Ito colorblind-safe palette
col_is655   <- "#0072B2"
col_is655_l <- "#8FB8D8"
col_spacer  <- "#BFBFBF"
col_seh     <- "#D55E00"
col_stop    <- "#000000"
col_tir     <- "#4A5A66"

# ---- element geometry (bp, element-local coordinates, single reading orientation) ----
is655_len  <- 1211
spacer_len <- 653
seh_len    <- 726

is655_start  <- 0
spacer_start <- is655_start + is655_len
seh_start    <- spacer_start + spacer_len
total_len    <- seh_start + seh_len

orf1 <- c(116, 397)
orf2 <- c(361, 1191)
orf1_stops <- c(161, 272)
orf2_stops <- c(1069)
overlap <- c(361, 397)
longest_orf <- c(442, 1071)   # 209 aa / 630 bp, nested in ORF2
tir_len <- 9

sd_start   <- spacer_start + 642; sd_end   <- spacer_start + 647   # AGGAGT
m10_start  <- spacer_start + 616; m10_end  <- spacer_start + 621   # TATAAT
m35_start  <- spacer_start + 594; m35_end  <- spacer_start + 599   # TTAAAA

seh_variant_pos <- seh_start + 365

# ---- helper: right-pointing gene arrow as a polygon ----
gene_arrow <- function(xmin, xmax, y, h, head_frac = 0.12, min_head = 8) {
  head_len <- min(max((xmax - xmin) * head_frac, min_head), (xmax - xmin) * 0.6)
  data.frame(
    x = c(xmin, xmax - head_len, xmax, xmax - head_len, xmin),
    y = c(y - h / 2, y - h / 2, y, y + h / 2, y + h / 2)
  )
}

# ---- rows (data frames) ----
rects <- data.frame(
  xmin = c(is655_start, spacer_start, seh_start),
  xmax = c(spacer_start, seh_start, seh_start + seh_len),
  y = 0, ymin = -0.42, ymax = 0.42,
  fill = c(col_is655, col_spacer, NA),
  label = c("IS655-like\n(1,211 bp)", "spacer\n(653 nt)", NA),
  stringsAsFactors = FALSE
)

seh_arrow <- gene_arrow(seh_start, seh_start + seh_len, y = 0, h = 0.84, head_frac = 0.08)

orf_y <- 1.30
orf_h <- 0.5
orf1_poly <- gene_arrow(orf1[1], orf1[2], orf_y, orf_h, head_frac = 0.18)
orf2_poly <- gene_arrow(orf2[1], orf2[2], orf_y, orf_h, head_frac = 0.10)
orf1_poly$orf <- "ORF1"
orf2_poly$orf <- "ORF2"
orf_polys <- rbind(orf1_poly, orf2_poly)

longest_y <- orf_y + orf_h / 2 + 0.30

tir_flag_h <- 0.16
tir_flag_top <- 0.42 + 0.05
tir_tri <- function(edge = c("left", "right")) {
  edge <- match.arg(edge)
  if (edge == "left") {
    x0 <- is655_start; x1 <- is655_start + tir_len
  } else {
    x0 <- spacer_start; x1 <- spacer_start - tir_len
  }
  data.frame(x = c(x0, x1, (x0 + x1) / 2), y = c(tir_flag_top + tir_flag_h, tir_flag_top + tir_flag_h, tir_flag_top))
}
tir_df <- rbind(
  cbind(tir_tri("left"), part = "L"),
  cbind(tir_tri("right"), part = "R")
)

stop_marks <- data.frame(x = c(orf1_stops, orf2_stops), y = orf_y)

overlap_rect <- data.frame(xmin = overlap[1], xmax = overlap[2], ymin = orf_y - orf_h / 2 - 0.06, ymax = orf_y + orf_h / 2 + 0.06)

p <- ggplot() +
  geom_rect(data = rects[1, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_is655, color = "grey25", linewidth = 0.3) +
  geom_rect(data = rects[2, ], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = col_spacer, color = "grey25", linewidth = 0.3) +
  geom_polygon(data = seh_arrow, aes(x = x, y = y), fill = col_seh, color = "grey25", linewidth = 0.3) +
  geom_polygon(data = tir_df, aes(x = x, y = y, group = part), fill = col_tir, color = NA) +
  geom_rect(data = overlap_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "grey20", linetype = "dashed", linewidth = 0.4) +
  geom_polygon(data = orf1_poly, aes(x = x, y = y), fill = col_is655_l, color = "grey15", linewidth = 0.3) +
  geom_polygon(data = orf2_poly, aes(x = x, y = y), fill = col_is655, color = "grey15", linewidth = 0.3) +
  geom_point(data = stop_marks, aes(x = x, y = y), shape = 4, size = 2.6, stroke = 1.1, color = col_stop) +
  annotate("segment", x = longest_orf[1], xend = longest_orf[2], y = longest_y, yend = longest_y,
           linewidth = 0.5, color = "grey20") +
  annotate("segment", x = c(longest_orf[1], longest_orf[2]), xend = c(longest_orf[1], longest_orf[2]),
           y = longest_y - 0.06, yend = longest_y + 0.06, linewidth = 0.5, color = "grey20") +
  annotate("text", x = mean(longest_orf), y = longest_y + 0.16,
           label = "Longest intact ORF: 209 aa (630 bp)", size = 2.7, family = font_family, color = "grey20") +
  annotate("segment", x = c(m35_start, m10_start, sd_start), xend = c(m35_start, m10_start, sd_start),
           y = -0.42, yend = -0.62, linewidth = 0.4, color = c("#7D8991", "#7D8991", "#333333")) +
  annotate("segment", x = mean(c(m35_start, sd_end)), xend = mean(c(m35_start, sd_end)),
           y = -0.62, yend = -0.95, linewidth = 0.4, color = "grey30") +
  annotate("text", x = mean(c(m35_start, sd_end)), y = -1.08,
           label = "sigma-A: TTAAAA / TATAAT   |   SD: AGGAGT (6 bp to ATG)",
           size = 2.55, family = font_family, color = "grey20") +
  annotate("point", x = seh_start + 2, y = 0, shape = 21, size = 1.8, fill = "white", color = "grey20") +
  annotate("point", x = seh_variant_pos, y = 0, shape = 8, size = 2.6, color = "black") +
  annotate("segment", x = seh_variant_pos, xend = seh_variant_pos, y = 0.42, yend = 0.60,
           linewidth = 0.35, color = "grey30") +
  annotate("text", x = seh_variant_pos, y = 0.66,
           label = "pos 365 / aa 122\n(GCA_001680935.1)", size = 2.5, family = font_family, color = "grey20",
           lineheight = 0.85, vjust = 0) +
  annotate("text", x = (is655_start + spacer_start) / 2, y = -0.62, label = "IS655-like (1,211 bp)",
           size = 2.9, family = font_family, fontface = "bold", color = "grey15") +
  annotate("text", x = (spacer_start + seh_start) / 2, y = -0.62, label = "spacer (653 nt)",
           size = 2.9, family = font_family, fontface = "bold", color = "grey15") +
  annotate("text", x = seh_start + seh_len / 2, y = -0.62, label = "seh (726 bp)",
           size = 2.9, family = font_family, fontface = "bold", color = "grey15") +
  annotate("text", x = is655_start + tir_len + 18, y = tir_flag_top + tir_flag_h + 0.02,
           label = "TIR 9 nt (ACATTTTTT)", size = 2.05, family = font_family, color = "grey30",
           hjust = 0, vjust = 0) +
  annotate("text", x = mean(orf2), y = orf_y + orf_h / 2 + 0.62,
           label = "DDN motif: D139, D157 retained; E192->N", size = 2.7, family = font_family, color = "grey15") +
  annotate("text", x = mean(orf1), y = orf_y + orf_h / 2 + 0.14, label = "ORF1", size = 2.3, family = font_family, color = "grey15", fontface = "bold") +
  annotate("text", x = mean(orf2), y = orf_y, label = "ORF2", size = 2.3, family = font_family, color = "white", fontface = "bold") +
  annotate("segment", x = mean(overlap), xend = mean(overlap), y = overlap_rect$ymin - 0.02, yend = 0.78,
           linewidth = 0.35, color = "grey30") +
  annotate("text", x = overlap[2] + 12, y = 0.74, label = "no frameshift motif",
           size = 2.05, family = font_family, color = "grey30", hjust = 0, vjust = 1) +
  annotate("segment", x = total_len - 500, xend = total_len, y = -1.55, yend = -1.55, linewidth = 0.6, color = "black") +
  annotate("text", x = total_len - 250, y = -1.72, label = "500 bp", size = 2.6, family = font_family) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-1.9, 2.35)) +
  theme_void(base_family = font_family) +
  theme(plot.margin = margin(10, 14, 6, 14))

svg_path <- file.path(out_dir, "manuscript_fig1_seh_is655_architecture.svg")
png_path <- file.path(out_dir, "manuscript_fig1_seh_is655_architecture.png")

ggsave(svg_path, p, device = svglite::svglite, width = 7.2, height = 3.6, units = "in")
ggsave(png_path, p, device = ragg::agg_png, width = 7.2, height = 3.6, units = "in", res = 300)

cat("Figure 1 SVG:", svg_path, "\n")
cat("Figure 1 PNG:", png_path, "\n")
