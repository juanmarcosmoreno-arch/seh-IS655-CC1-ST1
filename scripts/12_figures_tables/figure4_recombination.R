#!/usr/bin/env Rscript
# Figure 4: mosaic ancestry across the rlmH-cstA interval in ST2187 (cepa_18_barcode02).
# Ancestry-informative sites: results/seh_st1_bovine/st2187_breakpoints/ancestry_informative_sites.tsv
#   (cepa18_window_pos = position in the cepa_18_barcode02-normalized alignment window; ancestry_state =
#   ST1_like / ST2187_like, classified against 23EV612 and the ST2187 seh-negative reference)
# Structural breakpoint interval (recipient/cepa18 window coords):
#   results/seh_st1_bovine/st2187_breakpoints/st2187_breakpoint_coordinates.tsv
# 3SEQ breakpoint interval + module window: results/seh_st1_bovine/st2187_recombination_test/3seq_results.tsv

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

sites_file <- file.path(root, "results", "seh_st1_bovine", "st2187_breakpoints", "ancestry_informative_sites.tsv")
sites <- read.delim(sites_file, stringsAsFactors = FALSE)

col_st1   <- "#0072B2"
col_st2187 <- "#D55E00"
col_module <- "#F0E442"
col_struct <- "#E69F00"

# ---- sliding window (500 bp window, 50 bp step) proportion of ST1-like sites ----
win <- 500
step <- 50
pos_range <- range(sites$cepa18_window_pos)
centers <- seq(floor(pos_range[1] / step) * step, ceiling(pos_range[2] / step) * step, by = step)

roll <- lapply(centers, function(c) {
  in_win <- sites$cepa18_window_pos >= (c - win / 2) & sites$cepa18_window_pos <= (c + win / 2)
  n <- sum(in_win)
  if (n == 0) return(c(center = c, prop_st1 = NA_real_, n = 0))
  c(center = c, prop_st1 = mean(sites$ancestry_state[in_win] == "ST1_like"), n = n)
})
roll <- as.data.frame(do.call(rbind, roll))
roll <- roll[roll$n > 0, ]

# ---- breakpoints / module (cepa18 normalized window coordinates) ----
module_zone <- c(17424, 19432)
seq3_zones <- list(c(15459, 15479), c(20904, 20913))
struct_zones <- list(rlmH = c(14935, 15115), opposite = c(21065, 24615))

# ---- classification tally using the structural envelope, matches manuscript numbers exactly ----
inside <- sites$cepa18_window_pos >= struct_zones$rlmH[1] & sites$cepa18_window_pos <= struct_zones$opposite[2]
n_in_st1  <- sum(inside & sites$ancestry_state == "ST1_like")
n_in_neg  <- sum(inside & sites$ancestry_state == "ST2187_like")
n_out_st1 <- sum(!inside & sites$ancestry_state == "ST1_like")
n_out_neg <- sum(!inside & sites$ancestry_state == "ST2187_like")

x_lo <- max(0, pos_range[1] - 800)
x_hi <- pos_range[2] + 800

p <- ggplot() +
  annotate("rect", xmin = struct_zones$rlmH[1], xmax = struct_zones$rlmH[2], ymin = -0.08, ymax = 1.08,
           fill = col_struct, alpha = 0.30) +
  annotate("rect", xmin = struct_zones$opposite[1], xmax = struct_zones$opposite[2], ymin = -0.08, ymax = 1.08,
           fill = col_struct, alpha = 0.30) +
  annotate("rect", xmin = module_zone[1], xmax = module_zone[2], ymin = -0.08, ymax = 1.08,
           fill = col_module, alpha = 0.35) +
  annotate("rect", xmin = seq3_zones[[1]][1], xmax = seq3_zones[[1]][2], ymin = -0.08, ymax = 1.08,
           fill = col_st1, alpha = 0.85) +
  annotate("rect", xmin = seq3_zones[[2]][1], xmax = seq3_zones[[2]][2], ymin = -0.08, ymax = 1.08,
           fill = col_st1, alpha = 0.85) +
  geom_point(data = sites, aes(x = cepa18_window_pos, y = ifelse(ancestry_state == "ST1_like", 1, 0),
                                color = ancestry_state), size = 0.9, alpha = 0.55) +
  geom_line(data = roll, aes(x = center, y = prop_st1), color = "grey15", linewidth = 0.5) +
  scale_color_manual(values = c(ST1_like = col_st1, ST2187_like = col_st2187),
                      labels = c(ST1_like = "ST1-like site", ST2187_like = "ST2187-like site"),
                      name = NULL) +
  annotate("text", x = mean(module_zone), y = 1.28, label = "ST1-like module (IS655-like-spacer-seh)",
           size = 2.5, family = font_family, color = "grey20") +
  annotate("text", x = mean(struct_zones$rlmH), y = -0.22, label = "Structural\n(181 bp)",
           size = 2.15, family = font_family, color = "#8a5a00", lineheight = 0.85) +
  annotate("text", x = mean(struct_zones$opposite), y = -0.22, label = "Structural\n(3,551 bp)",
           size = 2.15, family = font_family, color = "#8a5a00", lineheight = 0.85) +
  annotate("text", x = mean(seq3_zones[[1]]), y = 1.14, label = "3SEQ", size = 2.15, family = font_family,
           color = "#00497a", hjust = 0.5) +
  annotate("text", x = mean(seq3_zones[[2]]), y = 1.14, label = "3SEQ", size = 2.15, family = font_family,
           color = "#00497a", hjust = 0.5) +
  annotate("text", x = mean(module_zone), y = 0.5,
           label = paste0("ST1-like (", n_in_st1, " sites)"), size = 2.6, family = font_family,
           color = "grey15", fontface = "bold") +
  annotate("text", x = x_lo + 1400, y = 0.5,
           label = paste0("ST2187-like\n(", n_out_neg, " sites)"), size = 2.5, family = font_family,
           color = "grey15", fontface = "bold", lineheight = 0.85) +
  scale_x_continuous(name = "Position in cepa_18_barcode02 normalized window (bp)",
                      labels = scales::comma) +
  scale_y_continuous(name = "Proportion ST1-like (500 bp sliding window)",
                      breaks = c(0, 0.5, 1), limits = c(-0.30, 1.42)) +
  coord_cartesian(xlim = c(x_lo, x_hi), clip = "off") +
  theme_minimal(base_family = font_family, base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(0.90, 0.55),
    legend.background = element_blank(),
    legend.text = element_text(size = 6.8, family = font_family),
    legend.key.size = unit(3, "mm"),
    axis.title = element_text(size = 8, family = font_family),
    axis.text = element_text(size = 7, family = font_family, color = "grey20"),
    plot.margin = margin(16, 10, 6, 6)
  )

svg_path <- file.path(out_dir, "manuscript_fig4_mosaic_ancestry.svg")
png_path <- file.path(out_dir, "manuscript_fig4_mosaic_ancestry.png")

ggsave(svg_path, p, device = svglite::svglite, width = 8.0, height = 4.2, units = "in")
ggsave(png_path, p, device = ragg::agg_png, width = 8.0, height = 4.2, units = "in", res = 300)

cat("Figure 4 SVG:", svg_path, "\n")
cat("Figure 4 PNG:", png_path, "\n")
cat(sprintf("Inside envelope: ST1-like=%d ST2187-like=%d | Outside: ST1-like=%d ST2187-like=%d\n",
            n_in_st1, n_in_neg, n_out_st1, n_out_neg))
