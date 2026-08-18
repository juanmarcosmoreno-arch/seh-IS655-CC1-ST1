#!/usr/bin/env Rscript
# Figure: CC1 phylogeny — module presence/absence on tips.
# v4: avoids gheatmap (ggplot2 4.0.0 incompatibility); uses color + shape only.

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)][1])
if (is.na(script_path)) script_path <- "scripts/140_cc1_phylogeny_figure.R"
root    <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
out_dir <- file.path(root, "results", "seh_st1_bovine", "cc1_phylogeny")

tree_file <- file.path(out_dir, "cc1_phylogeny_rooted_v2.nwk")
pres_file <- file.path(out_dir, "cc1_module_presence.tsv")

tree <- read.tree(tree_file)
tree$tip.label <- gsub("'", "", tree$tip.label)   # strip quotes added by dendropy
pres <- read.delim(pres_file, stringsAsFactors = FALSE)

# ---------- annotate ----------
pres <- pres %>%
  mutate(
    module_status = case_when(
      complete_module == "True"  & is655_like_present == "True"  ~ "Present",
      complete_module == "False" & is655_like_present == "True"  ~ "Partial (unknown)",
      TRUE                                                        ~ "Absent"
    ),
    cohort_shape_label = case_when(
      cohort == "ST772_historical" ~ "ST772 historical",
      cohort == "reference"        ~ "Reference",
      TRUE                         ~ "Jahnen 2025 (CC1/CC1660)"
    ),
    st_label = case_when(
      st %in% c("1", "81", "9300") ~ "ST1/ST81/ST9300",
      st == "1660"                  ~ "ST1660",
      st == "772"                   ~ "ST772",
      TRUE                          ~ "Unassigned"
    ),
    # label only the two reference genomes
    display_label = ifelse(sample_id %in% c("23EV612", "GCF_050908105.1"),
                           sample_id, "")
  )

# tip color = module status; shape = cohort; ST encoded via outline color
module_fill  <- c("Present" = "#1565C0",
                  "Absent"  = "#9E9E9E",
                  "Partial (unknown)" = "#FF8F00")
st_color     <- c("ST1/ST81/ST9300" = "#1565C0",
                  "ST1660"          = "#E65100",
                  "ST772"           = "#4E342E",
                  "Unassigned"      = "#9E9E9E")
cohort_shapes <- c("Jahnen 2025 (CC1/CC1660)" = 21,
                   "ST772 historical"          = 22,
                   "Reference"                 = 23)

# ---------- build plot ----------
p <- ggtree(tree, layout = "rectangular", ladderize = TRUE,
            size = 0.25, color = "gray50") %<+% pres

p <- p +
  # outer ring = ST group (stroke color)
  geom_tippoint(
    aes(fill = module_status, shape = cohort_shape_label, color = st_label),
    size = 2.5, stroke = 1.0
  ) +
  scale_fill_manual(
    values = module_fill, name = "Module status", na.translate = FALSE
  ) +
  scale_color_manual(
    values = st_color, name = "ST group", na.translate = FALSE
  ) +
  scale_shape_manual(
    values = cohort_shapes, name = "Dataset", na.translate = FALSE
  ) +
  geom_tiplab(
    aes(label = display_label),
    size = 2.2, hjust = -0.15, fontface = "italic", color = "black"
  ) +
  theme_tree2() +
  guides(
    fill  = guide_legend(
      override.aes = list(
        color = c("#9E9E9E", "#FF8F00", "#1565C0"),  # Absent / Partial / Present
        shape = 21, size = 3, stroke = 1
      )
    ),
    color = guide_legend(
      override.aes = list(fill = "gray80", shape = 21, size = 3, stroke = 1.5)
    ),
    shape = guide_legend(
      override.aes = list(fill = "gray70", color = "black", size = 3, stroke = 1)
    )
  ) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 7),
    legend.title    = element_text(size = 8, face = "bold"),
    plot.margin     = margin(4, 100, 4, 4),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank()
  ) +
  labs(
    title    = "CC1 core-genome phylogeny: seh+IS655-like module presence/absence",
    subtitle = sprintf(
      "Midpoint-rooted | %d tips | Parsimony: 1 gain, 0 losses (5 partials coded as unknown)",
      Ntip(tree)
    )
  )

# ---------- save ----------
ggsave(file.path(out_dir, "cc1_module_presence_tree_v2.pdf"),
       plot = p, width = 10, height = 16, dpi = 300)
ggsave(file.path(out_dir, "cc1_module_presence_tree_v2.png"),
       plot = p, width = 10, height = 16, dpi = 300)

message("Saved to: ", out_dir)
message(sprintf("Present: %d | Absent: %d | Partial (unknown): %d | Total tips: %d",
                sum(pres$module_status == "Present"),
                sum(pres$module_status == "Absent"),
                sum(pres$module_status == "Partial (unknown)"),
                Ntip(tree)))
