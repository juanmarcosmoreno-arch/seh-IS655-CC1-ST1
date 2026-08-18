#!/usr/bin/env Rscript
# Figure 2: joint ST1-ST2187 core-genome ML phylogeny.
# Tree: results/seh_st1_bovine/st1_st2187_phylogenomics/iqtree_gubbins_bootstrap/st1_st2187_recomb_filtered.treefile
#   (IQ-TREE, GTR+F+ASC, Gubbins-masked SKA alignment, 1000 UFBoot + 1000 SH-aLRT; node labels "SH-aLRT/UFBoot")
# Metadata: results/seh_st1_bovine/st1_st2187_phylogenomics/locus_state_matrix.tsv

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
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

tree_file <- file.path(root, "results", "seh_st1_bovine", "st1_st2187_phylogenomics",
                        "iqtree_gubbins_bootstrap", "st1_st2187_recomb_filtered.treefile")
matrix_file <- file.path(root, "results", "seh_st1_bovine", "st1_st2187_phylogenomics", "locus_state_matrix.tsv")

tree <- read.tree(tree_file)
mat <- read.delim(matrix_file, stringsAsFactors = FALSE)

n_tip <- Ntip(tree)

group_of <- function(st, seh_state) {
  dplyr::case_when(
    st == "ST1" ~ "ST1",
    st == "ST2187" & seh_state == "intact" ~ "ST2187 seh+",
    st == "ST2187" & seh_state == "absent" ~ "ST2187 seh-",
    TRUE ~ NA_character_
  )
}
mat$group <- group_of(mat$st, mat$seh_state)

tip_meta <- data.frame(node = match(mat$sample_id, tree$tip.label), group = mat$group, stringsAsFactors = FALSE)
tip_meta <- tip_meta[!is.na(tip_meta$node), ]

node_lab <- tree$node.label
node_ids <- (n_tip + 1):(n_tip + Nnode(tree))
support <- strsplit(node_lab, "/")
sh <- suppressWarnings(as.numeric(vapply(support, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))))
uf <- suppressWarnings(as.numeric(vapply(support, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))))
well_supported <- !is.na(sh) & !is.na(uf) & sh >= 80 & uf >= 95
support_class <- ifelse(is.na(sh) | is.na(uf), NA_character_,
                         ifelse(well_supported, "well supported", "not well supported"))
node_meta <- data.frame(
  node = node_ids,
  SHaLRT = sh,
  UFBoot = uf,
  well_supported = well_supported,
  support_class = support_class,
  stringsAsFactors = FALSE
)

meta <- dplyr::bind_rows(
  tip_meta %>% dplyr::mutate(SHaLRT = NA_real_, UFBoot = NA_real_, well_supported = NA, support_class = NA_character_),
  node_meta %>% dplyr::mutate(group = NA_character_)
)

st2187_tips <- mat$sample_id[mat$st == "ST2187"]
st2187_mrca <- getMRCA(tree, tip = st2187_tips)

pal <- c("ST1" = "#0072B2", "ST2187 seh+" = "#E69F00", "ST2187 seh-" = "#D55E00")
shape_pal <- c("well supported" = 16, "not well supported" = 1)

p <- ggtree(tree, size = 0.32, color = "grey30") %<+% meta

xmax <- max(p$data$x, na.rm = TRUE)

p <- p +
  geom_hilight(node = st2187_mrca, fill = "#F2E9DA", alpha = 0.7, extend = xmax * 0.12) +
  geom_tree(size = 0.32, color = "grey30") +
  geom_tippoint(aes(color = group), size = 1.5, show.legend = TRUE) +
  geom_tiplab(aes(color = group), size = 1.85, family = font_family, offset = xmax * 0.012, show.legend = FALSE) +
  geom_point2(aes(subset = !isTip, shape = support_class), size = 1.15, stroke = 0.35, color = "grey20", na.rm = TRUE) +
  geom_cladelab(node = st2187_mrca, label = "ST2187", offset = xmax * 0.30, offset.text = xmax * 0.02,
                barsize = 0.5, fontsize = 2.7, family = font_family, textcolor = "grey20", barcolor = "grey40") +
  scale_color_manual(values = pal, name = "Sequence type / seh status", na.translate = FALSE) +
  scale_shape_manual(values = shape_pal, name = "Node support\n(SH-aLRT>=80 & UFBoot>=95)", na.translate = FALSE) +
  geom_treescale(x = 0, y = -1.8, width = 0.0005, fontsize = 2.4, linesize = 0.4,
                  offset = 0.6, family = font_family) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::xlim(0, xmax * 1.42) +
  ggplot2::theme(
    legend.position = c(0.80, 0.14),
    legend.title = element_text(size = 6.6, family = font_family),
    legend.text = element_text(size = 6.2, family = font_family),
    legend.key.size = unit(3.0, "mm"),
    legend.background = element_blank(),
    plot.margin = margin(6, 10, 6, 6)
  )

svg_path <- file.path(out_dir, "manuscript_fig2_st1_st2187_phylogeny.svg")
png_path <- file.path(out_dir, "manuscript_fig2_st1_st2187_phylogeny.png")

ggsave(svg_path, p, device = svglite::svglite, width = 7.5, height = 9.5, units = "in")
ggsave(png_path, p, device = ragg::agg_png, width = 7.5, height = 9.5, units = "in", res = 300)

cat("Figure 2 SVG:", svg_path, "\n")
cat("Figure 2 PNG:", png_path, "\n")
