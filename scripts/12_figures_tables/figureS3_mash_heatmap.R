#!/usr/bin/env Rscript
# Figure S3: Mash distance heatmap for 53 ST1-ST2187 genomes.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

ROOT    <- normalizePath(file.path(dirname(
  sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
), ".."), mustWork = TRUE)
if (is.na(ROOT) || ROOT == "..") ROOT <- "/Volumes/JM/BioInf/GridION"

MASH    <- file.path(ROOT, "results/seh_st1_bovine/st1_st2187_phylogenomics/mash53_pairwise_distances.tsv")
META    <- file.path(ROOT, "results/seh_st1_bovine/phylogenomics_set/seh_st1_phylogenomics_set.tsv")
OUT     <- file.path(ROOT, "results/seh_st1_bovine/final_paper_figures")

# ---------- load data ----------
dist_df <- read.delim(MASH, stringsAsFactors = FALSE)
meta    <- read.delim(META, stringsAsFactors = FALSE)

# Build group label from metadata; fall back to accession-based guess
group_map <- meta %>%
  select(sample_id, tree_role, st, host, source) %>%
  mutate(group = case_when(
    grepl("ST2187", tree_role, ignore.case = TRUE) ~ "ST2187",
    st == "ST1" & grepl("Uruguay|local", tree_role, ignore.case = TRUE) ~ "Uruguay ST1 bovine",
    st == "ST1" & grepl("New Zealand|NZ", tree_role, ignore.case = TRUE) ~ "New Zealand ST1 bovine",
    st == "ST1" & grepl("China", tree_role, ignore.case = TRUE)          ~ "China ST1 bovine",
    st == "ST1" & grepl("Australia", tree_role, ignore.case = TRUE)      ~ "Australia ST1 bovine",
    st == "ST1" & grepl("non.bov|non_bov|human|clinical", tree_role, ignore.case = TRUE) ~ "Non-bovine ST1",
    st == "ST1" ~ "ST1",
    TRUE ~ tree_role
  ))

# All tips in the distance matrix
all_tips <- sort(unique(c(dist_df$sample_a, dist_df$sample_b)))

# Pivot to symmetric matrix
mat_wide <- dist_df %>%
  select(sample_a, sample_b, mash_distance) %>%
  bind_rows(rename(., sample_a = sample_b, sample_b = sample_a)) %>%
  distinct() %>%
  pivot_wider(names_from = sample_b, values_from = mash_distance, values_fill = 0) %>%
  tibble::column_to_rownames("sample_a")
mat_wide <- mat_wide[all_tips, all_tips]

# Cluster order
hc   <- hclust(as.dist(mat_wide), method = "ward.D2")
ord  <- hc$labels[hc$order]

# Long-format for ggplot, ordered
dist_long <- dist_df %>%
  bind_rows(rename(dist_df, sample_a = sample_b, sample_b = sample_a)) %>%
  distinct() %>%
  mutate(
    sample_a = factor(sample_a, levels = ord),
    sample_b = factor(sample_b, levels = rev(ord))
  )

# Group annotation strip (right side)
strip_df <- data.frame(sample_id = all_tips) %>%
  left_join(group_map, by = "sample_id") %>%
  mutate(
    group = ifelse(is.na(group), "ST2187", group),
    sample_id = factor(sample_id, levels = ord),
    x = 1
  )

group_colors <- c(
  "Uruguay ST1 bovine"     = "#1565C0",
  "New Zealand ST1 bovine" = "#29B6F6",
  "China ST1 bovine"       = "#EF5350",
  "Australia ST1 bovine"   = "#FF8F00",
  "Non-bovine ST1"         = "#66BB6A",
  "ST1"                    = "#78909C",
  "ST2187"                 = "#7B1FA2"
)

# ---------- heatmap ----------
p_heat <- ggplot(dist_long, aes(x = sample_a, y = sample_b, fill = mash_distance)) +
  geom_tile() +
  scale_fill_gradient(
    low = "#1A237E", high = "#FFECB3",
    name = "Mash distance",
    labels = scales::number_format(accuracy = 0.001)
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    plot.title   = element_text(size = 9, face = "bold"),
    plot.subtitle = element_text(size = 7, color = "gray40")
  ) +
  labs(
    title    = "Mash genomic distance — 53 ST1 and ST2187 genomes",
    subtitle = "Ward D2 clustering; darker = more similar"
  )

# Annotation strip (right)
p_strip <- ggplot(strip_df, aes(x = x, y = sample_id, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = group_colors, name = "Group", na.value = "#BDBDBD") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.y  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks   = element_blank(),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

# Compose with patchwork
library(patchwork)
p_final <- p_heat + p_strip +
  plot_layout(widths = c(10, 1), guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(OUT, "FigureS3_mash_heatmap_53genomes.pdf"),
       plot = p_final, width = 9, height = 8, dpi = 300)
ggsave(file.path(OUT, "FigureS3_mash_heatmap_53genomes.png"),
       plot = p_final, width = 9, height = 8, dpi = 300)

message("Saved FigureS3 to: ", OUT)
