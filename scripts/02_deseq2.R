# ==============================================================================
# 02_deseq2.R
# Analyse différentielle RNA-seq avec DESeq2
# Auteur: Roméo DJOMAN | AgroParisTech / Paris-Saclay
# ==============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
})

cat("[INFO] Chargement des données...\n")

# --- 1. Chargement des données ------------------------------------------------
counts <- read.csv("data/counts_matrix.csv", row.names = 1)
metadata <- read.csv("data/metadata.csv", row.names = 1)

# Vérification alignement
stopifnot(all(colnames(counts) == rownames(metadata)))
cat(sprintf("[INFO] %d gènes x %d échantillons\n", nrow(counts), ncol(counts)))

# --- 2. Construction objet DESeq2 --------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Filtrage minimum
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
cat(sprintf("[INFO] Gènes après filtrage : %d\n", nrow(dds)))

# --- 3. Analyse différentielle -----------------------------------------------
cat("[INFO] Exécution DESeq2...\n")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "control"),
               alpha = 0.05)

# Shrinkage LFC (apeglm)
resLFC <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")

# Résumé
cat("\n[INFO] Résumé DESeq2 :\n")
summary(res)

# --- 4. Export résultats -------------------------------------------------------
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(padj) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

write.csv(res_df, "results/DEG_results.csv", row.names = FALSE)
cat("[INFO] Résultats exportés : results/DEG_results.csv\n")

# --- 5. Volcano Plot ---------------------------------------------------------
cat("[INFO] Génération volcano plot...\n")

# Top gènes à annoter
top_genes <- res_df %>%
  filter(significance != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 15)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj + 1e-300),
                                    color = significance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_point(data = top_genes, size = 2, alpha = 0.9) +
  geom_label_repel(
    data = top_genes,
    aes(label = gene_id),
    size = 2.5, max.overlaps = 15, box.padding = 0.3,
    segment.color = "grey50", fontface = "italic"
  ) +
  scale_color_manual(
    values = c("Up" = "#E53935", "Down" = "#1565C0", "NS" = "#BDBDBD"),
    labels = c("Up" = "Surexprimé", "Down" = "Sous-exprimé", "NS" = "Non significatif")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "Volcano Plot — Analyse Différentielle (DESeq2)",
    subtitle = "Condition : Traité vs. Contrôle | seuils : |log2FC| > 1, padj < 0.05",
    x = "log2 Fold Change",
    y = "-log10(padj)",
    color = "Expression"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  annotate("text", x = 3.5, y = 2, label = sprintf("n Up = %d", sum(res_df$significance == "Up")),
           color = "#E53935", fontface = "bold", size = 3.5) +
  annotate("text", x = -3.5, y = 2, label = sprintf("n Down = %d", sum(res_df$significance == "Down")),
           color = "#1565C0", fontface = "bold", size = 3.5)

ggsave("figures/volcano_plot.png", volcano_plot, width = 9, height = 7, dpi = 200)
cat("[INFO] Figure sauvegardée : figures/volcano_plot.png\n")

# --- 6. Heatmap Top 50 -------------------------------------------------------
cat("[INFO] Génération heatmap...\n")

vsd <- vst(dds, blind = FALSE)
top50_genes <- res_df %>%
  filter(significance != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 50) %>%
  pull(gene_id)

mat <- assay(vsd)[top50_genes, ]
mat <- mat - rowMeans(mat)  # centrage

anno_col <- data.frame(Condition = metadata$condition, row.names = rownames(metadata))

png("figures/heatmap_top50.png", width = 1800, height = 2400, res = 200)
pheatmap(
  mat,
  annotation_col   = anno_col,
  annotation_colors = list(Condition = c(control = "#42A5F5", treated = "#EF5350")),
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  clustering_method = "ward.D2",
  color            = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  breaks           = seq(-2.5, 2.5, length.out = 101),
  fontsize_row     = 7,
  fontsize_col     = 10,
  main             = "Heatmap — Top 50 Gènes Différentiellement Exprimés (VST)",
  border_color     = NA
)
dev.off()
cat("[INFO] Figure sauvegardée : figures/heatmap_top50.png\n")

cat("\n[SUCCESS] Analyse DESeq2 terminée !\n")
