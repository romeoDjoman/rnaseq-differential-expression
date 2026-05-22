# ==============================================================================
# 03_enrichment.R
# Enrichissement fonctionnel GO et KEGG avec clusterProfiler
# Auteur: Roméo DJOMAN | AgroParisTech / Paris-Saclay
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(enrichplot)
})

cat("[INFO] Chargement des résultats DESeq2...\n")
res_df <- read.csv("results/DEG_results.csv")

# Gènes significatifs
deg_up <- res_df %>% filter(significance == "Up") %>% pull(gene_id)
deg_dn <- res_df %>% filter(significance == "Down") %>% pull(gene_id)
cat(sprintf("[INFO] %d gènes UP, %d gènes DOWN\n", length(deg_up), length(deg_dn)))

# --- Conversion des IDs Gene Symbol → Entrez ID ----------------------------
entrez_up <- bitr(deg_up, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)$ENTREZID
bg_genes   <- bitr(res_df$gene_id, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)$ENTREZID

# --- Enrichissement KEGG ---------------------------------------------------
cat("[INFO] Enrichissement KEGG...\n")
kegg_res <- enrichKEGG(
  gene         = entrez_up,
  universe     = bg_genes,
  organism     = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# Dotplot KEGG
if (nrow(as.data.frame(kegg_res)) > 0) {
  kegg_plot <- dotplot(kegg_res, showCategory = 15, title = "Enrichissement KEGG — Gènes surexprimés") +
    theme(
      plot.title   = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(size = 9)
    ) +
    scale_color_gradient(low = "#EF9A9A", high = "#B71C1C")
  ggsave("figures/kegg_enrichment.png", kegg_plot, width = 10, height = 7, dpi = 200)
  cat("[INFO] Figure KEGG sauvegardée : figures/kegg_enrichment.png\n")
} else {
  cat("[WARN] Aucune voie KEGG significative trouvée.\n")
}

# --- Enrichissement GO (Biological Process) --------------------------------
cat("[INFO] Enrichissement GO (BP)...\n")
go_res <- enrichGO(
  gene          = entrez_up,
  universe      = bg_genes,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

if (nrow(as.data.frame(go_res)) > 0) {
  go_plot <- barplot(go_res, showCategory = 20,
                     title = "Enrichissement GO (Biological Process)") +
    theme(plot.title = element_text(face = "bold"))
  ggsave("figures/go_enrichment.png", go_plot, width = 11, height = 8, dpi = 200)
  write.csv(as.data.frame(go_res), "results/GO_enrichment.csv", row.names = FALSE)
  cat("[INFO] Figure GO sauvegardée + résultats exportés.\n")
}

cat("\n[SUCCESS] Analyse d'enrichissement terminée !\n")
