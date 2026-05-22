# ==============================================================================
# utils.R — Fonctions utilitaires pour l'analyse RNA-seq
# Auteur: Roméo DJOMAN
# ==============================================================================

#' Génère un MAplot propre depuis les résultats DESeq2
#' @param res_df data.frame des résultats DESeq2
#' @param output_path chemin de sauvegarde
plot_maplot <- function(res_df, output_path = "figures/MA_plot.png") {
  library(ggplot2)
  res_df$sig <- ifelse(is.na(res_df$padj), "NS",
                  ifelse(res_df$padj < 0.05, "Significatif", "NS"))
  p <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("Significatif" = "#E53935", "NS" = "#BDBDBD")) +
    labs(title = "MA-Plot (moyenne vs. fold change)",
         x = "log10(Mean expression)", y = "log2 Fold Change", color = "") +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")
  ggsave(output_path, p, width = 8, height = 6, dpi = 180)
  invisible(p)
}

#' Résumé statistique rapide
#' @param res_df data.frame des résultats
summary_deg <- function(res_df) {
  cat(sprintf(
    "Gènes totaux : %d\nDEG (padj<0.05) : %d\n  > Surexprimés (FC>2) : %d\n  > Sous-exprimés (FC<0.5) : %d\n",
    nrow(res_df),
    sum(res_df$padj < 0.05, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE),
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
  ))
}
