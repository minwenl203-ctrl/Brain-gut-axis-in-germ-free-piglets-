###############################################################################
# GO BP enrichment for DEGs of selected cell types from the DEA results
#
# Reads all *_all_DEGs.csv files and, for the three neuronal cell types
# (GABAergic neurons / Glutamatergic neurons / Immature neurons), runs GO
# Biological Process enrichment on significant DEGs
# (|avg_log2FC| >= fc_cutoff, p_val < pval_cutoff), keeping the top 5 terms.
#
# Dependencies: clusterProfiler, org.Ss.eg.db, dplyr
###############################################################################

suppressMessages({
  library(clusterProfiler)
  library(org.Ss.eg.db)
  library(dplyr)
})

WORK_DIR    <- "./"
OUT_DIR     <- file.path(WORK_DIR, "GO_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

target_celltypes <- c("GABAergic neurons", "Glutamatergic neurons", "Immature neurons")
fc_cutoff   <- 0.25   # DEG thresholds (consistent with jjVolcano02.R)
pval_cutoff <- 0.05

GO_ONTOLOGY      <- "BP"   # Biological Process
GO_PVAL_CUTOFF   <- 0.05
GO_QVAL_CUTOFF   <- 0.2
GO_TOP_N         <- 5      # top N terms per cell type

csv_pattern <- '_all_DEGs\\.csv$'

run_go_enrichment <- function(gene_list, celltype_name, tissue_label) {
  if (length(gene_list) < 5) return(NULL)

  # SYMBOL -> ENTREZID
  entrez <- tryCatch(bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Ss.eg.db), error = function(e) NULL)
  if (is.null(entrez) || nrow(entrez) < 5) return(NULL)

  go_result <- tryCatch(enrichGO(gene = entrez$ENTREZID, OrgDb = org.Ss.eg.db, ont = GO_ONTOLOGY, pAdjustMethod = "BH", pvalueCutoff = GO_PVAL_CUTOFF, qvalueCutoff = GO_QVAL_CUTOFF, readable = TRUE), error = function(e) NULL)
  if (is.null(go_result) || nrow(go_result) == 0) return(NULL)

  go_df <- go_result@result %>%
    arrange(p.adjust) %>%
    slice_head(n = GO_TOP_N)

  go_df$CellType <- celltype_name
  go_df$Tissue   <- tissue_label
  return(go_df)
}

plot_files <- list.files(path = WORK_DIR, pattern = csv_pattern, full.names = FALSE)
all_go_results <- list()

for (csv_name in plot_files) {
  mydata <- read.csv(file.path(WORK_DIR, csv_name), stringsAsFactors = FALSE)
  if (nrow(mydata) == 0) next
  if (!"celltype" %in% colnames(mydata)) next  # skip files without a celltype column

  # Extract tissue label from filename (e.g. AMY, HIP)
  tissue_label <- gsub("^MMT_D30_vs_GF_D30_|_all_DEGs\\.csv$", "", csv_name)

  sig_data <- subset(mydata,
    !is.na(celltype) &
    celltype %in% target_celltypes &
    abs(avg_log2FC) >= fc_cutoff &
    p_val < pval_cutoff
  )
  if (nrow(sig_data) == 0) next

  for (ct in target_celltypes) {
    ct_genes <- unique(sig_data$gene[sig_data$celltype == ct])
    go_df <- run_go_enrichment(ct_genes, ct, tissue_label)
    if (!is.null(go_df)) all_go_results[[paste(csv_name, ct, sep = "_")]] <- go_df
  }
}

if (length(all_go_results) > 0) {
  go_all <- bind_rows(all_go_results)

  # Save top N results per tissue
  for (tissue in unique(go_all$Tissue)) {
    tissue_df <- go_all %>% filter(Tissue == tissue)
    write.csv(tissue_df, file.path(OUT_DIR, paste0("GO_BP_top5_", tissue, ".csv")), row.names = FALSE)
  }

  # Combined summary table
  write.csv(go_all, file.path(OUT_DIR, "GO_BP_top5_all_tissues.csv"), row.names = FALSE)
}
