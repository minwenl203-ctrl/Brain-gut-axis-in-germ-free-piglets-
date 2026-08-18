# Differential expression analysis: MMT_D30 vs GF_D30
# per cell type within each tissue (Seurat FindMarkers, Wilcoxon test)

library(Seurat)
library(dplyr)

output_dir <- "./"
scobj <- readRDS("Brain.rds")

stopifnot(all(c("subtype", "group", "Tissue") %in% colnames(scobj@meta.data)))

scobj$Tissue_label <- as.character(scobj$Tissue)
Tissues <- sort(unique(na.omit(scobj$Tissue_label)))

all_de_results <- list()
Tissue_de_results <- list()

for (current_Tissue in Tissues) {
  Tissue_obj <- subset(scobj, subset = Tissue_label == current_Tissue)
  subtypes <- sort(unique(na.omit(as.character(Tissue_obj$subtype))))
  ct_to_num <- setNames(seq_along(subtypes), subtypes)
  Tissue_result_list <- list()

  for (ct in subtypes) {
    subset_obj <- subset(Tissue_obj, subset = subtype == ct)
    valid_groups <- na.omit(unique(as.character(subset_obj$group)))
    if (!all(c("MMT_D30", "GF_D30") %in% valid_groups)) next

    Idents(subset_obj) <- "group"

    de_res <- tryCatch({
      FindMarkers(subset_obj, ident.1 = "MMT_D30", ident.2 = "GF_D30", test.use = "wilcox", min.pct = 0.25, logfc.threshold = 0.25)
    }, error = function(e) NULL)
    if (is.null(de_res) || nrow(de_res) == 0) next

    de_res$gene <- rownames(de_res)
    de_res$Tissue <- current_Tissue
    de_res$cluster <- ct_to_num[ct]
    de_res$celltype <- ct
    de_res$direction <- ifelse(de_res$avg_log2FC > 0, "up_in_MMT_D30", "up_in_GF_D30")
    de_res <- de_res[, c("gene", "Tissue", "cluster", "celltype", "direction", setdiff(colnames(de_res), c("gene", "Tissue", "cluster", "celltype", "direction"))), drop = FALSE]

    all_de_results[[paste(current_Tissue, ct, sep = "__")]] <- de_res
    Tissue_result_list[[ct]] <- de_res
  }

  if (length(Tissue_result_list) > 0) {
    Tissue_de_results[[current_Tissue]] <- bind_rows(Tissue_result_list)
  }
}

if (length(all_de_results) > 0) {
  all_de_combined <- bind_rows(all_de_results)
  all_de_combined$direction <- factor(all_de_combined$direction, levels = c("up_in_MMT_D30", "up_in_GF_D30"))
  write.csv(all_de_combined, file.path(output_dir, "MMT_D30_vs_GF_D30_all_markers.csv"), row.names = FALSE)

  for (current_Tissue in names(Tissue_de_results)) {
    Tissue_de <- Tissue_de_results[[current_Tissue]] %>% filter(p_val_adj < 0.05)
    if (nrow(Tissue_de) > 0) {
      filename <- paste0("MMT_D30_vs_GF_D30_", current_Tissue, "_all_DEGs.csv")
      write.csv(Tissue_de, file.path(output_dir, filename), row.names = FALSE)
    }
  }
}
