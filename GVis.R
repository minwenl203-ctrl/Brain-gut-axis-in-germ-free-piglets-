########################
# FindAllMarkers + ClusterGVis visualization + GO enrichment
########################

library(Seurat)
library(dplyr)
library(ggplot2)
suppressPackageStartupMessages(library(ClusterGVis))

if (!require("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

if (!require("jjAnno", quietly = TRUE)) devtools::install_github("junjunlab/jjAnno")
library(jjAnno)

input_rds <- "Brain_new.rds"  # adjust to your input file path
output_dir <- "./"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

scobj <- readRDS(input_rds)

stopifnot("group" %in% colnames(scobj@meta.data))

celltypes <- unique(scobj@meta.data$manual_celltype)

org.db.available <- FALSE
if (require("org.Ss.eg.db", quietly = TRUE)) {
  library(org.Ss.eg.db)
  org.db <- org.Ss.eg.db
  organism <- "ssc"
  org.db.available <- TRUE
}

for (ct in celltypes) {
  # Sanitize cell type name for folder name
  ct_safe <- gsub("[^a-zA-Z0-9_]", "_", ct)
  ct_dir <- file.path(output_dir, ct_safe)
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

  scobj_sub <- subset(scobj, subset = manual_celltype == ct)

  desired_order <- c("GF_D0", "GF_D30", "MMT_D30")
  present_groups <- intersect(desired_order, unique(scobj_sub@meta.data$group))
  scobj_sub@meta.data$group <- factor(scobj_sub@meta.data$group, levels = present_groups)
  Idents(scobj_sub) <- "group"

  n_groups <- length(unique(scobj_sub@meta.data$group))
  if (n_groups < 2) next

  # Step 1: FindAllMarkers
  markers.all <- tryCatch({
    Seurat::FindAllMarkers(scobj_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  }, error = function(e) NULL)
  if (is.null(markers.all) || nrow(markers.all) == 0) next
  write.csv(markers.all, file.path(ct_dir, "all_marker_genes.csv"), row.names = FALSE)

  # Top 5 markers per cluster by avg_log2FC
  markers <- markers.all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  write.csv(as.data.frame(markers), file.path(ct_dir, "top5_marker_genes.csv"), row.names = FALSE)

  # Step 2: Prepare visualization data and GO enrichment
  st.data1 <- prepareDataFromscRNA(object = scobj_sub, diffData = markers, showAverage = TRUE)

  enrich <- NULL
  enrich.available <- FALSE
  if (org.db.available) {
    enrich <- tryCatch({
      enrichCluster(object = st.data1, OrgDb = org.db, type = "BP", organism = organism, pvalueCutoff = 0.5, topn = 3, seed = 5201314)
    }, error = function(e) NULL)
    if (!is.null(enrich) && nrow(enrich) > 0) {
      write.csv(enrich, file.path(ct_dir, "go_enrichment_results.csv"), row.names = FALSE)
      enrich.available <- TRUE
    } else {
      enrich <- NULL
    }
  }

  # Step 3: Visualization
  markGenes <- unique(markers$gene)
  n_marker_genes <- length(markGenes)

  # Adaptive figure and font sizes
  hm_height <- max(6, n_marker_genes * 0.45 + 2)
  hm_width  <- max(8, n_groups * 2.5 + 6)
  line_width  <- max(8, n_marker_genes * 0.25 + 6)
  line_height <- max(5, n_marker_genes * 0.12 + 4)
  both_width  <- max(10, n_groups * 3 + 7)
  both_height <- if (enrich.available) max(6, n_marker_genes * 0.40 + 3) else max(5, n_marker_genes * 0.35 + 2)
  row_fontsize <- max(5, 10 - n_marker_genes * 0.2)
  col_fontsize <- max(6, 10 - n_groups * 0.5)

  # Line plot of cluster gene expression
  pdf(file.path(ct_dir, "cluster_gene_expression_line.pdf"), width = line_width, height = line_height)
  print(visCluster(object = st.data1, plot.type = "line"))
  dev.off()

  # Heatmap
  pdf(file.path(ct_dir, "cluster_gene_expression_heatmap.pdf"), width = hm_width, height = hm_height)
  print(visCluster(object = st.data1, plot.type = "heatmap", column_names_rot = 45, markGenes = markGenes, cluster.order = c(1:n_groups), row_names_gp = grid::gpar(fontsize = row_fontsize), column_names_gp = grid::gpar(fontsize = col_fontsize)))
  dev.off()

  # Combined plot: heatmap + GO annotation + line plot
  plot_both <- function(use_enrich = FALSE) {
    if (use_enrich) {
      n_enrich_terms <- nrow(enrich)
      go_colors <- rep(jjAnno::useMyCol("stallion", n = n_groups), each = 3)[1:n_enrich_terms]
      visCluster(object = st.data1, plot.type = "both", column_names_rot = 45, show_row_dend = FALSE, markGenes = markGenes, markGenes.side = "left", annoTerm.data = enrich, line.side = "left", cluster.order = c(1:n_groups), go.col = go_colors, add.bar = TRUE, row_names_gp = grid::gpar(fontsize = row_fontsize), column_names_gp = grid::gpar(fontsize = col_fontsize))
    } else {
      visCluster(object = st.data1, plot.type = "both", column_names_rot = 45, show_row_dend = FALSE, markGenes = markGenes, markGenes.side = "left", line.side = "left", cluster.order = c(1:n_groups), add.bar = TRUE, row_names_gp = grid::gpar(fontsize = row_fontsize), column_names_gp = grid::gpar(fontsize = col_fontsize))
    }
  }

  if (enrich.available) {
    # Fall back to a plot without GO annotation if enrichment plot fails
    result <- tryCatch({
      pdf(file.path(ct_dir, "cluster_heatmap_with_enrichment.pdf"), height = both_height, width = both_width, onefile = FALSE)
      print(plot_both(use_enrich = TRUE))
      dev.off()
      "enrich_ok"
    }, error = function(e) {
      if (length(dev.list()) > 0) dev.off()
      "fallback"
    })

    if (result == "fallback") {
      pdf(file.path(ct_dir, "cluster_heatmap_with_enrichment.pdf"), height = both_height, width = both_width, onefile = FALSE)
      print(plot_both(use_enrich = FALSE))
      dev.off()
    }
  } else {
    pdf(file.path(ct_dir, "cluster_heatmap_with_enrichment.pdf"), height = both_height, width = both_width, onefile = FALSE)
    print(plot_both(use_enrich = FALSE))
    dev.off()
  }
}
