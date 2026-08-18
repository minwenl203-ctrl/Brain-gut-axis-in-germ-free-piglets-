library(scMarkerViz)
library(ggplot2)
library(dplyr)

# Cell type colors (consistent with 07_plot.R)
celltype_cols <- c("Oligodendrocytes" = "#E3EDC0", "Microglia" = "#F2CA8D", "Astrocytes" = "#C1E6F3")

csv_pattern <- 'all_DEGs\\.csv$'
fc_cutoff <- 0.25    # |log2FC| threshold
pval_cutoff <- 0.05  # p-value threshold

work_dir <- './'
plot_files <- list.files(path = work_dir, pattern = csv_pattern, full.names = FALSE)

for (csv_name in plot_files) {
  csv_path <- file.path('./', csv_name)
  if (!file.exists(csv_path)) next

  mydata <- read.csv(csv_path)
  if (nrow(mydata) == 0) next
  if (!"celltype" %in% colnames(mydata)) next  # skip files without a celltype column

  plot_data <- subset(mydata, !is.na(celltype) & abs(avg_log2FC) >= fc_cutoff & p_val < pval_cutoff)
  plot_data <- subset(plot_data, celltype %in% names(celltype_cols))
  if (nrow(plot_data) == 0) next

  celltype_levels <- sort(unique(plot_data$celltype))

  # Map cell type names to colors; fill missing ones automatically
  tile_cols <- setNames(celltype_cols[celltype_levels], celltype_levels)
  na_cols <- is.na(tile_cols)
  if (any(na_cols)) tile_cols[na_cols] <- scales::hue_pal()(sum(na_cols))

  pdf_path <- sub('\\.csv$', '_gail_volcano_plot.pdf', csv_path)
  pdf(pdf_path, width = 8, height = 10)

  p <- marker_effect_plot(data = plot_data, group = "celltype", group_order = celltype_levels, group_label_size = 3, group_palette = tile_cols, effect_cutoff = fc_cutoff, significance_cutoff = pval_cutoff, significance_by = "p_value", show = "all", color_by = "direction", palette = c(down = '#835ac0ff', up = '#e1b85fff'), label = "top", label_n = 3, label_by = "effect", point_size = 0.5, point_alpha = 1, layout = "horizontal") +
    ggplot2::ylab("avg_log2FC")

  print(p)
  dev.off()
}
