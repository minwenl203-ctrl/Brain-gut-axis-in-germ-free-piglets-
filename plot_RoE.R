library(Startrac)
library(pheatmap)
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# Read Seurat object (must contain cell type annotation and group info)
sce <- readRDS("Brain_new.rds")  # adjust to your input file path
output_dir <- "./"

phe <- sce@meta.data

# Tissue*group interaction as cluster variable
phe$Tissue_group <- paste(phe$Tissue, phe$group, sep = "_")

group_order <- c("GF_D0", "GF_D30", "MMT_D30")
phe <- phe %>% filter(group %in% group_order)
phe$group <- factor(phe$group, levels = group_order)

# Keep selected cell types (glial subtypes)
subtype_keep <- c("Oligodendrocytes", "Microglia", "Astrocytes")
phe <- phe %>% filter(subtype %in% subtype_keep)
phe$subtype <- factor(phe$subtype, levels = subtype_keep)

Roe <- calTissueDist(phe, byPatient = FALSE, colname.patient = "orig.ident", colname.cluster = "Tissue_group", colname.tissue = "subtype", method = "chisq", min.rowSum = 0)

# Convert table to a plain numeric matrix
Roe <- unclass(Roe)

# Export Ro/e matrix as CSV
Roe_csv <- as.data.frame(Roe) %>% rownames_to_column("Tissue_Group")
write.csv(Roe_csv, file.path(output_dir, "RoE_table.csv"), row.names = FALSE, quote = FALSE)

# Symbol matrix for enrichment level display
symbol_mat <- matrix("", nrow = nrow(Roe), ncol = ncol(Roe), dimnames = dimnames(Roe))
symbol_mat[Roe > 2]               <- "+++"   # strong enrichment
symbol_mat[Roe > 1 & Roe <= 2]    <- "++"    # moderate enrichment
symbol_mat[Roe >= 0.5 & Roe <= 1] <- "+"     # mild enrichment / normal
symbol_mat[Roe > 0 & Roe < 0.5]   <- "+/-"   # mild depletion
symbol_mat[Roe == 0]              <- "-"     # absent

# Color breaks and palette (grey-blue gradient)
breaks <- c(0, 0.5, 1, 2, max(Roe))
colors <- c("#F1F1F1", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")

# Heatmap with enrichment symbols
pheatmap(Roe, scale = "none", display_numbers = symbol_mat, number_color = "black", fontsize_number = 10, color = colorRampPalette(colors)(100), breaks = seq(min(breaks), max(breaks), length.out = 101), cluster_rows = TRUE, cluster_cols = TRUE, border_color = "grey85", fontsize = 10, cellwidth = 30, cellheight = 25, filename = file.path(output_dir, "RoE_heatmap_symbol.pdf"))

# Bubble plot data: long format + cell counts
roe_long <- as.data.frame(Roe) %>% rownames_to_column("Tissue_Group") %>% pivot_longer(cols = -Tissue_Group, names_to = "CellType", values_to = "RoE")
cell_counts <- phe %>% count(Tissue_group, subtype, name = "Count") %>% dplyr::rename(Tissue_Group = Tissue_group, CellType = subtype)

plot_data <- roe_long %>%
  left_join(cell_counts, by = c("Tissue_Group", "CellType")) %>%
  mutate(
    Significance = case_when(
      RoE > 2 & Count > 50 ~ "***",
      RoE > 1.5 & Count > 30 ~ "**",
      RoE > 1.2 ~ "*",
      TRUE ~ ""
    ),
    log_count = log10(Count + 1)
  )

ggplot(plot_data, aes(x = Tissue_Group, y = CellType)) +
  geom_point(aes(size = log_count, fill = RoE), shape = 21, color = "black", stroke = 0.5) +
  geom_text(aes(label = Significance), color = "white", size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#738CC5", mid = "white", high = "#50BDBE", midpoint = 1, limits = c(0, max(plot_data$RoE)), name = "Ro/e") +
  scale_size_continuous(name = "Log10(Count+1)", range = c(3, 15)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_line(color = "grey90"), legend.position = "right") +
  labs(title = "Ro/e Enrichment with Cell Abundance (Tissue*Group vs CellType)", x = "Tissue*Group", y = "Cell Type")

ggsave(file.path(output_dir, "RoE_bubble_plot_glail.pdf"), width = 10, height = 5)

# Heatmap with numeric values
pheatmap(Roe, scale = "none", display_numbers = TRUE, number_color = "black", number_format = "%.2f", color = colorRampPalette(c("#738CC5", "white", "#50BDBE"))(100), cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 10, fontsize_row = 10, fontsize_col = 10, border_color = "grey80", cellwidth = 30, cellheight = 25, filename = file.path(output_dir, "RoE_heatmap.pdf"))
