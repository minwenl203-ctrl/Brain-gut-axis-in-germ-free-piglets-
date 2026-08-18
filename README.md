# scRNA-seq Analysis Pipeline

Analysis scripts for the single-cell RNA-seq study of brain tissues,
covering differential expression analysis, cell-type marker visualization,
GO enrichment, Ro/e tissue distribution analysis, and volcano plots.

## Repository Contents

| File | Description |
| --- | --- |
| `DEA.R` | Differential expression analysis (`MMT_D30` vs `GF_D30`) per tissue and cell type using Seurat `FindMarkers` (Wilcoxon test) |
| `plot_Volcano.R` | Volcano plots of significant DEGs per tissue using `scMarkerViz` |
| `go_enrichment.R` | GO Biological Process enrichment of DEGs from selected neuronal cell types |
| `GVis.R` | `FindAllMarkers` + `ClusterGVis` heatmap / line-plot visualization with GO annotation |
| `plot_RoE.sh` | SLURM job script for Ro/e tissue-distribution analysis (`STARTRAC`) with heatmaps and a bubble plot |
| `environment.yml` | Conda environment definition (name: `DEA`) |
| `install_packages.R` | Installs the GitHub-hosted R packages |

## Environment Setup

All analyses were performed in the conda environment `DEA`:

```bash
conda env create -f environment.yml
conda activate DEA
Rscript install_packages.R   # installs GitHub-hosted R packages
```

## Input Data

Place the following Seurat objects in the working directory
(adjust paths in the scripts if needed):

- `Brain.rds` — metadata columns required: `subtype`, `group`, `Tissue`
- `Brain_new.rds` — metadata columns required: `manual_celltype`, `group`

## Workflow

1. **DEA.R** — runs `FindMarkers` comparing `MMT_D30` vs `GF_D30` for every
   cell type within each tissue; outputs
   `MMT_D30_vs_GF_D30_<Tissue>_all_DEGs.csv`.
2. **plot_Volcano.R** — generates volcano plots (`*_gail_volcano_plot.pdf`)
   from the `*_all_DEGs.csv` files.
3. **go_enrichment.R** — GO BP enrichment of significant DEGs for
   GABAergic / Glutamatergic / Immature neurons; outputs `GO_results/`.
4. **GVis.R** — per-cell-type `FindAllMarkers`, `ClusterGVis` visualizations,
   and GO enrichment; outputs one folder per cell type.
5. **plot_RoE.sh** — Ro/e distribution analysis of glial subtypes
   (`Oligodendrocytes`, `Microglia`, `Astrocytes`); submit with
   `sbatch plot_RoE.sh` on SLURM clusters, or run the R block directly.

## Notes

- `scMarkerViz` (used by `plot_Volcano.R`) is not on CRAN or Bioconductor;
  it is installed from <https://github.com/DaoqiuWang/scMarkerViz> by
  `install_packages.R`.
- The scripts expect output files to be written to the working directory.
