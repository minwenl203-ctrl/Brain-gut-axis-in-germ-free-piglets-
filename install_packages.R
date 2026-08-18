# Install GitHub-hosted R packages required by this pipeline.
# Run after activating the conda environment:
#   conda activate DEA
#   Rscript install_packages.R

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("junjunlab/ClusterGVis")
remotes::install_github("junjunlab/jjAnno")
remotes::install_github("Japrin/STARTRAC")
remotes::install_github("DaoqiuWang/scMarkerViz")
