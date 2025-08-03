rm(list = ls())
pacman::p_load(Seurat, dplyr, data.table, ggplot2, cowplot, patchwork, ggrepel, hdf5r)
paad.data1 <- Read10X_h5(filename = here("1_data/TISCH_PAAD/PAAD_CRA001160_expression.h5"))
paad <- CreateSeuratObject(
  counts = paad.data1,
  project = "PAAD",
  min.cells = 3,
  min.features = 200
)
