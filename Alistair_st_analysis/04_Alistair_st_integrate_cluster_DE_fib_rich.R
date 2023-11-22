library(Seurat)
library(magrittr)
library(dplyr)

# this script performs cluster based DE expression analysis on fibroblast rich clusters only

dir = file.path("..", "analysis_output", "Alistair_st")

obj.integrated = readRDS(file.path(dir, "visium_alistair_integrate_diet.rds"))
obj.integrated@graphs = readRDS(file = file.path(dir, "visium_alistair_integrate_graphs.rds"))
obj.integrated@reductions = readRDS(file = file.path(dir, "visium_alistair_integrate_reductions.rds"))
metadata = readRDS(file = file.path(dir, "visium_alistair_integrate_cluster.rds"))
obj.integrated$seurat_clusters = metadata$seurat_clusters
obj.integrated$sample = metadata$sample

DE_df_list = list()

mal_rich_cluster_idx = c("0", "2", "7", "9", "10")
macro_rich_cluster_idx = c("4")
stroma_rich_cluster_idx = c("1")
fib_rich_cluster_idx = c("3", "5", "6", "8", "11")

for (c in fib_rich_cluster_idx) {
  c_barcodes = rownames(obj.integrated@meta.data)[obj.integrated@meta.data$seurat_clusters == c]
  others_barcodes = rownames(obj.integrated@meta.data)[obj.integrated@meta.data$seurat_clusters %in% fib_rich_cluster_idx[!fib_rich_cluster_idx %in% c]]
  
  DE_df = FindMarkers( # comparing each fibroblast-rich cluster to all other fibroblast-rich clusters 
    obj.integrated@assays$integrated,
    assay = "integrated",
    cells.1 = c_barcodes,
    cells.2 = others_barcodes,
    features = rownames(obj.integrated@assays$integrated)
  )
  
  if (dim(DE_df)[1] != 0) {
    DE_df_list[[c]] = DE_df %>%
      dplyr::filter(avg_log2FC  > 0, p_val_adj <= 0.05) %>%
      dplyr::arrange(desc(avg_log2FC), desc(pct.1 - pct.2))
  } else {
    DE_df_list[[c]] = data.frame()
  }
  
}

saveRDS(DE_df_list, file = file.path("intermediate_output", "integrate_fib_rich_cluster_DE.rds"))

DE_df_list = list()

for (c in fib_rich_cluster_idx) {
  
  c_barcodes = rownames(obj.integrated@meta.data)[obj.integrated@meta.data$seurat_clusters == c]
  others_barcodes = rownames(obj.integrated@meta.data)[obj.integrated@meta.data$seurat_clusters %in% c(mal_rich_cluster_idx, macro_rich_cluster_idx)]
  
  DE_df = FindMarkers( # comparing each fibroblast-rich cluster to all  non-fibroblast-rich clusters 
    obj.integrated@assays$integrated,
    assay = "integrated",
    cells.1 = c_barcodes,
    cells.2 = others_barcodes,
    features = rownames(obj.integrated@assays$integrated),
    min.pct = 0.5
  )
  
  if (dim(DE_df)[1] != 0) {
    DE_df_list[[c]] = DE_df %>%
      dplyr::filter(avg_log2FC  > 0, p_val_adj <= 0.05) %>%
      dplyr::arrange(desc(avg_log2FC), desc(pct.1 - pct.2))
  } else {
    DE_df_list[[c]] = data.frame()
  }
  
}

saveRDS(
  DE_df_list,
  file = file.path(
    "intermediate_output",
    "integrate_fib_vs_other_cluster_DE.rds"
  )
)
