library(Seurat)
library(magrittr)
library(dplyr)

# this script performs cluster based DE expression analysis 

dir = file.path("..", "analysis_output", "Alistair_st")

obj.integrated = readRDS(file.path(dir, "visium_alistair_integrate_diet.rds"))
obj.integrated@graphs = readRDS(file = file.path(dir, "visium_alistair_integrate_graphs.rds"))
obj.integrated@reductions = readRDS(file = file.path(dir, "visium_alistair_integrate_reductions.rds"))
metadata = readRDS(file = file.path(dir, "visium_alistair_integrate_cluster.rds"))
obj.integrated$seurat_clusters = metadata$seurat_clusters
obj.integrated$sample = metadata$sample

DE_df_list = list()

for (c in levels(obj.integrated$seurat_clusters)) {
  
  DE_df = FindMarkers(obj.integrated, 
                      assay = "integrated",
                      ident.1 = c,
                      group.by = "seurat_clusters",
                      ident.2 = NULL,
                      features = rownames(obj.integrated@assays$integrated),
                      logfc.threshold = 0.25,
                      min.pct = 0.25)
  
  if (dim(DE_df)[1] != 0) {
    DE_df_list[[c]] = DE_df %>%
      dplyr::filter(avg_log2FC  > 0, p_val_adj <= 0.05) %>%
      dplyr::arrange(desc(avg_log2FC), desc(pct.1 - pct.2))
  } else {
    DE_df_list[[c]] = data.frame()
  }
  
}

saveRDS(DE_df_list, file = file.path("intermediate_output", "integrate_cluster_DE.rds"))