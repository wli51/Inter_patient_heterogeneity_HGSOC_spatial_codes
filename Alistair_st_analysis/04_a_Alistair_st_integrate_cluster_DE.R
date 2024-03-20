library(Seurat)
library(magrittr)
library(dplyr)

### this script performs cluster based DE expression on integrated visium dataset with Seurat

### Prerequisites:  02_Alistair_st_integrate.R output

dir = file.path("..", "analysis_output", "Alistair_st")

obj.integrated = readRDS(file.path(dir, "visium_alistair_integrate_diet.rds"))

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
    DE_df_list[[c]] = DE_df
  } else {
    DE_df_list[[c]] = data.frame()
  }
  
}

saveRDS(DE_df_list, file = file.path("intermediate_output", "integrate_cluster_DE.rds"))