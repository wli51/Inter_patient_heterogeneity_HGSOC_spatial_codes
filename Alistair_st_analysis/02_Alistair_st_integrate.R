library(Seurat)
library(dplyr)
library(magrittr)

### This script loads prepossessed visium data and performs seurat anchor based integration 
# anchor with 3000 features chosen by SelectIntegrationFeatures

sample_abbrs = paste0("SP", 1:8)

# Load seurat objects 
obj.list = list()

for (i in 1:length(sample_abbrs)) {
  
  obj.list[[i]] = readRDS(file.path("..", "analysis_output", "Alistair_st", paste0("visium_", sample_abbrs[i], ".rds")))
  obj.list[[i]]@images = list() # having images will cause an error downstream in the integration workflow
  obj.list[[i]] = RenameCells(obj.list[[i]], add.cell.id = sample_abbrs[i])
  
}

common_gene_names = Reduce(intersect, sapply(obj.list, FUN = function(x) rownames(x@assays$SCT@scale.data), simplify = F))

for (i in 1:length(obj.list)) {
  svgs = obj.list[[i]]@assays$SCT@meta.features %>% filter(moransi.spatially.variable) %>% arrange(moransi.spatially.variable.rank)
}

int_features <- SelectIntegrationFeatures(object.list = obj.list, assay = rep("SCT", length(obj.list)), nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = obj.list, scale = FALSE, 
                                  reference = NULL, reduction = "cca", 
                                  anchor.features = int_features, assay = rep("SCT", length(obj.list)))
obj.integrated <- 
  IntegrateData(anchorset = anchors, dims = 1:30, 
                features.to.integrate	= common_gene_names # explicitly stating that we want the top svgs in the integrated dataset
  ) 

obj.integrated@meta.data$sample = sample

obj.integrated = ScaleData(obj.integrated, assay = "integrated", features = rownames(obj.integrated))
obj.integrated <- RunPCA(obj.integrated, verbose = FALSE, assay = "integrated")
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:30)
obj.integrated <- FindClusters(obj.integrated, resolution = 0.4)

saveRDS(obj.integrated@graphs, 
        file = file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_graphs.rds"))
saveRDS(obj.integrated@reductions, 
        file = file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_reductions.rds"))
saveRDS(obj.integrated@meta.data %>% dplyr::select(sample, barcode, seurat_clusters), 
        file = file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_cluster.rds"))

obj.integrated.diet = DietSeurat(obj.integrated, assays = "integrated", scale.data = FALSE)

saveRDS(obj.integrated.diet, file = file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_diet.rds"))


