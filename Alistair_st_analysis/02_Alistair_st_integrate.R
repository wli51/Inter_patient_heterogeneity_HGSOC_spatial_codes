library(Seurat)
library(dplyr)
library(magrittr)

### This script loads prepossessed visium data and performs seurat anchor based integration 
# anchor with 3000 features chosen by SelectIntegrationFeatures

### Prerequisites:  01_Alistair_st_preprocess.R output

# abbreviation for each sample to avoid having long names 
sample_abbrs = paste0("SP", 1:8)

# Load seurat objects 
obj.list = list()

for (i in 1:length(sample_abbrs)) {
  
  obj.list[[i]] = readRDS(file.path("..", "analysis_output", "Alistair_st", paste0("visium_", sample_abbrs[i], ".rds")))
  obj.list[[i]]@images = list() # having images will cause an error downstream in the integration workflow
  obj.list[[i]] = RenameCells(obj.list[[i]], add.cell.id = sample_abbrs[i])
  
}

common_gene_names = Reduce(intersect, sapply(obj.list, FUN = function(x) rownames(x@assays$SCT@data), simplify = F))

int_features <- SelectIntegrationFeatures(object.list = obj.list, 
                                          assay = rep("SCT", length(obj.list)), 
                                          nfeatures = 3000)

anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  assay = rep("SCT", length(obj.list)), 
                                  scale = TRUE, 
                                  reference = NULL, 
                                  reduction = "cca", 
                                  anchor.features = int_features)


obj.integrated <- 
  IntegrateData(anchorset = anchors, dims = 1:30, 
                features.to.integrate	= common_gene_names # explicitly stating that we want the top svgs in the integrated dataset
  ) 

rm(obj.list, anchors); gc()
obj.integrated.diet = DietSeurat(obj.integrated, assays = "integrated", layers = c("data"))

obj.integrated.diet = ScaleData(obj.integrated.diet, assay = "integrated", features = rownames(obj.integrated))
obj.integrated.diet <- RunPCA(obj.integrated.diet, verbose = FALSE, assay = "integrated")
obj.integrated.diet <- RunUMAP(obj.integrated.diet, dims = 1:30)
obj.integrated.diet <- FindNeighbors(obj.integrated.diet, dims = 1:30)
obj.integrated.diet <- FindClusters(obj.integrated.diet, resolution = 0.4)

obj.integrated.diet@assays$integrated@scale.data = matrix(nrow=0,ncol=0)

saveRDS(obj.integrated.diet, file = file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_diet.rds"))

metadata = obj.integrated.diet@meta.data

metadata %>% 
  dplyr::select(barcode, sample, seurat_clusters) %>%
  write.csv(file.path("intermediate_output", "seurat_cluster_metadata.csv"), 
            row.names=FALSE,
            quote = F)


