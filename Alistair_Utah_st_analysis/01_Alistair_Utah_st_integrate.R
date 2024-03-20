library(Seurat)
library(dplyr)
library(magrittr)

### This script loads prepossessed visium data and performs seurat anchor based integration 
# anchor with 3000 features chosen by SelectIntegrationFeatures

sources = c("Alistair", "Utah")
sample_abbrs_list = 
  list("Alistair" = paste0("SP", 1:8),
       "Utah" = c("2B", "3C", "4D"))
sample_abbrs = paste0("SP", 1:8)

hvgs = NULL

# Load seurat objects 
obj.list = list()

for (source in sources) {
  
  sample_abbrs = sample_abbrs_list[[source]]
  
  for (sample_abbr in sample_abbrs) {
    
    name = paste0(source, "_", sample_abbr)
    
    obj.list[[name]] = readRDS(file =
                                 file.path(
                                   "..",
                                   "analysis_output",
                                   paste0(source, "_st"),
                                   paste0("visium_", sample_abbr, ".rds")
                                 ))
    obj.list[[name]]@images = list() # having images will cause an error downstream in the integration workflow
    obj.list[[name]] = RenameCells(obj.list[[name]], add.cell.id = sample_abbr)
    hvgs = c(hvgs, VariableFeatures(obj.list[[name]]))
    
  }
  
}


common_gene_names = Reduce(intersect, sapply(obj.list, FUN = function(x) rownames(x@assays$SCT@data), simplify = F))
t = table(hvgs[hvgs %in% common_gene_names]) 
features2include = names(t)[t >= 2]
saveRDS(features2include, file = "var_genes.rds")

int_features <- SelectIntegrationFeatures(object.list = obj.list, 
                                          assay = rep("SCT", length(obj.list)), 
                                          nfeatures = 3000)

anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                  assay = rep("SCT", length(obj.list)), 
                                  scale = TRUE, 
                                  reference = NULL, 
                                  reduction = "cca", 
                                  anchor.features = int_features)

save(anchors, common_gene_names, file = "integrate_intermediate.RData")
load("integrate_intermediate.RData")
features2include = readRDS(file = "var_genes.rds")
int.genes = unique(union(anchors@anchor.features, features2include))

obj.integrated <- 
  IntegrateData(anchorset = anchors, dims = 1:30, 
                features.to.integrate	= int.genes # explicitly stating that we want the top svgs in the integrated dataset
  ) 

rm(obj.list, anchors); gc()
# obj.integrated.diet = DietSeurat(obj.integrated, assays = "integrated", layers = c("data"))

obj.integrated = ScaleData(obj.integrated, assay = "integrated", features = rownames(obj.integrated))
obj.integrated <- RunPCA(obj.integrated, verbose = FALSE, assay = "integrated")
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:30)
obj.integrated <- FindClusters(obj.integrated, resolution = 0.4)

obj.integrated@assays$integrated@scale.data = matrix(nrow=0,ncol=0)
obj.integrated@assays$SCT$scale.data = matrix(nrow=0,ncol=0)

saveRDS(obj.integrated, file = file.path("..", "analysis_output", "Alistair_Utah_st_integrate", "visium_alistair_utah_integrate_diet.rds"))


