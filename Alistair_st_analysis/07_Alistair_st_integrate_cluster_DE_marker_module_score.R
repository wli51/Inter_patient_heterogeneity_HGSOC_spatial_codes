library(Seurat)
library(magrittr)
library(dplyr)

dir = file.path("..", "analysis_output", "Alistair_st", )

# reads the DE genesets per spot on Alistair samples
genesets = readRDS(
  file = file.path(
    "intermediate_output",
    "integrate_cluster_DE_manual.rds"
  )
)

sample_abbrs = paste0("SP", 1:8)

obj.list = list()
for (i in 1:length(sample_abbrs)) {
  obj.list[[i]] = readRDS(file.path(dir, paste0("visium_", sample_abbrs[i], ".rds")))
}

for (i in 1:length(sample_abbrs)) {
  obj.list[[i]]@meta.data = obj.list[[i]]@meta.data %>% 
    dplyr::select(!starts_with("moudle."))
  
  obj.list[[i]] = AddModuleScore(
    obj.list[[i]],
    features = genesets,
    pool = NULL,
    nbin = 24,
    ctrl = 100,
    k = FALSE,
    assay = "SCT",
    name = "module.",
    seed = 1,
    search = FALSE
  )
  
  module_scores = obj.list[[i]]@meta.data %>% select(starts_with("module."))
  colnames(module_scores) = names(genesets) 
  saveRDS(module_scores, file = file.path("intermediate_output", paste0("integrate_DE_module_score_",sample_abbrs[i], ".rds")))
}

