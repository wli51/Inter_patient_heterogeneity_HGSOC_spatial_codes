library(Seurat)
library(dplyr)
library(magrittr)
library(R.utils)
library(readr)

### Preprocesses single cell RNA datasets 

### Prerequisites:  cellranger output  

dir <- file.path("..", "data", "Utah")
sample_ids <- c("16030X2", "16030X3", "16030X4")
sample_abbrs = sample_ids

for (i in 2:length(sample_ids)) {
  
  matrix_dir = file.path(dir, sample_ids[i], "filtered_feature_bc_matrix")
  
  counts = Read10X(
    data.dir = matrix_dir,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  
  sobj = CreateSeuratObject(
    counts,
    project = sample_abbrs[i],
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    meta.data = NULL
  )
  
  # some quality control 
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^RP[LS]")
  
  # cell cycle scoring 
  sobj <- subset(sobj, subset = nFeature_RNA > 500 & percent.mt <= 15)
  sobj <- NormalizeData(sobj)
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  ### SCT Normalization regressing out out cellcyle score and mt.gene percentage
  sobj <- SCTransform(sobj, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "percent.ribo", "nCount_RNA"), verbose = FALSE)
  
  sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj), assay = "SCT")
  sobj <- FindNeighbors(sobj, dims = 1:30)
  sobj <- RunUMAP(sobj, dims = 1:30, assay  = "SCT")
  sobj <- FindClusters(sobj, resolution = 0.4)
  
  saveRDS(sobj, file = file.path("..", "analysis_output", "Utah_sc", paste0("sobj_", sample_abbrs[i], ".rds")))
  rm(sobj)
  gc()
}