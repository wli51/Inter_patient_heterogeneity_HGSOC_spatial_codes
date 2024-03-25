library(Seurat)
library(dplyr)
library(magrittr)
library(R.utils)
library(readr)

### Preprocesses single cell RNA datasets 

### Prerequisites:  cellranger output  

dir = file.path("..", "data", "Alistair")
sample_abbrs = c("Y2", "Y3", "Y5", "MJ10", "MJ11")
sample_ids = paste0("GSM6506", as.character(105:109), "_counts_", sample_abbrs, ".txt.gz")

# for (i in 1:length(sample_ids)) {
  
  sample = sample_ids[i]
  abbr = sample_abbrs[i]
  counts <- read.table(file = gzfile(file.path(dir, sample)), sep = "\t")
  sobj = CreateSeuratObject(counts = counts, project = abbr, assay  = "RNA")
  
  # some quality control 
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^RP[LS]")
  
  # cell cycle scoring 
  sobj <- subset(sobj, subset = nFeature_RNA > 500 & percent.mt <= 15)
  sobj <- NormalizeData(sobj)
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  ### SCT Normalization regressing out out cellcyle score and mt.gene percentagei 
  sobj <- SCTransform(sobj, 
                      vars.to.regress = c("S.Score", "G2M.Score", "percent.mt", "percent.ribo", "nCount_RNA"), 
                      verbose = FALSE,
                      min_cells = 0
                      )
  
  sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj), assay = "SCT")
  sobj <- FindNeighbors(sobj, dims = 1:30)
  sobj <- RunUMAP(sobj, dims = 1:30, assay  = "SCT")
  sobj <- FindClusters(sobj, resolution = 0.4)
  
  saveRDS(sobj, file = file.path("..", "analysis_output", "Alistair_sc", paste0("sobj_", sample_abbrs[i], ".rds")))
  rm(counts, sobj);gc()
  
# }