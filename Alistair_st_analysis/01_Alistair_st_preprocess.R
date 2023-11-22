library(Matrix)
library(readr)
library(magrittr)
library(dplyr)
library(Seurat)

dir <- file.path("..", "data", "Alistair")
sample_folders = as.character(6506110:6506117)
sample_abbrs = paste0("SP", 1:8)


for (i in 1:length(sample_folders)) {
  
  sample = sample_folders[i]
  abbr = sample_abbrs[i]
  
  counts <- readMM(file.path(dir, sample, "filtered_feature_bc_matrix", "matrix.mtx.gz"))
  
  genes <- as.data.frame(read_tsv(file.path(dir, sample, "filtered_feature_bc_matrix", "features.tsv.gz"), col_names = FALSE))
  colnames(genes) = c("ENSembl", "symbol", "type")
  
  cell_ids <- as.data.frame(read_tsv(file.path(dir, sample, "filtered_feature_bc_matrix", "barcodes.tsv.gz"), col_names = FALSE))
  colnames(cell_ids) = c("cell_ids")
  
  counts = counts[!duplicated(genes$symbol), ]
  genes = genes[!duplicated(genes$symbol), ]
  
  colnames(counts) = cell_ids$cell_ids
  rownames(counts) = genes$symbol
  visium_obj = CreateSeuratObject(counts = counts, assay="Spatial")
  
  visium_img = Read10X_Image(
    file.path(dir, sample, "spatial"),
    image.name = "tissue_lowres_image.png",
    filter.matrix = TRUE
  )
  
  rm(counts, genes, cell_ids); gc()
  
  visium_obj@images$image = visium_img
  visium_obj@images$image@key = "lowres"
  
  visium_obj$barcode = colnames(visium_obj)
  visium_obj$sample = abbr
  
  # percent ribosomal and mitochondrial gene, no QC removal of spots based on these metrics are performed  
  visium_obj <- PercentageFeatureSet(visium_obj, "^MT-", col.name = "percent.mt")
  visium_obj <- PercentageFeatureSet(visium_obj, "^RP[SL]", col.name = "percent.ribo")
  
  # normalization, scaling and PCA
  visium_obj <- SCTransform(visium_obj, assay = "Spatial", verbose = FALSE, vars.to.regress = c("percent.mt", "percent.ribo", "nCount_Spatial"))
  visium_obj <- RunPCA(visium_obj, npcs = 30, assay = "SCT")
  
  # top 3000 hvgs identified, and non-mt and non-ribo hvgs are screened over to find svgs 
  visium_obj = FindVariableFeatures(object = visium_obj, selection.method = "vst", nfeatures= 3000)
  hvg = VariableFeatures(visium_obj)
  hvg = hvg[!grepl("^MT-", hvg) & !grepl("^RP[SL]", hvg)]
  visium_obj = FindSpatiallyVariableFeatures(visium_obj, assay = "SCT", features = hvg,
                                             selection.method = "moransi")
  saveRDS(visium_obj, file = file.path("..", "analysis_output", "Alistair_st", paste0("visium_", sample_abbrs[i],".rds")))
  
}