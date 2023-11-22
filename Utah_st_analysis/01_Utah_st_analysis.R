library(Seurat)
library(dplyr)
library(magrittr)

### Runs the analysis pipeline on Utah Visium samples. 
# Reads inputs from 10x visium ouput (h5 file, image, scale etc.) and saves intermediate object as rds file
# Preprocessing involves SCT normalization and regressing out read counts, percent.mt and percent.ribo genes per spot
# specifically, no QC based on read count or mt gene percentage is performed here
# Computes PCA but does no further clustering (as clusterings are to be generated from NMF of the integrated dataset) 
# Identifies spatially variable genes via the moran's I value 

# where visium data lives 
dir <- file.path("..", "data", "visium_JD24379", "Primary_data")

sample_folder_names <- c("1A-Bottom-Half-Tissue-set2_results", 
                         "2B-Top-Half-Tissue-Set2_results",
                         "3C-Bottom-Half-Tissue-Set1_results",
                         "4D-Top-Half-Tissue-Set1_results")

# abbreviation for each sample to avoid having long names 
sample_abbrs = c("1A", "2B", "3C", "4D")

for (i in 1:length(sample_abbrs)) {
  
  folder = sample_folder_names[i]
  abbr = sample_abbrs[i]
  
  visium_dir = file.path(dir, sample_folder_names[i])
  image_dir = file.path(visium_dir, "spatial")
  
  visium_img = Read10X_Image(image_dir,
                             image.name = "tissue_lowres_image.png",
                             filter.matrix = TRUE)
  
  visium_obj = Load10X_Spatial(
    visium_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "lowres",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = visium_img
  )
  
  visium_obj@meta.data$barcode = colnames(visium_obj)
  
  # percent ribosomal and mitochondrial gene, no QC removal of spots based on these metrics are performed 
  visium_obj <- PercentageFeatureSet(visium_obj, "^MT-", col.name = "percent.mt")
  visium_obj <- PercentageFeatureSet(visium_obj, "^RP[SL]", col.name = "percent.ribo")
  
  # normalization, scaling and PCA
  visium_obj <- SCTransform(visium_obj, assay = "Spatial", verbose = FALSE, vars.to.regress = c("percent.mt", "percent.ribo", "nCount_Spatial"))
  visium_obj <- ScaleData(visium_obj, assay = "SCT", features = rownames(visium_obj))
  visium_obj <- RunPCA(visium_obj, npcs = 30, assay = "SCT")
  
  # top 3000 hvgs identified, and non-mt and non-ribo hvgs are screened over to find svgs 
  visium_obj = FindVariableFeatures(object = visium_obj, selection.method = "vst", nfeatures= 3000)
  hvg = VariableFeatures(visium_obj)
  hvg = hvg[!grepl("^MT-", hvg) & !grepl("^RP[SL]", hvg)]
  visium_obj = FindSpatiallyVariableFeatures(visium_obj, assay = "SCT", features = hvg,
                                             selection.method = "moransi")
  
  saveRDS(visium_obj, file = file.path("..", "analysis_output", "Utah_st",paste0("visium_", abbr, ".rds")))
}