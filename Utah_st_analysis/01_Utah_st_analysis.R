library(Seurat)
library(dplyr)
library(magrittr)

### This script loads and preprocesses the Utah Visium samples, and then 
# performs seurat anchor based integration 
# anchor with 3000 features chosen by SelectIntegrationFeatures

### Prerequisites:  ** none

dir.create(file.path("..", "analysis_output", "Utah_st"))

# ----- visium spaceranger output dir ----- # 
dir <- file.path("..", "data", "visium_JD24379", "Primary_data")

sample_folder_names <- c("1A-Bottom-Half-Tissue-set2_results", 
                         "2B-Top-Half-Tissue-Set2_results",
                         "3C-Bottom-Half-Tissue-Set1_results",
                         "4D-Top-Half-Tissue-Set1_results")

# abbreviation for each sample to avoid having long names 
sample_abbrs = c("1A", "2B", "3C", "4D")

# ----- load and preprocessing ----- # 
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
  
  # normalization done with SCTransform, regressing out mito, ribo conc
  visium_obj <- SCTransform(visium_obj, assay = "Spatial", verbose = FALSE, vars.to.regress = c("percent.mt", "percent.ribo", "nCount_Spatial"))
  visium_obj <- ScaleData(visium_obj, assay = "SCT", features = rownames(visium_obj))
  visium_obj <- RunPCA(visium_obj, npcs = 30, assay = "SCT")
  
  # top 3000 hvgs identified, and non-mt and non-ribo hvgs are screened over to find svgs 
  visium_obj = FindVariableFeatures(object = visium_obj, selection.method = "vst", nfeatures= 3000)
  hvg = VariableFeatures(visium_obj)
  hvg = hvg[!grepl("^MT-", hvg) & !grepl("^RP[SL]", hvg)]
  visium_obj = FindSpatiallyVariableFeatures(visium_obj, assay = "SCT", features = hvg,
                                             selection.method = "moransi")
  
  saveRDS(visium_obj, file = file.path("..", "analysis_output", "Utah_st", paste0("visium_", abbr, ".rds")))
}

# extract top 1000 spatially variable features (pre-computed) per sample and take the union 
sv_features = NULL

for (i in 1:length(obj.list)) {
  svgs = obj.list[[i]]@assays$SCT@meta.features %>% 
    filter(moransi.spatially.variable) %>% 
    arrange(moransi.spatially.variable.rank)
  sv_features = union(sv_features, rownames(svgs)[1:1000])
}

int_features <- SelectIntegrationFeatures(object.list = obj.list, assay = rep("SCT", length(obj.list)), nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = obj.list, scale = FALSE, 
                                  reference = NULL, reduction = "rpca", 
                                  anchor.features = int_features, assay = rep("SCT", length(obj.list)))
obj.integrated <- 
  IntegrateData(anchorset = anchors, dims = 1:30, 
                features.to.integrate	= union(int_features, sv_features) # explicitly stating that we want the top svgs in the integrated dataset
  ) 

obj.integrated = ScaleData(obj.integrated, assay = "integrated", features = rownames(obj.integrated))
obj.integrated <- RunPCA(obj.integrated, verbose = FALSE, assay = "integrated")
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:30)
obj.integrated <- FindClusters(obj.integrated, resolution = 0.4)

saveRDS(obj.integrated, file = file.path("..", "analysis_output", "Utah_st", "visium_utah_integrate.rds"))
saveRDS(sv_features, file = file.path("..", "analysis_output", "Utah_st", "integrate_svg.rds")) # also want to save the union of spatially variable features for use later in NMF