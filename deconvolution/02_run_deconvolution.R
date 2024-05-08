library(Seurat)
library(dplyr)
library(magrittr)
library(readr)
library(spacexr)
library(Matrix)
library(CARD)

### This script performs de-convolution of the SRT samples using two packages: 
### 1) CARD 
### 2) RCTD

### Prerequisites:
# 1) output of all 01_run_deconvolution script 
# 2) output of all scripts under Utah_st_analysis and Alistair_st_analysis 
# 3) output of scVI integration (python script)

dir.create(file.path(".", 
                     "intermediate_output",
                     "CARD"))
dir.create(file.path(".", 
                     "intermediate_output",
                     "RCTD"))

# ----- load scRNA reference ----- #

count_mat_merge = readRDS(
  file = file.path(
    "intermediate_output",
    "count_mat.rds"
  )
)
meta_data_merge = readRDS(
  file = file.path(
    "intermediate_output",
    "ref_meta.rds"
  )
)

# ----- run RCTD deconvolution ----- #

cell_types= meta_data_merge$annotation
names(cell_types) = rownames(meta_data_merge)
reference <- Reference(count_mat_merge, cell_types)

sources = c("Alistair", "Utah")
sample_abbrs = list(
  "Alistair" = paste0("SP", 1:8),
  "Utah" = c("2B", "3C", "4D")
)
deconv.results.list = list()
deconv.props.list = list()

for (source in sources) {
  for (abbr in sample_abbrs[[source]]) {
    visium_obj = readRDS(file =
                           file.path(
                             "..",
                             "analysis_output",
                             paste0(source, "_st"),
                             paste0("visium_", abbr, ".rds")
                           ))
    
    st_counts = visium_obj@assays$SCT@counts
    st_coords = visium_obj@images[[1]]@coordinates %>%
      dplyr::mutate(x = imagerow, y = imagecol) %>%
      dplyr::select(x, y)
    st_coords = st_coords[colnames(st_counts),]
    colnames(st_coords) = c("X", "Y")
    sp_obj = SpatialRNA(st_coords, st_counts)
    rm(visium_obj)
    gc()
    
    myRCTD = create.RCTD(sp_obj, reference, max_cores = 8)
    myRCTD = run.RCTD(myRCTD, doublet_mode = 'full')
    
    deconv.results.list[[paste0(source, "_", abbr)]] <- myRCTD
    deconv.props.list[[paste0(source, "_", abbr)]] =
      normalize_weights(myRCTD@results$weights)
    
    rm(sp_obj, st_counts, st_coords, myRCTD)
    gc()
  }
}


saveRDS(
  deconv.results.list,
  file = file.path(
    "intermediate_output",
    "RCTD",
    "results.rds"
  )
)

saveRDS(
  deconv.props.list,
  file = file.path(
    "intermediate_output",
    "RCTD",
    "props.rds"
  )
)

# ----- run CARD deconvolution ----- #

sc_meta = meta_data_merge %>% 
  dplyr::mutate(cellID = row.names(.),
                cellType = annotation,
                sampleInfo = sample) %>%
  dplyr::select(cellID, cellType, sampleInfo)

sources = c("Alistair", "Utah")
sample_abbrs = list(
  "Alistair" = paste0("SP", 1:8),
  "Utah" = c("2B", "3C", "4D")
)

deconv.props.list = list()
for (source in sources) {
  for (sample_abbr in sample_abbrs[[source]]) {
    visium_obj = readRDS(file =
                           file.path(
                             "..",
                             "analysis_output",
                             paste0(source, "_st"),
                             paste0("visium_", sample_abbr, ".rds")
                           ))
    
    st_counts = visium_obj@assays$SCT@counts
    st_location = visium_obj@images[[1]]@coordinates %>%
      dplyr::mutate(x = imagerow, y = imagecol) %>%
      dplyr::select(x, y)
    st_location = st_location[colnames(st_counts),]
    rm(visium_obj)
    gc()
    
    CARD_obj = createCARDObject(
      sc_count = count_mat_merge,
      sc_meta = sc_meta,
      spatial_count = st_counts,
      spatial_location = st_location,
      ct.varname = "cellType",
      ct.select = unique(sc_meta$cellType),
      sample.varname = "sampleInfo",
      minCountGene = 100,
      minCountSpot = 10
    )
    
    CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
    deconv.props.list[[paste0(source, "_", sample_abbr)]] = CARD_obj@Proportion_CARD
    rm(CARD_obj, st_counts)
    gc()
  }
}

saveRDS(deconv.props.list,
        file = file.path("intermediate_output",
                         "CARD",
                         "props.rds"))