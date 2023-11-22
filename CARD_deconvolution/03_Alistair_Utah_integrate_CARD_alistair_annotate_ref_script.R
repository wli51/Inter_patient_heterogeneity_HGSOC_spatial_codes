library(CARD)
library(magrittr)
library(dplyr)

sc_counts = readRDS(
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_reference_count_mat_alistair_ref.rds"
  )
)
sc_meta = readRDS(file = file.path(
  "intermediate_output",
  "CARD_reference",
  "CARD_alistair_reference_meta.rds"
))
sc_meta = as.data.frame(sc_meta)
rownames(sc_meta) = sc_meta$cellID

sample_abbrs = c(c("2B", "3C", "4D"), paste0("SP", 1:8))

for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  visium_obj = readRDS(file =
                         file.path(
                           "..",
                           "analysis_output",
                           paste0("visium_", sample_abbr, ".rds")
                         ))
  
  st_counts = visium_obj@assays$SCT@counts
  
  st_location = visium_obj@images[[1]]@coordinates %>%
    dplyr::mutate(x = imagerow, y = imagecol) %>%
    dplyr::select(x, y)
  
  st_location = st_location[colnames(st_counts), ]
  
  rm(visium_obj)
  gc()
  
  CARD_obj = createCARDObject(
    sc_count = sc_counts,
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
  saveRDS(
    CARD_obj@Proportion_CARD,
    file = file.path(
      "intermediate_output",
      "CARD_samplewise_output",
      paste0("CARD_alistair_ref_output_", sample_abbr, ".rds")
    )
  )
  
  rm(CARD_obj, st_counts)
  gc()
}

