### Generates reference count matrix and metadata needed by CARD pipeline 

### Prerequisites:  ** scVI celltype annotation 
###                 ** Alistair reference annotation 
###                 ** post-QC scRNA seurat objects   

library(Seurat)
library(magrittr)
library(dplyr)

# scVI export annotation results, with cell barcode column, sample identifier column and celltype columns
reference_annotation = read.csv(file.path("..", "data", "annotate_meta_scvi.csv"), header=TRUE)
reference_annotation$cell_id = gsub("-", "_", reference_annotation$cell_id)

# load all scRNA seurat objects, where SCT counts are to be extracted from 
sample_abbrs = c("Y2", "Y3", "Y5", "MJ10", "MJ11", "16030X2", "16030X3", "16030X4")
obj.list = list()
for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  obj.list[[sample_abbr]] = readRDS(file = file.path(
    "..",
    "analysis_output",
    paste0("sobj_annotate_", sample_abbr, ".rds")
  ))
  obj.list[[sample_abbr]]$sample = sample_abbr
}

gene_names = sapply(
  obj.list,
  FUN = function(x)
    rownames(x@assays$SCT@counts)
)
common_gene_names = Reduce(intersect, gene_names)

# merging count matrix form all scRNA samples
count_mat_merge = NULL
for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  count_mat_merge = cbind(count_mat_merge, obj.list[[sample_abbr]]@assays$SCT@counts[common_gene_names, ])
}
colnames(count_mat_merge) = gsub("-", "_", colnames(count_mat_merge))

rm(obj.list)
gc()

# matching up referene metadata with merged count matrix 
match_idx = match(colnames(count_mat_merge), reference_annotation$cell_id)

ref_meta = reference_annotation[match_idx,] %>% 
  dplyr::mutate(cellID = cell_id, 
                cellType = marker_annotation, 
                sampleInfo = orig.ident) %>% 
  dplyr::select(cellID, cellType, sampleInfo)
head(ref_meta)

# save reference count matrix and metadata
saveRDS(
  count_mat_merge,
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_reference_count_mat_full.rds"
  )
)
saveRDS(
  ref_meta,
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_scvi_reference_meta_full.rds"
  )
)

# subset matrix and metadata to generate reference data only with alistair scRNA samples
count_mat_merge_alistair = count_mat_merge[,!grepl("^16030", colnames(count_mat_merge))]
match_idx = match(colnames(count_mat_merge_alistair), ref_meta$cellID)

saveRDS(
  count_mat_merge_alistair,
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_reference_count_mat_alistair.rds"
  )
)
saveRDS(
  ref_meta_alistair,
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_scvi_reference_meta_alistair.rds"
  )
)

reference_annotation = read_csv(file.path("data", "Alistair_reference_annotation.csv"))
reference_annotation$cell_rename_id = paste0(reference_annotation[["Sample"]], "_", gsub("-", "_", reference_annotation[["Cell ID"]]))
barcodes_intersect = intersect(colnames(count_mat_merge), reference_annotation$cell_rename_id)
count_mat_merge_alistair_ref = count_mat_merge[,barcodes_intersect]

ref_meta_alistair_ref = 
  reference_annotation[match(colnames(count_mat_merge_alistair_ref), reference_annotation$cell_rename_id),] %>% 
  dplyr::mutate(cellID = cell_rename_id, 
                cellType = `Fine-grain annotations`, 
                sampleInfo = Sample) %>% 
  dplyr::select(cellID, cellType, sampleInfo)

saveRDS(
  count_mat_merge_alistair_ref,
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_reference_count_mat_alistair_ref.rds"
  )
)
saveRDS(
  ref_meta_alistair_ref,
  file = file.path(
    "intermediate_output",
    "CARD_reference",
    "CARD_alistair_reference_meta.rds"
  )
)
