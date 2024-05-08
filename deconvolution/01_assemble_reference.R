library(Seurat)
library(dplyr)
library(magrittr)
library(readr)

### This script performs de-convolution of the SRT samples using two packages: 
### 1) CARD 
### 2) RCTD

### Prerequisites: 
# 1) raw cellranger output file (filtered_feature_bc_matrix) of Utah and Alistair samples under ../data/Alistair and ../data/Utah
# 2) output of scVI integration (python script)

dir.create(file.path(".", "intermediate_output"))

# ----- load scRNA integration results and assemble scRNA reference ----- #

# load single cell RNA objects and keep only the count matrix
obj.list = list()
dir = file.path("..", "data")
source = "Utah"
sample_abbrs = c("16030X2",
                 "16030X3",
                 "16030X4",
                 "19459X1",
                 "19595X1",
                 "19833X1",
                 "19833X2")
for (abbr in sample_abbrs) {
  matrix_dir = file.path(dir, source, abbr, "filtered_feature_bc_matrix")
  
  counts = Read10X(
    data.dir = matrix_dir,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  
  obj.list[[paste0(source, "-", abbr)]] = CreateSeuratObject(
    counts,
    project = paste0(source, "-", abbr),
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    meta.data = NULL
  )
  
  obj.list[[paste0(source, "-", abbr)]]@meta.data$source = source
  obj.list[[paste0(source, "-", abbr)]]@meta.data$sample = abbr
  obj.list[[paste0(source, "-", abbr)]]@meta.data$cell_id = paste0(abbr,
                                                                   "_",
                                                                   rownames(obj.list[[paste0(source, "-", abbr)]]@meta.data))
}

source = "Alistair"
sample_abbrs = c("Y2", "Y3", "Y5", "MJ10", "MJ11")
sample_fnames = paste0("GSM6506",
                       as.character(105:109),
                       "_counts_",
                       sample_abbrs,
                       ".txt.gz")
for (i in 1:length(sample_abbrs)) {
  sample_dir = sample_fnames[i]
  abbr = sample_abbrs[i]
  counts = read.table(file = gzfile(file.path(dir, source, sample_dir)), sep = "\t")
  obj.list[[paste0(source, "-", abbr)]] = CreateSeuratObject(
    counts = counts,
    project = paste0(source, "-", abbr),
    assay  = "RNA"
  )
  
  obj.list[[paste0(source, "-", abbr)]]@meta.data$source = source
  obj.list[[paste0(source, "-", abbr)]]@meta.data$sample = abbr
  obj.list[[paste0(source, "-", abbr)]]@meta.data$cell_id = paste0(abbr,
                                                                   "_",
                                                                   rownames(obj.list[[paste0(source, "-", abbr)]]@meta.data))
}


gene_names = sapply(obj.list, FUN = function(x) rownames(x@assays$RNA$counts))
common_gene_names = Reduce(intersect, gene_names)

count_mat_merge = NULL
meta_data_merge = data.frame()
for (i in 1:length(obj.list)) {
  count_mat = obj.list[[i]]@assays$RNA$counts[common_gene_names,]
  barcodes = colnames(count_mat)
  barcodes_rename = obj.list[[i]]@meta.data$cell_id
  colnames(count_mat) = barcodes_rename
  count_mat_merge = cbind(count_mat_merge, count_mat)
  
  meta_data = data.frame(sample = obj.list[[i]]@meta.data$sample, 
                         barcode = barcodes, 
                         row.names = barcodes_rename)
  meta_data_merge = rbind(meta_data_merge, meta_data)
}
colnames(count_mat_merge) = gsub("-", "_", colnames(count_mat_merge))
rownames(meta_data_merge) = gsub("-", "_",rownames(meta_data_merge))

# ----- load scVI integration and annotation results, generate metadata with annotation corresponding to the merged count matrix ----- #

annotation_dict = read_csv(file.path("..", "data", "scVI_integration_output", "annotation_dict.csv"),
                           col_names = F)
annotate_metadata = read_csv(file.path("..", "data", "scVI_integration_output", "annotate_metadata.csv"))
annotate_metadata$cell_id = paste0(annotate_metadata$sample, "_", annotate_metadata$barcode)
annotate_metadata$cell_id = gsub("-", "_", annotate_metadata$cell_id)
annotate_metadata$marker_annotation = factor(annotate_metadata$marker_annotation, 
                                             levels = unique(annotation_dict$X2))

mask = rownames(meta_data_merge) %in% annotate_metadata$cell_id
count_mat_merge = count_mat_merge[,mask]
meta_data_merge = meta_data_merge[mask,]
match_idx = match(rownames(meta_data_merge), annotate_metadata$cell_id)
meta_data_merge$annotation = annotate_metadata$marker_annotation[match_idx]

saveRDS(
  count_mat_merge,
  file = file.path(
    "intermediate_output",
    "count_mat.rds"
  )
)
saveRDS(
  meta_data_merge,
  file = file.path(
    "intermediate_output",
    "ref_meta.rds"
  )
)