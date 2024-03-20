library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ccfindR)

### This script performs NMF cluster based differential gene expression analysis on integrated Utah st datasets 

### Prerequisites:  ** 01_Utah_st_analysis.R outputs(preprocessing output)
###                 ** 02_Utah_st_NMF.R outputs(NMF clustering)

dir.create(file.path("intermediate_output", "DE"))

# ----- load NMF clustering info ----- # 
vb_nmf <- readRDS(file.path("intermediate_output", "vb_nmf_old.rds"))
# extract NMF cluster id
cluster_ids = cluster_id(vb_nmf, rank=optimal_rank(vb_nmf)$ropt)
barcodes = sapply(strsplit(colnames(vb_nmf), "_"), "[[", 2) # spot names from integration step are of the format [sample id]_[raw barcode]
sample = sapply(strsplit(colnames(vb_nmf), "_"), "[[", 1) # get the sample id per barcode 

# ----- load integrated object ----- # 
visium_dir = file.path("..", "analysis_output", "Utah_st")
obj.integrated = readRDS(file = file.path(visium_dir, "visium_utah_integrate.rds"))

match_idx = match(rownames(obj.integrated@meta.data), names(cluster_ids))
any(is.na(match_idx))

obj.integrated$NMF_cluster = cluster_ids[match_idx]
obj.integrated$NMF_cluster = factor(obj.integrated$NMF_cluster)

# ----- DE ----- # 
de_df_list_raw = list()
for (c in levels(obj.integrated$NMF_cluster)) {
  de_markers <-
    FindMarkers(
      obj.integrated,
      ident.1 = c,
      # program cluster id
      group.by = "NMF_cluster",
      logfc.threshold = 0
    )
  
  de_df_list_raw[[c]] = de_markers
}

saveRDS(de_df_list_raw, file = file.path("intermediate_output", "DE", "integrate_NMF_de_df_list.rds"))

filter_de_df <- function(df) {
  df %>%
    dplyr::filter(p_val_adj <= 0.05,
                  pct.1 >= 0.75,
                  avg_log2FC > 0.5,
                  !grepl("^MT-", rownames(.))) %>%
    dplyr::arrange(p_val_adj, desc(avg_log2FC))
}

de_df_list_filter = sapply(de_df_list_raw, FUN = filter_de_df, simplify = F)

saveRDS(de_df_list_filter, file = file.path("intermediate_output", "DE", "integrate_NMF_de_df_list_filter.rds"))

