library(Seurat)
library(magrittr)
library(dplyr)
library(ccfindR)
library(nebula)

### This script performs NMF cluster based differential gene expression analysis using nebula

### Prerequisites:  
## 1) output from 01_Utah_st_analysis.R (preprocessing)
## 2) output from 02_Utah_st_NMF.R outputs(NMF clustering)

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
obj.integrated$sample = sapply(strsplit(rownames(obj.integrated@meta.data), "_"), "[[", 1)
obj.integrated$subject = "PF10E4D"
match_ordering = match(rownames(obj.integrated@meta.data), names(cluster_ids))
any(is.na(match_ordering))
obj.integrated@meta.data$NMF_cluster = cluster_ids[match_ordering]
obj.integrated$NMF_cluster = factor(as.character(obj.integrated$NMF_cluster))

# ----- nebula DE ----- # 
neb_seurat <- scToNeb(obj = obj.integrated, 
                      assay = "SCT", 
                      id = "sample", 
                      pred = c("NMF_cluster"), 
                      offset="nCount_Spatial")
neb_seurat$count = neb_seurat$count[rownames(obj.integrated@assays$integrated@data),]

pred_clusts = c(1, 2, 3, 6, 7, 8, 9, 10)
des_mat = sapply(pred_clusts, FUN = function(x) ifelse(neb_seurat$pred$NMF_cluster == as.character(x), 1, 0))
des_mat = cbind(1, des_mat)
colnames(des_mat) = c("(Intercept)", paste0("NMF_cluster", pred_clusts))
rownames(des_mat) = rownames(neb_seurat$pred)

re.utah.nmf = nebula(neb_seurat$count,
                     neb_seurat$id,
                     pred = des_mat,
                     offset = neb_seurat$offset)

cols_lfc = re.utah.nmf$summary %>% 
  dplyr::select(starts_with("logFC_")) %>% 
  colnames()

cols_pval = re.utah.nmf$summary %>% 
  dplyr::select(starts_with("p_")) %>% 
  colnames()

predictors = sapply(strsplit(cols_lfc, "_"), FUN = function(x) x[[length(x)]])

row_lfc_summary_df = 
  data.frame(t(apply(re.utah.nmf$summary[cols_lfc],1,FUN = function(x) order(x, decreasing = T)[1:2], simplify = T)),
             t(apply(re.utah.nmf$summary[cols_lfc],1,FUN = function(x) x[order(x, decreasing = T)][1:2], simplify = T))
  )
colnames(row_lfc_summary_df) = c("max_cluster", "sec_max_cluster", "lfc_max_cluster", "lfc_sec_max_cluster")

row_pval_summary_df = NULL 
for (i in 1:dim(row_lfc_summary_df)[1]) {
  row_pval_summary_df = rbind(row_pval_summary_df,
                              c(re.utah.nmf$summary[i, cols_pval[row_lfc_summary_df[i, 1]]], 
                                re.utah.nmf$summary[i, cols_pval[row_lfc_summary_df[i, 2]]]))
}
colnames(row_pval_summary_df) = c("pval_max_cluster", "pval_sec_max_cluster")

row_lfc_summary_df$max_cluster = predictors[row_lfc_summary_df$max_cluster]
row_lfc_summary_df$sec_max_cluster = predictors[row_lfc_summary_df$sec_max_cluster]

nmf_de_genes_df = cbind(re.utah.nmf$summary[c("gene_id", "gene")], 
                        re.utah.nmf$overdispersion,
                        row_lfc_summary_df, 
                        row_pval_summary_df
) %>% 
  dplyr::filter(pval_max_cluster <= 0.05,
                !grepl("^MT-", gene)) %>% 
  group_by(max_cluster) %>%
  dplyr::arrange(max_cluster, pval_max_cluster) 

nebula_features = split(nmf_de_genes_df$gene, nmf_de_genes_df$max_cluster)
nebula_top_features = sapply(nebula_features, FUN =function(x) x[1:min(15, length(x))], simplify = F)

saveRDS(nmf_de_genes_df, file = file.path("intermediate_output", "DE", "nebula_de_df.rds"))
saveRDS(re.utah.nmf, file = file.path("intermediate_output", "DE", "nebula_utah_nmf.rds"))
