library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(nebula)
library(corrplot)
library(caret)

### this script performs cluster based DE expression on integrated visium dataset with nebula

### Prerequisites:  02_Alistair_st_integrate.R output

obj.integrated = readRDS(file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_diet.rds"))

CARD_summary_df = readRDS(file = file.path("..", "analysis_output", "CARD", "CARD_deconv_summary_scvi_reference.rds"))
match_idx = match(rownames(obj.integrated@meta.data), CARD_summary_df$cell_id)
CARD_summary_df = CARD_summary_df[match_idx,] %>% dplyr::select(-c("source", "sample", "cell_id"))
obj.integrated@meta.data = cbind(obj.integrated@meta.data, CARD_summary_df)
celltypes = colnames(CARD_summary_df)
metadata = obj.integrated@meta.data

sample_abbrs = paste0("SP", 1:8)

# Load seurat objects 
obj.list = list()
for (sample_abbr in sample_abbrs) {
  obj.list[[sample_abbr]] = readRDS(file.path("..", "analysis_output", "Alistair_st", paste0("visium_", sample_abbr, ".rds")))
}
common_gene_names = Reduce(intersect, sapply(obj.list, FUN = function(x) rownames(x@assays$SCT$counts), simplify = F))

count_mat_merge = NULL
cell_ids = NULL
for (sample_abbr in sample_abbrs) {
  cmat = obj.list[[sample_abbr]]@assays$SCT$counts[common_gene_names,]
  barcodes = colnames(cmat)
  count_mat_merge = cbind(count_mat_merge, cmat)
  cell_ids = c(cell_ids, paste0(sample_abbr, "_", barcodes))
}
colnames(count_mat_merge) = cell_ids
match_idx = match(colnames(obj.integrated), colnames(count_mat_merge))
sct_count_assay = CreateAssay5Object(counts = count_mat_merge[,match_idx])
obj.integrated[["SCT"]] = sct_count_assay

# ----- Nebula object setup ----- # 
neb_seurat <- scToNeb(obj = obj.integrated, 
                      assay = "SCT", 
                      id = "sample", 
                      pred = c("seurat_clusters", celltypes), 
                      offset="nCount_SCT")
neb_seurat$count = neb_seurat$count[rownames(obj.integrated@assays$integrated@data),]

# Nebula design matrix construction
celltypes_format = gsub("[\\(\\);_ ]", "\\.", celltypes) # format celltype names to comply with perper dataframe column naming
ref_cluster = "1"
des_mat_cluster = model.matrix( ~ seurat_clusters + 0, data = neb_seurat$pred)
des_mat_cluster = des_mat_cluster[, colnames(des_mat_cluster)[colnames(des_mat_cluster) != paste0("seurat_clusters", ref_cluster)]]
des_mat_prop = model.matrix(formula(paste0("~ ", paste(celltypes_format, collapse = " + "), " + 0")), data = neb_seurat$pred)
des_mat = cbind(`(Intercept)` = 1, des_mat_cluster, des_mat_prop)
df = des_mat[,2:dim(des_mat)[2]]
linear_dep_idx2remove = findCorrelation(cor(df), cutoff = 0.7) # remove celltype proportions that are colinear with other predictors
des_mat = des_mat[,-(linear_dep_idx2remove+1)]

re.alistair.cluster.prop = nebula(neb_seurat$count,
                                  neb_seurat$id,
                                  pred = des_mat,
                                  offset = neb_seurat$offset)

if (!dir.exists(file.path("intermediate_output", "nebula"))) {
  dir.create(file.path("intermediate_output", "nebula"))
}

saveRDS(re.alistair.cluster.prop, file = file.path("intermediate_output", "nebula", "nebula_alistair_integrate_cluster_prop_obj.rds"))
des_mat_scale_prop = des_mat
des_mat_scale_prop[, colnames(des_mat_scale_prop)[colnames(des_mat_scale_prop) %in% celltypes_format]] = (scale(des_mat_scale_prop[, colnames(des_mat_scale_prop)[colnames(des_mat_scale_prop) %in% celltypes_format]]))
re.alistair.cluster.prop.scale = nebula(neb_seurat$count,
                                        neb_seurat$id,
                                        pred = des_mat_scale_prop,
                                        offset = neb_seurat$offset)
saveRDS(re.alistair.cluster.prop.scale, file = file.path("intermediate_output", "nebula","nebula_alistair_integrate_cluster_prop_scale_obj.rds"))

# ----- find DEGs significantly associated with cluster ----- # 
de_cluster = re.alistair.cluster.prop$summary %>% 
  dplyr::select(colnames(.)[grepl("seurat_clusters", colnames(.))], 
                colnames(.)[grepl("Intercept", colnames(.))],
                "gene", "gene_id")
clusters = gsub("^logFC_seurat_clusters", "", colnames(de_cluster)[grepl("logFC_seurat_clusters", colnames(de_cluster))])
de_cluster_pval = de_cluster %>% dplyr::select(starts_with("p_"))
de_cluster_pval_adj = as.data.frame(apply(
  de_cluster_pval,
  MARGIN = 2,
  FUN = function(x)
    p.adjust(x, p.adjust.methods[7])
))

colnames(de_cluster_pval_adj) = gsub("^p_", "", colnames(de_cluster_pval_adj))
de_cluster_lfc = de_cluster %>% dplyr::select(starts_with("logFC_"))
colnames(de_cluster_lfc) = gsub("^logFC_", "", colnames(de_cluster_lfc))
de_cluster_lfc = de_cluster_lfc[,colnames(de_cluster_pval_adj)]

de_cluster_se = de_cluster %>% dplyr::select(starts_with("se_"))
colnames(de_cluster_se) = gsub("^se_", "", colnames(de_cluster_se))
de_cluster_se = de_cluster_se[,colnames(de_cluster_pval_adj)]

sig_pos_table = matrix(FALSE, nrow = nrow(de_cluster), ncol = length(clusters))
rownames(sig_pos_table) = de_cluster$gene
colnames(sig_pos_table) = clusters
lfc_table = matrix(NaN, nrow = nrow(de_cluster), ncol = length(clusters))
rownames(lfc_table) = de_cluster$gene
colnames(lfc_table) = clusters
pval_table = lfc_table
ci_lower_table = lfc_table
ci_upper_table = lfc_table

# iterate over all nebula genes, for each gene, compare the nebula DE clusters 
# coefficients that are indicative of significant positive log fold changes,
# if a gene is positively up-regulated for more than 1 cluster, and the cluster
# coefficient with the highest logFC has non-overlapping CI with others, then 
# the gene is deemed DEG for that highest logFC cluster, otherwise, the gene is 
# considered to be up-regualted in all clusters with positive, significant 
# nebula coefficinet
for (i in 1:nrow(de_cluster)) {
  de_cluster_row = de_cluster[i,]
  row_gene = de_cluster_row[, "gene"]
  
  pvals = de_cluster_pval_adj[i,] %>% dplyr::select(-c("(Intercept)"))
  ses = de_cluster_se[i, ] %>% dplyr::select(-c("(Intercept)"))
  lfcs = de_cluster_lfc[i,] %>% dplyr::select(-c("(Intercept)"))
  cis = rbind(lfcs - 2 * ses,
              lfcs + 2 * ses)
  
  pval_table[row_gene,] = as.numeric(pvals)
  lfc_table[row_gene,] = as.numeric(lfcs)
  ci_lower_table[row_gene,] = as.numeric(cis[1,])
  ci_upper_table[row_gene,] = as.numeric(cis[2,])
  
  is_sig = pvals <= 0.05
  is_pos = cis[1,] > 0 
  
  pos_sig_mask = is_sig & is_pos
  pos_sig_clusters = gsub("seurat_clusters", "", names(pvals)[pos_sig_mask])
  
  lfcs = lfcs[pos_sig_mask]
  cis = cis[, pos_sig_mask]
  
  if (length(pos_sig_clusters) == 1) {
    sig_pos_table[row_gene, pos_sig_clusters] = TRUE
    next
    
  } else if (length(pos_sig_clusters) == 0) {
    next
  }
  
  non_overlapping_ci = sapply(
    cis[1, ],
    FUN = function(x)
      sum(x > cis[2, ]) > 0,
    simplify = T
  )
  
  if (any(non_overlapping_ci)) {
    pos_sig_clusters = pos_sig_clusters[non_overlapping_ci]
  } 
  
  sig_pos_table[row_gene, pos_sig_clusters] = TRUE
}
# filter for nebula DEGs only found to be significantly up-regulated in exactly 
# one cluster
num_pos_sig_cluster = apply(sig_pos_table, MARGIN = 1, FUN = sum, simplify = T)
gene_one_pos_sig = names(num_pos_sig_cluster)[num_pos_sig_cluster == 1]
sig_pos_cluster_index = apply(sig_pos_table[gene_one_pos_sig,], MARGIN = 1, FUN = function(x) which(x), simplify = T)
sig_pos_cluster = colnames(sig_pos_table)[sig_pos_cluster_index]
nebula_cluster_degs = split(names(sig_pos_cluster_index), sig_pos_cluster)
nebula_cluster_degs_filter = list()
for (c in names(nebula_cluster_degs)) {
  degs = nebula_cluster_degs[[c]]
  lfcs = NULL
  pvals = NULL
  for (g in degs) {
    lfcs = c(lfcs, lfc_table[g, c])
    pvals = c(pvals, pval_table[g, c])
  }
  names(lfcs) = degs
  lfcs = lfcs[lfcs >= 0.25]
  names(pvals) = degs
  pvals = pvals[order(pvals, decreasing = F)]
  
  nebula_cluster_degs_filter[[c]] = lfcs[names(pvals)[names(pvals) %in% names(lfcs)]]
}
saveRDS(nebula_cluster_degs_filter, file =file.path("intermediate_output", "nebula", "nebula_alistair_integrate_top_cluster_DEGs.rds") )


# ----- find DEGs significantly associated with deconvolution proportion ----- # 
de_ct = re.alistair.cluster.prop$summary %>% 
  dplyr::select(colnames(.)[!grepl("seurat_clusters", colnames(.))] & 
                  colnames(.)[!grepl("Intercept", colnames(.))],
                "gene", "gene_id")
de_ct_pval = de_ct %>% dplyr::select(starts_with("p_"))
de_ct_pval_adj = as.data.frame(apply(
  de_ct_pval,
  MARGIN = 2,
  FUN = function(x)
    p.adjust(x, "fdr")
))
colnames(de_ct_pval_adj) = gsub("^p_", "", colnames(de_ct_pval_adj))
cts = colnames(de_ct_pval_adj)
de_ct_lfc = de_ct %>% dplyr::select(starts_with("logFC_"))
colnames(de_ct_lfc) = gsub("^logFC_", "", colnames(de_ct_lfc))
de_ct_lfc = de_ct_lfc[,cts]
de_ct_se = de_ct %>% dplyr::select(starts_with("se_"))
colnames(de_ct_se) = gsub("^se_", "", colnames(de_ct_se))
de_ct_se = de_ct_se[,cts]

sig_pos_table_ct = matrix(FALSE, nrow = nrow(de_ct), ncol = length(cts))
rownames(sig_pos_table_ct) = de_ct$gene
colnames(sig_pos_table_ct) = cts
lfc_table_ct = matrix(NaN, nrow = nrow(de_ct), ncol = length(cts))
rownames(lfc_table_ct) = de_ct$gene
colnames(lfc_table_ct) = cts
pval_table_ct = lfc_table_ct
ci_lower_table_ct = lfc_table_ct
ci_upper_table_ct = lfc_table_ct

for (i in 1:nrow(de_ct)) {
  de_ct_row = de_ct[i,]
  row_gene = de_ct_row[, "gene"]
  
  pvals = de_ct_pval_adj[i,]
  ses = de_ct_se[i,]
  lfcs = de_ct_lfc[i,]
  cis = rbind(lfcs - 2 * ses,
              lfcs + 2 * ses)
  
  pval_table_ct[row_gene,] = as.numeric(pvals)
  lfc_table_ct[row_gene,] = as.numeric(lfcs)
  ci_lower_table_ct[row_gene,] = as.numeric(cis[1,])
  ci_upper_table_ct[row_gene,] = as.numeric(cis[2,])
  
  is_sig = pvals <= 0.05
  is_pos = cis[1,] > 0 
  
  pos_sig_mask = is_sig & is_pos
  pos_sig_cts = names(pvals)[pos_sig_mask]
  
  lfcs = lfcs[pos_sig_mask]
  cis = cis[, pos_sig_mask]
  
  if (length(pos_sig_cts) == 1) {
    sig_pos_table_ct[row_gene, pos_sig_cts] = TRUE
    next
    
  } else if (length(pos_sig_cts) == 0) {
    next
  }
  
  non_overlapping_ci = sapply(
    cis[1, ],
    FUN = function(x)
      sum(x > cis[2, ]) > 0,
    simplify = T
  )
  
  if (any(non_overlapping_ci)) {
    pos_sig_cts = pos_sig_cts[non_overlapping_ci]
  } 
  
  sig_pos_table_ct[row_gene, pos_sig_cts] = TRUE
}

num_pos_sig_cts = apply(sig_pos_table_ct, MARGIN = 1, FUN = sum, simplify = T)
gene_one_pos_sig = names(num_pos_sig_cts)[num_pos_sig_cts == 1]
sig_pos_ct_index = apply(sig_pos_table_ct[gene_one_pos_sig,], MARGIN = 1, FUN = function(x) which(x), simplify = T)
sig_pos_ct = colnames(sig_pos_table_ct)[sig_pos_ct_index]
nebula_ct_degs = split(names(sig_pos_ct_index), sig_pos_ct)
nebula_ct_degs_filter = list()
for (c in names(nebula_ct_degs)) {
  degs = nebula_ct_degs[[c]]
  lfcs = NULL
  pvals = NULL
  for (g in degs) {
    lfcs = c(lfcs, lfc_table_ct[g, c])
    pvals = c(pvals, pval_table_ct[g, c])
  }
  names(lfcs) = degs
  lfcs = lfcs[lfcs >= 0.25]
  names(pvals) = degs
  pvals = pvals[order(pvals, decreasing = F)]
  
  nebula_ct_degs_filter[[c]] = lfcs[names(pvals)[names(pvals) %in% names(lfcs)]]
}
saveRDS(nebula_ct_degs_filter, file = file.path("intermediate_output", "nebula", "nebula_alistair_integrate_top_ct_DEGs.rds") )

# ----- identify DEGs only associated with cluster and not deconvolution proportions, export ----- # 

ct_deg_table = pval_table_ct <= 0.05 & lfc_table_ct >= 0.25
ct_degs = apply(ct_deg_table, MARGIN = 1, FUN = sum)
ct_degs = names(ct_degs[ct_degs >= 1])
nebula_cluster_degs_ct_independent = sapply(nebula_cluster_degs_filter, FUN = function(x) x[!names(x) %in% ct_degs])
nebula_cluster_degs_ct_dependent = sapply(nebula_cluster_degs_filter, FUN = function(x) x[names(x) %in% ct_degs])
nebula_cluster_degs_ct_independent = nebula_cluster_degs_ct_independent[as.character(0:11)[as.character(0:11) %in% names(nebula_cluster_degs_ct_independent)]]
nebula_cluster_degs_ct_dependent = nebula_cluster_degs_ct_dependent[as.character(0:11)[as.character(0:11) %in% names(nebula_cluster_degs_ct_dependent)]]
saveRDS(nebula_cluster_degs_ct_independent, file =file.path("intermediate_output", "nebula", "nebula_alistair_integrate_top_cluster_DEGs_ct_ind.rds") )
saveRDS(nebula_cluster_degs_ct_dependent, file =file.path("intermediate_output", "nebula", "nebula_alistair_integrate_top_cluster_DEGs_ct_dep.rds") )

# ----- export all DEGs associated with cluster ----- # 

all_deg_table = de_cluster_pval_adj <= 0.05 & de_cluster_lfc >= 0.5
nebula_cluster_all_deg= apply(all_deg_table, MARGIN = 2, FUN = function(x) de_cluster$gene[x], simplify = F)
names(nebula_cluster_all_deg) = gsub("seurat_clusters", "", names(nebula_cluster_all_deg))
nebula_cluster_all_deg = nebula_cluster_all_deg[sapply(nebula_cluster_all_deg, length, simplify = ) > 0]

all_ct_deg_table = de_ct_pval_adj <= 0.05 & de_ct_lfc >= 0.25
all_ct_degs = de_cluster$gene[apply(all_ct_deg_table, MARGIN = 1, FUN=function(x) any(x), simplify = T)]

nebula_cluster_all_deg_ct_independent = sapply(nebula_cluster_all_deg, FUN = function(x) x[!x %in% all_ct_degs])

saveRDS(nebula_cluster_all_deg, file = file.path("intermediate_output", "nebula", "nebula_alistair_integrate_all_DEGs.rds") )
saveRDS(nebula_cluster_all_deg_ct_independent, file = file.path("intermediate_output", "nebula", "nebula_alistair_integrate_all_DEGs_ct_independent.rds") )


