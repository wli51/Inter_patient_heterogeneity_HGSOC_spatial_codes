library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(nebula)
library(corrplot)
library(caret)
library(readr)

### this script performs cluster based DE expression on integrated visium dataset with nebula

if (!dir.exists(file.path("intermediate_output", "nebula"))) {
  dir.create(file.path("intermediate_output", "nebula"))
}
if (!dir.exists(file.path("intermediate_output", "nebula", "obj"))) {
  dir.create(file.path("intermediate_output", "nebula", "obj"))
}

if (!dir.exists(file.path("intermediate_output", "nebula", "DEGs"))) {
  dir.create(file.path("intermediate_output", "nebula", "DEGs"))
}

### Prerequisites:  02_Alistair_st_integrate.R output

# ----- Load integrated visium dataset ----- # 
obj.integrated = readRDS(file.path("..", "analysis_output", "Alistair_st", "visium_alistair_integrate_diet.rds"))

# ----- Load deconvolution results for use as Nebula covariates ----- # 
annotation_dict = read_csv(file.path("..",
                                     "data", 
                                     "scVI_integration_output", 
                                     "annotation_dict.csv"),
                           col_names = F)

deconv.props.list.RCTD = readRDS(file = file.path("..", 
                                                  "deconvolution", 
                                                  "intermediate_output", 
                                                  "RCTD",
                                                  "props.rds"))

deconv.props.list.CARD = readRDS(file = file.path("..", 
                                                  "deconvolution", 
                                                  "intermediate_output", 
                                                  "CARD",
                                                  "props.rds"))
deconv.props.bind.RCTD = NULL
deconv.props.bind.CARD = NULL
for (x in names(deconv.props.list.RCTD)) {
  source = strsplit(x, "_")[[1]][1]
  abbr = strsplit(x, "_")[[1]][2] 
  
  df = tibble::as_tibble(deconv.props.list.RCTD[[x]]) %>%
    dplyr::mutate(spot_id = paste0(abbr, "_", row.names(deconv.props.list.RCTD[[x]])),
                  sample = abbr)
  
  deconv.props.bind.RCTD = 
    rbind(deconv.props.bind.RCTD,
          df)
  
  df = tibble::as_tibble(deconv.props.list.CARD[[x]]) %>%
    dplyr::mutate(spot_id = paste0(abbr, "_", row.names(deconv.props.list.CARD[[x]])),
                  sample = abbr)
  
  deconv.props.bind.CARD = 
    rbind(deconv.props.bind.CARD,
          df)
}

# Load indiviudal seurat objects corresponding to the samples integrated in the integrated dataset and extract SCT processed counts
sample_abbrs = paste0("SP", 1:8)
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
rm(obj.list)
gc()

# assign merged SCT count matrix as an assay
match_idx = match(colnames(obj.integrated), colnames(count_mat_merge))
sct_count_assay = CreateAssay5Object(counts = count_mat_merge[,match_idx])
obj.integrated[["SCT"]] = sct_count_assay
rm(count_mat_merge, sct_count_assay)
gc()

# ----- Run Nebula pipeline for all combinations of deconvolution method proportion and cluster ----- # 
deconv.types = c("RCTD", "CARD")
deconv.props.bind.list = list(
  "RCTD" = deconv.props.bind.RCTD,
  "CARD" = deconv.props.bind.CARD
)

for (deconv.type in deconv.types) {
  
  deconv.props.bind = deconv.props.bind.list[[deconv.type]]
  
  match_idx = match(rownames(obj.integrated@meta.data), deconv.props.bind$spot_id)
  deconv.props.bind = deconv.props.bind[match_idx,] %>% dplyr::select(-c("spot_id", "sample"))
  obj.integrated@meta.data = obj.integrated@meta.data[,!colnames(obj.integrated@meta.data) %in% unique(annotation_dict$X2)]
  obj.integrated@meta.data = cbind(obj.integrated@meta.data, deconv.props.bind)
  celltypes = colnames(deconv.props.bind)
  metadata = obj.integrated@meta.data
  
  # ----- Nebula object setup ----- # 
  neb_seurat <- scToNeb(obj = obj.integrated, 
                        assay = "SCT", 
                        id = "sample", 
                        pred = c("seurat_clusters", celltypes, "sample"), 
                        offset="nCount_SCT")
  
  neb_seurat$count = neb_seurat$count[rownames(obj.integrated@assays$integrated@data),]
  
  # Nebula design matrix construction
  
  # construct the seurat cluster one-hot encoding design matrix
  # manually define seurat cluster 1 as reference cluster 
  ref_cluster = "1"
  levels_reorder = levels(neb_seurat$pred[,"seurat_clusters"])
  levels_reorder = c(ref_cluster, levels_reorder[levels_reorder != ref_cluster])
  neb_seurat$pred[,"seurat_clusters"] = factor(neb_seurat$pred[,"seurat_clusters"], levels = levels_reorder)
  des_mat_cluster = model.matrix( ~ seurat_clusters, data = neb_seurat$pred)
  
  # construct the deconvolution proportion design matrix
  celltypes_format = gsub("[\\(\\);_ ]", "\\.", celltypes) # format celltype names to comply with R dataframe column naming requirement
  des_mat_prop = model.matrix(formula(paste0("~ ", paste(celltypes_format, collapse = " + "), " + 0")), data = neb_seurat$pred)
  des_mat = cbind(des_mat_cluster, des_mat_prop)
  linear_dep_idx2remove = findCorrelation(cor(des_mat[,2:dim(des_mat)[2]]), cutoff = 0.7) # remove celltype proportions that are colinear with other predictors
  des_mat = des_mat[,-(linear_dep_idx2remove+1)]
  # run Nebula with deconvolution proportions as covariates
  
  re = nebula(neb_seurat$count,
              neb_seurat$id,
              pred = des_mat,
              offset = neb_seurat$offset)
  
  saveRDS(re,
          file = file.path("intermediate_output",
                           "nebula",
                           "obj",
                           paste0("nebula_", deconv.type, "_corrected_obj.rds")))
  
}

# ----- Run Nebula pipeline for all combinations of deconvolution method proportion and sample ----- # 

deconv.type = deconv.types[[2]]

deconv.props.bind = deconv.props.bind.list[[deconv.type]]

match_idx = match(rownames(obj.integrated@meta.data), deconv.props.bind$spot_id)
deconv.props.bind = deconv.props.bind[match_idx,] %>% dplyr::select(-c("spot_id", "sample"))
obj.integrated@meta.data = obj.integrated@meta.data[,!colnames(obj.integrated@meta.data) %in% unique(annotation_dict$X2)]
obj.integrated@meta.data = cbind(obj.integrated@meta.data, deconv.props.bind)
celltypes = colnames(deconv.props.bind)
metadata = obj.integrated@meta.data

# ----- Nebula object setup ----- # 
neb_seurat <- scToNeb(obj = obj.integrated, 
                      assay = "SCT", 
                      id = "sample", 
                      pred = c("seurat_clusters", celltypes, "sample"), 
                      offset="nCount_SCT")

# Nebula design matrix construction

# construct the sample id one-hot encoding design matrix
des_mat_sample = model.matrix( ~ sample, data = neb_seurat$pred)

# construct the deconvolution proportion design matrix
celltypes_format = gsub("[\\(\\);_ ]", "\\.", celltypes) # format celltype names to comply with R dataframe column naming requirement
des_mat_prop = model.matrix(formula(paste0("~ ", paste(celltypes_format, collapse = " + "), " + 0")), data = neb_seurat$pred)
des_mat = cbind(des_mat_sample, des_mat_prop)
linear_dep_idx2remove = findCorrelation(cor(des_mat[,2:dim(des_mat)[2]]), cutoff = 0.7) # remove celltype proportions that are colinear with other predictors
des_mat = des_mat[,-(linear_dep_idx2remove+1)]

re = nebula(neb_seurat$count,
            neb_seurat$id,
            pred = des_mat,
            offset = neb_seurat$offset)

saveRDS(re,
        file = file.path("intermediate_output",
                         "nebula",
                         "obj",
                         paste0("nebula_sample_", deconv.type, "_corrected_obj.rds")))

# ----- Run Nebula without including the deconvolution proportions ----- # 
re = nebula(neb_seurat$count,
            neb_seurat$id,
            pred = des_mat_cluster,
            offset = neb_seurat$offset)

saveRDS(re,
        file = file.path("intermediate_output",
                         "nebula",
                         "obj",
                         paste0("nebula_", "uncorrected_obj.rds")))

# ----- helper function to find DEGs significantly associated with cluster ----- # 
filter_DEGs <- function(nebula.obj,
                        covariate_prefix,
                        lfc_thresh = 0.1,
                        factor_levels = NULL) {
  
  if (covariate_prefix == "") {
    if (is.null(factor_levels)) {
      return(NULL)
    }
    
    columns = colnames(nebula.obj$summary)[sapply(colnames(nebula.obj$summary), FUN = function(y) 
      any(sapply(factor_levels, FUN = function(x) grepl(x, y), simplify = T)),
      simplify = T)]
    
    de_df = nebula.obj$summary %>% 
      dplyr::select(all_of(columns),
                    colnames(.)[grepl("Intercept", colnames(.))],
                    "gene", "gene_id"
      )
    
    factor_levels = gsub("^logFC_", "", columns[grepl("^logFC_", columns)])
  } else {
    de_df = nebula.obj$summary %>%
      dplyr::select(colnames(.)[grepl(covariate_prefix, colnames(.))],
                    colnames(.)[grepl("Intercept", colnames(.))],
                    "gene", "gene_id")
    if (is.null(factor_levels)) {
      factor_levels = gsub(paste0("^logFC_", covariate_prefix), "",
                           colnames(de_df)[grepl(paste0("^logFC_", covariate_prefix), colnames(de_df))])
      
    } else {
      return(NULL)
    }
  }
  
  de_df_pval = de_df %>% dplyr::select(starts_with("p_"))
  de_df_pval_adj = as.data.frame(apply(
    de_df_pval,
    MARGIN = 2,
    FUN = function(x)
      p.adjust(x, p.adjust.methods[7])
  ))
  
  colnames(de_df_pval_adj) = gsub("^p_", "", colnames(de_df_pval_adj))
  de_df_lfc = de_df %>% dplyr::select(starts_with("logFC_"))
  colnames(de_df_lfc) = gsub("^logFC_", "", colnames(de_df_lfc))
  de_df_lfc = de_df_lfc[, colnames(de_df_pval_adj)]
  
  de_df_se = de_df %>% dplyr::select(starts_with("se_"))
  colnames(de_df_se) = gsub("^se_", "", colnames(de_df_se))
  de_df_se = de_df_se[, colnames(de_df_pval_adj)]
  
  sig_pos_table = matrix(FALSE,
                         nrow = nrow(de_df),
                         ncol = length(factor_levels))
  rownames(sig_pos_table) = de_df$gene
  colnames(sig_pos_table) = factor_levels
  lfc_table = matrix(NaN,
                     nrow = nrow(de_df),
                     ncol = length(factor_levels))
  rownames(lfc_table) = de_df$gene
  colnames(lfc_table) = factor_levels
  pval_table = lfc_table
  ci_lower_table = lfc_table
  ci_upper_table = lfc_table
  
  # iterate over all nebula genes, for each gene, compare the nebula DE factor_levels
  # coefficients that are indicative of significant positive log fold changes,
  # if a gene is positively up-regulated for more than 1 cluster, and the cluster
  # coefficient with the highest logFC has non-overlapping CI with others, then
  # the gene is deemed DEG for that highest logFC cluster, otherwise, the gene is
  # considered to be up-regualted in all factor_levels with positive, significant
  # nebula coefficinet
  for (i in 1:nrow(de_df)) {
    de_df_row = de_df[i, ]
    row_gene = de_df_row[, "gene"]
    
    pvals = de_df_pval_adj[i, ] %>% dplyr::select(-c("(Intercept)"))
    ses = de_df_se[i,] %>% dplyr::select(-c("(Intercept)"))
    lfcs = de_df_lfc[i, ] %>% dplyr::select(-c("(Intercept)"))
    cis = rbind(lfcs - 2 * ses,
                lfcs + 2 * ses)
    
    pval_table[row_gene, ] = as.numeric(pvals)
    lfc_table[row_gene, ] = as.numeric(lfcs)
    ci_lower_table[row_gene, ] = as.numeric(cis[1, ])
    ci_upper_table[row_gene, ] = as.numeric(cis[2, ])
    
    is_sig = pvals <= 0.05
    is_pos = cis[1, ] > 0
    
    pos_sig_mask = is_sig & is_pos
    pos_sig_factor_levels = gsub(covariate_prefix, "", names(pvals)[pos_sig_mask])
    
    lfcs = lfcs[pos_sig_mask]
    cis = cis[, pos_sig_mask]
    
    if (length(pos_sig_factor_levels) == 1) {
      sig_pos_table[row_gene, pos_sig_factor_levels] = TRUE
      next
      
    } else if (length(pos_sig_factor_levels) == 0) {
      next
    }
    
    non_overlapping_ci = sapply(
      cis[1,],
      FUN = function(x)
        sum(x > cis[2,]) > 0,
      simplify = T
    )
    
    if (any(non_overlapping_ci)) {
      pos_sig_factor_levels = pos_sig_factor_levels[non_overlapping_ci]
    }
    
    sig_pos_table[row_gene, pos_sig_factor_levels] = TRUE
  }
  
  # filter for nebula DEGs only found to be significantly up-regulated for one factor level
  num_pos_sig = apply(sig_pos_table,
                      MARGIN = 1,
                      FUN = sum,
                      simplify = T)
  gene_one_pos_sig = names(num_pos_sig)[num_pos_sig == 1]
  sig_pos_index = apply(
    sig_pos_table[gene_one_pos_sig, ],
    MARGIN = 1,
    FUN = function(x)
      which(x),
    simplify = T
  )
  sig_pos = colnames(sig_pos_table)[sig_pos_index]
  nebula_degs = split(names(sig_pos_index), sig_pos)
  
  nebula_degs_filter = list()
  for (c in names(nebula_degs)) {
    degs = nebula_degs[[c]]
    lfcs = NULL
    pvals = NULL
    for (g in degs) {
      lfcs = c(lfcs, lfc_table[g, c])
      pvals = c(pvals, pval_table[g, c])
    }
    names(lfcs) = degs
    lfcs = lfcs[lfcs >= lfc_thresh]
    names(pvals) = degs
    pvals = pvals[order(pvals, decreasing = F)]
    
    nebula_degs_filter[[c]] = lfcs[names(pvals)[names(pvals) %in% names(lfcs)]]
  }
  
  list("degs" = nebula_degs,
       "degs_filter" = nebula_degs_filter,
       "lfc" = lfc_table,
       "pval" = pval_table,
       "sig" = sig_pos_table
  )
}

# ----- Select DEGs ----- # 

for (fname in c(
  "nebula_RCTD_corrected",
  "nebula_CARD_corrected",
  "nebula_uncorrected"
)) {
  re = readRDS(file = file.path(
    "intermediate_output",
    "nebula",
    "obj",
    paste0(fname, "_obj.rds")
  ))
  # TODO
  DEGs.cluster = filter_DEGs(re, covariate_prefix = "seurat_clusters")
  # DEGs.cluster = filter_DEGs(re, covariate_prefix = "seurat_clusters")
  
  saveRDS(DEGs.cluster,
          file = file.path(
            "intermediate_output",
            "nebula",
            "DEGs",
            paste0(fname, "_DEGs.rds")
          ))
  
  if (grepl("_corrected", fname)) {
    DEGs.ct = filter_DEGs(re, covariate_prefix = "", factor_levels = celltypes_format)
    
    saveRDS(DEGs.ct,
            file = file.path(
              "intermediate_output",
              "nebula",
              "DEGs",
              paste0(fname, "_ct_DEGs.rds")
            ))
  }
}

# map deconvolution information 
match_idx = match(rownames(obj.integrated@meta.data), deconv_summary_df$spot_id)
deconv_summary_df = deconv_summary_df[match_idx,] %>% dplyr::select(-c("spot_id", "sample"))
obj.integrated@meta.data = cbind(obj.integrated@meta.data, deconv_summary_df)
celltypes = colnames(deconv_summary_df)
metadata = obj.integrated@meta.data

# ----- Nebula object setup ----- # 
neb_seurat <- scToNeb(obj = obj.integrated, 
                      assay = "SCT", 
                      id = "sample", 
                      pred = c("seurat_clusters", celltypes, "sample"), 
                      offset="nCount_SCT")

neb_seurat$count = neb_seurat$count[rownames(obj.integrated@assays$integrated@data),]

# Nebula design matrix construction

ref_cluster = "1"
levels_reorder = levels(neb_seurat$pred[,"seurat_clusters"])
levels_reorder = c(ref_cluster, levels_reorder[levels_reorder != ref_cluster])
neb_seurat$pred[,"seurat_clusters"] = factor(neb_seurat$pred[,"seurat_clusters"], levels = levels_reorder)
des_mat_cluster = model.matrix( ~ seurat_clusters, data = neb_seurat$pred)

celltypes_format = gsub("[\\(\\);_ ]", "\\.", celltypes) # format celltype names to comply with perper dataframe column naming
des_mat_prop = model.matrix(formula(paste0("~ ", paste(celltypes_format, collapse = " + "), " + 0")), data = neb_seurat$pred)
des_mat = cbind(des_mat_cluster, des_mat_prop)
# df = des_mat[,2:dim(des_mat)[2]]
linear_dep_idx2remove = findCorrelation(cor(des_mat[,2:dim(des_mat)[2]]), cutoff = 0.7) # remove celltype proportions that are colinear with other predictors
des_mat = des_mat[,-(linear_dep_idx2remove+1)]

re.alistair.cluster = nebula(neb_seurat$count,
                             neb_seurat$id,
                             pred = des_mat_cluster,
                             offset = neb_seurat$offset)

saveRDS(re.alistair.cluster, file = file.path("intermediate_output", "nebula", "nebula_alistair_integrate_cluster_obj.rds"))


# helper function to find DEGs significantly associated with cluster
filter_DEGs <- function(nebula.obj,
                        covariate_prefix,
                        lfc_thresh = 0.1) {
  
  de_df = nebula.obj$summary %>% 
    dplyr::select(colnames(.)[grepl(covariate_prefix, colnames(.))], 
                  colnames(.)[grepl("Intercept", colnames(.))],
                  "gene", "gene_id")
  factor_levels = gsub(paste0("^logFC_", covariate_prefix), "", 
                       colnames(de_df)[grepl(paste0("^logFC_", covariate_prefix), colnames(de_df))])
  de_df_pval = de_df %>% dplyr::select(starts_with("p_"))
  de_df_pval_adj = as.data.frame(apply(
    de_df_pval,
    MARGIN = 2,
    FUN = function(x)
      p.adjust(x, p.adjust.methods[7])
  ))
  
  colnames(de_df_pval_adj) = gsub("^p_", "", colnames(de_df_pval_adj))
  de_df_lfc = de_df %>% dplyr::select(starts_with("logFC_"))
  colnames(de_df_lfc) = gsub("^logFC_", "", colnames(de_df_lfc))
  de_df_lfc = de_df_lfc[, colnames(de_df_pval_adj)]
  
  de_df_se = de_df %>% dplyr::select(starts_with("se_"))
  colnames(de_df_se) = gsub("^se_", "", colnames(de_df_se))
  de_df_se = de_df_se[, colnames(de_df_pval_adj)]
  
  sig_pos_table = matrix(FALSE,
                         nrow = nrow(de_df),
                         ncol = length(factor_levels))
  rownames(sig_pos_table) = de_df$gene
  colnames(sig_pos_table) = factor_levels
  lfc_table = matrix(NaN,
                     nrow = nrow(de_df),
                     ncol = length(factor_levels))
  rownames(lfc_table) = de_df$gene
  colnames(lfc_table) = factor_levels
  pval_table = lfc_table
  ci_lower_table = lfc_table
  ci_upper_table = lfc_table
  
  # iterate over all nebula genes, for each gene, compare the nebula DE factor_levels
  # coefficients that are indicative of significant positive log fold changes,
  # if a gene is positively up-regulated for more than 1 cluster, and the cluster
  # coefficient with the highest logFC has non-overlapping CI with others, then
  # the gene is deemed DEG for that highest logFC cluster, otherwise, the gene is
  # considered to be up-regualted in all factor_levels with positive, significant
  # nebula coefficinet
  for (i in 1:nrow(de_df)) {
    de_df_row = de_df[i, ]
    row_gene = de_df_row[, "gene"]
    
    pvals = de_df_pval_adj[i, ] %>% dplyr::select(-c("(Intercept)"))
    ses = de_df_se[i,] %>% dplyr::select(-c("(Intercept)"))
    lfcs = de_df_lfc[i, ] %>% dplyr::select(-c("(Intercept)"))
    cis = rbind(lfcs - 2 * ses,
                lfcs + 2 * ses)
    
    pval_table[row_gene, ] = as.numeric(pvals)
    lfc_table[row_gene, ] = as.numeric(lfcs)
    ci_lower_table[row_gene, ] = as.numeric(cis[1, ])
    ci_upper_table[row_gene, ] = as.numeric(cis[2, ])
    
    is_sig = pvals <= 0.05
    is_pos = cis[1, ] > 0
    
    pos_sig_mask = is_sig & is_pos
    pos_sig_factor_levels = gsub(covariate_prefix, "", names(pvals)[pos_sig_mask])
    
    lfcs = lfcs[pos_sig_mask]
    cis = cis[, pos_sig_mask]
    
    if (length(pos_sig_factor_levels) == 1) {
      sig_pos_table[row_gene, pos_sig_factor_levels] = TRUE
      next
      
    } else if (length(pos_sig_factor_levels) == 0) {
      next
    }
    
    non_overlapping_ci = sapply(
      cis[1,],
      FUN = function(x)
        sum(x > cis[2,]) > 0,
      simplify = T
    )
    
    if (any(non_overlapping_ci)) {
      pos_sig_factor_levels = pos_sig_factor_levels[non_overlapping_ci]
    }
    
    sig_pos_table[row_gene, pos_sig_factor_levels] = TRUE
  }
  
  # filter for nebula DEGs only found to be significantly up-regulated for one factor level
  num_pos_sig = apply(sig_pos_table,
                      MARGIN = 1,
                      FUN = sum,
                      simplify = T)
  gene_one_pos_sig = names(num_pos_sig)[num_pos_sig == 1]
  sig_pos_index = apply(
    sig_pos_table[gene_one_pos_sig, ],
    MARGIN = 1,
    FUN = function(x)
      which(x),
    simplify = T
  )
  sig_pos = colnames(sig_pos_table)[sig_pos_index]
  nebula_degs = split(names(sig_pos_index), sig_pos)
  
  nebula_degs_filter = list()
  for (c in names(nebula_degs)) {
    degs = nebula_degs[[c]]
    lfcs = NULL
    pvals = NULL
    for (g in degs) {
      lfcs = c(lfcs, lfc_table[g, c])
      pvals = c(pvals, pval_table[g, c])
    }
    names(lfcs) = degs
    lfcs = lfcs[lfcs >= lfc_thresh]
    names(pvals) = degs
    pvals = pvals[order(pvals, decreasing = F)]
    
    nebula_degs_filter[[c]] = lfcs[names(pvals)[names(pvals) %in% names(lfcs)]]
  }
  
  list("degs" = nebula_degs,
       "degs_filter" = nebula_degs_filter,
       "lfc" = lfc_table,
       "pval" = pval_table,
       "sig" = sig_pos_table
  )
}

# ----- Run Nebula without including the deconvolution proportions ----- # 
nebula_filter = filter_DEGs(re.alistair.cluster, 
                            covariate_prefix = "seurat_clusters")

saveRDS(nebula_filter, file = file.path("intermediate_output", 
                                        "nebula", 
                                        "nebula_alistair_uncorrected_DEGs.rds"))
