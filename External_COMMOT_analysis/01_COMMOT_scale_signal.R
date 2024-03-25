library(Seurat)
library(magrittr)
library(dplyr)
library(readr)
library(stringr)

### Load individual Visium samples and corresponding clustering assignment from integrated analysis 

### Prerequisites:  ** exported outgoing signal summary csv from COMMOT

metadata = readRDS(file = file.path("..", "analysis_output", "visium_alistair_integrate_cluster.rds"))
sample_abbrs = paste0("SP", 1:8)
obj.list = list()

for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  barcodes_sample = metadata[metadata$sample == sample_abbr,]$barcode
  cluster_sample = metadata[metadata$sample == sample_abbr,]$seurat_clusters
  
  obj.list[[i]] = readRDS(file.path("..", "analysis_output",  paste0("visium_", sample_abbr,".rds")))
  match_idx = match(rownames(obj.list[[i]]@meta.data), barcodes_sample)
  any(is.na(match_idx))
  obj.list[[i]]$integrate_cluster = cluster_sample[match_idx]
}


### COMMOT ligand receptor analysis output (generated per-sample)
commot_sum_sender_lr_list = list()
commot_sum_sender_pathway_list = list()
commot_info_list = list()

for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  x =  read_csv(file.path("..", "data", "COMMOT", "Alistair", sample_abbr, "commot_sender.csv"))
  info = read_csv(file.path("..", "data", "COMMOT", "Alistair", sample_abbr, "commot_info.csv"))
  
  pathway_summary = as.data.frame(x[,paste0("s-", unique(info$pathway))])
  rownames(pathway_summary) = x[[1]]
  lr_summary = as.data.frame(x) %>% dplyr::select(-paste0("s-", unique(info$pathway)), -`...1`)
  rownames(lr_summary) = x[[1]]
  
  # pathways = colnames(x)
  # x = data.frame(x)
  # rownames(x) = x[[1]]
  # x = x[,str_count(pathways, "-") == 2]
  commot_sum_sender_lr_list[[i]] = as.data.frame(lr_summary)
  colnames(commot_sum_sender_lr_list[[i]]) = gsub("s-", "", colnames(commot_sum_sender_lr_list[[i]]))
  colnames(commot_sum_sender_lr_list[[i]]) = gsub("-", "\\.", colnames(commot_sum_sender_lr_list[[i]]))
  
  commot_sum_sender_pathway_list[[i]] = as.data.frame(pathway_summary)
  colnames(commot_sum_sender_pathway_list[[i]]) = gsub("s-", "", colnames(commot_sum_sender_pathway_list[[i]]))
  colnames(commot_sum_sender_pathway_list[[i]]) = gsub("-", "\\.", colnames(commot_sum_sender_pathway_list[[i]]))
  
  commot_info_list[[i]] = info
}

df_lr_info = dplyr::bind_rows(commot_info_list) %>%
  dplyr::select(ligand, receptor, pathway) %>%
  dplyr::distinct(ligand, receptor, pathway)

saveRDS(df_lr_info, 
        file = file.path("intermediate_output", "COMMOT_lr_info.rds"))

commot_sum_sender_lr_list_scale = list()
commot_sum_sender_pathway_list_scale = list()
for (i in 1:length(sample_abbrs)) {
  commot_sum_sender_lr_list_scale[[i]] = commot_sum_sender_lr_list[[i]]/commot_sum_sender_lr_list[[i]][,'total.total'] * 100
  commot_sum_sender_lr_list_scale[[i]][commot_sum_sender_lr_list[[i]][,'total.total'] == 0,] = 0
  
  commot_sum_sender_pathway_list_scale[[i]] = commot_sum_sender_pathway_list[[i]]/commot_sum_sender_lr_list[[i]][,'total.total'] * 100
  commot_sum_sender_pathway_list_scale[[i]][commot_sum_sender_lr_list[[i]][,'total.total'] == 0,] = 0
}

# Concatenate results from all samples 
sample_factor_level = c("SP1", "SP4", "SP7", "SP2", "SP3", "SP8", "SP5", "SP6")
sample2subtype = setNames(c("DIF", "MES", "MES", "DIF", "IMR", "IMR", "DIF", "MES"),
                          paste0("SP", 1:8))

df_all_sample_lr = NULL
for (i in 1:length(sample_abbrs)) {
  df = commot_sum_sender_lr_list_scale[[i]] %>%
    dplyr::select(-c("total.total")) %>%
    dplyr::mutate(
      barcode = row.names(.),
      sample = sample_abbrs[i],
      total.signaling = commot_sum_sender_lr_list[[i]]$total.total
    ) %>%
    tidyr::pivot_longer(cols = -c("barcode", "sample", "total.signaling")) %>%
    dplyr::mutate(sample = sample_abbrs[i])
  df_all_sample_lr = rbind(df_all_sample_lr,
                           df)
}
df_all_sample_lr$subtype = sample2subtype[df_all_sample_lr$sample]
saveRDS(df_all_sample_lr, 
        file = file.path("intermediate_output", "COMMOT_scaled_lr.rds"))

df_all_sample_pathway = NULL
for (i in 1:length(sample_abbrs)) {
  df = commot_sum_sender_pathway_list_scale[[i]] %>%
    dplyr::mutate(
      barcode = row.names(.),
      sample = sample_abbrs[i],
      total.signaling = commot_sum_sender_lr_list[[i]]$total.total
    ) %>%
    tidyr::pivot_longer(cols = -c("barcode", "sample", "total.signaling")) %>%
    dplyr::mutate(sample = sample_abbrs[i])
  df_all_sample_pathway = rbind(df_all_sample_pathway,
                           df)
}
df_all_sample_pathway$subtype = sample2subtype[df_all_sample_pathway$sample]
saveRDS(df_all_sample_pathway, 
        file = file.path("intermediate_output", "COMMOT_scaled_pathway.rds"))


# Test for within subtype differentially regulated signaling
subtype2indices = list()
subtype2indices[["DIF"]] = c(1, 4, 7)
subtype2indices[["IMR"]] = c(5, 6)
subtype2indices[["MES"]] = c(2, 3, 8)

within_subtype_test_results_df = list()

for (subtype in names(subtype2indices)) {
  df_subtype_lr = NULL
  for (i in subtype2indices[[subtype]]) {
    df = commot_sum_sender_lr_list_scale[[i]] %>%
      dplyr::select(-c("total.total")) %>%
      dplyr::mutate(
        barcode = row.names(.),
        sample = sample_abbrs[i],
        total.signaling = commot_sum_sender_lr_list[[i]]$total.total
      ) %>%
      tidyr::pivot_longer(cols = -c("barcode", "sample", "total.signaling")) %>%
      dplyr::mutate(sample = sample_abbrs[i])
    df_subtype_lr = rbind(df_subtype_lr, df)
  }
  
  within_subtype_test_results_df[[subtype]] = NULL
  
  for (lr in unique(df_subtype_lr$name)) {
    
    xy = df_subtype_lr %>%
      dplyr::filter(total.signaling >= 5,
                    value >= 1,
                    name == lr) %>%
      dplyr::select(sample, name, value) %>%
      dplyr::mutate(sample = factor(sample))
    
    if (length(unique(xy$sample)) < length(subtype2indices[[subtype]])) {
      within_subtype_test_results_df[[subtype]] =
        rbind(within_subtype_test_results_df[[subtype]],
              data.frame(lr = lr, chisq = NA, pval = NA))
      next
    }
    results = kruskal.test(x = xy$value,
                           g = xy$sample,
                           na.action = "na.omit")
    
    within_subtype_test_results_df[[subtype]] =
      rbind(
        within_subtype_test_results_df[[subtype]],
        data.frame(
          lr = lr,
          chisq = results$statistic,
          pval = results$p.value
        )
      )
  }
  
  within_subtype_test_results_df[[subtype]]$pval.adj = p.adjust(within_subtype_test_results_df[[subtype]]$pval, method = "fdr")
}
# filter for significant tests
within_subtype_test_results_df_filter =  sapply(
  within_subtype_test_results_df,
  FUN = 
    function(x)
      x %>% dplyr::filter(complete.cases(.),
                          pval.adj <= 0.05) %>%
    dplyr::arrange(pval.adj),
  simplify = F
) 
saveRDS(within_subtype_test_results_df_filter, 
        file = file.path("intermediate_output", "COMMOT_within_subtype_lr.rds"))

### Test for between subtype regulated signaling 
between_subtype_test_results_df = NULL
for (lr in unique(df_all_sample_lr$name)) {
  
  # only test in spots which total signal is above 5 
  xy = df_all_sample_lr %>%
    dplyr::filter(total.signaling >= 5,
                  value >= 1,
                  name == lr) %>%
    dplyr::select(sample, name, value) %>%
    dplyr::mutate(sample = factor(sample))
  
  if (length(unique(xy$sample)) < length(sample_abbrs)) {
    between_subtype_test_results_df = 
      rbind(between_subtype_test_results_df,
            data.frame(lr = lr, chisq = NA, pval = NA))
    next
  }
  results = kruskal.test(x = xy$value,
                         g = xy$sample,
                         na.action = "na.omit")
  
  between_subtype_test_results_df = 
    rbind(between_subtype_test_results_df,
          data.frame(lr = lr, chisq = results$statistic, pval = results$p.value))
}

between_subtype_test_results_df$pval.adj = p.adjust(between_subtype_test_results_df$pval, method = "fdr")
saveRDS(between_subtype_test_results_df, 
        file = file.path("intermediate_output", "COMMOT_between_subtype_lr.rds"))


