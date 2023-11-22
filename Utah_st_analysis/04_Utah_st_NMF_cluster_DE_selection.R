

nmf_de_df_lists = readRDS(file = file.path("intermediate_output", "NMF_programs_de_df_lists.rds"))
nmf_de_df_lists = sapply(nmf_de_df_lists, FUN = function(x) x[!sapply(x, is.null, simplify = T)])
nmf_de_df_lists = sapply(nmf_de_df_lists, FUN = function(de_df_list) sapply(de_df_list, FUN = function(x) x[x$avg_log2FC > 0.25,], simplify = F), simplify = F)

cluster_present_in_n_samples = sapply(nmf_de_df_lists, length)
sample_specific_clusters = as.character(which(cluster_present_in_n_samples == 1))

nmf_de_intersect = sapply(nmf_de_df_lists, FUN = function(x) Reduce(intersect, sapply(x, rownames, simplify = F)))
nmf_de_intersect_unspecific = unique(unlist(nmf_de_intersect)[duplicated(unlist(nmf_de_intersect))])

unspecific_gene_assign = NULL

for (g in nmf_de_intersect_unspecific) {
  
  min_lfc = NULL 
  cluster_idx = NULL
  
  for (i in 1:length(nmf_de_intersect)) {
    if (g %in% nmf_de_intersect[[i]]) {
      de_df_list = nmf_de_df_lists[[i]]
      min_lfc = c(min_lfc,
                  min(sapply(nmf_de_df_lists[[i]], FUN = function(x) x[g, ][["avg_log2FC"]], simplify = T)))
      cluster_idx = c(cluster_idx, i)
    }
  }
  
  lfc_diff = max(min_lfc) - min_lfc 
  if (sum(lfc_diff > 0.1) == (length(min_lfc) - 1)) {
    unspecific_gene_assign = c(unspecific_gene_assign, cluster_idx[which(min_lfc == max(min_lfc))])
  } else {
    unspecific_gene_assign = c(unspecific_gene_assign, -1)
  }
  
}

nmf_de_intersect_manual = sapply(nmf_de_intersect, FUN = function(x) x[! x %in% nmf_de_intersect_unspecific])

for (i in 1:length(unspecific_gene_assign)) {
  if (unspecific_gene_assign[i] != -1) {
    nmf_de_intersect_manual[[unspecific_gene_assign[i]]] = c(nmf_de_intersect_manual[[unspecific_gene_assign[i]]], 
                                                             nmf_de_intersect_unspecific[i])
  }
}
names(nmf_de_intersect_manual) = as.character(1:10)
nmf_de_intersect_manual = nmf_de_intersect_manual[sapply(nmf_de_intersect_manual, length, simplify = T) > 0]

nmf_de_intersect_manual_lfc = list()
for (c in names(nmf_de_intersect_manual)) {
  lfc = NULL
  for (g in nmf_de_intersect_manual[[c]]) {
    lfc = c(lfc, 
            min(sapply(nmf_de_df_lists[[as.numeric(c)]], FUN = function(x) x[g, ][["avg_log2FC"]], simplify = T)))
  }
  nmf_de_intersect_manual_lfc[[c]] = setNames(lfc, nmf_de_intersect_manual[[c]])
}
nmf_de_intersect_manual_lfc = sapply(nmf_de_intersect_manual_lfc, FUN = function(x) x[order(x, decreasing = T)], simplify = F)
nmf_de_intersect_manual_lfc = readRDS(nmf_de_intersect_manual_lfc, file = file.path("intermediate_output", "NMF_intersect_program_lfc.rds"))
nmf_de_intersect_manual = sapply(nmf_de_intersect_manual_lfc, FUN = function(x) names(x))
nmf_de_intersect_manual = nmf_de_intersect_manual[sapply(nmf_de_intersect_manual, length, simplify = T) >= 5]

nmf_de_intersect_de_type = rep("consensus", length(nmf_de_intersect_manual))
nmf_de_intersect_de_type[names(nmf_de_intersect_manual) %in% sample_specific_clusters] = "specific"

saveRDS(nmf_de_intersect_manual_lfc, file = file.path("intermediate_output", "NMF_intersect_program_lfc.rds"))
saveRDS(nmf_de_intersect_manual, file = file.path("intermediate_output", "NMF_intersect_program.rds"))
saveRDS(nmf_de_intersect_de_type, file = file.path("intermediate_output", "NMF_intersect_type.rds"))
