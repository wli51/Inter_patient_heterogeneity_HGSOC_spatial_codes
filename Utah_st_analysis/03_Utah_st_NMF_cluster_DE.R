library(Seurat)
library(ggplot2)
library(dplyr)
library(magrittr)
library(ccfindR)

### This script takes the output of NMF and select/refine top genes for each program

# nmf output (ccfindR object)
vb_nmf <- readRDS(file.path("intermediate_output", "vb_nmf.rds"))

# extract the basis and coefficent matrices
basis_mat = basis(vb_nmf)[ranks(vb_nmf) == optimal_rank(vb_nmf)$ropt][[1]]
coef_mat = coeff(vb_nmf)[ranks(vb_nmf) == optimal_rank(vb_nmf)$ropt][[1]]

# extract NMF cluster id
cluster_ids = cluster_id(vb_nmf, rank=optimal_rank(vb_nmf)$ropt)
barcodes = sapply(strsplit(colnames(vb_nmf), "_"), "[[", 2) # spot names from integration step are of the format [sample id]_[raw barcode]
sample = sapply(strsplit(colnames(vb_nmf), "_"), "[[", 1) # get the sample id per barcode 

save(cluster_ids, barcodes, sample, file = file.path("intermediate_output", "vb_nmf_cluster_info.RData"))

# ----- Selection of genes for each program ----- # 

# extract genes for each program
gene_programs = list() # a list of the same length as the number of programs (NMF rank)
col_max_index = max.col(basis_mat)

for (i in 1:dim(basis_mat)[2]) {
  # every gene will be assigned to the program where it scores the highest in basis matrix
  scores = basis_mat[col_max_index == i , i]
  # within each program, order the genes by the scores (the ordering does not impact any analysis after)
  scores = scores[order(scores, decreasing = TRUE)]
  gene_programs[[i]] = names(scores)
}

# ----- Refinement of genes for each program ----- # 

# The visium data are needed for refinement
sample_abbrs = c("2B", "3C", "4D")

obj.list = list()
for (i in 1:length(sample_abbrs)) {
  obj.list[[i]] = readRDS(file.path("..", "analysis_output", paste0("visium_", sample_abbrs[i], ".rds")))
}

# append the NMF cluster assignment to visium data 
for (i in 1:length(sample_abbrs)) {
  filter = sample == sample_abbrs[i]
  clusters = cluster_ids[filter]
  
  match_ordering = match(colnames(obj.list[[i]]), barcodes[filter])
  
  obj.list[[i]]@meta.data$NMF_cluster = clusters[match_ordering]
}

min.cells.group = 50
nmf_de_df_lists = list()
nmf_de_df_raw_lists = list()

# iterate over the programs, performing differential expression of program genes-
# comparing spots from the same program cluster with all other spots  
for (i in 1:length(gene_programs)) {
  
  de_df_list = list()
  de_df_list_raw = list()
  
  for (j in 1:length(obj.list)) {
    
    if (sum(obj.list[[j]]@meta.data$NMF_cluster == i) >= min.cells.group) {
      
      de_markers <-
        FindMarkers(
          obj.list[[j]],
          ident.1 = i, # program cluster id
          group.by = "NMF_cluster",
          logfc.threshold = 0,
          min.cells.group = min.cells.group
        )
      
      de_df_list_raw[[j]] = de_markers
      
      de_markers %>% dplyr::filter(p_val_adj <= 0.05, 
                                   avg_log2FC >= 0.25,
                            !grepl("^RP[LS]", rownames(de_markers)),
                            !grepl("^MT-", rownames(de_markers)))
      
      
      de_df_list[[j]] = de_markers
    }
  }
  
  nmf_de_df_lists[[i]] = de_df_list
  nmf_de_df_raw_lists[[i]] = de_df_list
  
}

saveRDS(nmf_de_df_lists, file = file.path("intermediate_output", "NMF_programs_de_df_lists.rds"))
saveRDS(nmf_de_df_raw_lists, file = file.path("intermediate_output", "NMF_programs_de_df_raw_lists.rds"))
