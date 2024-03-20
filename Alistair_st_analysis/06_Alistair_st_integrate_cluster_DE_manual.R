library(Seurat)
library(magrittr)
library(dplyr)


### this script involves manual curation of integrate cluster DE results

### Prerequisites:  04_Alistair_st_integrate_cluster_DE_fib_rich.R

dir = file.path("..", "analysis_output", "Alistair_st")

malignant_cluster_idx = c("0", "2", "7", "9", "10")
MHC_cluster_idx = c("4")
storma_cluster_idx = c("1")
ecm_cluster_idx = c("1", "3", "5", "6", "8", "11")
cluster_order = c(malignant_cluster_idx, MHC_cluster_idx, storma_cluster_idx, ecm_cluster_idx)

sample_abbrs = paste0("SP", 1:8)
obj.list = list()

for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  obj.list[[i]] = readRDS(file.path(dir,  paste0("visium_", sample_abbr,".rds")))
}
common_gene_names = Reduce(intersect, sapply(obj.list, FUN = function(x) rownames(x@assays[["SCT"]])))

# All Major histocompatiblity complex genes 
hla_genes = common_gene_names[grep("^HLA-[ABCEF]", common_gene_names)]
hla_d_genes = common_gene_names[grep("^HLA-D[PQRMO]", common_gene_names)]
hla_features = c(hla_genes[order(hla_genes)][1:4], hla_d_genes[order(hla_d_genes)][1:6])

# All MT genes
mt_genes = common_gene_names[grep("^MT-", common_gene_names)]
mt_features = mt_genes[order(mt_genes)][1:10]

# All Ribo genes 
rp_s_genes = common_gene_names[grep("^RPS", common_gene_names)]
rp_l_genes = common_gene_names[grep("^RPL", common_gene_names)]
rp_features = c(rp_s_genes[order(rp_s_genes)][1:5], rp_l_genes[order(rp_l_genes)][1:5])

# All immunoglobulin coding genes
ig_genes= common_gene_names[grep("^IG[KLH]", common_gene_names)]
ig_features = ig_genes[order(ig_genes)]

# isolate common fibroblast DE genes 
de_genes_fibro_vs_list = readRDS(file = file.path("intermediate_output", "integrate_fib_vs_other_cluster_DE.rds"))
fib_common_features = Reduce(intersect, sapply(de_genes_fibro_vs_list, FUN = rownames, simplify = F))
fib_common_features = fib_common_features[! fib_common_features %in% c(hla_genes, hla_d_genes, mt_genes, rp_s_genes, rp_s_genes, ig_genes)]

# Remove all the genesets defined above from the standard cluster based DE output 
DE_df_list = readRDS(file = file.path("intermediate_output", "integrate_cluster_DE.rds"))
DE_genes = sapply(DE_df_list, rownames)
DE_genes = DE_genes[cluster_order[cluster_order %in% names(DE_genes)]]
DE_genes_filter = DE_genes
DE_genes_filter = sapply(
  DE_genes_filter,
  FUN = function(x)
    x[!x %in% c(
      hla_genes,
      hla_d_genes,
      mt_genes,
      rp_s_genes,
      rp_s_genes,
      ig_genes,
      fib_common_features
    )]
)

# remove additional cluster genes 
unspecific_de_genes = unlist(DE_genes)[duplicated(unlist(DE_genes))]

DE_genes_filter = sapply(
  DE_genes_filter,
  FUN = function(x)
    x[!x %in% unspecific_de_genes]
)

DE_genes_filter[["MHC"]] = hla_features
DE_genes_filter[["IG"]] = ig_features
DE_genes_filter[["Fibroblast common DEGs"]] = fib_common_features
DE_genes_filter[["MT"]] = mt_features
DE_genes_filter[["Ribo"]] = rp_features

DE_genes_filter = DE_genes_filter[sapply(DE_genes_filter, length) >= 5]


saveRDS(
  DE_genes_filter,
  file = file.path(
    "intermediate_output",
    "integrate_cluster_DE_manual.rds"
  )
)

