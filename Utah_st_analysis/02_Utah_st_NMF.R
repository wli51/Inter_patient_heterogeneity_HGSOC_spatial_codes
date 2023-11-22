library(Seurat)
library(dplyr)
library(ccfindR)
library(S4Vectors)

sample_abbrs = c("1A", "2B", "3C", "4D")

obj.list = list()
for (i in 1:length(sample_abbrs)) {
  
  obj.list[[i]] = readRDS(file.path("..", "analysis_output", "Utah_st", paste0("visium_", sample_abbrs[i], ".rds")))
  obj.list[[i]] = RenameCells(obj.list[[i]], add.cell.id = sample_abbrs[i])
  
}
common_gene_names = Reduce(intersect, sapply(obj.list, FUN = function(x) rownames(x@assays$SCT@counts)))
svgs = Reduce(intersect, sapply(obj.list, FUN = function(x) x@assays$SCT@meta.features %>% 
                              dplyr::filter(moransi.spatially.variable) %>% 
                              dplyr::arrange(moransi.spatially.variable.rank) %>% 
                              top_n(1000, wt = -moransi.spatially.variable.rank) %>%
                              row.names(.)))

mat = NULL
meta.data = NULL
for (i in 1:length(sample_abbrs)) {
  mat = rbind(mat,
              obj.list[[i]]@assays$SCT@counts[svgs,])
  meta.data = rbind(meta.data,
                    obj.list[[i]]@meta.data)
}

sc <- scNMFSet(count = mat)
colData(sc) = DataFrame(meta.data)

sc <- vb_factorize(sc, ranks = seq(8,16), nrun = 5, verbose = 1, Tol = 1e-4, progress.bar = FALSE)

saveRDS(sc, file = file.path("intermediate_output", "vb_nmf.rds"))