library(Seurat)
library(ggplot2)
library(pheatmap)
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(AnnotationDbi)

### This script runs cluster profiler function on nebula DE results 

### Prerequisites:  04_Alistair_st_integrate_cluster_nebula_DE.R

if (!dir.exists(file.path("intermediate_output", "GO"))) {
  dir.create(file.path("intermediate_output", "GO"))
}

organism = 'org.Hs.eg.db'

for (fname in c("nebula_RCTD_corrected",
                "nebula_CARD_corrected",
                "nebula_uncorrected"
)) {
  DEGs = readRDS(file = file.path(
    "intermediate_output",
    "nebula",
    "DEGs",
    paste0(fname, "_DEGs.rds")
  ))
  
  # ----- analysis of significant DEGs from the uncorrected analysis ----- #
  
  deg_entrez = list()
  for (c in names(DEGs$degs_filter)) {
    symbols = names(DEGs$degs_filter[[c]])
    if (length(symbols) == 0) {
      next
    }
    id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
    mapping = id_map[which(!is.na(id_map))]
    map_idx = match(symbols, names(mapping))
    mapped = mapping[map_idx]
    mapped = mapped[!is.na(map_idx)]
    deg_entrez[[c]] = mapped
  }
  
  xx.ego.BP <-
    compareCluster(
      deg_entrez,
      fun = "enrichGO",
      ont = "BP",
      OrgDb = organism,
      pvalueCutoff = 0.05
    )
  xx.ego.BP <- pairwise_termsim(xx.ego.BP)
  saveRDS(
    xx.ego.BP,
    file = file.path(
      "intermediate_output",
      "GO",
      paste0("GO_", fname, ".rds")
    )
  )
  
  if (grepl("_corrected", fname)) {
    
    DEGs.ct = readRDS(file = file.path(
      "intermediate_output",
      "nebula",
      "DEGs",
      paste0(fname, "_ct_DEGs.rds")
    ))
    
    deg_entrez = list()
    for (c in names(DEGs.ct$degs_filter)) {
      symbols = names(DEGs.ct$degs_filter[[c]])
      if (length(symbols) == 0) {
        next
      }
      id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
      mapping = id_map[which(!is.na(id_map))]
      map_idx = match(symbols, names(mapping))
      mapped = mapping[map_idx]
      mapped = mapped[!is.na(map_idx)]
      deg_entrez[[c]] = mapped
    }
    
    xx.ego.BP <-
      compareCluster(
        deg_entrez,
        fun = "enrichGO",
        ont = "BP",
        OrgDb = organism,
        pvalueCutoff = 0.05
      )
    xx.ego.BP <- pairwise_termsim(xx.ego.BP)
    saveRDS(
      xx.ego.BP,
      file = file.path(
        "intermediate_output",
        "GO",
        paste0("GO_ct_", fname, ".rds")
      )
    )
    
    deg_cluster_only = list()
    ct_degs = Reduce(union, sapply(DEGs.ct$degs_filter, FUN = function(x) names(x)))
    for (n in names(DEGs$degs_filter)) {
      degs = DEGs$degs_filter[[n]] 
      deg_cluster_only[[n]] = degs[!names(degs) %in% ct_degs]
    }
    
    deg_entrez = list()
    for (c in names(deg_cluster_only)) {
      symbols = names(deg_cluster_only[[c]])
      if (length(symbols) == 0) {
        next
      }
      id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
      mapping = id_map[which(!is.na(id_map))]
      map_idx = match(symbols, names(mapping))
      mapped = mapping[map_idx]
      mapped = mapped[!is.na(map_idx)]
      deg_entrez[[c]] = mapped
    }
    
    xx.ego.BP <-
      compareCluster(
        deg_entrez,
        fun = "enrichGO",
        ont = "BP",
        OrgDb = organism,
        pvalueCutoff = 0.05
      )
    xx.ego.BP <- pairwise_termsim(xx.ego.BP)
    saveRDS(
      xx.ego.BP,
      file = file.path(
        "intermediate_output",
        "GO",
        paste0("GO_cluster_only_", fname, ".rds")
      )
    )
  }
}
