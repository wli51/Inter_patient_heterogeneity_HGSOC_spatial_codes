library(Seurat)
library(ggplot2)
library(pheatmap)
library('org.Hs.eg.db')
library(DOSE)
library(enrichplot)
library(clusterProfiler)

### This script runs cluster profiler function on nebula DE results 

### Prerequisites:  04_b_Alistair_st_integrate_cluster_nebula_DE.R

# list of cluster based seurat differential expression analysis outputs, 1 for each cluster 
nebula_cluster_degs_filter = readRDS(file =file.path("intermediate_output", "nebula", "nebula_alistair_integrate_cluster_DEGs.rds") )
nebula_ct_degs_filter = readRDS(file = file.path("intermediate_output", "nebula", "nebula_alistair_integrate_ct_DEGs.rds") )
nebula_cluster_degs_ct_independent = readRDS(file =file.path("intermediate_output", "nebula", "nebula_alistair_integrate_cluster_DEGs_ct_ind.rds") )
nebula_cluster_degs_ct_dependent = readRDS(file =file.path("intermediate_output", "nebula", "nebula_alistair_integrate_cluster_DEGs_ct_dep.rds") )

nebula_cluster_deg_ct_independent_entrez = list()
for (c in names(nebula_cluster_degs_ct_independent)) {
  symbols = names(nebula_cluster_degs_ct_independent[[c]])
  if (length(symbols) == 0) {
    next
  }
  id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  mapping = id_map[which(!is.na(id_map))]
  map_idx = match(symbols, names(mapping))
  mapped = mapping[map_idx]
  mapped = mapped[!is.na(map_idx)]
  nebula_cluster_deg_ct_independent_entrez[[c]] = mapped
}

organism = 'org.Hs.eg.db'
xx.ego.BP <-
  compareCluster(
    nebula_cluster_deg_ct_independent_entrez,
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
    "nebula",
    "nebula_deg_ct_independent_ego_bp.rds"
  )
)

# ----- analysis of DEGs only associated with cluster ----- # 

nebula_cluster_deg_ct_dependent_entrez = list()
for (c in names(nebula_cluster_degs_ct_dependent)) {
  symbols = names(nebula_cluster_degs_ct_dependent[[c]])
  if (length(symbols) == 0) {
    next
  }
  id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  mapping = id_map[which(!is.na(id_map))]
  map_idx = match(symbols, names(mapping))
  mapped = mapping[map_idx]
  mapped = mapped[!is.na(map_idx)]
  nebula_cluster_deg_ct_dependent_entrez[[c]] = mapped
}

xx.ego.BP <-
  compareCluster(
    nebula_cluster_deg_ct_dependent_entrez,
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
    "nebula",
    "nebula_deg_ct_dependent_ego_bp.rds"
  )
)

# ----- analysis of all significant DEGs ----- # 

nebula_cluster_deg_entrez = list()
for (c in names(nebula_cluster_degs_filter)) {
  symbols = names(nebula_cluster_degs_filter[[c]])
  if (length(symbols) == 0) {
    next
  }
  id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  mapping = id_map[which(!is.na(id_map))]
  map_idx = match(symbols, names(mapping))
  mapped = mapping[map_idx]
  mapped = mapped[!is.na(map_idx)]
  nebula_cluster_deg_entrez[[c]] = mapped
}

xx.ego.BP <-
  compareCluster(
    nebula_cluster_deg_entrez,
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
    "nebula",
    "nebula_all_deg_ego_bp.rds"
  )
)

# ----- analysis of DEGs assocaited with deconvolution proportions ----- # 

nebula_ct_deg_entrez = list()
for (c in names(nebula_ct_degs_filter)) {
  symbols = names(nebula_ct_degs_filter[[c]])
  if (length(symbols) == 0) {
    next
  }
  id_map = AnnotationDbi::mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  mapping = id_map[which(!is.na(id_map))]
  map_idx = match(symbols, names(mapping))
  mapped = mapping[map_idx]
  mapped = mapped[!is.na(map_idx)]
  nebula_ct_deg_entrez[[c]] = mapped
}

xx.ego.BP <-
  compareCluster(
    nebula_ct_deg_entrez,
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
    "nebula",
    "nebula_ct_deg_ego_bp.rds"
  )
)
