nmf_de_intersect_manual = readRDS(file = file.path("intermediate_output", "NMF_intersect_program_annotate.rds"))

for (i in 1:length(sample_abbrs)) {
  sample_abbr = sample_abbrs[i]
  
  sobj = readRDS(file.path("..", "analysis_output", paste0("visium_", sample_abbr,".rds")))
  sobj = FindSpatiallyVariableFeatures(sobj, assay = "SCT", slot = "scale.data", features = unlist(nmf_de_intersect_manual), selection.method = "moransi")
  saveRDS(sobj@assays$SCT@meta.features[unlist(nmf_de_intersect_manual),], file = file.path("intermediate_output", paste0(sample_abbr, "_NMF_cluster_DE_moransi.rds")))
}