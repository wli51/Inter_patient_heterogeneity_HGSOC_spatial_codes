### This script involves manual annoatation of clustering of the integrated 
# dataset informed by analysis of preliminary results

### Prerequisites:  02_Alistair_st_integrate.R output

fibroblast_cluster_idx = c("3", "5", "6", "8", "11")
macrophage_cluster_idx = c("4")
malignant_cluster_idx = c("0", "2", "7", "9", "10")
stroma_cluster_idx = c("1")

cluster_order = c(malignant_cluster_idx, macrophage_cluster_idx, stroma_cluster_idx, fibroblast_cluster_idx)

save(fibroblast_cluster_idx, 
     macrophage_cluster_idx,
     malignant_cluster_idx,
     stroma_cluster_idx,
     cluster_order,
     file = file.path("intermediate_output", "integration_annotation.RData"))

