nmf_de_intersect_manual = readRDS(file = file.path("intermediate_output", "NMF_intersect_program.rds"))
names(nmf_de_intersect_manual) = c("WFDC2/CLU/SLPI",
                                   "DAPL1",
                                   "IGKC/IGHG",
                                   "VEGFA",
                                   "SPP1/FTL",
                                   "MALAT1",
                                   "COL1A1/VIM/FN1",
                                   "MDK")
saveRDS(nmf_de_intersect_manual, file = file.path("intermediate_output", "NMF_intersect_program_annotate.rds"))