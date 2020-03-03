####################################################################################
####################################################################################
# split seurat objects by batch 
all.pbmcs.list <- SplitObject(all_pbmcs, split.by = "batch")


####################################################################################
####################################################################################
# normalize all batches and find top 2000 features
for (i in 1:length(all.pbmcs.list)) {
  all.pbmcs.list[[i]] <- NormalizeData(all.pbmcs.list[[i]], verbose = FALSE)
  all.pbmcs.list[[i]] <- FindVariableFeatures(all.pbmcs.list[[i]], 
                                              selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

reference.list <- all.pbmcs.list[c("batch1", "batch2", "batch3", "batch4", 
                                   "batch5", "batch6", "batch7")]
all.pbmcs.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated.pbmcs <- IntegrateData(anchorset = all.pbmcs.anchors, dims = 1:30)
saveRDS(all.pbmcs.anchors, "all.pbmc.anchors.rds")
saveRDS(integrated.pbmcs, "integrated.pbmcs.rds")