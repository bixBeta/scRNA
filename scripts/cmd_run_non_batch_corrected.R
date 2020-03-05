library("Seurat")
library("ggplot2")

# load seurat object
s2.0.cmd <- readRDS("RDS/s2.0.Filtered.SeuratObject.rds")

# fast analysis 
options(future.globals.maxSize = 9000 * 1024^2)

# run sctransform
s2.0.cmd <- SCTransform(s2.0.cmd, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)

# Perform dimensionality reduction by PCA and UMAP embedding
# These are now standard steps in the Seurat workflow for visualization and clustering

s2.0.cmd <- RunPCA(s2.0.cmd, verbose = FALSE)
s2.0.cmd <- RunUMAP(s2.0.cmd, dims = 1:30, verbose = FALSE)
s2.0.cmd <- FindNeighbors(s2.0.cmd, dims = 1:30, verbose = FALSE)
s2.0.cmd <- FindClusters(s2.0.cmd, verbose = FALSE)

png(filename = "plots/cmd_version_non_batch_corrected_clustering.png", width = 2400, height = 1200)
DimPlot(s2.0.cmd, label = TRUE, group.by = "batch")
dev.off()

saveRDS(s2.0.cmd, file = "RDS/cmd.version.s2.0.cmd.non.batch.corrected.scTransformed.rds")