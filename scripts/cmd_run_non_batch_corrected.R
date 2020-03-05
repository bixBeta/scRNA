library("Seurat")
library("ggplot2")

# load seurat object
s2.0.cmd <- readRDS("RDS/s2.0.Filtered.SeuratObject.rds")

# fast analysis 
options(future.globals.maxSize = 9000 * 1024^2)

# run sctransform
s2.0 <- SCTransform(s2.0, vars.to.regress = c("percent.mt", "nCount_RNA"), verbose = T)

# Perform dimensionality reduction by PCA and UMAP embedding
# These are now standard steps in the Seurat workflow for visualization and clustering

s2.0 <- RunPCA(s2.0, verbose = FALSE)
s2.0 <- RunUMAP(s2.0, dims = 1:30, verbose = FALSE)
s2.0 <- FindNeighbors(s2.0, dims = 1:30, verbose = FALSE)
s2.0 <- FindClusters(s2.0, verbose = FALSE)

png(filename = "plots/cmd_version_non_batch_corrected_clustering.png", width = 4000, height = 2000)
DimPlot(s2.0, label = TRUE) + NoLegend()
dev.off()

saveRDS(s2.0.cmd, file = "RDS/cmd.version.s2.0.non.batch.corrected.scTransformed.rds")l