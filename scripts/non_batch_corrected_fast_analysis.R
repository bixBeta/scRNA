####################################################################################
####################################################################################
# load seurat object
all_pbmcs <- readRDS("MEGA_MECFS_CORRECTED_META_DATA_merged_seurat_object.rds")


####################################################################################
####################################################################################
# store mitochondrial percentage in object meta data
all_pbmcs <- PercentageFeatureSet(all_pbmcs, pattern = "^MT-", col.name = "percent.mt")

# plot raw distributions
VlnPlot(all_pbmcs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1 <- FeatureScatter(all_pbmcs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_pbmcs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

####################################################################################
####################################################################################
# filter using distributions 
s2.0 <- subset(all_pbmcs, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
plot3 <- FeatureScatter(s2.0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(s2.0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))

# save filtered object 
# use this for all analysis
saveRDS(s2.0, "s2.0.Filtered.SeuratObject.rds")

####################################################################################
####################################################################################
# fast analysis 

# run sctransform
s2.0 <- SCTransform(s2.0, vars.to.regress = c("percent.mt", "ncount_RNA"), verbose = FALSE)

# Perform dimensionality reduction by PCA and UMAP embedding

# These are now standard steps in the Seurat workflow for visualization and clustering
s2.0 <- RunPCA(s2.0, verbose = FALSE)
s2.0 <- RunUMAP(s2.0, dims = 1:30, verbose = FALSE)

s2.0 <- FindNeighbors(s2.0, dims = 1:30, verbose = FALSE)
s2.0 <- FindClusters(s2.0, verbose = FALSE)
DimPlot(s2.0, label = TRUE) + NoLegend()