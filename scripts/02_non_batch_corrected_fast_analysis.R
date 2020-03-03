library("Seurat")
library("dplyr")
library("future")
library("tibble")
plan("multiprocess", workers = 4)
plan()
####################################################################################
####################################################################################
all_pbmcs = readRDS("MEGA_MECFS_merged_seurat_object.rds")

# extract phenodata from all_pbmcs
raw.pheno.data <- all_pbmcs@meta.data
table(raw.pheno.data$orig.ident == "2250")

# correct 2250 to 2350
raw.pheno.data = raw.pheno.data %>%
  rownames_to_column(., "rownames") %>%
  mutate(orig.ident = replace(orig.ident, which(orig.ident == "2250"), "2350")) %>%
  column_to_rownames("rownames")

# check 
table(raw.pheno.data$orig.ident == "2250") 
table(raw.pheno.data$orig.ident == "2350")

write.table(raw.pheno.data, "metadata/raw.phenoData.txt", sep = "\t", quote = F, col.names = NA)

####################################################################################
####################################################################################
# combine Liz's and raw metadata and update the Seurat Object

lizs.pheno <- read.table("metadata/liz_phenoData.txt", header = T)
lizs.pheno <- 
  lizs.pheno %>%
  mutate_each(funs = (as.character), library_id)


raw.pheno.data$cell.barcode = rownames(raw.pheno.data)
pheno1 = left_join(raw.pheno.data, lizs.pheno, by = c("orig.ident" = "library_id"))

pheno1 = pheno1 %>% select(4,everything())
rownames(pheno1) = pheno1$cell.barcode
pheno1 = pheno1 %>% select(-1)

str(pheno1)

all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$batch, col.name = "batch")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$ENID, col.name = "ENID")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$Phenotype, col.name = "Phenotype")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$CPET.Day, col.name = "CPET.Day")


saveRDS(all_pbmcs, file = "MEGA_MECFS_CORRECTED_META_DATA_merged_seurat_object.rds")



####################################################################################
####################################################################################
# load seurat object
all_pbmcs <- readRDS("RDS/MEGA_MECFS_CORRECTED_META_DATA_merged_seurat_object.rds")


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


#***********************************************************************************
# save filtered object 
# use this for all analysis
saveRDS(s2.0, "s2.0.Filtered.SeuratObject.rds")


####################################################################################
####################################################################################
# Add number of genes per UMI for each cell to metadata
s2.0$log10GenesPerUMI <- log10(s2.0$nFeature_RNA) / log10(s2.0$nCount_RNA)

# Compute percent mito ratio
s2.0$mitoRatio <- PercentageFeatureSet(object = s2.0, pattern = "^MT-")
s2.0$mitoRatio <- s2.0@meta.data$mitoRatio / 100


metadata <- s2.0@meta.data


# Visualize the number of cell counts per sample
library(ggplot2)
library(viridis)
library(plotly)
library(forcats)

ggplotly(metadata %>% 
           ggplot(aes(x=forcats::fct_infreq(orig.ident), fill=batch)) + 
           geom_bar() + scale_fill_manual(values=c(viridis(13))) +
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
           theme(plot.title = element_text(hjust=0.5, face="bold")) +
           ggtitle("NCells"))





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