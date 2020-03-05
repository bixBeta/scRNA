library("Seurat")
library("dplyr")
library("future")
library("tibble")
plan("multiprocess", workers = 4)
plan()
####################################################################################
####################################################################################
all_pbmcs = readRDS("RDS/MEGA_MECFS_merged_seurat_object.rds")

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
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$orig.ident, col.name = "orig.ident")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$batch, col.name = "batch")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$ENID, col.name = "ENID")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$Phenotype, col.name = "Phenotype")
all_pbmcs <- AddMetaData(object = all_pbmcs, metadata = pheno1$CPET.Day, col.name = "CPET.Day")


saveRDS(all_pbmcs, file = "RDS/MEGA_MECFS_CORRECTED_META_DATA_merged_seurat_object.rds")



####################################################################################
####################################################################################
# load seurat object
all_pbmcs <- readRDS("RDS/MEGA_MECFS_CORRECTED_META_DATA_merged_seurat_object.rds")


####################################################################################
####################################################################################
# store mitochondrial percentage in object meta data
all_pbmcs <- PercentageFeatureSet(all_pbmcs, pattern = "^MT-", col.name = "percent.mt")

# apparently my method of just updating the metadata did not work :(  
all_pbmcs <- RenameIdents(all_pbmcs, `2250` = "2350")

# plot raw distributions
png(filename = "plots/raw_distributions.png", width = 4200, height = 2200)
VlnPlot(all_pbmcs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

png(filename = "plots/featureScatter_Plots.png", width = 1920, height = 1080)
plot1 <- FeatureScatter(all_pbmcs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all_pbmcs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

####################################################################################
####################################################################################
# filter using distributions 
s2.0 <- subset(all_pbmcs, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)

png(filename = "plots/trimmed_raw_distributions.png", width = 4200, height = 2200)
VlnPlot(s2.0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

png(filename = "plots/filtered_featureScatter_Plots.png", width = 1920, height = 1080)
plot3 <- FeatureScatter(s2.0, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(s2.0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot3, plot4))
dev.off()

#***********************************************************************************
# save filtered object 
# use this for all analysis
saveRDS(s2.0, "RDS/s2.0.Filtered.SeuratObject.rds")
#***********************************************************************************

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


metadata$Phenoday <- as.factor(paste0(metadata$Phenotype, metadata$CPET.Day))


png(filename = "plots/number_of_cells_per_phenoday.png", width = 800, height = 600)
metadata %>% 
           ggplot(aes(x=forcats::fct_infreq(Phenoday), fill=Phenoday)) + 
           geom_bar() + scale_fill_manual(values=c(viridis(4))) +
           theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
           theme(plot.title = element_text(hjust=0.5, face="bold")) +
           ggtitle("NCells")
dev.off()


# Visualize the number UMIs/transcripts per cell
png(filename = "plots/number_of_transcripts_per_cell.png", width = 1000, height = 800)
metadata %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()



# Visualize the distribution of genes detected per cell/sample via histogram
metadata %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell/sample via boxplot
png(filename = "plots/number_of_features_per_sample.png", width = 2600, height = 1200)
metadata %>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")
dev.off()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Phenoday)


raw.counts.s2.0 <- GetAssayData(object = s2.0, slot = "counts")

####################################################################################
####################################################################################
# normalization and cell cycle scoring
library(RCurl)
library(cowplot)

seurat_phase <- NormalizeData(s2.0)

# Load cell cycle markers
load("RDS/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)                                

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
png(filename = "plots/cell_cycle_effect.png", width = 1920, height = 1080)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
dev.off()
####################################################################################
####################################################################################
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
png(filename = "plots/non_batch_corrected_clustering.png", width = 1920, height = 1080)
DimPlot(s2.0, label = TRUE) + NoLegend()
dev.off()

