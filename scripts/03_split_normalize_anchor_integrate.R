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