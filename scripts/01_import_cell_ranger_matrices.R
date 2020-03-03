library(Seurat)
library(dplyr)

cellranger.matrix.path <- read.table(file = "metadata/all_samples_file_paths_scRNA.txt", sep = "\t", header = T)

####################################################################################
####################################################################################
# read in filtered matrices for all samples 

for (i in 1:nrow(cellranger.matrix.path)) {
  assign(x = paste0("RD10x_",cellranger.matrix.path$sample_id[i]), 
         Read10X(data.dir = as.character(cellranger.matrix.path$filtered_matrix[i])), 
         envir = .GlobalEnv)
}

####################################################################################
####################################################################################
# create seurat objects from imported matrices

ten_x_pattern <- objects(pattern = "RD10x_")

for (i in 1:length(ten_x_pattern)) {
  assign(x = paste0("SObj_", strsplit(ten_x_pattern[i], "_")[[1]][2]), 
         value = CreateSeuratObject(counts = get(ten_x_pattern[i]), 
                                    project = strsplit(ten_x_pattern[i], "_")[[1]][2] , 
                                    min.cells = 3, min.features = 200),
         envir = .GlobalEnv)
}

####################################################################################
####################################################################################
# merge all seurat objects into one main object

sobj_pattern <- noquote(objects(pattern = "SObj_"))
all_pbmcs <- merge(get(sobj_pattern[1]), 
                   y = mget(sobj_pattern[2:length(sobj_pattern)]))

saveRDS(all_pbmcs, file = "MEGA_MECFS_merged_seurat_object.rds")
