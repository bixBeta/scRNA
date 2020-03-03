library(sleepwalk)

sleepwalk( s1.2.rename$umap@cell.embeddings, s1.2.rename$pca@cell.embeddings )


Idents(s1.2.rename)

sce.s1.2.rename = as.SingleCellExperiment(s1.2.rename)
sce.s1.2.rename@assays$data$counts
s1.2.rename@assays$RNA@



new_cell_data_set(expression_data =sce.s1.2.rename@assays$data$counts, cell_metadata = s1.2.rename@meta.data )


pct  = s1.2.rename@reductions$pca@stdev / sum(s1.2.rename@reductions$pca@stdev) * 100

# Calculate cumulative percents for each PC
cum <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]

co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.

co2
# Minimum of the two calculation
pcs <- min(co1, co2) # change to any other number

pcs
