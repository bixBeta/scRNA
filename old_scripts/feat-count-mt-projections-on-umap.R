# Determine metrics to plot present in seurat_control@meta.data
metrics <-  c("nCount_RNA", "nFeature_RNA","percent.mt")

# Extract the UMAP coordinates for each cell and include information about the metrics to plot
qc_data <- FetchData(s1.2.rename, 
                     vars = c(metrics, "seurat_clusters", "UMAP_1", "UMAP_2"))

umap_label <- FetchData(s1.2.rename, 
                        vars = c("seurat_clusters", "UMAP_1", "UMAP_2"))  %>%
  as.data.frame() %>% 
  group_by(orig.ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))


# Plot a UMAP plot for each metric
map(metrics, function(qc){
  ggplot(qc_data,
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=qc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=orig.ident, x, y)) +
    ggtitle(qc)
}) %>%
  plot_grid(plotlist = .)


