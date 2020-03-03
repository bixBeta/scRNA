s1.2.RE <- RenameIdents(s1.2, 
                        `0` = "0-T",
                        `2` = "2-T",
                        `4` = "4-T",
                        `5` = "5-T",
                        `11` = "11-T",
                        `21` = "21-T",
                        
                        `7` = "7-B",
                        `9` = "9-B",
                        
                        `6` = "6-CD8A",
                        `8` = "8-CD8A",
                        `13` = "13-CD8A",
                        `16` = "16-CD8A",
                        
                        `3` = "3-NK",
                        `12` = "12-NK",
                        `18` = "18-NK",
                        
                        
                        `17` = "17-Platelet",
                        `20` = "20-Platelet",
                        
                        `1` = "1-Mono",
                        `10` = "10-Mono",
                        `14` = "14-Mono",
                        `15` = "15-Mono",
                        `19` = "19-Mono",
                        `22` = "22-Mono",
                        `23` = "23-Mono",
                        `24` = "24-Mono"
)
#phenoday.cluster.averages.RE <- AverageExpression(s1.2.RE, return.seurat = TRUE, add.ident = "PhenoDay")
DoHeatmap(phenoday.cluster.averages.RE, features = goi, size = 5 ,draw.lines =T, slot = "scale.data")
DoHeatmap(phenoday.cluster.averages.RE, features = as.character(rownames(s1.2.RE)), size = 5 ,draw.lines =T, slot = "scale.data")

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#------------------------------------------------------------------------------------------------------------------------------------------------
# Cluster Specific Markers 

all.markers.s1.2.RE <- FindAllMarkers(s1.2.RE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


y <- list()

y.ft <- all.markers.s1.2.RE %>%
            filter(p_val_adj < 0.05)

write.table(y.ft, "Filtered.all.markers.s1.2.RE.txt", sep = "\t", quote = F, col.names = NA)


for (i in 1:length(levels(y.ft$cluster))) {
  f1 <- y.ft %>%
    filter(cluster == levels(y.ft$cluster)[i]) 
  y[[i]] <- f1$gene
  names(y)[i] <- levels(y.ft$cluster)[i]
}

upset(fromList(y), order.by = "freq",point.size = 3, 
      text.scale = 3, sets =levels(y.ft$cluster), keep.order = T)





