# IL2RB (Il2rb)
# SCART1 (Cd163l1)
# IL18RAP (Il18rap)
# STC2 (Stc2)
# AQP9 (Aqp9)
# GAS7 (Gas7)
# GZMM (Gzmm)
# TOX
# PDCD1
# CD274

goi <- c("IL2RB",
"SCART1",
"IL18RAP",
"AQP9",
"GAS7",
"GZMM",
"TOX",
"PDCD1",
"CD274")


#patient.cluster.averages <- AverageExpression(s1.2, return.seurat = TRUE, add.ident = "ENID")
DoHeatmap(cluster.averages, features = goi, size = 5 ,draw.lines =T, slot = "scale.data")
VlnPlot(cluster.averages, features = goi,pt.size = 0)


for (i in 0:24) {
  assign(x = paste0("cluster", i , ".s1.2.Markers"),
         value = FindMarkers(s1.2, ident.1 = "D1-case", 
                            group.by = "PhenoDay", subset.ident = i, 
                            only.pos = F, logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}


FeaturePlot(s1.2,goi,pt.size = 1)




s1.2.rename <- RenameIdents(s1.2, 
                            `0` = "T-Cells",
                            `2` = "T-Cells",
                            `3` = "NK-Cells",
                            `4` = "T-Cells",
                            `5` = "T-Cells",
                            `6` = "CD8A",
                            `7` = "B-Cells",
                            `8` = "CD8A",
                            `9` = "B-Cells",
                            `11` = "T-Cells",
                            `12` = "NK-Cells",
                            `13` = "CD8A",
                            `16` = "CD8A",
                            `17` = "Platelet",
                            `18` = "NK-Cells",
                            `20` = "Platelet",
                            `21` = "T-Cells",
                            `1` = "Monocytes",
                            `10` = "Monocytes",
                            `14` = "Monocytes",
                            `15` = "Monocytes",
                            `19` = "Monocytes",
                            `22` = "Monocytes",
                            `23` = "Monocytes",
                            `24` = "Monocytes"
                            
)

#phenoday.cluster.averages.renamed <- AverageExpression(s1.2.rename, return.seurat = TRUE, add.ident = "PhenoDay")

DoHeatmap(patient.cluster.averages, features = goi, size = 5 ,draw.lines =T, slot = "scale.data")
DoHeatmap(phenoday.cluster.averages.renamed, features = goi, size = 5 ,draw.lines =T, slot = "scale.data")


# Number of Cells Per PhenoDay Per Cluster 
x <- as.data.frame(table(s1.2.rename@active.ident, s1.2.rename@meta.data$PhenoDay))
colnames(x) <- c("clusterIdent", "PhenoDay", "Freq")
library(dplyr)
x <- x %>% arrange(Freq)

ggplot(data = x, mapping = aes(x = clusterIdent, y = Freq)) +
  geom_line(color="Black", alpha = 0.3 )+ 
  geom_point(aes(color=PhenoDay),alpha = 0.5) + theme_minimal() + 
  scale_color_brewer(palette = "Set1") +
  xlab("Cluster Ident")+ylab("Cell Counts")+
  ggtitle(" Number of cells per cluster per phenoDay")






