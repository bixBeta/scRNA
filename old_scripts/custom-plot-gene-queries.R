Idents(s1.2.RE) %>% levels()
nk3.subset = subset(s1.2.RE, idents = "3-NK")
nk3.subset$PhenoDay

nk3.subset.cavg  =  AverageExpression(nk3.subset, return.seurat = TRUE, add.ident = "cellID" )

nk3.pheno = nk3.subset@meta.data

nk3.subset <- AddMetaData(object = nk3.subset, metadata = rownames(nk3.subset@meta.data), col.name = "cellID")

DoHeatmap(nk3.subset.cavg, features = goi2, size = 5 ,draw.lines =T, slot = "scale.data")

goi2 = c("TOX", "PDCD1", "CD274")

DoHeatmap(phenoday.cluster.averages.RE, features = goi2, size = 5 ,draw.lines =T, slot = "scale.data")

# geneExpression <- function(gene){
# 
#   queryExpression <- s1.2.RE@assays$RNA@counts[rownames(s1.2.RE@assays$RNA@counts) == gene,]
#   query <- as.data.frame(t(queryExpression))
#   query$UMI <- rownames(query)
#   jnd <- left_join(localMeta, query, by = "UMI")
#   jnd <- jnd %>%
#     group_by(PhenoDay, seurat_clusters) %>%
#     summarise_each(funs = (mean), !!gene)
#   #print(jnd)
# }
# gQuery <- as.data.frame(geneExpression("TOX"))
# 

phenoday.cluster.averages.RE@meta.data

qAvgs = as.data.frame(phenoday.cluster.averages.RE@assays$RNA@data)
q2plot = qAvgs[rownames(qAvgs) %in% goi2,]

library(reshape)
q2plot.m = melt(q2plot)

tempslot = as.character(q2plot.m$variable)

m <- vector("integer", length(tempslot))
mn <- vector("integer", length(tempslot))
for (i in 1:length(tempslot)) {

  mn[i] = strsplit(x =tempslot[i], split = "_")[[1]][1]
  
}

q2plot.m$PhenoDay = m
q2plot.m$Clust = mn

library(viridis)
library(ggplot2)
library(ggbeeswarm)

q2ploted = ggplot(q2plot.m, aes(Clust, value, fill=PhenoDay)) +
                geom_boxplot()+
                geom_quasirandom(aes(color= PhenoDay),alpha = 0.75,dodge.width=0.75)+
                theme_bw() + scale_color_manual(values=c(viridis(4)))+
                scale_fill_manual(values=c(viridis(4)))
q2ploted = q2ploted + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("filename.png", plot = q2ploted,
       scale = 1,
       dpi = 200)

write.table(q2plot.m, "q2plot.m.txt", sep = "\t", quote = F, col.names = NA)
write.table(q2plot, "q2plot.txt", sep = "\t", quote = F)



DimPlot(s1.2.rename)
