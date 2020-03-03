old.so = readRDS("/workdir/singleCellData/R_session_for_SEURAT/commandLine/filter-2-normalized_scTransformed_pbmc.rds")
DimPlot(old.so, reduction = "umap", group.by = "batch")

# all markers distibutions for GSEA
gsea.all.markers = FindAllMarkers(s1.2.RE, only.pos = F, min.pct = 0.2, logfc.threshold = 0, verbose = T, assay = "RNA")


for (i in 1:length(levels(Idents(s1.2.RE)))) {
  print(levels(Idents(s1.2.RE))[i])
  assign(x = paste0(levels(Idents(s1.2.RE))[i], ".gsea.case.v.ctrl"),
         value = FindMarkers(s1.2.RE, ident.1 = "case", group.by = "Phenotype",
                             subset.ident = levels(Idents(s1.2.RE))[i], 
                             only.pos = F,logfc.threshold = 0),
                            
         envir =.GlobalEnv  )
}


gsea.case.v.ctrl <- ls(pattern = ".gsea.case.v.ctrl")

for (i in 1:length(gsea.case.v.ctrl)) {
  write.table(mget(gsea.case.v.ctrl[i]), 
              file = paste0(gsea.case.v.ctrl[i], ".txt"), 
              sep = "\t", quote = F, col.names = NA)
}

