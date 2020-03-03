db 
select(as.data.frame(get(db[1])), contains("p_val_adj")) %>% filter(p_val_adj < 0.05)

library(tidyverse) 

for (i in 1:length(db)) {
  print(db[1])
  assign(x = paste0("filtered.0.05.",db[i]),
        value = 
          as.data.frame(get(db[i])) %>%
          rownames_to_column("gene")  %>%
          filter( p_val_adj < 0.05) %>%
          column_to_rownames("gene"),
        
        envir = .GlobalEnv)
  
}

dbf <- ls(pattern = "filtered.0.05")
for (i in 1:length(dbf)) {
  write.table(mget(dbf[i]), file = paste0(dbf[i], ".txt"), sep = "\t", quote = F, col.names = NA)
  
}

avg.expression.per.cluster.hSC.base = AverageExpression(s1.2.RE, return.seurat = TRUE)
avg.expression.per.cluster.hSC.base.df  = as.data.frame(avg.expression.per.cluster.hSC.base@assays$RNA[])

write.table(avg.expression.per.cluster.hSC.base.df, file = "avg.expression.per.cluster.hSC.base.df.txt", col.names = NA,  quote = F, sep = "\t")


scaled.avg.expression.per.cluster.hSC.base = AverageExpression(s1.2.RE, use.scale = T)
