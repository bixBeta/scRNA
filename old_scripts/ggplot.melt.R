melted = readRDS("/workdir/fa286/sandBox/ggplot/melted/melted.yx1.Rds")
head(melted)
summary(melted$value)
table(melted$value)
ggplot(data = melted, aes(x = value)) + geom_histogram(binwidth = 1)



# comp.RM1.lv = readRDS("lv-distances/comp.RM1.lv.Rds")
# comp.RM2.lv = readRDS("lv-distances/comp.RM2.lv.Rds")
# comp.YX1.lv = readRDS("lv-distances/comp.YX1.lv.Rds")
# comp.YX2.lv = readRDS("lv-distances/comp.YX2.lv.Rds")

load("lv-distances/lv.matrices.w.colnamesWLandRownamesWL.Rdata")



library(reshape)
melted.YX1 =  melt(comp.YX1.lv)
colnames(melted.YX1) = c("row.YX1", "col.WL", "value")
melted.YX2 =  melt(comp.YX2.lv)
colnames(melted.YX2) = c("row.YX2", "col.WL", "value")

melted.RM1 =  melt(comp.RM1.lv)
colnames(melted.RM1) = c("row.RM1", "col.WL", "value")
melted.RM2 =  melt(comp.RM2.lv)
colnames(melted.RM2) = c("row.RM2", "col.WL", "value")


head(melted.RM2)

summary(melted.YX1$value)
summary(melted.YX2$value)
summary(melted.RM1$value)
summary(melted.RM2$value)

df.yx1 = as.data.frame(table(melted.YX1$value))
colnames(df.yx1) = c("YX1.LV", "YX1.Freq")

df.yx2 = as.data.frame(table(melted.YX2$value))
colnames(df.yx2) = c("YX2.LV", "YX2.Freq")

df.rm1 = as.data.frame(table(melted.RM1$value))
colnames(df.rm1) = c("RM1.LV", "RM1.Freq")

df.rm2 = as.data.frame(table(melted.RM2$value))
colnames(df.rm2) = c("RM2.LV", "RM2.Freq")

lv.distribution = cbind(df.yx1,df.yx2,df.rm1,df.rm2)
mlv = melt(lv.distribution)


write.table(lv.distribution, "lv.distribution.txt", sep = "\t", col.names = T)

g1 = ggplot(data = melted.YX1, aes(x = value)) + geom_histogram(binwidth = 1) # Sturges approx ≈ 0.56, i am using 1 
g2 = ggplot(data = melted.YX2, aes(x = value)) + geom_histogram(binwidth = 1) 
g3 = ggplot(data = melted.RM1, aes(x = value)) + geom_histogram(binwidth = 1) 
g4 = ggplot(data = melted.RM2, aes(x = value)) + geom_histogram(binwidth = 1) 

