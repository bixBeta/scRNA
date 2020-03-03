#------------------------------------------------------------------------------------------------------------------------------------------------
# DGEs

#------------------------------------------------------------------------------------------------------------------------------------------------
# Within Cluster PhenoDay Changes
for (i in 1:length(levels(Idents(s1.2.rename)))) {
  print(levels(Idents(s1.2.rename))[i])
  assign(x = paste0(levels(Idents(s1.2.rename))[i], ".s1.2.PhenoDay"),
         value = FindMarkers(s1.2.rename, ident.1 = "D1-case", 
                             group.by = "PhenoDay", subset.ident = levels(Idents(s1.2.rename))[i], 
                             only.pos = F),
         envir =.GlobalEnv  )
}

phenoDay <- list()
for (i in 1:length(ls(pattern="s1.2.PhenoDay"))) { 
  phenoDay[[i]] <- mget(ls(pattern="s1.2.PhenoDay")[i])
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster CaseD1 v. CaseD2
for (i in 1:length(levels(Idents(s1.2.rename)))) {
  print(levels(Idents(s1.2.rename))[i])
  assign(x = paste0(levels(Idents(s1.2.rename))[i], ".s1.2.CaseD1.v.CaseD2"),
         value = FindMarkers(s1.2.rename, ident.1 = "D1-case", 
                             ident.2 = "D2-case", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.rename))[i], 
                             only.pos = F),
         envir =.GlobalEnv  )
}

caseD1D2 <- list()
for (i in 1:length(ls(pattern="CaseD1.v.CaseD2"))) { 
  caseD1D2[[i]] <- mget(ls(pattern="CaseD1.v.CaseD2")[i])
}


  #------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster ControlD1 v. ControlD2
for (i in 1:length(levels(Idents(s1.2.rename)))) {
  print(levels(Idents(s1.2.rename))[i])
  assign(x = paste0(levels(Idents(s1.2.rename))[i], ".s1.2.CtrlD1.v.CtrlD2"),
         value = FindMarkers(s1.2.rename, ident.1 = "D1-control", 
                             ident.2 = "D2-control", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.rename))[i], 
                             only.pos = F,logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}

ctrlD1D2 <- list()
for (i in 1:length(ls(pattern="CtrlD1.v.CtrlD2"))) { 
  ctrlD1D2[[i]] <- mget(ls(pattern="CtrlD1.v.CtrlD2")[i])
}

#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster Case v. Control in D1
for (i in 1:length(levels(Idents(s1.2.rename)))) {
  print(levels(Idents(s1.2.rename))[i])
  assign(x = paste0(levels(Idents(s1.2.rename))[i], ".s1.2.CaseD1.v.ControlD1"),
         value = FindMarkers(s1.2.rename, ident.1 = "D1-case", 
                             ident.2 = "D1-control", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.rename))[i], 
                             only.pos = F,logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}


caseD1ctrlD1 <- list()
for (i in 1:length(ls(pattern=".s1.2.CaseD1.v.ControlD1"))) { 
  caseD1ctrlD1[[i]] <- mget(ls(pattern=".s1.2.CaseD1.v.ControlD1")[i])
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster Case v. Control in D2
for (i in 1:length(levels(Idents(s1.2.rename)))) {
  print(levels(Idents(s1.2.rename))[i])
  assign(x = paste0(levels(Idents(s1.2.rename))[i], ".s1.2.CaseD2.v.CtrlD2"),
         value = FindMarkers(s1.2.rename, ident.1 = "D2-case", 
                             ident.2 = "D2-control", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.rename))[i], 
                             only.pos = F,logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}

caseD2ctrlD2 <- list()
for (i in 1:length(ls(pattern=".s1.2.CaseD2.v.CtrlD2"))) { 
  caseD2ctrlD2[[i]] <- mget(ls(pattern=".s1.2.CaseD2.v.CtrlD2")[i])
}

#------------------------------------------------------------------------------------------------------------------------------------------------
p1 <- objects(pattern = ".s1.2.DE.Mar")
p2 <- objects(pattern = ".s1.2.Marker")

for (i in 1:length(p1)) {
  print(p1[i])
  print(goi %in% rownames(p1[i]))
  
}


da <- ls(pattern = "s1.2.C")
for (i in 1:length(da)) {
  write.table(mget(da[i]), file = paste0(da[i], ".txt"), sep = "\t", quote = F, col.names = NA)
  
}

#------------------------------------------------------------------------------------------------------------------------------------------------
