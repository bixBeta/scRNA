levels(Idents(s1.2.RE))
saveRDS(s1.2.RE, file = "s1.2.RE_hyphenated_w_superCluster.RData")

#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster CaseD1 v. CaseD2
for (i in 1:length(levels(Idents(s1.2.RE)))) {
  print(levels(Idents(s1.2.RE))[i])
  assign(x = paste0(levels(Idents(s1.2.RE))[i], ".s1.2.caseD1.v.caseD2"),
         value = FindMarkers(s1.2.RE, ident.1 = "D1-case", 
                             ident.2 = "D2-case", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.RE))[i], 
                             only.pos = F, logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}

# hyphenated super cluster (s1.2.RE)
caseD1D2_hSC <- list()
for (i in 1:length(ls(pattern="caseD1.v.caseD2"))) { 
  caseD1D2_hSC[[i]] <- mget(ls(pattern="caseD1.v.caseD2")[i])
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster ControlD1 v. ControlD2
for (i in 1:length(levels(Idents(s1.2.RE)))) {
  print(levels(Idents(s1.2.RE))[i])
  assign(x = paste0(levels(Idents(s1.2.RE))[i], ".s1.2.ctrlD1.v.ctrlD2"),
         value = FindMarkers(s1.2.RE, ident.1 = "D1-control", 
                             ident.2 = "D2-control", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.RE))[i], 
                             only.pos = F,logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}

ctrlD1D2_hSC <- list()
for (i in 1:length(ls(pattern="ctrlD1.v.ctrlD2"))) { 
  ctrlD1D2_hSC[[i]] <- mget(ls(pattern="ctrlD1.v.ctrlD2")[i])
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster Case v. Control in D1
for (i in 1:length(levels(Idents(s1.2.RE)))) {
  print(levels(Idents(s1.2.RE))[i])
  assign(x = paste0(levels(Idents(s1.2.RE))[i], ".s1.2.caseD1.v.controlD1"),
         value = FindMarkers(s1.2.RE, ident.1 = "D1-case", 
                             ident.2 = "D1-control", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.RE))[i], 
                             only.pos = F,logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}

caseD1ctrlD1_hSC <- list()
for (i in 1:length(ls(pattern=".s1.2.caseD1.v.controlD1"))) { 
  caseD1ctrlD1_hSC[[i]] <- mget(ls(pattern=".s1.2.caseD1.v.controlD1")[i])
}


#------------------------------------------------------------------------------------------------------------------------------------------------
# Within CLuster Case v. Control in D2
for (i in 1:length(levels(Idents(s1.2.RE)))) {
  print(levels(Idents(s1.2.RE))[i])
  assign(x = paste0(levels(Idents(s1.2.RE))[i], ".s1.2.caseD2.v.ctrlD2"),
         value = FindMarkers(s1.2.RE, ident.1 = "D2-case", 
                             ident.2 = "D2-control", group.by = "PhenoDay",
                             subset.ident = levels(Idents(s1.2.RE))[i], 
                             only.pos = F,logfc.threshold = 0.1),
         envir =.GlobalEnv  )
}

caseD2ctrlD2_hSC <- list()
for (i in 1:length(ls(pattern=".s1.2.caseD2.v.ctrlD2"))) { 
  caseD2ctrlD2_hSC[[i]] <- mget(ls(pattern=".s1.2.caseD2.v.ctrlD2")[i])
}


db <- ls(pattern = "s1.2.c")
for (i in 1:length(db)) {
  write.table(mget(db[i]), file = paste0(db[i], ".txt"), sep = "\t", quote = F, col.names = NA)
  
}

