annot = as.data.frame(data.table::fread("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/pfams_annotation.tsv", header=T, check.names = F,fill = TRUE))
annot = annot[!duplicated(annot$'#'), ]
