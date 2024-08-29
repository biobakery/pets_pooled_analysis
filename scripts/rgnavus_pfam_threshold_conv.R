pfam_table = as.data.frame(readr::read_tsv("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/g__Blautia.s__Ruminococcus_gnavus.tsv"))
rownames(pfam_table) = pfam_table$"# Gene Family"
pfam_table$"# Gene Family" = NULL

pfam_table[ pfam_table<0.0000001 ] <- 0
pfam_table[ pfam_table>=0.0000001 ] <- 1
write.table(pfam_table, file='/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/g__Blautia.s__Ruminococcus_gnavus_binary.tsv', quote=FALSE, sep='\t')

errors = load('/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/errors.RData')

