#This code was written by Tobyn Branck and provides the AMR analysis and Figure 6 visualization

library("RNOmni")
library("ggplot2")
library("grid")
library("vegan")
library(tidyverse)
library(lme4)
library(factoextra)
library(ade4)
library(FactoMineR)
library("plyr")
library("dbplyr")
library(scales)
library(ggbreak)
library(grid)
library(gtable)
library(gridExtra)
library(scatterplot3d)
library(reshape2)
library(gtools)

#set the file path for input & output files
input_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs")
output_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs/AMR_section")

### import metadata
metadata = as.data.frame(data.table::fread(file.path(input_dir,"metadata.csv"), header=T, check.names = F))
rownames(metadata) = metadata$sample_id_metaphlan

### import AMR genes families (uniref90s that mapped to CARD database) ###
AMR_genes = as.data.frame(data.table::fread(file.path(input_dir,"gene_fams_relab_final9.9_unstratified_CARD_headers_rename.tsv"), header=T, check.names = F))
rownames(AMR_genes) = AMR_genes$'# Gene Family'
AMR_genes$'# Gene Family' = NULL
colnames(AMR_genes) = gsub("_Abundance-RPKs","",colnames(AMR_genes))

dog_samps = rownames(metadata)[metadata$species=="dog"]
cat_samps = rownames(metadata)[metadata$species=="cat"]
human_samps = rownames(metadata)[metadata$species=="human"]
hmp1_II_samps = rownames(metadata)[metadata$study_readable=="HMP1-II"]
madagascar_samps = rownames(metadata)[metadata$study_readable=="Madagascar"]

# remove madagascar samples
metadata = metadata[!metadata$study_readable=="Madagascar",]
AMR_genes = AMR_genes[,colnames(AMR_genes) %in% rownames(metadata)]

dog_AMR = AMR_genes[,colnames(AMR_genes) %in% dog_samps]
cat_AMR = AMR_genes[,colnames(AMR_genes) %in% cat_samps]
hmp_AMR = AMR_genes[,colnames(AMR_genes) %in% hmp1_II_samps]

#get prevalence value for "must be in at least 3 samples" for each host
dog_prev = 3/dim(dog_AMR)[2]
cat_prev = 3/dim(cat_AMR)[2]
hmp_prev = 3/dim(hmp_AMR)[2]

#sort genes by highest prevalence
dog_AMR$prev = rowSums(dog_AMR > 0)/dim(dog_AMR)[2]
dog_AMR = dog_AMR[order(dog_AMR$prev,decreasing=TRUE),]
dog_AMR = dog_AMR[dog_AMR$prev>dog_prev,]

cat_AMR$prev = rowSums(cat_AMR > 0)/dim(cat_AMR)[2]
cat_AMR = cat_AMR[order(cat_AMR$prev,decreasing=TRUE),]
cat_AMR = cat_AMR[cat_AMR$prev>cat_prev,]

hmp_AMR$prev = rowSums(hmp_AMR > 0)/dim(hmp_AMR)[2]
hmp_AMR = hmp_AMR[order(hmp_AMR$prev,decreasing=TRUE),]
hmp_AMR = hmp_AMR[hmp_AMR$prev>hmp_prev,]

amr_genes_keep = union(rownames(cat_AMR),rownames(dog_AMR))
amr_genes_keep = union(amr_genes_keep,rownames(hmp_AMR))

AMR_genes_filtered = AMR_genes[rownames(AMR_genes) %in% amr_genes_keep,]

AMR_genes_filtered = as.data.frame(t(AMR_genes_filtered))
AMR_genes_filtered = AMR_genes_filtered[rowSums(AMR_genes_filtered) != 0, ]
metadata = metadata[rownames(metadata) %in% rownames(AMR_genes_filtered),]
metadata = metadata[mixedorder(rownames(metadata)),]
AMR_genes_filtered = AMR_genes_filtered[mixedorder(rownames(AMR_genes_filtered)),]

metadata_for_amr_adonis = metadata[rownames(metadata) %in% rownames(AMR_genes_filtered),]
AMR_genes_filtered = AMR_genes_filtered[gtools::mixedorder(rownames(AMR_genes_filtered)),]
amr_adonis = adonis(AMR_genes_filtered ~ metadata_for_amr_adonis$species, data = metadata_for_amr_adonis, method = "bray",na.rm=TRUE)

AMR_genes_filtered_maas = as.data.frame(t(AMR_genes_filtered))

write.csv(AMR_genes_filtered_maas,file.path(output_dir,"S13_ARGs_abundance_table.csv"))

library(Maaslin2)
Maaslin2(AMR_genes_filtered_maas,metadata_for_amr_adonis,output = file.path(output_dir,"maaslin_pets_HostsStrat_hmpRef_CARD_AMR_Noprevfilter"),fixed_effects = c("species"),reference=c("species,human"), random_effects = c("study_readable"),min_prevalence = 0)
Maaslin2(AMR_genes_filtered_maas,metadata_for_amr_adonis,output = file.path(output_dir,"maaslin_pets_HostsStrat_dogRef_CARD_AMR_Noprevfilter"),fixed_effects = c("species"),reference=c("species,dog"), random_effects = c("study_readable"),min_prevalence = 0)
Maaslin2(AMR_genes_filtered_maas,metadata_for_amr_adonis,output = file.path(output_dir,"maaslin_pets_HostsStrat_catRef_CARD_AMR_Noprevfilter"),fixed_effects = c("species"),reference=c("species,cat"), random_effects = c("study_readable"),min_prevalence = 0)

## PcoA
# Bray-Curtis dissimilarity
all = data.frame(AMR_genes_filtered, metadata)
amr_dist = vegdist(AMR_genes_filtered, method = "bray")

weights = c(case_when(metadata$species == "human" ~ 1/238,
                      metadata$species == "cat" ~ 1/367,
                      metadata$species == "dog" ~ 1/2272)) # weighted by sample size

pcoa = wcmdscale(amr_dist, k = 2, eig = TRUE, w = weights)

positions = pcoa$points
colnames(positions) = c("pcoa1","pcoa2")

100*pcoa$eig / sum(pcoa$eig)
percent_var_explained = metagMisc::eig_perc(pcoa$eig, positive = T, plot = F)
percent_var_explained = format(round(percent_var_explained[1:2],digits = 2),nsmall=1,trim=TRUE)
labs = c(glue::glue("PCoA1 [{percent_var_explained[1]}%]"),
         glue::glue("PCoA2 [{percent_var_explained[2]}%]"))
positions = as.data.frame(positions)

metadata_ord = metadata[rownames(metadata) %in% rownames(positions),]
metadata_ord$host_stratified = NA
metadata_ord$host_stratified[metadata_ord$study_readable == 'HMP1-II'] = 'HMP1-II'
metadata_ord$host_stratified[metadata_ord$species == 'dog'] = 'dog'
metadata_ord$host_stratified[metadata_ord$species == 'cat'] = 'cat'

positions$host = metadata_ord$species[rownames(positions) %in% rownames(metadata_ord)]
positions$host2 = metadata_ord$host_stratified[rownames(positions) %in% rownames(metadata_ord)]

# calculate weighted average scores of genes driving the cluster of PCoA
wascores = data.frame(wascores(positions[1:2], AMR_genes_filtered))

# maaslin results
uniref_catRef = as.data.frame(data.table::fread(file.path(output_dir,"maaslin_pets_HostsStrat_catRef_CARD_AMR_Noprevfilter/significant_results.tsv"), header=T, check.names = F)) %>% mutate(reference = "cat")
uniref_dogRef = as.data.frame(data.table::fread(file.path(output_dir,"maaslin_pets_HostsStrat_dogRef_CARD_AMR_Noprevfilter/significant_results.tsv"), header=T, check.names = F)) %>% mutate(reference = "dog")
uniref_hmpRef = as.data.frame(data.table::fread(file.path(output_dir,"maaslin_pets_HostsStrat_hmpRef_CARD_AMR_Noprevfilter/significant_results.tsv"), header=T, check.names = F)) %>% mutate(reference = "HMP1-II")

# top 4 coefficients genes among the significant results in each host species
uniref_catRef_dog_vs_cat = uniref_catRef[order(abs(uniref_catRef$coef), decreasing = TRUE), ] %>% filter(value == "dog")
uniref_hmpRef_dog_vs_human = uniref_hmpRef[order(abs(uniref_hmpRef$coef), decreasing = TRUE), ] %>% filter(value == "dog")
uniref_hmpRef_cat_vs_human = uniref_hmpRef[order(abs(uniref_hmpRef$coef), decreasing = TRUE), ] %>% filter(value == "cat")

maaslin_sig_list = union(uniref_catRef_dog_vs_cat[1:4, 1], uniref_hmpRef_dog_vs_human[1:4, 1]) %>% 
  union(uniref_hmpRef_cat_vs_human[1:4, 1])

# reformat the maaslin genes annotation
maaslin_sig_list = sub("\\.\\.", ": ", maaslin_sig_list)
maaslin_sig_list = gsub("\\.", " ", maaslin_sig_list)
maaslin_sig_list = gsub("\\:.*", "", maaslin_sig_list)

# grep lines containing masslin_sig_list
matching_rows = sapply(rownames(wascores), function(row_name) { 
  any(sapply(maaslin_sig_list, function(pattern) {
    grepl(pattern, row_name)
  }))
})

# subset the dataframe based on the rows that matched
data_wa_sub = wascores[matching_rows, ]
data_wa_sub2 = wascores[c('UniRef90_A0A2R4H2Z9: rRNA methylase',
                          'UniRef90_D6PYN7: DfrA14 (Fragment)',
                          'UniRef90_A0A127SV09: Chloramphenicol acetyltransferase',
                          'UniRef90_A0A1Y1CU85: Aminoglycoside N(6&apos )-acetyltransferase',
                          'UniRef90_A0A174XA02: Beta-lactamase',
                          'UniRef90_P0A0C2: Bifunctional AAC/APH'), ] # selected by visualization


data_wa_sub3 = rbind(data_wa_sub, data_wa_sub2) # combined 
rownames(data_wa_sub3) = gsub("UniRef90_", "", rownames(data_wa_sub3))

cairo_pdf("CARDAMR_ord_Allhosts.pdf", width = 8, height = 6)

ggplot() + 
  geom_point(data = positions %>% as_tibble(rownames = "samples",eig=TRUE), 
             aes(x = pcoa1, y = pcoa2, color = host), size = 4, alpha = 0.5, show.legend = F) + 
  geom_text(data = data_wa_sub3, aes(x = pcoa1, y = pcoa2, label = rownames(data_wa_sub3)), 
            check_overlap = F, vjust = 0, hjust = 0.5, size = 3, color = "black") +
  theme_classic() + 
  theme(legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.position = "right",
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12, face = "bold")) +
  scale_color_manual(values=c("indianred","#E69F00","#56B4E9")) + 
  labs(x = labs[1],y = labs[2], title = NULL)

dev.off()

## heatmap
# top 20 most prevalent per host
top_cat = rownames(cat_AMR)[1:15]
top_dog = rownames(dog_AMR)[1:15]
top_hmp = rownames(hmp_AMR)[1:15]

top_union = union(top_cat,top_dog)
top_union = union(top_union,top_hmp)

# top 15 coefficients genes among significant maaslin results
maaslin_sig8_list = union(unique(uniref_catRef_dog_vs_cat$feature)[1:15], 
                          unique(uniref_hmpRef_dog_vs_human$feature)[1:15]) %>% 
  union(unique(uniref_hmpRef_dog_vs_human$feature)[1:15])

# reformat the maaslin genes annotation
maaslin_sig8_list = sub("\\.\\.", ": ", maaslin_sig8_list)
maaslin_sig8_list = gsub("\\.", " ", maaslin_sig8_list)
maaslin_sig8_list = gsub("\\:.*", "", maaslin_sig8_list)

# grep lines containing masslin_sig8_list
matching_rows = sapply(rownames(AMR_genes), function(row_name) {
  any(sapply(maaslin_sig8_list, function(pattern) {
    grepl(pattern, row_name)
  }))
})

most_prev = AMR_genes[matching_rows, ]

min = min(most_prev[most_prev>0])/10 # set 0 to be 1/10 of the minimum relative abundance = 1.79949e-09
most_prev[most_prev == 0] = min
most_prev = log10(most_prev)

most_prev = most_prev[,colnames(most_prev) %in% rownames(metadata)] # samples matched with metadata
rownames(most_prev) = gsub("UniRef90_", "", rownames(most_prev))
colnames(most_prev) = NULL # matrix

metadata_hm = metadata %>% select(c("Study ID" = study_readable, "Host Species" = species)) # metadata

# color setting
col = list(`Host Species`= c("cat" = "indianred","dog" = "#E69F00","human" = "#56B4E9"),
           `Study ID` = c('HMP1-II' = "#56B4E9", 'Cross-sectional Study2' = "#b7aaba", 'Cross-sectional Study3' = "#d44ab4", 
                          'Deusch et al. (2014)' = "#639CCF", 'Allaway et al.' = "#31378c", 'Deusch et al. (2015)' = "#6d81ae", 
                          'DietInt Study4' = "#798C7E", 'DietInt Study5' = "#b05c63", 'Cross-sectional Study1' = "#FFFFB3", 
                          'DietInt Study1' = "#d69d68", 'DietInt Study3' = "#7eda86", 'DietInt Study2' = "#fbc086",
                          'Young et al.' = "#74d0d1", 'Ateba et al.' = "#e0a8f5", 'Liu et al.' = "#d2bed7", 'Tanprasertsuk et al.' = "#c7a679",
                          'Ma et al.' = "#612141", 'Xu et al.' = "#ae4109", 'Alessandri et al.' = "#b435ee", 'Wang et al.' = "#e9fdc3",
                          'Maldonado-Contreras et al.' = "#c7f4fc",'Xu et al.' ="#ae4109",'Coelho et al.'="#4D643A",'Yarlagadda et al.'="#B5E6D4"))


# plot
library(ComplexHeatmap)
library(circlize)
cairo_pdf("CARDAMR_heatmap_top30_maaslin_byhost.pdf",width = 14,height = 7)
pheatmap(mat = as.matrix(most_prev),
         name = "Log10 Relative Abundance",
         show_colnames = FALSE,
         color = colorRamp2(c(-8.744851,-1.890752), c("white","#333366")),
         border_color = NA, 
         annotation_col = metadata_hm,
         annotation_colors = col,
         #cutree_cols = length(metadata_hm$`Host Species`),
         column_split = metadata_hm$`Host Species`,
         cluster_cols = TRUE,
         legend = T,
         annotation_legend = FALSE)
dev.off()

###bring in CARD annotations
card_annotations = as.data.frame(data.table::fread(file.path(input_dir,"aro_index.tsv"), header=T, check.names = F))

card_uniprot_align = as.data.frame(data.table::fread(file.path(input_dir,"output.tsv"), header=F, check.names = F))
names(card_uniprot_align) = c("CARD_ID","uniref_ID","query_seq_length", "qstart","qend","slen","sstart","send","alignment_length","bscore",
                              "eval","percent_identity","qcovhsp","scovhsp")
card_uniprot_align$uniref_ID = gsub("\\|.*","",card_uniprot_align$uniref_ID)
card_uniprot_align$query_cov = ((card_uniprot_align$qend - card_uniprot_align$qstart + 1) / card_uniprot_align$query_seq_length)*100
card_uniprot_align$sub_cov = ((card_uniprot_align$send - card_uniprot_align$sstart + 1) / card_uniprot_align$slen)*100
card_uniprot_align = card_uniprot_align[card_uniprot_align$percent_identity>=90 & card_uniprot_align$query_cov>=80 & card_uniprot_align$sub_cov>=80,]
card_uniprot_align = card_uniprot_align %>% distinct %>% group_by(CARD_ID) %>% top_n(1, percent_identity)
card_uniprot_align = card_uniprot_align %>% distinct %>% group_by(CARD_ID) %>% top_n(1, query_cov)
card_uniprot_align = card_uniprot_align %>% distinct %>% group_by(CARD_ID) %>% top_n(1, sub_cov)
card_uniprot_align$CARD_ID = str_split_fixed(card_uniprot_align$CARD_ID, '\\|', 4)[,3]

AMR_genes_card_annotation = AMR_genes
rownames(AMR_genes_card_annotation) = gsub("\\:.*","",rownames(AMR_genes_card_annotation))

AMR_genes_card_annotation_CARD_ID = AMR_genes_card_annotation
AMR_genes_card_annotation_CARD_ID$CARD_ID = card_uniprot_align$CARD_ID[match(rownames(AMR_genes_card_annotation), card_uniprot_align$uniref_ID)]
AMR_genes_card_annotation_CARD_ID$drug_class = card_annotations$`Drug Class`[match(AMR_genes_card_annotation_CARD_ID$CARD_ID, card_annotations$`ARO Accession`)]
AMR_genes_card_annotation_CARD_ID$resistance_mech = card_annotations$`Resistance Mechanism`[match(AMR_genes_card_annotation_CARD_ID$CARD_ID, card_annotations$`ARO Accession`)]
AMR_genes_card_annotation_CARD_ID$uniref90 = rownames(AMR_genes_card_annotation_CARD_ID)

AMR_genes_card_annotation_CARD_ID_drugStrat = AMR_genes_card_annotation_CARD_ID %>% separate_rows(drug_class,sep = ";")

abx_df = melt(AMR_genes_card_annotation_CARD_ID_drugStrat, id='drug_class')
abx_df$species = metadata$species[match(abx_df$variable,rownames(metadata))]
abx_df$value = as.numeric(abx_df$value)
abx_df = aggregate(value ~ drug_class + variable + species, data=abx_df, sum)

abx_df_log = abx_df
abx_df_log[abx_df_log==0] = 0.0000000001
abx_df_log$value = log10(abx_df_log$value)
abx_df_log$value = as.numeric(abx_df_log$value)

#"un-melt" abx_df and run maaslin and then annotated figure
abx_df_pivot = dcast(data = abx_df,formula = drug_class~variable,fun.aggregate = sum,value.var = "value")
rownames(abx_df_pivot) = abx_df_pivot$drug_class
abx_df_pivot$drug_class = NULL

# abx --> drug classes
AMR_genes_card_annotation_CARD_ID_drugStrat$ABX_Class = AMR_genes_card_annotation_CARD_ID_drugStrat$drug_class
AMR_genes_card_annotation_CARD_ID_drugStrat_Action = AMR_genes_card_annotation_CARD_ID_drugStrat %>% 
  mutate(ABX_Class = case_when(drug_class == "tetracycline antibiotic" ~ "Tetracyclines",
                                 drug_class == "aminocoumarin antibiotic"  ~ "Aminocoumarins",
                                 drug_class == "aminoglycoside antibiotic" ~ "Aminoglycosides",
                                 drug_class == "carbapenem" ~ "B-lactams",
                                 drug_class == "cephalosporin" ~ "B-lactams",
                                 drug_class == "cephamycin" ~ "B-lactams",
                                 drug_class == "diaminopyrimidine antibiotic" ~ "Diaminopyrimidines",
                                 drug_class == "disinfecting agents and antiseptics" ~ "Disinfecting agents and antiseptics",
                                 drug_class == "fluoroquinolone antibiotic" ~ "Quinolones",
                                 drug_class == "fusidane antibiotic" ~ "Fusidanes",
                                 drug_class == "glycopeptide antibiotic" ~ "Glycopeptides",
                                 drug_class == "glycylcycline" ~ "Tetracyclines",
                                 drug_class == "lincosamide antibiotic" ~ "Lincosamides",
                                 drug_class == "macrolide antibiotic" ~ "Macrolides",
                                 drug_class == "monobactam" ~ "B-lactams",
                                 drug_class == "mupirocin-like antibiotic" ~ "Pleuromutilins",
                                 drug_class == "nitrofuran antibiotic" ~ "Nitrofurans",
                                 drug_class == "nitroimidazole antibiotic" ~ "Nitroimidazoles",
                                 drug_class == "nucleoside antibiotic" ~ "Nucleosides",
                                 drug_class == "oxazolidinone antibiotic" ~ "Oxazolidinones",
                                 drug_class == "penam" ~ "B-lactams",
                                 drug_class == "penem" ~ "B-lactams",
                                 drug_class == "peptide antibiotic" ~ "Peptide antibiotics",
                                 drug_class == "phenicol antibiotic" ~ "Phenicols",
                                 drug_class == "phosphonic acid antibiotic" ~ "Phosphonic acid antibiotics",
                                 drug_class == "pleuromutilin antibiotic" ~ "Pleuromutilins",
                                 drug_class == "rifamycin antibiotic" ~ "Ansamycins",
                                 drug_class == "streptogramin A antibiotic" ~ "Streptogramins",
                                 drug_class == "streptogramin antibiotic" ~ "Streptogramins",
                                 drug_class == "streptogramin B antibiotic" ~ "Streptogramins",
                                 drug_class == "sulfonamide antibiotic" ~ "Sulfonamides",
                                 drug_class == "sulfone antibiotic" ~ "Sulfones"))

CARDABX_df = melt(AMR_genes_card_annotation_CARD_ID_drugStrat_Action, id='ABX_Class')
CARDABX_df$species = metadata$species[match(CARDABX_df$variable,rownames(metadata))]
CARDABX_df$value = as.numeric(CARDABX_df$value)
CARDABX_df = aggregate(value ~ ABX_Class + variable + species, data=CARDABX_df, sum)
write.csv(CARDABX_df,file.path(output_dir,"ABX_general_classes.csv"), row.names = FALSE)

CARDABX_df_log = CARDABX_df
CARDABX_df_log[CARDABX_df_log==0] = min(CARDABX_df_log$value[CARDABX_df_log$value>0])/10
CARDABX_df_log$value = log10(CARDABX_df_log$value)
CARDABX_df_log$value = as.numeric(CARDABX_df_log$value)

top_abx = c('Tetracyclines','Lincosamides','Macrolides','B-lactams','Aminoglycosides','Streptogramins','Quinolones')

####################################################################
### #Figure 6B boxplots and supplemental figure@@@ full boxplots ###
####################################################################
proportion_zeros = CARDABX_df %>%
  dplyr::group_by(ABX_Class,species) %>% 
  dplyr::summarise(prop_0 = sum(value==0)/dplyr::n())

proportion_zeroes = ggplot(proportion_zeros, aes(x=reorder(factor(ABX_Class),-prop_0), y=prop_0, fill = factor(species)))+
  geom_bar(position="dodge", stat="identity",width = 0.60)+
  theme( legend.position = "none" ) +
  coord_flip() +
  scale_fill_manual(values=c(cat="indianred",dog="#E69F00",human="#56B4E9")) +
  ylab("Proportion of zeroes") +
  xlab("") +
  theme_light() +
  theme(legend.position="none") +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(limits = c(0,1),breaks = c(0.0,0.5,1.0))

cairo_pdf(file.path(output_dir,"ABX_proportion_zeros.pdf"),width = 2,height = 20)
proportion_zeroes
dev.off()

CARDABX_df_log_new = CARDABX_df
CARDABX_df_log_new = CARDABX_df_log_new[CARDABX_df_log_new$value!=0,]

CARDABX_df_log_new$value = log10(CARDABX_df_log_new$value)
CARDABX_df_log_new$value = as.numeric(CARDABX_df_log_new$value)

levels = c("Tetracyclines", "Lincosamides", "Macrolides","B-lactams","Aminoglycosides",
           "Streptogramins","Phenicols","Diaminopyrimidines","Nucleosides","Peptide antibiotics",
           "Ansamycins","Quinolones","Disinfecting agents and antiseptics","Pleuromutilins",
           "Aminocoumarins","Nitroimidazoles","Phosphonic acid antibiotics","Glycopeptides",
           "Sulfones","Sulfonamides","Oxazolidinones","Nitrofurans","Fusidanes")
levels = rev(levels)

CARDABX_df_log_new$ABX_Class <- factor(CARDABX_df_log_new$ABX_Class, levels=levels)

ABX_revised = ggplot(CARDABX_df_log_new, aes(x=ABX_Class, y=value, color = factor(species)))+
  geom_boxplot(lwd=1) +
  theme( legend.position = "none" ) +
  coord_flip() +
  scale_color_manual(values=c(cat="indianred",dog="#E69F00",human="#56B4E9")) +
  ylab("log10(relative abundance)") +
  xlab("") +
  theme_light() +
  theme(legend.position="none") +
  theme(text = element_text(size = 20))

cairo_pdf(file.path(output_dir,"ABX_target_generalClasses_Nozeros.pdf"),width = 6,height = 10)
ABX_revised
dev.off()

cairo_pdf(file.path(output_dir,"ABX_revised_fig.pdf"),width =16,height = 18)
ggpubr::ggarrange(proportion_zeroes, ABX_revised + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank() ), nrow = 1,widths = c(0.75, 1))
dev.off()

#revised truncated
proportion_zeros_trun = proportion_zeros[proportion_zeros$ABX_Class %in% top_abx,]
CARDABX_df_log_new_trun = CARDABX_df_log_new[CARDABX_df_log_new$ABX_Class %in% top_abx,]

proportion_zeroes_trun = ggplot(proportion_zeros_trun, aes(x=reorder(factor(ABX_Class),-prop_0), y=prop_0, fill = factor(species)))+
  geom_bar(position="dodge", stat="identity",width = 0.60)+
  theme( legend.position = "none" ) +
  coord_flip() +
  scale_fill_manual(values=c(cat="indianred",dog="#E69F00",human="#56B4E9")) +
  ylab("Proportion of zeroes") +
  xlab("") +
  theme_light() +
  theme(legend.position="none") +
  theme(text = element_text(size = 20)) +
  scale_y_continuous(limits = c(0,1),breaks = c(0.0,0.5,1.0))

ABX_revised_trun = ggplot(CARDABX_df_log_new_trun, aes(x=ABX_Class, y=value, color = factor(species)))+
  geom_boxplot(lwd=1) +
  theme( legend.position = "none" ) +
  coord_flip() +
  scale_color_manual(values=c(cat="indianred",dog="#E69F00",human="#56B4E9")) +
  ylab("log10(relative abundance)") +
  xlab("") +
  theme_light() +
  theme(legend.position="none") +
  theme(text = element_text(size = 20))

cairo_pdf(file.path(output_dir,"ABX_revised_truncated_fig.pdf"),width =8,height = 8)
ggpubr::ggarrange(proportion_zeroes_trun, ABX_revised_trun + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank() ), nrow = 1,widths = c(0.75, 1))
dev.off()

CARDABX_df_pivot = dcast(data = CARDABX_df,formula = ABX_Class~variable,fun.aggregate = sum,value.var = "value")
rownames(CARDABX_df_pivot) = CARDABX_df_pivot$ABX_Class
CARDABX_df_pivot$ABX_Class = NULL

write.csv(CARDABX_df_pivot,file.path(output_dir,"S14_ARGs_ABXclasses_abundance_table.csv"))

#maaslin for abx targets
Maaslin2(CARDABX_df_pivot,metadata_for_amr_adonis,output = file.path(output_dir,"maaslin_pets_HostsStrat_hmpRef_ABXtargets_superclasses"),fixed_effects = c("species"),reference=c("species,human"), random_effects = c("study_readable"))
Maaslin2(CARDABX_df_pivot,metadata_for_amr_adonis,output = file.path(output_dir,"maaslin_pets_HostsStrat_dogRef_ABXtargets_superclasses"),fixed_effects = c("species"),reference=c("species,dog"), random_effects = c("study_readable"))
Maaslin2(CARDABX_df_pivot,metadata_for_amr_adonis,output = file.path(output_dir,"maaslin_pets_HostsStrat_catRef_ABXtargets_superclasses"),fixed_effects = c("species"),reference=c("species,cat"), random_effects = c("study_readable"))

#individual AMR gene boxplots
maaslin_df = as.data.frame(t(AMR_genes_filtered_maas))
maaslin_df$samples = rownames(maaslin_df)
maaslin_df = melt(maaslin_df)
maaslin_df$species = metadata$species[match(maaslin_df$sample,rownames(metadata))]

###### Individual gene distributions (Supplemental Figure 16) #########
Q8G3I2_df = maaslin_df[maaslin_df$variable == 'UniRef90_Q8G3I2: Isoleucine--tRNA ligase',]

cairo_pdf(file.path(output_dir,"UniRef90_Q8G3I2: Isoleucine--tRNA ligase"),width = 6,height = 6)
ggplot(Q8G3I2_df, mapping= aes(x = variable, y = value))+
  geom_boxplot(aes(color = species), outlier.shape = NA ) +
  geom_point(aes(color= species), alpha = 0.5, size=4,
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_color_manual(values=c(cat="indianred",dog="#E69F00",human="#56B4E9")) +
  xlab("") +
  ylab("Relative abundance") +
  ggtitle("UniRef90_Q8G3I2: Isoleucine--tRNA ligase") +
  theme_classic() +
  theme(legend.position="none") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
dev.off()

A0A086ZQ37_df = maaslin_df[maaslin_df$variable == 'UniRef90_A0A086ZQ37: DNA-directed RNA polymerase subunit beta',]

cairo_pdf(file.path(output_dir,"UniRef90_A0A086ZQ37: DNA-directed RNA polymerase subunit beta"),width = 6,height = 6)
ggplot(A0A086ZQ37_df, mapping= aes(x = variable, y = value))+
  geom_boxplot(aes(color = species), outlier.shape = NA ) +
  geom_point(aes(color= species), alpha = 0.5, size=4,
             position = position_jitterdodge(jitter.width = 0.1)) +
  scale_color_manual(values=c(cat="indianred",dog="#E69F00",human="#56B4E9")) +
  xlab("") +
  ylab("Relative abundance") +
  ggtitle("UniRef90_A0A086ZQ37: DNA-directed RNA polymerase subunit beta") +
  theme_classic() +
  theme(legend.position="none") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
dev.off()