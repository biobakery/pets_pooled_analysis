#This function takes in pairwise anpan results from the gene model and filters for significant genes with effect sizes > 1
#function also finds enriched genes for cats, dogs, and humans. Note: in the anpan pairwise tests, host is the outcome. 
#in the dog vs. cats pairwise test, dogs are given the outcome value 1, and cats are given the outcome value 0. In the
#humans vs. dogs pairwise test, dogs = 1, humans = 0. In the humans vs. cats pairwise test, cats = 1, humans = 0.
#species-specific total gene values (used to normalize by every 1000 genes) were derived from the species' pangenomes

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

#set the file path for input & output files
input_dir_cat_dog = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/cat_dog_pfam")
input_dir_human_dog = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/human_dog_pfam")
input_dir_human_cat = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/human_cat_pfam")
output_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs")

###################################
#### R. gnavus heatmap section ####
###################################

##Note the same code below was used to make the main figure anpan heatmap and the supplemental anpan heatmap
##the main figure was generated from the union of the top 20 differential genes for each host pair, and the supplemental figure
##used the union of the top 50 differential genes per host pair
cat_dog_Rgnavus_heat = as.data.frame(data.table::fread(file.path(input_dir_cat_dog,"cat_dog_filtered_Ruminococcus_gnavus.tsv"), header=T, check.names = F))
cat_dog_Rgnavus_heat$species[cat_dog_Rgnavus_heat$species==TRUE] = 'dog'
cat_dog_Rgnavus_heat$species[cat_dog_Rgnavus_heat$species==FALSE] = 'cat'

human_dog_Rgnavus_heat = as.data.frame(data.table::fread(file.path(input_dir_human_dog,"human_dog_filtered_Ruminococcus_gnavus.tsv"), header=T, check.names = F))
human_dog_Rgnavus_heat$species[human_dog_Rgnavus_heat$species==TRUE] = 'dog'
human_dog_Rgnavus_heat$species[human_dog_Rgnavus_heat$species==FALSE] = 'human'

human_cat_Rgnavus_heat = as.data.frame(data.table::fread(file.path(input_dir_human_cat,"human_cat_filtered_Ruminococcus_gnavus.tsv"), header=T, check.names = F))
human_cat_Rgnavus_heat$species[human_cat_Rgnavus_heat$species==TRUE] = 'cat'
human_cat_Rgnavus_heat$species[human_cat_Rgnavus_heat$species==FALSE] = 'human'

### Bring in gene annotation file (I grepped for the UniRef90 annotations from an annotated gene family table on the cluster that was created for the samples in this study and annotated by Humann)
pfam_annot = as.data.frame(data.table::fread("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs/map_pfam_name.txt", header=T, check.names = F))

#cat_dog,human_dog, and human_cat genes with largest effect size
cat_dog_rgnavus = as.data.frame(data.table::fread(file.path(input_dir_cat_dog,"cat_dog_Ruminococcus_gnavus_gene_terms.tsv"), header=T, check.names = F))
human_dog_rgnavus = as.data.frame(data.table::fread(file.path(input_dir_human_dog,"human_dog_Ruminococcus_gnavus_gene_terms.tsv"), header=T, check.names = F))
human_cat_rgnavus = as.data.frame(data.table::fread(file.path(input_dir_human_cat,"human_cat_Ruminococcus_gnavus_gene_terms.tsv"), header=T, check.names = F))

comb_rgnavus = bind_rows(cat_dog_rgnavus,human_dog_rgnavus,human_cat_rgnavus)
comb_rgnavus$fdr_adjusted = p.adjust(comb_rgnavus$p.value, method="BH")
cat_dog_rgnavus = comb_rgnavus[1:252,1:7]
human_dog_rgnavus = comb_rgnavus[253:408,1:7]
human_cat_rgnavus = comb_rgnavus[409:905,1:7]

cat_dog_rgnavus = cat_dog_rgnavus[cat_dog_rgnavus$fdr_adjusted<0.10,]
cat_dog_rgnavus = cat_dog_rgnavus[abs(cat_dog_rgnavus$estimate)>=1,]

human_dog_rgnavus = human_dog_rgnavus[human_dog_rgnavus$fdr_adjusted<0.10,]
human_dog_rgnavus = human_dog_rgnavus[abs(human_dog_rgnavus$estimate)>=1,]

human_cat_rgnavus = human_cat_rgnavus[human_cat_rgnavus$fdr_adjusted<0.10,]
human_cat_rgnavus = human_cat_rgnavus[abs(human_cat_rgnavus$estimate)>=1,]


cat_dog_rgnavus_for_vis = cat_dog_rgnavus[order(abs(cat_dog_rgnavus$estimate), decreasing = TRUE),]
cat_dog_rgnavus_topgenes = cat_dog_rgnavus_for_vis$gene[1:20]

human_cat_rgnavus_for_vis = human_cat_rgnavus[order(abs(human_cat_rgnavus$estimate), decreasing = TRUE),]
human_cat_rgnavus_topgenes = human_cat_rgnavus_for_vis$gene[1:20]

human_dog_rgnavus_for_vis = human_dog_rgnavus[order(abs(human_dog_rgnavus$estimate), decreasing = TRUE),]
human_dog_rgnavus_topgenes = human_dog_rgnavus_for_vis$gene[1:20]

top_genes = union(cat_dog_rgnavus_topgenes,human_dog_rgnavus_topgenes)
top_genes = union(top_genes, human_cat_rgnavus_topgenes)
write.table(top_genes,file=file.path(output_dir,"top50ea_differential_PFAMS_anpan.tsv"),row.names=FALSE,col.names = FALSE,sep="\t",quote = FALSE)

heat_df_combined = bind_rows(cat_dog_Rgnavus_heat,human_dog_Rgnavus_heat,human_cat_Rgnavus_heat)
heat_df_combined$housing = NULL
#note: when the dataframes for each pairwise test were combined, samples were duplicated in the new df. I noticed that for some samples,
#the duplicates had an "NA" for a given gene in one of the duplicates, that means that that gene was not interrogated in anpan for one of the
#pairwise tests. how could a gene be "there" in a sample but considered not there in the same sample when being tested in a different pairwise 
# test...turns out anpan also filters genes (drops genes) if they are consistently there across samples in both conditions and therefore an NA was
#assigned. The below line combines the information for the sample duplicates from multiple pairwise tests (note a given sample, actually host, was 
#used in up to 2 pairwise tests). the data used for the heatmaps in anpan is gene profiles for samples in the form of TRUE/FALSE for present/absent 
#of genes that passed anpan's filtering criteria. the below line gives 1 if TRUE and 0 if FALSE or NA. after the summarise function, genes with >1
#are "present" and 0 = "absent". a value of two means TRUE from both pairwise tests, 1 means TRUE from one pairwise test and that gene was dropped in
#the other pairwise test because it was too consistent over all samples of both conditions, 0 means FALSE or NA in one or both pairwise tests or not present
heat_df_combined = heat_df_combined %>% 
  dplyr::group_by(sample_id, species) %>% 
  dplyr::summarise(across(starts_with("PF"), ~sum(., na.rm = TRUE)))

heat_df_combined[heat_df_combined==2] <- 1

#remove genes that are not in 'top_genes'
heat_df_combined_top_genes = heat_df_combined[,colnames(heat_df_combined) %in% top_genes]
df_for_vis = cbind(heat_df_combined[,c('sample_id','species')],heat_df_combined_top_genes)
df_for_vis_meta = df_for_vis[,c('sample_id','species')]

df_for_vis_meta = as.data.frame(df_for_vis_meta)
rownames(df_for_vis_meta) = df_for_vis_meta$sample_id
df_for_vis_meta$sample_id = NULL

df_for_vis = as.data.frame(df_for_vis)
rownames(df_for_vis) = df_for_vis$sample_id
df_for_vis$sample_id = NULL
df_for_vis$species = NULL
df_for_vis = as.data.frame(t(df_for_vis))

#annotate genes using gene annotations file
pfam_annot$annotation = paste0(pfam_annot$gene,": ",pfam_annot$annotation)
rownames(df_for_vis) <- pfam_annot$annotation[match(rownames(df_for_vis),pfam_annot$gene)]
#rownames(df_for_vis) = gsub("UniRef90_","",rownames(df_for_vis))

library(ComplexHeatmap)
library(circlize) 
colnames(df_for_vis) = NULL
col=list(species=c("cat"="indianred","dog"="#E69F00","human"="#56B4E9"))
ha = HeatmapAnnotation(`species`=df_for_vis_meta$species,col=col)

cairo_pdf(file.path(output_dir,"Rgnavus_Anpan_heatmap_top20_pfam.pdf"),width = 14,height = 18)
Heatmap(as.matrix(df_for_vis),name="df_for_vis_meta",top_annotation = ha,col=colorRamp2(c(0,1), c("white","#333366")),show_heatmap_legend = FALSE)
dev.off()