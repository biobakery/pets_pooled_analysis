#This function that takes in pairwise anpan results from gene model and filters for significant genes with effect sizes > 1
#funtion also finds enriched genes for cats, dogs, and humans. Note: in the anpan pairwise tests, host is the outcome. 
#in the dog vs. cats pairwise test, dogs are given the outcome value 1, and cats are given the outcome value 0. In the
#humans vs. dogs pairwise test, dogs = 1, humans = 0. In the humans vs. cats pairwise test, cats = 1, humans = 0.
#species-specific total gene values (used to normalize by every 1000 genes) were derived from the species pangenomes

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

process_data <- function(file_paths, total_genes) {
  dfs <- lapply(file_paths, function(path) {
    data <- as.data.frame(data.table::fread(path, header = TRUE, check.names = FALSE))
    return(data)
  })
  
  combined_data <- bind_rows(dfs[[1]], dfs[[2]], dfs[[3]])
  combined_data$fdr_adjusted <- p.adjust(combined_data$p.value, method = "BH")
  
  cat_dog_data <- combined_data[1:dim(dfs[[1]])[1], 1:7]
  human_dog_data <- combined_data[(dim(dfs[[1]])[1] + 1):(dim(dfs[[1]])[1] + dim(dfs[[2]])[1]), 1:7]
  human_cat_data <- combined_data[(dim(dfs[[1]])[1] + dim(dfs[[2]])[1] + 1):(dim(dfs[[1]])[1] + dim(dfs[[2]])[1] + dim(dfs[[3]])[1]), 1:7]
  
  cat_dog_data <- cat_dog_data[cat_dog_data$fdr_adjusted < 0.10, ]
  cat_dog_data <- cat_dog_data[abs(cat_dog_data$estimate) >= 1, ]
  
  human_dog_data <- human_dog_data[human_dog_data$fdr_adjusted < 0.10, ]
  human_dog_data <- human_dog_data[abs(human_dog_data$estimate) >= 1, ]
  
  human_cat_data <- human_cat_data[human_cat_data$fdr_adjusted < 0.10, ]
  human_cat_data <- human_cat_data[abs(human_cat_data$estimate) >= 1, ]
  
  cat_dog_Dog_enriched <- cat_dog_data[cat_dog_data$estimate > 0, ]
  human_dog_Dog_enriched <- human_dog_data[human_dog_data$estimate > 0, ]
  gnavus_dog_enriched <- union(cat_dog_Dog_enriched$gene, human_dog_Dog_enriched$gene)
  
  cat_dog_Cat_enriched <- cat_dog_data[cat_dog_data$estimate < 0, ]
  human_cat_Cat_enriched <- human_cat_data[human_cat_data$estimate > 0, ]
  gnavus_cat_enriched <- union(cat_dog_Cat_enriched$gene, human_cat_Cat_enriched$gene)
  
  human_cat_Human_enriched <- human_cat_data[human_cat_data$estimate < 0, ]
  human_dog_Human_enriched <- human_dog_data[human_dog_data$estimate < 0, ]
  gnavus_human_enriched <- union(human_cat_Human_enriched$gene, human_dog_Human_enriched$gene)
  
  result <- c(
    dim(cat_dog_data)[1], dim(cat_dog_Dog_enriched)[1], dim(cat_dog_Cat_enriched)[1],
    dim(human_dog_data)[1], dim(human_dog_Dog_enriched)[1], dim(human_dog_Human_enriched)[1],
    dim(human_cat_data)[1], dim(human_cat_Cat_enriched)[1], dim(human_cat_Human_enriched)[1],
    length(gnavus_dog_enriched), length(gnavus_cat_enriched), length(gnavus_human_enriched),
    total_genes
  )
  
  return(result)
}

#### R. gnavus ####
# File paths for R.gnavus
file_paths_rgnavus <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Ruminococcus_gnavus_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Ruminococcus_gnavus_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Ruminococcus_gnavus_gene_terms.tsv"
)
# Total genes for R.gnavus
total_genes_rgnavus <- 9844
# Call the function for R.gnavus
R.gnavus <- process_data(file_paths_rgnavus, total_genes_rgnavus)
print(R.gnavus)

#### Blautia wexlerae ####
file_paths_bwexlerae <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Blautia_wexlerae_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Blautia_wexlerae_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Blautia_wexlerae_gene_terms.tsv"
)
# Total genes for Blautia wexlerae
total_genes_bwexlerae <- 5413
# Call the function for Blautia wexlerae
B.wexlerae <- process_data(file_paths_bwexlerae, total_genes_bwexlerae)
print(B.wexlerae)

#### P. copri ####
file_paths_pcopri <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Prevotella_copri_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Prevotella_copri_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Prevotella_copri_gene_terms.tsv"
)
# Total genes for pcopri
total_genes_pcopri <- 11818
# Call the function for pcopri
P.copri <- process_data(file_paths_pcopri, total_genes_pcopri)
print(P.copri)

#### A. hadrus ####
file_paths_ahadrus <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Anaerostipes_hadrus_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Anaerostipes_hadrus_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Anaerostipes_hadrus_gene_terms.tsv"
)
# Total genes for ahadrus
total_genes_ahadrus <- 7437
# Call the function for ahadrus
A.hadrus <- process_data(file_paths_ahadrus, total_genes_ahadrus)
print(A.hadrus)

#### B. stercoris ####
file_paths_bstercoris <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Bacteroides_stercoris_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Bacteroides_stercoris_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Bacteroides_stercoris_gene_terms.tsv"
)
# Total genes for bstercoris
total_genes_bstercoris <- 7149
# Call the function for bstercoris
B.stercoris <- process_data(file_paths_bstercoris, total_genes_bstercoris)
print(B.stercoris)

#### B. uniformis ####
file_paths_buniformis <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Bacteroides_uniformis_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Bacteroides_uniformis_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Bacteroides_uniformis_gene_terms.tsv"
)
# Total genes for buniformis
total_genes_buniformis <- 16710
# Call the function for buniformis
B.uniformis <- process_data(file_paths_buniformis, total_genes_buniformis)
print(B.uniformis)

#### B.vulgatus ####
file_paths_bvulgatus <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Bacteroides_vulgatus_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Bacteroides_vulgatus_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Bacteroides_vulgatus_gene_terms.tsv"
)
# Total genes for bvulgatus
total_genes_bvulgatus <- 17726
# Call the function for bvulgatus
B.vulgatus <- process_data(file_paths_bvulgatus, total_genes_bvulgatus)
print(B.vulgatus)

#### B. longum ####
file_paths_blongum <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Bifidobacterium_longum_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Bifidobacterium_longum_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Bifidobacterium_longum_gene_terms.tsv"
)
# Total genes for blongum
total_genes_blongum <- 12638
# Call the function for blongum
B.longum <- process_data(file_paths_blongum, total_genes_blongum)
print(B.longum)

#### Collinsella aerofaciens ####
file_paths_caerofaciens <- c(
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Collinsella_aerofaciens_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Collinsella_aerofaciens_gene_terms.tsv",
  "/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Collinsella_aerofaciens_gene_terms.tsv"
)
# Total genes for caerofaciens
total_genes_caerofaciens <- 3974
# Call the function for caerofaciens
C.aerofaciens <- process_data(file_paths_caerofaciens, total_genes_caerofaciens)
print(C.aerofaciens)

###The following microbes were not present in all hosts and make applying the above function a little tricky, 
###so each has its own code for finding the enriched genes in the respective hosts they are found in
### Megasphaera elsdensii
cat_dog_Mels = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Megasphaera_elsdenii_gene_terms.tsv", header=T, check.names = F))

cat_dog_Mels = cat_dog_Mels[cat_dog_Mels$q_bug_wise<0.10,]
cat_dog_Mels = cat_dog_Mels[abs(cat_dog_Mels$estimate)>=1,]

Mels_cat_enriched = cat_dog_Mels[cat_dog_Mels$estimate<0,]
Mels_dog_enriched = cat_dog_Mels[cat_dog_Mels$estimate>0,]

M.elsdensii_total_genes = 3663
M.elsdensii = c(dim(cat_dog_Mels)[1],"","","","","","","","",dim(Mels_dog_enriched)[1],dim(Mels_cat_enriched)[1],"",M.elsdensii_total_genes)

###Bifido pullorum
cat_dog_Bpullorum = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Bifidobacterium_pullorum_gene_terms.tsv", header=T, check.names = F))

cat_dog_Bpullorum = cat_dog_Bpullorum[cat_dog_Bpullorum$q_bug_wise<0.10,]
cat_dog_Bpullorum = cat_dog_Bpullorum[abs(cat_dog_Bpullorum$estimate)>=1,]

Bpullorum_cat_enriched = cat_dog_Bpullorum[cat_dog_Bpullorum$estimate<0,]
Bpullorum_dog_enriched = cat_dog_Bpullorum[cat_dog_Bpullorum$estimate>0,]

B.pullorum_total_genes = 1675
B.pullorum = c(dim(cat_dog_Bpullorum)[1],"","","","","","","","",dim(Bpullorum_dog_enriched)[1],dim(Bpullorum_cat_enriched)[1],"",B.pullorum_total_genes)

###Lactobacillus acidophilus
cat_dog_Lacido = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Lactobacillus_acidophilus_gene_terms.tsv", header=T, check.names = F))

cat_dog_Lacido = cat_dog_Lacido[cat_dog_Lacido$q_bug_wise<0.10,]
cat_dog_Lacido = cat_dog_Lacido[abs(cat_dog_Lacido$estimate)>=1,]

Lacido_cat_enriched = cat_dog_Lacido[cat_dog_Lacido$estimate<0,]
Lacido_dog_enriched = cat_dog_Lacido[cat_dog_Lacido$estimate>0,]

L.acido_total_genes = 4007
L.acido = c(dim(cat_dog_Lacido)[1],"","","","","","","","",dim(Lacido_dog_enriched)[1],dim(Lacido_cat_enriched)[1],"",L.acido_total_genes)

#Acidaminococcus_intestini_gene_terms.tsv.gz
human_cat_Aintestini = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Acidaminococcus_intestini_gene_terms.tsv", header=T, check.names = F))

human_cat_Aintestini = human_cat_Aintestini[human_cat_Aintestini$q_bug_wise<0.10,]
human_cat_Aintestini = human_cat_Aintestini[abs(human_cat_Aintestini$estimate)>=1,]

Aintestini_human_enriched = human_cat_Aintestini[human_cat_Aintestini$estimate<0,]
Aintestini_cat_enriched = human_cat_Aintestini[human_cat_Aintestini$estimate>0,]

A.intestini_total_genes = 2459
A.intestini = c("","","","","","",dim(human_cat_Aintestini)[1],"","","",dim(Aintestini_cat_enriched)[1],dim(Aintestini_human_enriched)[1],A.intestini_total_genes)

#Gemmiger_formicilis_gene_terms.tsv.gz
human_cat_Gform = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Gemmiger_formicilis_gene_terms.tsv", header=T, check.names = F))

human_cat_Gform = human_cat_Gform[human_cat_Gform$q_bug_wise<0.10,]
human_cat_Gform = human_cat_Gform[abs(human_cat_Gform$estimate)>=1,]

Gform_human_enriched = human_cat_Gform[human_cat_Gform$estimate<0,]
Gform_cat_enriched = human_cat_Gform[human_cat_Gform$estimate>0,]

G.form_total_genes = 2745
G.form = c("","","","","","",dim(human_cat_Gform)[1],"","","",dim(Gform_cat_enriched)[1],dim(Gform_human_enriched)[1],G.form_total_genes)

#Blautia obeum
human_cat_Bobeum = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Blautia_obeum_gene_terms.tsv", header=T, check.names = F))

human_cat_Bobeum = human_cat_Bobeum[human_cat_Bobeum$q_bug_wise<0.10,]
human_cat_Bobeum = human_cat_Bobeum[abs(human_cat_Bobeum$estimate)>=1,]

Bobeum_human_enriched = human_cat_Bobeum[human_cat_Bobeum$estimate<0,]
Bobeum_cat_enriched = human_cat_Bobeum[human_cat_Bobeum$estimate>0,]

B.obeum_total_genes = 17634
B.obeum = c("","","","","","",dim(human_cat_Bobeum)[1],"","","",dim(Bobeum_cat_enriched)[1],dim(Bobeum_human_enriched)[1],B.obeum_total_genes)

anpan_gene_diffs_df = data.frame(R.gnavus,P.copri,B.wexlerae,A.hadrus,B.longum,B.stercoris,B.uniformis,B.vulgatus,C.aerofaciens,M.elsdensii,B.pullorum,L.acido,A.intestini,G.form,B.obeum)
rownames(anpan_gene_diffs_df) = c("cat_dog","cat_dog_DOGenriched","cat_dog_CATenriched","human_dog","human_dog_DOGenriched","human_dog_HUMenriched","human_cat","human_cat_CATenriched","human_cat_HUMenriched","dog_enriched","cat_enriched","human_enriched","number_genes")
anpan_gene_diffs_df = as.data.frame(t(anpan_gene_diffs_df))
anpan_bugs = rownames(anpan_gene_diffs_df)
anpan_gene_diffs_df = as.data.frame(sapply(anpan_gene_diffs_df,as.numeric))
anpan_gene_diffs_df$correction_value = (1000/anpan_gene_diffs_df$number_genes)
anpan_gene_diffs_df[,1:12] = anpan_gene_diffs_df[,1:12]*anpan_gene_diffs_df$correction_value

rownames(anpan_gene_diffs_df) = anpan_bugs

# this dataframe is the basis for the anpan summary figure
anpan_gene_diffs_df[,1:12] = round(anpan_gene_diffs_df[,1:12],digits = 0)

###################################
#### R. gnavus heatmap section ####
###################################

##Note the same code below was used to make the main figure anpan heatmap and the supplemental anpan heatmap
##the main figure was generated from the union of the top 20 differential genes for each host pair, and the supplemental figure
##used the union of the top 50 differential genes per host pair
cat_dog_Rgnavus_heat = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/filtered_Ruminococcus_gnavus_catdog.tsv", header=T, check.names = F))
cat_dog_Rgnavus_heat$species[cat_dog_Rgnavus_heat$species==TRUE] = 'dog'
cat_dog_Rgnavus_heat$species[cat_dog_Rgnavus_heat$species==FALSE] = 'cat'

human_dog_Rgnavus_heat = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/filtered_Ruminococcus_gnavus_humandog.tsv", header=T, check.names = F))
human_dog_Rgnavus_heat$species[human_dog_Rgnavus_heat$species==TRUE] = 'dog'
human_dog_Rgnavus_heat$species[human_dog_Rgnavus_heat$species==FALSE] = 'human'

human_cat_Rgnavus_heat = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/filtered_Ruminococcus_gnavus_humancat.tsv", header=T, check.names = F))
human_cat_Rgnavus_heat$species[human_cat_Rgnavus_heat$species==TRUE] = 'cat'
human_cat_Rgnavus_heat$species[human_cat_Rgnavus_heat$species==FALSE] = 'human'

### Bring in gene annotation file (I grepped for the UniRef90 annotations from an annotated gene family on the cluster, created for the samples in this study and annotated by Humann)
gene_annot = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/gene_annotation_pets.tsv", header=T, check.names = F))

#cat_dog,human_dog, and human_cat genes with largest effect size
cat_dog_rgnavus = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/cat_dog/Ruminococcus_gnavus_gene_terms.tsv", header=T, check.names = F))
human_dog_rgnavus = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_dog/Ruminococcus_gnavus_gene_terms.tsv", header=T, check.names = F))
human_cat_rgnavus = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/anpan_exploratory/human_cat/Ruminococcus_gnavus_gene_terms.tsv", header=T, check.names = F))

comb_rgnavus = bind_rows(cat_dog_rgnavus,human_dog_rgnavus,human_cat_rgnavus)
comb_rgnavus$fdr_adjusted = p.adjust(comb_rgnavus$p.value, method="BH")
cat_dog_rgnavus = comb_rgnavus[1:4639,1:7]
human_dog_rgnavus = comb_rgnavus[4640:6353,1:7]
human_cat_rgnavus = comb_rgnavus[6354:11401,1:7]

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
write.table(top_genes,file='/Users/tobynbranck/Documents/pets_meta/top50ea_differential_genes_anpan.tsv',row.names=FALSE,col.names = FALSE,sep="\t",quote = FALSE)

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
  dplyr::summarise(across(starts_with("Uniref"), ~sum(., na.rm = TRUE)))

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
gene_annot$annotation = paste0(gene_annot$gene,": ",gene_annot$annotation)
rownames(df_for_vis) <- gene_annot$annotation[match(rownames(df_for_vis),gene_annot$gene)]
rownames(df_for_vis) = gsub("UniRef90_","",rownames(df_for_vis))

library(ComplexHeatmap)
library(circlize) 
colnames(df_for_vis) = NULL
col=list(species=c("cat"="indianred","dog"="#E69F00","human"="#56B4E9"))
ha = HeatmapAnnotation(`species`=df_for_vis_meta$species,col=col)

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/Rgnavus_Anpan_heatmap_top20.pdf",width = 16,height = 12)
Heatmap(as.matrix(df_for_vis),name="df_for_vis_meta",top_annotation = ha,col=colorRamp2(c(0,1), c("white","#333366")),show_heatmap_legend = FALSE)
dev.off()

############################################################
#####glycosyl hydrolases - potential supplement heatmap#####
############################################################
###the glycosyl hydrolases supplemental figure in the manuscript was built by finding the top 500 differential genes per pair
glycosyl = df_for_vis %>% filter(grepl('glyc|Glyc',rownames(df_for_vis)))
glycosyl = glycosyl %>% filter(!grepl('D4L2K7|A0A2N5PSH9|A7B3Y2|A0A2N5PD97|A7B3K7',rownames(glycosyl)))

colnames(glycosyl) = NULL
col=list(species=c("cat"="indianred","dog"="#E69F00","human"="#56B4E9"))
ha = HeatmapAnnotation(`species`=df_for_vis_meta$species,col=col)

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/Rgnavus_Anpan_heatmap_Glycosyl.pdf",width = 18,height = 10)
Heatmap(as.matrix(glycosyl),name="df_for_vis_meta",top_annotation = ha,col=colorRamp2(c(0,1), c("white","#333366")),show_heatmap_legend = FALSE)
dev.off()
