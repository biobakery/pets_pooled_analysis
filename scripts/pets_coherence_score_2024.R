#This code is adapted by Tobyn Branck from Dr. Ali Rahnavard's coherent_score.R written for his work in "Epidemiological associations with genomic variation in SARS-CoV-2"
# calculates a coherence score for each SGB which measures the host-specificity / niche-specificity of SGB subclades
library("ggplot2")
library("grid")
library(tidyverse)
library(factoextra)
library(ade4)
library(FactoMineR)
library("plyr")
library("dbplyr")
library(ape)
library(usedist)
library(vegan)

#set the file path for input & output files
input_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs")
output_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs")
source("/home/tobynb/shared-folder/tbranck/keep_SGB.R")

cluster_metadata = "host"

#define cluster score function
suppervised_sil_score <- function(D , metadata, cluster_metadata, cluster_name1, cluster_name2){
  cluster_a <- colnames(D)[metadata[colnames(D),cluster_metadata] == cluster_name1]
  if (is.na(cluster_name2))
    cluster_b <- colnames(D)[metadata[colnames(D),cluster_metadata] != cluster_name1]
  else
    cluster_b <- colnames(D)[metadata[colnames(D),cluster_metadata] == cluster_name2]
  diag(D) <- NA
  sil_scores <- vector("numeric", length(cluster_a))
  for (i in 1:length(cluster_a)){
    if (length(cluster_a) <5 | length(cluster_b) <5)
      #a = 0.0 #this was the original way to skip clusters (hosts) with fewer than 5 strains (samples)
      sil_scores[i] = NaN
    else{
      a <- mean(as.numeric(D[cluster_a[i], cluster_a[-i]]), na.rm = TRUE)
      b <- mean(as.numeric(D[cluster_a[i], cluster_b]), na.rm = TRUE)
      s <- (b-a)/max(a,b)
      sil_scores[i] <- s
    }
  }
  #print (sil_scores)
  cluster_sil_score <- mean(sil_scores, na.rm = TRUE)
  return (cluster_sil_score)
}

metaphlan_table = as.data.frame(data.table::fread(file.path(output_dir,"pets_metaphlan_v4_formatted_SGBtree_Jun23tax.csv"), header=T, check.names = F))
metaphlan_table$V1 = gsub("SP_","t_",metaphlan_table$V1)
rownames(metaphlan_table) = metaphlan_table$V1
metaphlan_table$V1 = NULL
metaphlan_table = keep_SGB(metaphlan_table)
metaphlan_table$taxa = rownames(metaphlan_table)
taxa_df = as.data.frame(metaphlan_table$taxa)
taxa_df$SGB = gsub(".*SGB","",taxa_df$`metaphlan_table$taxa`)
taxa_df$SGB = paste0("SGB",taxa_df$SGB,sep="")

metadata = as.data.frame(data.table::fread(file.path(input_dir,"metadata.csv"), header=T, check.names = F))
madagascar_samp_mapping = as.data.frame(data.table::fread(file.path(input_dir,"madagascar_SraRunTable.txt"), header=T, check.names = F))
metadata$sample_id_metaphlan[metadata$study_readable=='Madagascar'] = madagascar_samp_mapping$Run[match(metadata$sample_id_metaphlan,madagascar_samp_mapping$`Sample Name`)][239:350]
metadata_full = metadata
metadata = metadata[,c('sample_id_metaphlan','species')]
names(metadata) = c('sample','host')
#metadata$sample = paste0(metadata$sample,'_bowtie2')

rownames(metadata) = metadata$sample
metadata$sample = NULL
#rownames(metadata) = gsub("_bowtie2","",rownames(metadata))

dir<-file.path(input_dir,"distmats")
#get a list of all files with distmat in the name in your directory
files<-list.files(path=dir, pattern='.distmat', full.names = TRUE)
clusters = c("human","cat","dog")

scores_dataframe = data.frame(taxa = c(""),
                              human_sil_score = c(""),
                              cat_sil_score = c(""),
                              dog_sil_score = c(""))

# loops through Kimura-2 distance matrices (1 per SGB) and calculates the coherence score (defined in the function above). The coherence score is calculated for each
# host species (if enough of a host's samples carry the organism - see Methods).
# All of the SGBs' host-specific coherence scores are stored in a data frame and then an averaged coherence score is calculated for each SGB = final coherence score.
for (file in files){
  name <- file
  distance = as.data.frame(data.table::fread(file, header=F, check.names = F))
  distance$V1 = NULL
  distance[,ncol(distance)] = gsub("_bowtie2","",distance[,ncol(distance)])
  distance[,ncol(distance)] = sub(" .*", "", distance[,ncol(distance)])
  rownames(distance) = distance[,ncol(distance)]
  colnames(distance) = distance[,ncol(distance)]
  distance <- distance[1:(length(distance)-2)]
  
  metadata_bug_specific = metadata[rownames(metadata) %in% rownames(distance),]
  num_dog_samps = length(metadata_bug_specific[metadata_bug_specific=="dog"])
  num_cat_samps = length(metadata_bug_specific[metadata_bug_specific=="cat"])
  num_hum_samps = length(metadata_bug_specific[metadata_bug_specific=="human"])
  host_numbers_list = c(num_dog_samps,num_cat_samps,num_hum_samps)
  if (sum(host_numbers_list >= 5) < 2) { #skip bugs that are specific to a host (only in one host) (skips if less than 2 hosts do not have at least 5 samples with this bug in it)
    next
  }
  print (name)
  bugname = gsub(file.path(input_dir,"distmats/t__"),"",name)
  bugname = gsub("_group","",bugname)
  bugname = gsub(".distmat","",bugname)
  taxa_name = taxa_df$`metaphlan_table$taxa`[taxa_df$SGB == bugname]
  print(taxa_name)
  
  individual_bugs_scores = list()
  for (cl in clusters){
    temp_sil <- suppervised_sil_score(distance, metadata,
                                      cluster_metadata=cluster_metadata,
                                      cluster_name1 = cl, cluster_name2 = NA)
    individual_bugs_scores = append(individual_bugs_scores,temp_sil)
  }
  individual_bugs_scores = unlist(individual_bugs_scores)
  indiv_bug_df = data.frame(taxa = taxa_name,
                            human_sil_score = individual_bugs_scores[1],
                            cat_sil_score = individual_bugs_scores[2],
                            dog_sil_score = individual_bugs_scores[3])
  scores_dataframe = rbind(scores_dataframe,indiv_bug_df)
  
  #generates pcoa of distances
  distance[lower.tri(distance)] <- t(distance)[lower.tri(distance)]
  distance = distance[, which(colMeans(!is.na(distance)) > 0.5)]
  distance = distance[which(rowMeans(!is.na(distance)) > 0.5), ]
  if (any(is.na(distance))){
    next
  }
  e.sir.pcoa <- cmdscale( distance, eig = T )
  variance <- head(eigenvals(e.sir.pcoa)/sum(eigenvals(e.sir.pcoa)))
  x_variance <- as.integer(variance[1]*100)
  y_variance <- as.integer(variance[2]*100)
  e.sir.scores <- as.data.frame( e.sir.pcoa$points )

  e.sir.scores.meta <- merge( e.sir.scores, metadata, by = 'row.names' )
  rownames( e.sir.scores.meta ) <- e.sir.scores.meta[,1]
  e.sir.scores.meta[,1] <- NULL

  colnames(e.sir.scores.meta) <- c( "PCo1", "PCo2", "host" )
  percent_var_explained = metagMisc::eig_perc(e.sir.pcoa$eig, positive = T, plot = T)

  #e.sir.scores.meta_test = e.sir.scores.meta
  e.sir.scores.meta$newnames = e.sir.scores.meta$host
  e.sir.scores.meta$newnames = metadata_full$study_readable[match(rownames(e.sir.scores.meta),metadata_full$sample_id_metaphlan)] 
  e.sir.scores.meta$host[e.sir.scores.meta$newnames=='HMP1-II'] = "HMP1-II"
  e.sir.scores.meta$host[e.sir.scores.meta$newnames=='Madagascar'] = "Madagascar"
  e.sir.scores.meta$newnames = NULL
  
  cairo_pdf(sprintf(file.path(output_dir,"distmat_PCoAs/%s_phyloDist_PCoA.pdf"),taxa_name))
  print(ggplot( e.sir.scores.meta, aes(PCo1, PCo2, color=host) ) +
    geom_point(size = 4, alpha = 0.6) + theme_classic() +
    theme(axis.line.x = element_line(colour = 'black', linewidth=0.75, linetype='solid'),
          axis.line.y = element_line(colour = 'black', linewidth=0.75, linetype='solid'),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    xlab(paste("PCo1 (",percent_var_explained[1],"% variance explained)")) + ylab(paste("PCo2 (",percent_var_explained[2],"% variance explained)")) +
    scale_color_manual(values=c(cat='indianred',dog='#E69F00',"HMP1-II" = "#56B4E9","Madagascar"="#0f3b50")) +
    ggtitle(sprintf("%s",taxa_name)))
  dev.off()
}

write.csv(scores_dataframe,file.path(output_dir,"all_silhouette_scores_per_bug.csv"))

# Reading in the coherence scores for the SGBs (if want to skip the above calculation)
scores_dataframe = as.data.frame(read.csv(file.path(output_dir,"all_silhouette_scores_per_bug.csv")))
scores_dataframe$X = NULL

scores_dataframe = scores_dataframe %>% mutate_all(~ifelse(is.nan(.), NA, .))
scores_dataframe[scores_dataframe<0] = 0
# Average across the host-specific coherence scores for each SGB to get final coherence scores (1 per SGB)
scores_dataframe$mean_score = apply(scores_dataframe[,2:4], 1, mean, na.rm=TRUE)

# Formatting data for visualization of the coherence score heatmap and barplot (Figure 4 [highest 25 and lowest 25 coherence scores], and Supplemental figures [extended versions of Figure 4])
scores_dataframe = scores_dataframe[order(scores_dataframe$mean_score, decreasing = TRUE), ]  
scores_dataframe = scores_dataframe[2:nrow(scores_dataframe)-1,]

scores_melted = reshape2::melt(scores_dataframe, id=c('taxa','mean_score'))
scores_melted$value = as.numeric(scores_melted$value)
scores_melted$value=round(scores_melted$value,2)

###getting dataframe of SGBs and their respective phyla####
metaphlan_table = as.data.frame(data.table::fread(file.path(output_dir,"pets_metaphlan_v4_formatted_SGBtree_Jun23tax.csv"), header=T, check.names = F))
metaphlan_table$V1 = gsub("SP_","t_",metaphlan_table$V1)
rownames(metaphlan_table) = metaphlan_table$V1
metaphlan_table$V1 = NULL
phyla = rownames(metaphlan_table)

metaphlan_table = keep_SGB(metaphlan_table)
metaphlan_table$taxa = rownames(metaphlan_table)
taxa_df = as.data.frame(metaphlan_table$taxa)
taxa_df$SGB = gsub(".*SGB","",taxa_df$`metaphlan_table$taxa`)
taxa_df$SGB = paste0("SGB",taxa_df$SGB,sep="")

taxa_df$phyla = phyla
taxa_df$phyla = gsub("|c__.*", "", taxa_df$phyla)
taxa_df$phyla = gsub("k__Bacteria\\|p__", "",taxa_df$phyla)
taxa_df$phyla = gsub("\\|", "",taxa_df$phyla)

taxa_df = taxa_df %>%
  mutate(phyla_color = case_when(
    phyla == "Actinobacteria" ~ "#FF69B4",
    phyla == "Deinococcus_Thermus" ~ "#800080",
    phyla == "Bacteroidota" ~ "#008080",
    phyla == "Bacteria_unclassified" ~ "#FEF3AB",
    phyla == "Candidatus_Melainabacteria" ~ "#CDAD00",
    phyla == "Candidatus_Saccharibacteria" ~ "#FFFF99",
    phyla == "Candidatus_Thermoplasmatota" ~ "#F5FBAF",
    phyla == "Elusimicrobia" ~ "#E55748",
    phyla == "Euryarchaeota" ~ "#009ACD",
    phyla == "Firmicutes" ~ "#9BCD9B",
    phyla == "Fusobacteria" ~ "#DCF199",
    phyla == "Lentisphaerae" ~ "#660066",
    phyla == "Proteobacteria" ~ "#CDC1C5",
    phyla == "Spirochaetes" ~ "#808080",
    phyla == "Synergistetes" ~ "#BEAED4",
    phyla == "Tenericutes" ~ "#4DA7B0",
    phyla == "Thaumarchaeota" ~ "#B0171F",
    phyla == "Verrucomicrobia" ~ "#FDB164"
  ))

scores_melted$phyla_color <- taxa_df$phyla_color[match(scores_melted$taxa, taxa_df$`metaphlan_table$taxa`)]
scores_melted$phyla <- taxa_df$phyla[match(scores_melted$taxa, taxa_df$`metaphlan_table$taxa`)]

library(ggtext)
scores_melted$taxa <- paste0("<span style=\"color: ", scores_melted$phyla_color, "\">", scores_melted$taxa, "</span>")

cairo_pdf(file.path(output_dir,"silhouette_score_heatmap_mean.pdf"),height=15, width = 7)
ggplot(scores_melted,aes(variable,reorder(taxa,mean_score),fill=value,label=value)) + 
  geom_tile(color="white") +
  #geom_text(label=scores_melted$value,size=2.5) +
  coord_equal(ratio=.4) +
  viridis::scale_fill_viridis(discrete=FALSE,direction=-1) +
  xlab("") + 
  ylab("") +
  theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  theme(axis.text.y = ggtext::element_markdown(size = 10)) +
  geom_text(aes(color = after_scale(prismatic::best_contrast(fill, c("white", "black")))),size = 2.5)
dev.off()

scores_dataframe$phyla_color <- taxa_df$phyla_color[match(scores_dataframe$taxa, taxa_df$`metaphlan_table$taxa`)]
scores_dataframe$phyla <- taxa_df$phyla[match(scores_dataframe$taxa, taxa_df$`metaphlan_table$taxa`)]
scores_dataframe1 = scores_dataframe
scores_dataframe$taxa <- paste0("<span style=\"color: ", scores_dataframe$phyla_color, "\">", scores_dataframe$taxa, "</span>")
write.csv(scores_dataframe,file.path(output_dir,"coherence_scores_df.csv"))
cairo_pdf(file.path(output_dir,"silhouette_score_barplot_mean.pdf"),height=15, width = 7)
ggplot(scores_dataframe,aes(x=reorder(taxa,mean_score),y=mean_score,fill=mean_score))+
  geom_bar(stat="identity") + 
  coord_flip() +
  xlab("") +
  ylab("Coherence score") +
  #scale_fill_distiller(palette = "Blues",direction = 1) +
  viridis::scale_fill_viridis(direction = -1) +
  theme_classic() +
  theme(axis.text.y = ggtext::element_markdown(size = 10))
dev.off()

# Truncating the extended version of the data to visualize a summary of coherence scores in main text (Figure 4)
scores_barplot_truncated = rbind(scores_dataframe1[1:25,],scores_dataframe1[100:nrow(scores_dataframe1),])
scores_barplot_truncated$category = NA
scores_barplot_truncated$category[1:25] = "Top 25"
scores_barplot_truncated$category[26:50] = "Bottom 25"
scores_barplot_truncated$category = factor(scores_barplot_truncated$category,levels=c("Top 25","Bottom 25"))

scores_barplot_truncated$phyla_color <- taxa_df$phyla_color[match(scores_barplot_truncated$taxa, taxa_df$`metaphlan_table$taxa`)]
scores_barplot_truncated$phyla <- taxa_df$phyla[match(scores_barplot_truncated$taxa, taxa_df$`metaphlan_table$taxa`)]
scores_barplot_truncated$taxa <- paste0("<span style=\"color: ", scores_barplot_truncated$phyla_color, "\">", scores_barplot_truncated$taxa, "</span>")

cairo_pdf(file.path(output_dir,"silhouette_score_barplot_truncated.pdf"),height=10, width = 7.5)
ggplot(scores_barplot_truncated,aes(x=reorder(taxa,mean_score),y=mean_score,fill=mean_score))+
  geom_bar(stat="identity") + 
  xlab("") +
  ylab("Coherence score") +
  viridis::scale_fill_viridis(direction = -1) +
  theme_classic() +
  facet_grid(category~.,scales="free") +
  coord_flip() +
  theme(axis.text.y = ggtext::element_markdown(size = 12))
dev.off()
