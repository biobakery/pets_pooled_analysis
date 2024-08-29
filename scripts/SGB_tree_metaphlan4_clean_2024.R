# This script takes in as input 
# - the formatted metaphlan 4 relative abundance tables from companion animals and humans (HMP1-II and Madagascar cohorts combined) - formatted using the Pets_pooled_analysis_clean_2024.R script
# - the formatted metadata file generated from metadata.R
# - 
# Author: Tobyn Branck

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

#set the file path for input & output files
input_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs")
output_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs")

###################
####DATA IMPORT####
###################
# Import pets data
pets_MAGs = as.data.frame(data.table::fread(file.path(output_dir,"pets_MAGs_formatted_for_SGBtree_Jun23tax.csv"), header = T, check.names = F))
pets_metaphlan = as.data.frame(data.table::fread(file.path(output_dir,"pets_metaphlan_v4_formatted_SGBtree_Jun23tax.csv"), header=T, check.names = F))

# Import human (HMP1-II and Madagascar) files
human_metaphlan = as.data.frame(data.table::fread(file.path(output_dir,"ALLhuman_metaphlan4_SGBs_formatted_SGBtree_Jun23.csv"), header=T, check.names = F))

# Import metadata
metadata = as.data.frame(data.table::fread(file.path(input_dir,"metadata.csv"), header=T, check.names = F))

# Import SGBs_info.tsv, the output from the assembly pipeline; extracting the novel SGBs
pets_novel_sgbs = as.data.frame(data.table::fread(file.path(input_dir,"SGB_info.tsv"), header=T, check.names = F))
novel_clusters_genome_reps = pets_novel_sgbs$genome
pets_novel_sgbs = pets_novel_sgbs %>% separate_rows(cluster_members,sep=",")
pets_novel_sgbs$cluster_members = gsub(".fa","",pets_novel_sgbs$cluster_members)

# Import file that maps the MAG to SGB in Oct22 database
mapping_df = as.data.frame(data.table::fread(file.path(input_dir,"sequences_oct22_curtis.tsv"), header=T, check.names = F))
mapping_df$genome_id = gsub("Curtis__","",mapping_df$genome_id)
mapping_df = mapping_df[mapping_df$genome_id %in% pets_novel_sgbs$cluster_members,]
# novel_SGBs is a list of SGBs that the genomes from our novel clusters were assigned to in the new Oct22 database
novel_SGBs = unique(mapping_df$sgb_id)

#Import soil sample files
soil_metaphlan = as.data.frame(data.table::fread(file.path(input_dir,"soil_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
soil_metaphlan = soil_metaphlan %>% filter(grepl("t__SGB", soil_metaphlan$`# taxonomy`))
rownames(soil_metaphlan) = soil_metaphlan$`# taxonomy`
soil_metaphlan$`# taxonomy` = NULL
rownames(soil_metaphlan) = gsub(".*t__SGB","",rownames(soil_metaphlan))
rownames(soil_metaphlan) = gsub("_group","",rownames(soil_metaphlan))
soil_metaphlan$counts = rowSums(soil_metaphlan>.00001)
soil_SGBs = rownames(soil_metaphlan)

#renaming headers and rownames
metadata$sample_id = NULL
names(metadata)[1] = "sample"
names(pets_metaphlan)[1] = "taxonomy"
names(human_metaphlan)[1] = "taxonomy"
rownames(pets_metaphlan) = pets_metaphlan$taxonomy
rownames(human_metaphlan) = human_metaphlan$taxonomy

###############################################################################
####Prepping taxonomic files for tree annotation (presence/absence of SGBs)####
###############################################################################
human_df = human_metaphlan
human_df$taxonomy = NULL
pets_df = pets_metaphlan
pets_df$taxonomy = NULL
#pets_df = pets_df[1:1483,]

combined_data_frame = merge(pets_df,human_df,by="row.names",all=TRUE)
combined_data_frame[is.na(combined_data_frame)] = 0
rownames(combined_data_frame) = combined_data_frame$Row.names
combined_data_frame$Row.names = NULL

human_filt = as.data.frame(t(human_df))
human_filt$host = "human"
human_filt$sample = rownames(human_filt)
human_filt = reshape2::melt(human_filt,id=c("sample","host"))
pet_filt = as.data.frame(t(pets_df))
pet_filt$host = metadata$species[match(rownames(pet_filt),metadata$sample)]
pet_filt$sample = rownames(pet_filt)
pet_filt = reshape2::melt(pet_filt,id=c("sample","host"))

# combine long format for all hosts: bug abundance by sample 
combined_test = rbind(pet_filt,human_filt)

#count how many samples (by host) an SGB has > 0.00001 relative abundance
combined_prev_abund_filtering = aggregate(value ~ host + variable, combined_test, function(x) sum(x > 0.00001, na.rm = TRUE))

names(combined_prev_abund_filtering)[3] = 'number_greater_0.00001'
num_dog = count(metadata$species=='dog')[2,2]
num_cat = count(metadata$species=='cat')[2,2]
num_human = count(metadata$species=='human')[2,2]
combined_prev_abund_filtering$prevalence = ifelse(combined_prev_abund_filtering$host=="dog",combined_prev_abund_filtering$number_greater_0.00001/num_dog , ifelse(combined_prev_abund_filtering$host=="cat",combined_prev_abund_filtering$number_greater_0.00001/num_cat, 
                                                                                                                                                               combined_prev_abund_filtering$number_greater_0.00001/num_human))
combined_prev_abund_filtering = combined_prev_abund_filtering %>% dplyr::group_by(variable) %>% 
  dplyr::mutate(prev_dog = ifelse(host=='dog',prevalence,NA)) %>% 
  tidyr::fill(prev_dog,.direction="updown") %>% 
  ungroup()

combined_prev_abund_filtering = combined_prev_abund_filtering %>% dplyr::group_by(variable) %>% 
  dplyr::mutate(prev_cat = ifelse((host=='cat'),prevalence,NA)) %>% 
  tidyr::fill(prev_cat, .direction="updown") %>% 
  ungroup()

combined_prev_abund_filtering = combined_prev_abund_filtering %>% dplyr::group_by(variable) %>% 
  dplyr::mutate(prev_human = ifelse((host=='human'),prevalence,NA)) %>% 
  tidyr::fill(prev_human, .direction="updown") %>% 
  ungroup()

combined_prev_abund_filtering[is.na(combined_prev_abund_filtering)] = 0

combined_prev_abund_filtering_filtered = combined_prev_abund_filtering[combined_prev_abund_filtering$number_greater_0.00001 >=3, ] 
length(unique(combined_prev_abund_filtering_filtered$variable))

combined_SGBs = unique(combined_prev_abund_filtering_filtered$variable)
length(grep("SGB",combined_SGBs))
combined_SGBs = str_split_fixed(combined_SGBs,"\\|",8)[,8]
combined_SGBs = gsub("SP__","SGB_",combined_SGBs)

#### create a dataframe that will have tree annotation information ####
tree_metadata = combined_prev_abund_filtering_filtered
tree_metadata = tree_metadata %>% dplyr::group_by(variable) %>% dplyr::mutate(dog = any(host=='dog'))
tree_metadata = tree_metadata %>% dplyr::group_by(variable) %>% dplyr::mutate(cat = any(host=='cat'))
tree_metadata = tree_metadata %>% dplyr::group_by(variable) %>% dplyr::mutate(human = any(host=='human'))
tree_metadata = tree_metadata[!duplicated(tree_metadata$variable),] # when SGBs are found in more than 1 host, the current df has a row for each host it appears in. Since the dog, cat, human cols are already set up in this df, we can remove duplicate occurrences of SGBs (keep distinct rows)
tree_metadata = tree_metadata %>% dplyr::group_by(variable) %>% dplyr::mutate(Shared_CA = if(dog==TRUE & cat==TRUE & human==FALSE) "TRUE" else "FALSE")
tree_metadata = tree_metadata %>% dplyr::group_by(variable) %>% dplyr::mutate(Shared_All_Hosts = if(dog==TRUE & cat==TRUE & human==TRUE) "TRUE" else "FALSE")
tree_metadata$knownness = as.numeric(grepl('_SGB', tree_metadata$variable, ignore.case=T))
tree_metadata$knownness[tree_metadata$knownness==1] = "uSGB"
tree_metadata$knownness[tree_metadata$knownness==0] = "kSGB"

tree_metadata$SGB = gsub(".*SP__","", tree_metadata$variable)
tree_metadata$SGB = gsub("_group","",tree_metadata$SGB)
novel_SGBs = gsub("SGB","",novel_SGBs)
tree_metadata$knownness[tree_metadata$SGB %in% novel_SGBs] = "novel"

tree_metadata$phylum = gsub(".*p__","",tree_metadata$variable)
tree_metadata$phylum = gsub("\\|.*","",tree_metadata$phylum)

tree_metadata = as.data.frame(tree_metadata %>% distinct(variable, .keep_all=TRUE))
rownames(tree_metadata) = gsub(".*\\|", "", tree_metadata$variable)
rownames(tree_metadata) = gsub("SP__","",rownames(tree_metadata))
rownames(tree_metadata) = gsub("_group","",rownames(tree_metadata))

#adding in soil information -- note: 10 SGBs are present in soil if we don't filter for 10^-5 in @ least 3 samples (only 2 remain if we do)
tree_metadata$soil = NA
tree_metadata$soil[rownames(tree_metadata) %in% soil_SGBs] = "TRUE"
tree_metadata$soil[is.na(tree_metadata$soil)] = "FALSE"

library(RColorBrewer)
library(ggtree)
library(ggnewscale)
library(ggstar)
library(ggtreeExtra)
n <- 20
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors <- c("#9ACD32", "#EE6A50", "#87CEFA", "#FFC125", "#D15FEE", "#8DEEEE", "#800000",
            "#006400", "#800080", "#66C2A5", "#BEAED4", "#FFFF99","#808080", "#B0171F", "#191970", "#7B68EE",
            "#00CD00", "Black")

###combined_SGBs = what should be included in the tree
combined_SGBs_formatted = gsub("SGB_","",combined_SGBs)
combined_SGBs_formatted = gsub("_group","",combined_SGBs_formatted)

tree_metadata$tips = rownames(tree_metadata)
tree_metadata = tree_metadata[,c("tips","host","variable","number_greater_0.00001","prevalence","prev_dog","prev_cat","prev_human","dog","cat","human","Shared_CA","Shared_All_Hosts","knownness","phylum","soil")]
tree_metadata$novel[tree_metadata$knownness=="novel"] = "novel"
tree_metadata$novel[is.na(tree_metadata$novel)] = "non-novel"

#######################################
##### graphlan annotation section #####
#######################################
tree_metadata = tree_metadata %>%
  mutate(clade_marker_color = case_when(
    phylum == "Actinobacteria" ~ "#FF69B4",
    phylum == "Deinococcus_Thermus" ~ "#800080",
    phylum == "Bacteroidota" ~ "#008080",
    phylum == "Bacteria_unclassified" ~ "#FEF3AB",
    phylum == "Candidatus_Melainabacteria" ~ "#CDAD00",
    phylum == "Candidatus_Saccharibacteria" ~ "#FFFF99",
    phylum == "Candidatus_Thermoplasmatota" ~ "#F5FBAF",
    phylum == "Elusimicrobia" ~ "#E55748",
    phylum == "Euryarchaeota" ~ "#009ACD",
    phylum == "Firmicutes" ~ "#9BCD9B",
    phylum == "Fusobacteria" ~ "#DCF199",
    phylum == "Lentisphaerae" ~ "#660066",
    phylum == "Proteobacteria" ~ "#CDC1C5",
    phylum == "Spirochaetes" ~ "#808080",
    phylum == "Synergistetes" ~ "#BEAED4",
    phylum == "Tenericutes" ~ "#4DA7B0",
    phylum == "Thaumarchaeota" ~ "#B0171F",
    phylum == "Verrucomicrobia" ~ "#FDB164"
  ))

tree_metadata = tree_metadata %>%
  mutate(clade_marker_shape =  "o")

tree_metadata = tree_metadata %>%
  mutate(clade_marker_size = "25")

tree_metadata = tree_metadata %>%
  mutate(ring_color_1 = case_when(
    knownness == "kSGB" ~ "#FFFFFF",
    knownness == "uSGB" ~ "#6F6F6F",
    knownness == "novel" ~ "#6F6F6F"))

tree_metadata = tree_metadata %>%
  mutate(ring_color_2 = case_when(
    human == "FALSE" ~ "white",
    prev_human > 0 & prev_human <= 0.25 ~ "#229DE2",
    prev_human > 0.25 & prev_human <= 0.5 ~ "#177AB2",
    prev_human > 0.5 & prev_human <= 0.75 ~ "#10577E",
    prev_human > 0.75 & prev_human <= 1 ~ "#072739"))

tree_metadata = tree_metadata %>%
  mutate(ring_color_3 = case_when(
    dog == "FALSE" ~ "white",
    prev_dog > 0 & prev_dog <= 0.25 ~ "#FFB50E",
    prev_dog > 0.25 & prev_dog <= 0.5 ~ "#E69f00",
    prev_dog > 0.5 & prev_dog <= 0.75 ~ "#AB7600",
    prev_dog > 0.75 & prev_dog <= 1 ~ "#704E00"))

tree_metadata = tree_metadata %>%
  mutate(ring_color_4 = case_when(
    cat == "FALSE" ~ "white",
    prev_cat > 0 & prev_cat <= 0.25 ~ "#CD5C5C",
    prev_cat > 0.25 & prev_cat <= 0.5 ~ "#B63838",
    prev_cat > 0.5 & prev_cat <= 0.75 ~ "#892A2A",
    prev_cat > 0.75 & prev_cat <= 1 ~ "#5C1C1C"))

tree_metadata = tree_metadata %>%
  mutate(ring_color_5 = case_when(
    Shared_CA == "TRUE" ~ "#003333",
    Shared_CA == "FALSE" ~ "white"))

tree_metadata = tree_metadata %>%
  mutate(ring_color_6 = case_when(
    Shared_All_Hosts == "TRUE" ~ "#003333",
    Shared_All_Hosts == "FALSE" ~ "white"))

tree_metadata = tree_metadata %>%
  mutate(ring_shape_7 = case_when(
    novel == "novel" ~ "v",
    novel == "non-novel" ~ ""))

tree_metadata = tree_metadata %>%
  mutate(ring_color_7 = case_when(
    novel == "novel" ~ "black",
    novel == "non-novel" ~ "white"))

for(i in 1:6) {
  new <- rep(0.5, nrow(tree_metadata))
  tree_metadata[ , ncol(tree_metadata) + 1] <- new
  colnames(tree_metadata)[ncol(tree_metadata)] <- paste0("ring_width_", i)
}

annotation_file_combined = tree_metadata[,c("tips","clade_marker_color","clade_marker_shape","clade_marker_size","ring_color_1",
                                            "ring_color_2","ring_color_3","ring_color_4","ring_color_5","ring_color_6","ring_color_7",
                                            "ring_shape_7","ring_width_1","ring_width_2","ring_width_3","ring_width_4",
                                            "ring_width_5","ring_width_6")]
annotation_file_combined = reshape2::melt(annotation_file_combined,id=c("tips"))
global_combined <- t(
  data.frame(
    start_rotation = c(270, "", ""),
    total_plotted_degrees = c(330, "", ""),
    branch_bracket_depth =	c(0, "", ""),
    clade_marker_size	= c(0, "", ""),
    clade_marker_edge_color = c("grey", "", ""),
    ring_label = c(1, "kSGB/uSGB", ""),
    ring_label = c(2, "Human", ""),
    ring_label = c(3, "Dog", ""),
    ring_label = c(4, "Cat", ""),
    ring_label = c(5, "Companion Animals", ""),
    ring_label = c(6, "All hosts", ""),
    ring_label = c(7, "Novelty", ""),
    ring_label_font_size = c(1, 10, ""),
    ring_label_font_size = c(2, 10, ""),
    ring_label_font_size = c(3, 10, ""),
    ring_label_font_size = c(4, 10, ""),
    ring_label_font_size = c(5, 10, ""),
    ring_label_font_size = c(6, 10, ""),
    ring_label_font_size = c(7, 10, ""),
    ring_internal_separator_thickness = c(1, .5, ""),
    ring_external_separator_thickness = c(1, .5, ""),
    ring_external_separator_thickness = c(2, .5, ""),
    ring_external_separator_thickness = c(3, .5, ""),
    ring_external_separator_thickness = c(4, .5, ""),
    ring_external_separator_thickness = c(5, .5, ""),
    ring_external_separator_thickness = c(6, .5, ""),
    ring_separator_color = c(1, "#626262", ""),
    ring_separator_color = c(2, "#626262", ""),
    ring_separator_color = c(3, "#626262", ""),
    ring_separator_color = c(4, "#626262", ""),
    ring_separator_color = c(5, "#626262", ""),
    ring_separator_color = c(6, "#626262", ""),
    Actinobacteria = c("clade_marker_color","#FF69B4",""),
    Deinococcus_Thermus = c("clade_marker_color","#800080",""),
    Bacteroidota = c("clade_marker_color","#008080",""),
    Bacteria_unclassified = c("clade_marker_color","#FEF3AB",""),
    Candidatus_Melainabacteria = c("clade_marker_color","#CDAD00",""),
    Candidatus_Saccharibacteria = c("clade_marker_color","#FFFF99",""),
    Candidatus_Thermoplasmatota = c("clade_marker_color","#F5FBAF",""),
    Elusimicrobia = c("clade_marker_color","#E55748",""),
    Euryarchaeota = c("clade_marker_color","#009ACD",""),
    Firmicutes = c("clade_marker_color","#9BCD9B",""),
    Fusobacteria = c("clade_marker_color","#DCF199",""),
    Lentisphaerae = c("clade_marker_color","#660066",""),
    Proteobacteria = c("clade_marker_color","#CDC1C5",""),
    Spirochaetes = c("clade_marker_color","#808080",""),
    Synergistetes = c("clade_marker_color","#BEAED4",""),
    Tenericutes = c("clade_marker_color","#4DA7B0",""),
    Thaumarchaeota = c("clade_marker_color","#B0171F",""),
    Verrucomicrobia = c("clade_marker_color","#FDB164","")
  )
) %>%
  as.data.frame() %>%
  rownames_to_column("tips") %>%
  setNames(c("tips", "option", "level", "value")) %>%
  mutate(tips = gsub("\\..*", "", tips))

global_combined$value = NULL
names(global_combined) = c("tips","variable","value")
annotation_file_combined = rbind(global_combined,annotation_file_combined)

sub_for_format = annotation_file_combined[7011:nrow(annotation_file_combined),] %>% 
  extract(variable, into = c("variable", "ring_number"), "(.*)_([^_]+)$")
sub_for_format = sub_for_format %>% setNames(c("tips","variable","value","value2"))
annotation_file_combined$value2 = NA
annotation_file_combined[7011:nrow(annotation_file_combined),] = sub_for_format


write.table(annotation_file_combined, file = file.path(output_dir,"annotation_file_metaphlanv4Oct22_7.11.txt"), sep = "\t",
            row.names = FALSE,col.names = FALSE,na = "",quote=FALSE)

tree = read.tree(file.path(input_dir,"Oct22.nwk"))
metaphlan_tree = castor::get_subtree_with_tips(tree,only_tips = tree_metadata$tips)$subtree
ape::write.tree(metaphlan_tree,file=file.path(output_dir,"subsetted_tree_Oct22.txt"))

#library(rbiom)
test = merge(x = pets_metaphlan, y = human_metaphlan, by = "taxonomy",
      all = TRUE)
rownames(test) = test$taxonomy
test$taxonomy = NULL
test[is.na(test)] <- 0
rownames(test) = gsub(".*SP__","",rownames(test))

library(phyloseq)
#OTU = otu_table(test, taxa_are_rows = TRUE)
#physq = phyloseq(OTU,metaphlan_tree)
#unweighted_unifrac = UniFrac(physq, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
#weighted_unifrac = UniFrac(physq, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
#saveRDS(unweighted_unifrac, file = file.path(output_dir,"unweighted_unifrac_distance.rds"))
#saveRDS(weighted_unifrac, file = file.path(output_dir,"weighted_unifrac_distance.rds"))

weights = c(case_when(metadata$species == "human" ~ 1/350,
                      metadata$species == "cat" ~ 1/367,
                      metadata$species == "dog" ~ 1/2273))

#reading in weighted and unweighted unifrac data (since it takes hours to produce, saved as .rds)
unweighted_unifrac = readRDS("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs/unweighted_unifrac_distance.rds")
weighted_unifrac = readRDS("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs/weighted_unifrac_distance.rds")

unweight_ordination = wcmdscale(unweighted_unifrac, k=2,eig=TRUE,w = weights)
percent_var_explained_uw = metagMisc::eig_perc(unweight_ordination$eig, positive = T, plot = F)
percent_var_explained_uw = format(round(percent_var_explained_uw[1:2],digits = 2),nsmall=1,trim=TRUE)
labs_uw = c(glue::glue("PCoA1 [{percent_var_explained_uw[1]}%]"),
         glue::glue("PCoA2 [{percent_var_explained_uw[2]}%]"))

unweight_ordination <- data.frame(unweight_ordination$points)

weighted_ordination = wcmdscale(weighted_unifrac, k=2,eig=TRUE,w = weights)
percent_var_explained_w = metagMisc::eig_perc(weighted_ordination$eig, positive = T, plot = F)
percent_var_explained_w = format(round(percent_var_explained_w[1:2],digits = 2),nsmall=1,trim=TRUE)
labs_w = c(glue::glue("PCoA1 [{percent_var_explained_w[1]}%]"),
            glue::glue("PCoA2 [{percent_var_explained_w[2]}%]"))

weighted_ordination <- data.frame(weighted_ordination$points)

write.csv(as.data.frame(unweight_ordination),file = file.path(output_dir,"unweighted_ordination_weighted.csv"))
write.csv(as.data.frame(weighted_ordination),file = file.path(output_dir,"weighted_ordination_weighted.csv"))

#unweight_ordination <- as.data.frame(unweight_ordination)
metadata_ord = metadata
metadata_ord$species[metadata_ord$study_readable=="Madagascar"] = "Madagascar"

colnames(unweight_ordination) <- c("PCo1", "PCo2")
unweight_ordination$Species <- factor(metadata_ord$species) #add group of interest, mine was Morphospecies in the data frame cal_fem_data2
colnames(weighted_ordination) <- c("PCo1", "PCo2")
weighted_ordination$Species <- factor(metadata_ord$species) #add group of interest, mine was Morphospecies in the data frame cal_fem_data2
#set the order of samples for plotting purposes
reference = c("dog","cat","human")
unweight_ordination = unweight_ordination[order(factor(unweight_ordination$Species, levels = reference)),]
weighted_ordination = weighted_ordination[order(factor(weighted_ordination$Species, levels = reference)),]

cairo_pdf(file.path(output_dir,"unweighted_ordination_freq_weighted.pdf"))
ggplot(unweight_ordination, aes(x = PCo1, y = PCo2, color = Species)) + 
  geom_point(size = 3,alpha=.4) +
  xlab("PCo1") +
  ylab("PCo2") + 
  ggtitle("unweighted unifrac") +
  scale_color_manual(values = c("indianred","#E69F00","#56B4E9","#0f3b50")) +
  theme_classic() +
  labs(x = labs_uw[1],y = labs_uw[2], title = NULL)
dev.off()

cairo_pdf(file.path(output_dir,"weighted_ordination_freq_weighted.pdf"))
ggplot(weighted_ordination, aes(x = PCo1, y = PCo2, color = Species)) + 
  geom_point(size = 3,alpha=.4) +
  xlab("PCo1") +
  ylab("PCo2") + 
  ggtitle("weighted unifrac") +
  scale_color_manual(values = c("indianred","#E69F00","#56B4E9","#0f3b50")) +
  theme_classic() +
  labs(x = labs_w[1],y = labs_w[2], title = NULL)
dev.off()

###################################################################################
######### Figure 2: novel SGBs tree with surrounding kSGBs + common gut bugs#######
###################################################################################

##### obtaining known neighbors for the 9 SGBs that survived metaphlan filtering and were therefore (the ones that didn't survive choco aren't placed in the tree)
##### included in the metaphlan databaasae
test = tree_metadata[tree_metadata$knownness=='novel',]

tree_metadata_kSGB_novel = tree_metadata[tree_metadata$knownness=='novel' | tree_metadata$knownness=='kSGB',]
our_tree = castor::get_subtree_with_tips(tree,only_tips = tree_metadata_kSGB_novel$tips)$subtree
labels = our_tree$tip.label

tips_for_subset = list()
tips_next = list()
tips_prev = list() # dump the prev and next SGBs in here to understand the tree
for (item in test$tips){ #for the nine novel SGBs that were present in metaphlan, get neighboring SGBs
  id = match(item,labels)
  id_prev = id-1
  id_next = id+1
  prev_item = labels[id_prev]
  next_item = labels[id_next]
  tips_prev[[length(tips_prev) +1]] = prev_item
  tips_next[[length(tips_next) +1]] = next_item
  tips_for_subset[[length(tips_for_subset) +1]] = item
  tips_for_subset[[length(tips_for_subset) +1]] = prev_item
  tips_for_subset[[length(tips_for_subset) +1]] = next_item
}

closest_knowns = append(unique(unlist(tips_prev)),unique(unlist(tips_next)))

#####get list of genomes corresponding to common gut bugs + known neighbors to add to phylophlan run of novel genomes
#####common gut bugs were manually selected to add context to the novel genomes and represent several phyla.
#####FYI the representative genomes for the novel SGBs were the representative genomes provided in the "genome" column
#####of the SGB_info.tsv generated by Will's assembly pipeline.
mapping = as.data.frame(data.table::fread(file.path(input_dir,"sequences_oct22_curtis.tsv"), header=T, check.names = F))
mapping$genome_id = gsub("Curtis__","",mapping$genome_id)
common_bugs = c("10068","14535","1815","17244","4540","2318","4285","9283","1830","1836","4584","15318")
common_bugs_and_knowns = append(common_bugs,closest_knowns)
common_bugs_and_knowns = paste0("SGB",common_bugs_and_knowns)
write.table(common_bugs_and_knowns,file="common_gut_bugs_and_knowns.txt",quote=FALSE,row.names = FALSE,col.names=FALSE)

###import tree constructed from phylophlan for below visualization
common_knowns_novel_tree = read.tree("/Users/tobynbranck/Documents/pets_meta/RAxML_bestTree.Novel_tree_run_5.5.23.tre")
common_bugs_and_knowns_trun = gsub("SGB","",common_bugs_and_knowns)

### The SGB_info.tsv was generated during assembly and provides the information for which MAGs were clustered into the 
### novel SGBs
pets_novel_sgbs = as.data.frame(data.table::fread(file.path(input_dir,"SGB_info.tsv"), header=T, check.names = F))
pets_novel_sgbs$genome = gsub(".fa","",pets_novel_sgbs$genome)

novel_clusters_genome_reps = pets_novel_sgbs$genome
#for annotation purposes, must replace a few of the novel genome names since some of the genomes weren't transferred to 
#nicola's group for some reason, and thus are not in the genome--> SGB file. I'm replacing them with a different member
#belonging to the same novel SGB cluster, which is ok, since we're representing the SGBs in the tree, rather than highlighting the genome itself.
novel_clusters_genome_reps = gsub("reERR318679.bin.30","reERR318688.bin.84",novel_clusters_genome_reps)
novel_clusters_genome_reps = gsub("reERR4183454.bin.17","ERR4183468.bin.13",novel_clusters_genome_reps)
novel_clusters_genome_reps = gsub("reERR4183462.bin.30","ERR4183491.bin.20",novel_clusters_genome_reps)
novel_clusters_genome_reps = gsub("reERR878254.bin.25","reERR318685.bin.81",novel_clusters_genome_reps)

mapping_df = as.data.frame(data.table::fread(file.path(input_dir,"sequences_oct22_curtis.tsv"), header=T, check.names = F))
mapping_df$genome_id = gsub("Curtis__","",mapping_df$genome_id)

mapping_df = mapping_df[mapping_df$genome_id %in% novel_clusters_genome_reps,]

mapping_df$genome_id = gsub("reERR318688.bin.84","reERR318679.bin.30",mapping_df$genome_id)
mapping_df$genome_id = gsub("ERR4183468.bin.13","reERR4183454.bin.17",mapping_df$genome_id)
mapping_df$genome_id = gsub("ERR4183491.bin.20","reERR4183462.bin.30",mapping_df$genome_id)
mapping_df$genome_id = gsub("reERR318685.bin.81","reERR878254.bin.25",mapping_df$genome_id)

test_names = ifelse(match(common_knowns_novel_tree$tip.label,mapping_df$genome_id),mapping_df$sgb_id[match(common_knowns_novel_tree$tip.label,mapping_df$genome_id)],common_knowns_novel_tree$tip.label)
test_names[is.na(test_names)] <- common_knowns_novel_tree$tip.label[is.na(test_names)]
common_knowns_novel_tree$tip.label = test_names

#### Note: removing tips from tree: c("SRR8638156.bin.9","SGB105987__M1683799196__Curtis__reERR318702.bin.82","SGB106059__M1131443709__Curtis__reERR318703.bin.138")
#SRR8638156.bin.9 because it didn't end up in the list of genomes provided to Nicola's group (reason for 19 novel SGBs instead of 20), and the other two because they were included as "closest knowns
#to two of the 9 novel SGBs that were identified in metaphlan but are also actually 2 of the 10 novel SGBs that we found that were NOT identified in metaphlan - so their 
#corresponding SGBs are already represented in the tree
common_knowns_novel_tree = ape::drop.tip(common_knowns_novel_tree, c("SRR8638156.bin.9","SGB105987__M1683799196__Curtis__reERR318702.bin.82","SGB106059__M1131443709__Curtis__reERR318703.bin.138"))

common_knowns_novel_tree$tip.label = gsub('_.*','',common_knowns_novel_tree$tip.label)
common_knowns_novel_tree$tip.label = gsub("SGB","",common_knowns_novel_tree$tip.label)

test_novel_metadata = tree_metadata[tree_metadata$tips %in% common_knowns_novel_tree$tip.label,]
novel_tips_not_in_metaph = setdiff(common_knowns_novel_tree$tip.label,test_novel_metadata$tips)

rownames(test_novel_metadata) = test_novel_metadata$variable
rownames(test_novel_metadata) = gsub("SP__","t__",rownames(test_novel_metadata))
source("/Users/tobynbranck/Documents/keep_SGB.R")
test_novel_metadata = keep_SGB(test_novel_metadata)
test_novel_metadata$abbrev_tax = rownames(test_novel_metadata)

test_novel_metadata$tips = paste0("SGB",test_novel_metadata$tips)
common_knowns_novel_tree$tip.label = paste0("SGB",common_knowns_novel_tree$tip.label)

excluded_SGB_info = as.data.frame(read.csv(file.path(input_dir,"SGB_taxa.csv"), header=F, check.names = F))
names(excluded_SGB_info) = c("tips","taxa")
rownames(excluded_SGB_info) = excluded_SGB_info$taxa
excluded_SGB_info = keep_SGB(excluded_SGB_info)
rownames(excluded_SGB_info) = gsub("SGBSGB","SGB",rownames(excluded_SGB_info))
excluded_SGB_info$phylum = gsub(".*p__","",excluded_SGB_info$taxa)
excluded_SGB_info$phylum = gsub("\\|.*","",excluded_SGB_info$phylum)

#now need to add some kind of annotation to test_novel_metadata for the 10 novel SGBs that were not 
#found in metaphlan (therefore is not included in the tree_metadata file)
test_novel_metadata[nrow(test_novel_metadata) + 10,] = NA
test_novel_metadata$tips[38:47] = paste0("SGB",excluded_SGB_info$tips)
test_novel_metadata$knownness[38:47] = "novel"
test_novel_metadata$novel[38:47] = "novel"
test_novel_metadata$phylum[38:47] = excluded_SGB_info$phylum
test_novel_metadata$abbrev_tax[38:47] = rownames(excluded_SGB_info)
rownames(test_novel_metadata)[38:47] = test_novel_metadata$abbrev_tax[38:47]


common_knowns_novel_tree$tip.label = test_novel_metadata$abbrev_tax[match(common_knowns_novel_tree$tip.label,test_novel_metadata$tips)]

test_novel_metadata$SGB = test_novel_metadata$tips
test_novel_metadata$tips = rownames(test_novel_metadata)

new_novel_tree = ggtree(common_knowns_novel_tree,layout="rectangular") %<+% test_novel_metadata + 
  geom_tippoint(aes(shape=novel,color=phylum),size=4.5,alpha=1,stroke=.1) + geom_tiplab(aes(label=label,color=novel),hjust=-.05,size=3) +
  scale_fill_manual(values=c("#FF69B4","#008080","#9BCD9B","#CDC1C5"),labels = c("Human","Dog","Firmicutes","Cat","Proteobacteria"),name="") + 
  scale_color_manual(values=c("#FF69B4","#008080","#9BCD9B","gray48","black","#CDC1C5","#808080")) +
  xlim(0, .82) + 
  theme_tree2() +
  xlab("evolutionary distance") +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=11)) +
  theme(axis.text.x   = element_text(family = "Arial"),
        axis.title.x  = element_text(family = "Arial"),
        axis.title.y  = element_text(family = "Arial"))
  

novel_melt = test_novel_metadata[,c("tips","prev_dog","prev_cat","prev_human")]
novel_melt = as.data.frame(reshape2::melt(novel_melt))

#order based on tree order 
novel_melt$tips <- factor(novel_melt$tips, levels=get_taxa_name(new_novel_tree))
novel_melt$tips = fct_rev(novel_melt$tips)

testplot = ggplot(novel_melt[which(novel_melt$value>0 | is.na(novel_melt$value)),], aes(x = value, y= tips,color=variable,group=tips)) +
  #geom_line(color="grey") +
  geom_point(size = 4.5,alpha=.8)  +
  cowplot::theme_minimal_hgrid() +
  ylab("") + xlab("prevalence") +
  scale_color_manual(values=c("#E69F00","indianred","#56B4E9"),labels=c('Dog', 'Cat','Human')) +
  scale_y_discrete(expand = c(0,.5)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.title= element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=11),
        legend.position = "none") + 
  theme(axis.text.x   = element_text(family = "Arial"),
        axis.title.x  = element_text(family = "Arial"),
        axis.title.y  = element_text(family = "Arial"))

cairo_pdf("novel_tree_prevalence_complete.pdf",width=16,height = 11)
cowplot::plot_grid(new_novel_tree,testplot,label_size=12,align = 'h',rel_widths = c(2,1.25))
dev.off()