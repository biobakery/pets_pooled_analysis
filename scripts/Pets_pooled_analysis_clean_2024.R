# This script takes in as input the metaphlan 4 relative abundance tables from companion animal samples, HMP1-II, and Madagascar cohorts
# and the formatted metadata file generated from metadata.R.
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
library(scatterplot3d)
library(reshape2)

#### Set the file path for input & output files
input_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/inputs")
output_dir = file.path("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/outputs")

#### create taxonomic mapping for Oct22 and Jun23
Oct22_assign = as.data.frame(data.table::fread(file.path(input_dir,"SGB.Oct22.txt"), header=T, check.names = F))
Oct22_assign = Oct22_assign[Oct22_assign$`# Label` == 'SGB',]
Oct22_assign$SGB = paste(Oct22_assign$`# Label`,Oct22_assign$ID,sep='')
Oct22_assign = Oct22_assign[,c('SGB','Assigned taxonomy')]
names(Oct22_assign)[2] = 'Oct22_taxonomy'

Jun23_assign = as.data.frame(data.table::fread(file.path(input_dir,"SGB.Jun23.txt"), header=T, check.names = F))
Jun23_assign = Jun23_assign[Jun23_assign$`# Label` == 'SGB',]
Jun23_assign$SGB = paste(Jun23_assign$`# Label`,Jun23_assign$ID,sep='')
Jun23_assign = Jun23_assign[,c('SGB','Assigned taxonomy')]
names(Jun23_assign)[2] = 'Jun23_taxonomy'

map = merge(Oct22_assign,Jun23_assign,by="SGB") 
write.csv(map,file.path(output_dir,"mapping_Oct22_Jun23.tsv"),row.names = FALSE)

#### Import the companion animal metaphlan 4 taxonomic profiles relative abundance
bulk_samps_metaphlan4 = as.data.frame(data.table::fread(file.path(input_dir,"pets_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
names(bulk_samps_metaphlan4) = gsub("_taxonomic_profile","",names(bulk_samps_metaphlan4))
rownames(bulk_samps_metaphlan4) = bulk_samps_metaphlan4[,1]

#### As the Yarlagadda et al. and Coelho et al. datasets were added post-processing of the original samples, we are importing and merging here with the original. Out of the 107 samples in the Yarlagadda et al. data, the authors provided duplicated samples (n=18), which we are removing to avoid redundancy.
yar_metaphlan = as.data.frame(data.table::fread(file.path(input_dir,"Yarlagadda_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
names(yar_metaphlan) = gsub("_taxonomic_profile","",names(yar_metaphlan))
rownames(yar_metaphlan) = yar_metaphlan[,1]

#yar duplicates check analysis - this analysis shows that the duplicates are correlated with each other and justifies dropping the duplicates (except for SAADFHOND1, which seems to have been contaminated with a different sample, so we are dropping it completely)
yar_dup_check = as.data.frame(data.table::fread("/home/tobynb/shared-folder/tbranck/Pets_pooled_analysis/yar_dup_check.txt"), header=T, check.names = F)
yar_check_df = yar_metaphlan[,colnames(yar_metaphlan) %in% yar_dup_check$sample]
yar_check_df_t__ = yar_check_df %>% filter(grepl("t__SGB", rownames(yar_check_df)))
yar_check_df_t__ = yar_check_df_t__/100
colnames(yar_check_df_t__) = yar_dup_check$Individual[match(colnames(yar_check_df_t__),yar_dup_check$sample)]
col<- colorRampPalette(c("blue", "white", "red"))(20)
corr = cor(yar_check_df_t__,method="spearman")
pheatmap::pheatmap(corr) #this showed us that this study's duplicates are very similar to each other, so we are dropping duplicates (dropping duplicates corresponding to 'SAADFHOND1', since these duplicates look nothing like eachother and one is probably contaminated, but we don't have information on which that would be)

#bring in list of the duplicates and dropping the sample SRR14842322 (animal 'SAADFHOND1') 
yar_duplicates = as.data.frame(data.table::fread(file.path(input_dir,"PRJNA738608_duplicate.txt"), header=F, check.names = F))
yar_metaphlan = yar_metaphlan[,-which(names(yar_metaphlan) %in% yar_duplicates$V1)]
yar_metaphlan$SRR14842322 = NULL

coelho_metaphlan = as.data.frame(data.table::fread(file.path(input_dir,"Coelho_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
names(coelho_metaphlan) = gsub("_taxonomic_profile","",names(coelho_metaphlan))
rownames(coelho_metaphlan) = coelho_metaphlan[,1]

### Merge all animal taxonomic profiles
temp_join = merge(x = bulk_samps_metaphlan4, y = yar_metaphlan, all = TRUE)
metaphlan4 = merge(x = temp_join, y = coelho_metaphlan, all = TRUE)

rownames(metaphlan4) = metaphlan4$'# taxonomy'
metaphlan4[is.na(metaphlan4)] <- 0

#### The keep_SGB function renames the taxonomy to the lowest known taxonomic level and preserves the SGB number (function written by Jacob Nearing)
keep_SGB <-function(bug4){
  temp <- bug4[grepl("t__", rownames(bug4)),]
  SGB_name <- rownames(temp) 
  #get suffix
  SGB_name_suffix <- gsub(".*t__", "", SGB_name)
  #get prefix
  SGB_name_prefix <- sapply(SGB_name, function(x) str_extract(x, "^(.*?)(?=(\\|g__GGB|\\|o__OFGB|\\|c__CFGB|\\|p__PFGB|\\|t__))"))
  #fix prefix to the last known assignment
  SGB_name_prefix <- gsub('.*[kpcofgst]__', "", SGB_name_prefix)
  #combined
  SGB_comb <- paste(SGB_name_prefix, paste0("SGB",SGB_name_suffix,sep=""), sep="_")
  rownames(temp) <- SGB_comb
  return(temp)
}

# Since the metaphlan4 table is stratified by taxonomic level, this next step greps for the SGBs only, which are noted by the "t__" prefix
# replaces the Oct22 taxonomic assignments with those from Jun23, puts on a 0-1 scale (default is 0-100), and renormalizes (since we removed the unclassified row)
metaphlan4_t__ = metaphlan4 %>% filter(grepl("t__SGB", rownames(metaphlan4)))
metaphlan4_t__$'# taxonomy' = gsub('_group','',metaphlan4_t__$'# taxonomy') #new
metaphlan4_t__$SGB = paste('SGB',gsub('.*\\|t__SGB','',metaphlan4_t__$'# taxonomy'),sep='') # new
metaphlan4_t__$'# taxonomy' = map$Jun23_taxonomy[match(metaphlan4_t__$SGB,map$SGB)] # new
#where SGB == 'SGB97147', keep Oct22 name because this SGB is not in the Jun23 db
metaphlan4_t__$'# taxonomy'[metaphlan4_t__$SGB == 'SGB97147'] = rownames(metaphlan4_t__)[metaphlan4_t__$SGB == 'SGB97147'] #new
#where SGB == 'SGB102330', keep Oct22 name because this SGB is not in the Jun23 db
metaphlan4_t__$'# taxonomy'[metaphlan4_t__$SGB == 'SGB102330'] = rownames(metaphlan4_t__)[metaphlan4_t__$SGB == 'SGB102330'] #new
metaphlan4_t__$SGB = NULL #new
rownames(metaphlan4_t__) = metaphlan4_t__$'# taxonomy' #new

metaphlan4_t__[,2:dim(metaphlan4_t__)[2]] = metaphlan4_t__[,2:dim(metaphlan4_t__)[2]]/100
metaphlan4_t__ = cbind(metaphlan4_t__$'# taxonomy',sweep(metaphlan4_t__[,2:dim(metaphlan4_t__)[2]],2,colSums(metaphlan4_t__[,2:dim(metaphlan4_t__)[2]]),'/'))

#### Bring in metadata formatted using the metadata.R script
metadata_complete = as.data.frame(data.table::fread(file.path(input_dir,"metadata.csv"), header=T, check.names = F))
rownames(metadata_complete) = metadata_complete$sample_id_metaphlan

######mean and standard deviation of weights and ages######
cat_meta = metadata_complete[metadata_complete$species=='cat',]
dog_meta = metadata_complete[metadata_complete$species=='dog',]
#ages
cat_age_mean = mean(cat_meta$age)
cat_age_stnddev = sd(cat_meta$age)
dog_age_mean = mean(dog_meta$age,na.rm=TRUE)
dog_age_stnddev = sd(dog_meta$age,na.rm=TRUE)
#weights
cat_weight_mean = mean(cat_meta$weight,na.rm=TRUE)
cat_weight_sd = sd(cat_meta$weight,na.rm=TRUE)
dog_weight_mean = mean(dog_meta$weight,na.rm=TRUE)
dog_weight_sd = sd(dog_meta$weight,na.rm=TRUE)
############################################################

# Substituting t__SGB for SP__ (as in "subspecies") so that I can grep for things that are unknown versus known SGBs
names(metaphlan4_t__)[1] = 'Taxonomy'
metaphlan4_t__$Taxonomy = gsub("t__SGB","SP__",metaphlan4_t__$Taxonomy)
#adding the unknown row (and putting the unclassified on 0-1 scale)
metaphlan4_t__[nrow(metaphlan4_t__) + 1,] = metaphlan4[c("UNCLASSIFIED"),]
metaphlan4_t__[nrow(metaphlan4_t__),][2:ncol(metaphlan4_t__)] = metaphlan4_t__[nrow(metaphlan4_t__),][2:ncol(metaphlan4_t__)]/100

# Bringing in MAGs for MAG quality plots
MAGs = as.data.frame(data.table::fread(file.path(input_dir,"checkm_qa_and_n50.tsv"), header=T, check.names = F))
MAGs$sample = sub("\\..*", "", MAGs$bin_id)

# Remove duplicated samples from the MAGs data; see explanation of duplicates below
yar_duplicates = as.list(yar_duplicates$V1)
MAGs <- MAGs[ ! MAGs$sample %in% yar_duplicates, ]

# Matching MAGs' corresponding samples to samples in metadata
MAGs$study = metadata_complete$study_readable[match(MAGs$sample,metadata_complete$sample_id_metaphlan)]
MAGs = MAGs %>% mutate(sample = ifelse((!study %in% c("DietInt_Study2", "DietInt_Study4")), 
                                   gsub("_L00[1-9]", "", sample), sample))

# Keep only samples that exist in metadata (more samples were assembled than what we are using in the analysis - i.e., duplicated samples) 
MAGs = MAGs[MAGs$sample %in% metadata_complete$sample_id_metaphlan,]

# Checks for and removes samples (columns) if 100% unknown or 0% relative abundance
all_unknown_list = names(metaphlan4_t__[,metaphlan4_t__[nrow(metaphlan4_t__),]>=1])
metaphlan4_t__ = metaphlan4_t__[, !colnames(metaphlan4_t__) %in% all_unknown_list]
MAGs <- MAGs[ ! MAGs$sample %in% all_unknown_list, ]
write.csv(MAGs,file.path(output_dir,"pets_MAGs_formatted_for_SGBtree_Jun23tax.csv"))

# Matching with metadata: the corrected metadata provided by Hill's does not contain duplicates
# of the samples that were provided by Hill's and processed by us (duplicates have "_b" in sample name)
# for instance, the samples "DC-029_b_S85" and "DC-029_S85" were originally provided by Hills and
# processed by us, but the new metadata provided by Hill's does not contain "DC-029_b_S85" and we
# are removing it from the dataset. Also removed here are the samples that were part of a sequencing
# depth experiment, i.e., samples were deep sequenced, 6M, and 3M reads. we are removing those that 
# were sequenced at lesser depths (6M and 3M) 
taxa = metaphlan4_t__$Taxonomy
metaphlan4_t__ = metaphlan4_t__[,colnames(metaphlan4_t__) %in% rownames(metadata_complete)]
rownames(metaphlan4_t__) = taxa
write.csv(head(metaphlan4_t__,-1),file.path(output_dir,"pets_metaphlan_v4_formatted_SGBtree_Jun23tax.csv")) #saving formatted taxonomic profiles without the UNCLASSIFIED column

### Getting raw data table (metaphlan output that has taxonomy stratified) that has only the samples we want for the Supplement
bulk = as.data.frame(data.table::fread(file.path(input_dir,"pets_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
yarlagadda = as.data.frame(data.table::fread(file.path(input_dir,"Yarlagadda_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
coelho = as.data.frame(data.table::fread(file.path(input_dir,"Coelho_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
temp_join_supp = merge(x = bulk, y = yarlagadda, all = TRUE)
metaphlan4_supp = merge(x = temp_join_supp, y = coelho, all = TRUE)
rownames(metaphlan4_supp) = metaphlan4_supp$'# taxonomy'
metaphlan4_supp$'# taxonomy' = NULL
metaphlan4_supp[is.na(metaphlan4_supp)] <- 0
names(metaphlan4_supp) = gsub("_taxonomic_profile","",names(metaphlan4_supp))
metaphlan4_supp = metaphlan4_supp[,colnames(metaphlan4_supp) %in% colnames(metaphlan4_t__)]
write.csv(metaphlan4_supp,file.path(output_dir,"rawPets_Metaphlan4_FilteredSamps.csv"))

# Get number MAGs categorized into each quality bin
table(MAGs$quality)

## number of low, medium, and high quality MAGs ##
checkm_pal <- c(
  "high_quality"="midnightblue",
  "medium_quality"="skyblue",
  "low_quality"="grey"
)
checkm_scatter = ggplot(MAGs, aes(x=contamination, y=completeness)) + 
  geom_point(aes(fill=quality), pch=21, colour="black") +
  geom_hline(yintercept=c(50, 90), linetype="dashed", color = "red", linewidth=1) +
  geom_vline(xintercept=c(5, 10), linetype="dashed", color = "red", linewidth=1) +
  theme_classic() +
  theme(
    axis.title = element_text(size=12.5, face="bold"),
    axis.text = element_text(size=12.5),
    legend.title = element_text(size=10, face="bold"),
    legend.text = element_text(size=10),
    legend.justification = c(1, 0),
    legend.position = "inside",
    legend.position.inside = c(.98,.2),
    legend.background = element_rect(colour ="black", linetype = "solid", linewidth=0.5)) +
  labs(x="Contamination (sqrt (%))", y="Completeness (%)", fill="MAG quality") +
  scale_x_sqrt(breaks=c(1,5,10,100)) +
  scale_y_continuous(breaks=c(50,90,100)) +
  scale_fill_manual(values = checkm_pal) +
  guides(fill = guide_legend(override.aes = list(size=2.5)))

cairo_pdf(file.path(output_dir,"MAG_quality_scatter_Jun23tax.pdf"),width=6,height=5)
checkm_scatter
dev.off()

## Number of low, medium, and high quality MAGs per host species
metadata_complete$sample = rownames(metadata_complete)
test = MAGs[,c('sample','quality')]
new = merge(test,metadata_complete,by="sample")
bar_df = as.data.frame(new %>% dplyr::group_by(species,quality) %>% tally())

bar_df$quality <- factor(bar_df$quality, levels=c('high_quality', 'medium_quality', 'low_quality'))

checkm_bar = ggplot(bar_df, aes(x=species, y=n, fill=quality, label=comma(n))) +
  geom_bar(stat="identity", position="dodge", colour="black") +
  geom_text(position=position_dodge(width=0.9), vjust=-0.25) +
  theme_classic() +
  theme(
    axis.title = element_text(size=12.5, face="bold"),
    axis.text = element_text(size=12.5),
    legend.position = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),label=comma) +
  scale_fill_manual(values = checkm_pal) +
  labs(x="Species", y="Number of MAGs")

cairo_pdf(file.path(output_dir,"MAG_quality_BAR_Jun23tax.pdf"),width=6,height=5)
checkm_bar
dev.off()

#pre-filtering, how many unknown vs. known SGBs?
sum(!grepl("_SGB",rownames(metaphlan4_t__)))
sum(grepl("_SGB", rownames(metaphlan4_t__)))

#### Proportion of uSGBs per sample #####
#first, count number of rows greater than zero (subtracting 1 so that we do not count the UNKNOWN column)
total_greater_zero = colSums(metaphlan4_t__ > 0) - 1
#add a column that denotes whether uSGB or kSGB
metaphlan4_t__ = metaphlan4_t__ %>% mutate(SGB_novelty = case_when(grepl("_SGB", rownames(metaphlan4_t__)) ~ "uSGB"))
metaphlan4_t__$SGB_novelty <- metaphlan4_t__$SGB_novelty %>% replace_na('kSGB')
table(metaphlan4_t__$SGB_novelty)
#for each column, count how many occurrences of greater than 0 and uSGB
num_uSGB_lst = c()
for (i in 1:ncol(metaphlan4_t__[,-1])) {
  sum_per_sample = sum(metaphlan4_t__[,i] > 0 & metaphlan4_t__$SGB_novelty == 'uSGB', na.rm=TRUE)
  num_uSGB_lst = append(num_uSGB_lst,sum_per_sample)
}

pdf(file.path(output_dir,"Total_#SGBs_per_sample.pdf"))
hist(total_greater_zero,breaks = 30)
dev.off()

total_greater_zero_df = as.data.frame(total_greater_zero)
total_greater_zero_df$species = metadata_complete$species[match(row.names(total_greater_zero_df), row.names(metadata_complete))]

proportion_uSGB = num_uSGB_lst/total_greater_zero
proportion_uSGB[ncol(metaphlan4_t__)] = NA
metaphlan4_t__ = rbind(metaphlan4_t__,proportion_uSGB)

proportion_uSGB_df = as.data.frame(proportion_uSGB)
metadata_complete$proportion <- proportion_uSGB_df$proportion_uSGB[match(row.names(metadata_complete), row.names(proportion_uSGB_df))]
study_color_vec = c("indianred","#E69F00","#E69F00","#E69F00","#E69F00","indianred","#E69F00","indianred","indianred","#E69F00",
                    "#E69F00","#E69F00","#E69F00","#E69F00","#E69F00","indianred","#E69F00","#E69F00","#E69F00","darkorchid3")

rownames(metaphlan4_t__) = gsub("SP_","t_",rownames(metaphlan4_t__))

#### Next section sets up data to view the number of shared vs. unique SGB
df_for_overlap = metaphlan4_t__
df_for_overlap = head(df_for_overlap,-2)
df_for_overlap$SGB_novelty = NULL

df_for_overlap = as.data.frame(t(df_for_overlap))
df_for_overlap$host = metadata_complete$species[match(rownames(df_for_overlap),metadata_complete$sample)]
df_for_overlap$sample = rownames(df_for_overlap)

test = melt(df_for_overlap,id=c("sample","host"))
test_for_human = test

#### Importing human host metaphlan v4 output ####
human_metaphlan4 = as.data.frame(data.table::fread(file.path(input_dir,"hmp1_II_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))
madagascar_metaphlan4 = as.data.frame(data.table::fread(file.path(input_dir,"madagascar_metaphlan_taxonomic_profiles.tsv"), header=T, check.names = F))

# Map Madagascar sample names from SRA metadata to match sample names used in our metadata
names(madagascar_metaphlan4) = gsub("_taxonomic_profile","",names(madagascar_metaphlan4))
madagascar_samp_mapping = as.data.frame(data.table::fread(file.path(input_dir,"madagascar_SraRunTable.txt"), header=T, check.names = F))
colnames(madagascar_metaphlan4) = madagascar_samp_mapping$`Library Name`[match(colnames(madagascar_metaphlan4),madagascar_samp_mapping$Run)]
colnames(madagascar_metaphlan4)[1] = '# taxonomy'

# Merge both HMP1-II and Madagascar dfs into one "human" df
human_metaphlan4 = merge(human_metaphlan4, madagascar_metaphlan4, by = "# taxonomy", all = TRUE)
human_metaphlan4[is.na(human_metaphlan4)] = 0

names(human_metaphlan4) = gsub("_taxonomic_profile","",names(human_metaphlan4))
rownames(human_metaphlan4) = human_metaphlan4[,1]

human_metaphlan4_t__ = human_metaphlan4 %>% filter(grepl("t__SGB", rownames(human_metaphlan4)))

#### Remapping human taxonomic profiles Oct22 --> Jun23
human_metaphlan4_t__$'# taxonomy' = gsub('_group','',human_metaphlan4_t__$'# taxonomy') #new
human_metaphlan4_t__$SGB = paste('SGB',gsub('.*\\|t__SGB','',human_metaphlan4_t__$'# taxonomy'),sep='') # new
human_metaphlan4_t__$'# taxonomy' = map$Jun23_taxonomy[match(human_metaphlan4_t__$SGB,map$SGB)] # new
#where SGB == 'SGB97147', keep Oct22 name because this SGB is not in the Jun23 db
human_metaphlan4_t__$'# taxonomy'[human_metaphlan4_t__$SGB == 'SGB97147'] = rownames(human_metaphlan4_t__)[human_metaphlan4_t__$SGB == 'SGB97147'] #new
#where SGB == 'SGB102330', keep Oct22 name because this SGB is not in the Jun23 db
human_metaphlan4_t__$'# taxonomy'[human_metaphlan4_t__$SGB == 'SGB102330'] = rownames(human_metaphlan4_t__)[human_metaphlan4_t__$SGB == 'SGB102330'] #new
human_metaphlan4_t__$SGB = NULL #new

human_metaphlan4_t__$`# taxonomy` = gsub("t__SGB","SP__",human_metaphlan4_t__$`# taxonomy`)
rownames(human_metaphlan4_t__) = human_metaphlan4_t__$`# taxonomy`

human_metaphlan4_t__ = human_metaphlan4_t__[,names(human_metaphlan4_t__) %in% metadata_complete$sample]
human_metaphlan4_t__ = human_metaphlan4_t__/100
human_metaphlan4_t__ = sweep(human_metaphlan4_t__,2,colSums(human_metaphlan4_t__),'/')

write.csv(human_metaphlan4_t__,file.path(output_dir,"ALLhuman_metaphlan4_SGBs_formatted_SGBtree_Jun23.csv"))

rownames(human_metaphlan4_t__) = gsub("SP_","t_",rownames(human_metaphlan4_t__))

set_SGB_mapping = union(rownames(human_metaphlan4_t__),rownames(metaphlan4_t__))
set_SGB_mapping = as.data.frame(set_SGB_mapping)
names(set_SGB_mapping)[1] = 'full_tax'
set_SGB_mapping$SGB_knownness = ifelse(grepl("_SGB",set_SGB_mapping$full_tax),"kSGB","uSGB")
set_SGB_mapping$SGB_ID = gsub(".*t__","",set_SGB_mapping$full_tax)
set_SGB_mapping$SGB_ID = paste0("SGB",set_SGB_mapping$SGB_ID)
rownames(set_SGB_mapping) = set_SGB_mapping$full_tax
set_SGB_mapping = keep_SGB(set_SGB_mapping)
set_SGB_mapping$abbreviated_tax = rownames(set_SGB_mapping)

for_uSGB_proportion_human = human_metaphlan4_t__
for_uSGB_proportion_pets = metaphlan4_t__
for_uSGB_proportion_human$SGB = ifelse(grepl('_SGB',rownames(for_uSGB_proportion_human)),"uSGB","kSGB")
for_uSGB_proportion_pets$SGB = ifelse(grepl('_SGB',rownames(for_uSGB_proportion_pets)),"uSGB","kSGB")

#### Abundance-weighted version of proportion of uSGB per sample per host####
for_uSGB_proportion_human_weighted = for_uSGB_proportion_human
for_uSGB_proportion_pets_weighted = head(for_uSGB_proportion_pets,-2)
for_uSGB_proportion_pets_weighted$SGB_novelty = NULL

for_uSGB_proportion_pets_weighted[ for_uSGB_proportion_pets_weighted<0.00001 ] <- 0
for_uSGB_proportion_human_weighted[ for_uSGB_proportion_human_weighted<0.00001 ] <- 0

last_column = ncol(for_uSGB_proportion_pets_weighted)-1
last_column_human = ncol(for_uSGB_proportion_human_weighted)-1

pets_weighted_uSGB_prop = colSums(for_uSGB_proportion_pets_weighted[,1:last_column][for_uSGB_proportion_pets_weighted$SGB=='uSGB',])
human_weighted_uSGB_prop = colSums(for_uSGB_proportion_human_weighted[,1:last_column_human][for_uSGB_proportion_human_weighted$SGB=='uSGB',])

pets_weighted_uSGB_prop = as.data.frame(pets_weighted_uSGB_prop)
human_weighted_uSGB_prop = as.data.frame(human_weighted_uSGB_prop)

names(human_weighted_uSGB_prop)[1] = "prop_uSGB"
names(pets_weighted_uSGB_prop)[1] = "prop_uSGB"

#setting up two columns, one for host species (cat, dog, human), and one for host species with human studies stratified
pets_weighted_uSGB_prop$host1 = metadata_complete$species[rownames(metadata_complete) %in% rownames(pets_weighted_uSGB_prop)]
pets_weighted_uSGB_prop$host2 = metadata_complete$species[rownames(metadata_complete) %in% rownames(pets_weighted_uSGB_prop)]

human_weighted_uSGB_prop$host1 = "human"
human_weighted_uSGB_prop$host2 = metadata_complete$study_readable[metadata_complete$sample %in% rownames(human_weighted_uSGB_prop)]

#### Plot for proportions of uSGBs per sample, abundance weighted ####
proportion_combined = rbind(human_weighted_uSGB_prop,pets_weighted_uSGB_prop)
proportion_combined$host1 <- factor(proportion_combined$host1, levels=c("human","cat","dog"))
proportion_combined$host2 <- factor(proportion_combined$host2, levels=c("HMP1-II","Madagascar","cat","dog"))

prop_plot <- ggplot(proportion_combined, aes(x=host2, y=prop_uSGB,col=host2,alpha=host2)) + 
  geom_violin(linewidth=3) + coord_flip()

cairo_pdf("proportion_uSGBs_alldata_abund_weighted_Jun23tax.pdf",width=5,height=5)
prop_plot + geom_point(aes(col=host2),size=3,position="jitter",stroke=1) +
  theme_classic() + 
  xlab("") +
  ylab("proportion of uSGBs") + 
  scale_color_manual(values = c("#56B4E9","#0f3b50","indianred","#E69F00")) +
  theme(legend.title= element_blank(),
        legend.position = "top",
        legend.key.size = unit(.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(fill = c("#56B4E9","#0f3b50","indianred","#E69F00")))) +
  scale_alpha_manual(name = "host2", values = c(.2, .2, .2, .1)) +
  theme(legend.position = "none")
dev.off()
########################################################################

human_overlap = human_metaphlan4_t__
human_overlap = as.data.frame(t(human_overlap))
human_overlap$host = "human"
human_overlap$sample = rownames(human_overlap)
human_test = melt(human_overlap,id=c("sample","host"))

combined_test = rbind(test_for_human,human_test)

# Counting how many samples (by host) an SGB is present in at > 0.00001 relative abundance
combined_prev_abund_filtering = aggregate(value ~ host + variable, combined_test, function(x) sum(x > 0.00001, na.rm = TRUE))

names(combined_prev_abund_filtering)[3] = 'number_greater_0.00001'
num_dog = count(metadata_complete$species=='dog')[2,2]
num_cat = count(metadata_complete$species=='cat')[2,2]
num_human = count(metadata_complete$species=='human')[2,2]

combined_prev_abund_filtering$prevalence = ifelse(combined_prev_abund_filtering$host=="dog",combined_prev_abund_filtering$number_greater_0.00001/num_dog , ifelse(combined_prev_abund_filtering$host=="cat",combined_prev_abund_filtering$number_greater_0.00001/num_cat, 
                                                                                                                                                               combined_prev_abund_filtering$number_greater_0.00001/num_human))
combined_prev_abund_filtering_filtered = combined_prev_abund_filtering[combined_prev_abund_filtering$number_greater_0.00001 >=3, ] 

combined_prev_abund_filtering_filtered_dcast = dcast(combined_prev_abund_filtering_filtered, variable ~ host,value.var='number_greater_0.00001')
combined_prev_abund_filtering_filtered_dcast$variable = as.character(combined_prev_abund_filtering_filtered_dcast$variable)
write.csv(combined_prev_abund_filtering_filtered_dcast,file.path(output_dir,"SGB_Prevalence_By_Host_Jun23tax.csv"),row.names = FALSE)

#### Helicobacter table (Supplemental table) ####
helicob = combined_prev_abund_filtering_filtered_dcast[grep('Helicobacter',combined_prev_abund_filtering_filtered_dcast$variable),]
helicob$cat = helicob$cat/num_cat
helicob$dog = helicob$dog/num_dog
helicob$human = helicob$human/num_human
helicob = helicob %>% mutate(across(where(is.numeric), ~ round(., 3)))
write.csv(helicob,file.path(output_dir,"helicobacter_df_Jun23tax.csv"),row.names = FALSE)
##################################################

#### Campylobacter table (Supplemental table) ####
camp = combined_prev_abund_filtering_filtered_dcast[grep('Campylobacter_',combined_prev_abund_filtering_filtered_dcast$variable),]
camp$cat = camp$cat/num_cat
camp$dog = camp$dog/num_dog
camp$human = camp$human/num_human
camp = camp %>% mutate(across(where(is.numeric), ~ round(., 3)))
write.csv(camp,file.path(output_dir,"campylobacter_df.csv"),row.names = FALSE)
##################################################

#### Venn diagram ################################
cat_list = combined_prev_abund_filtering_filtered_dcast$variable[!is.na(combined_prev_abund_filtering_filtered_dcast$cat)]
dog_list = combined_prev_abund_filtering_filtered_dcast$variable[!is.na(combined_prev_abund_filtering_filtered_dcast$dog)]
human_list = combined_prev_abund_filtering_filtered_dcast$variable[!is.na(combined_prev_abund_filtering_filtered_dcast$human)]

library(ggVennDiagram)
library(ggplot2)
library(sf)
cat_dog = intersect(cat_list,dog_list)
cat_human = intersect(cat_list,human_list)
dog_human = intersect(dog_list,human_list)
cat_dog_human = intersect(cat_dog,human_list)

x = list(cat=cat_list,dog=dog_list,human=human_list)

library(ggvenn)
cairo_pdf(file.path(output_dir,"v4_VENN_0.00001_atLeast3samples.pdf"),width = 5,height=5)
ggvenn(x, fill_color = c("indianred", "#E69F00", "#56B4E9"), stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE, text_size = 5)
dev.off()
###################################################

#### Figure 3B and Supplemental plots for universally shared, companion animals only, host-specific; note - MEDIAN ABUNDANCE WHEN PRESENT
venn <- ggVennDiagram::Venn(x)
venn_data = process_data(venn)
test = as.data.frame(venn_region(venn_data))
shared_all = unlist(test$item[7])
shared_all_df = combined_test[combined_test$variable %in% shared_all,]
compan_animal_only = unlist(test$item[4])
compan_animal_only_df = combined_test[combined_test$variable %in% compan_animal_only,]
human_dog = unlist(test$item[6])
human_dog_df = combined_test[combined_test$variable %in% human_dog,]
human_cat = unlist(test$item[5])
human_cat_df = combined_test[combined_test$variable %in% human_cat,]
dog_ = unlist(test$item[2])
cat_ = unlist(test$item[1])
human_ = unlist(test$item[3])
uniques = c(dog_,cat_,human_)
unique_df = combined_test[combined_test$variable %in% uniques,]

### shared boxplots
shared_all_df$value[shared_all_df$value==0] = NA
shared_all_grouped = shared_all_df %>%
  dplyr::group_by(host, variable) %>%
  dplyr::summarise(median_abundance = median(value, na.rm=TRUE))
shared_all_grouped = as.data.frame(shared_all_grouped)

shared_all_grouped_trunc = shared_all_grouped %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(max_median = max(median_abundance))
shared_all_grouped_trunc = as.data.frame(shared_all_grouped_trunc)
shared_all_grouped_trunc = shared_all_grouped_trunc[order(shared_all_grouped_trunc$max_median,decreasing=TRUE),]

shared_all_grouped_trunc = shared_all_grouped_trunc %>% mutate(SGB_novelty = case_when(grepl("_SGB", shared_all_grouped_trunc$variable) ~ "uSGB"))
shared_all_grouped_trunc$SGB_novelty <- shared_all_grouped_trunc$SGB_novelty %>% replace_na('kSGB')

shared_all_grouped_trunc = shared_all_grouped_trunc %>% group_by(SGB_novelty) %>% slice(1:8) %>% ungroup()
shared_all_grouped_trunc = as.data.frame(shared_all_grouped_trunc)
rownames(shared_all_grouped_trunc) = shared_all_grouped_trunc$variable
shared_all_grouped_trunc = keep_SGB(shared_all_grouped_trunc)
shared_all_grouped_trunc$abbrev_tax = rownames(shared_all_grouped_trunc)

shared_all_grouped_trunc20 = shared_all_grouped[shared_all_grouped$variable %in% shared_all_grouped_trunc$variable,] 
shared_all_grouped_trunc20$abbrev_tax = shared_all_grouped_trunc$abbrev_tax[match(shared_all_grouped_trunc20$variable,shared_all_grouped_trunc$variable)]
shared_all_grouped_trunc20$SGB_novelty = shared_all_grouped_trunc$SGB_novelty[match(shared_all_grouped_trunc20$variable,shared_all_grouped_trunc$variable)]
shared_all_grouped_trunc20$abbrev_tax = gsub("_group","",shared_all_grouped_trunc20$abbrev_tax)
shared_all_grouped_trunc20$abbrev_tax = gsub("_"," ",shared_all_grouped_trunc20$abbrev_tax)

### companion animal shared boxplots
compan_animal_only_df$value[compan_animal_only_df$value==0] = NA
compan_animal_only_df = compan_animal_only_df[!compan_animal_only_df$host == "human",]
compan_animal_only_grouped = compan_animal_only_df %>%
  dplyr::group_by(host, variable) %>%
  dplyr::summarise(median_abundance = median(value, na.rm=TRUE))
compan_animal_only_grouped = as.data.frame(compan_animal_only_grouped)

compan_animal_only_grouped_trunc = compan_animal_only_grouped %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(max_median = max(median_abundance))
compan_animal_only_grouped_trunc = as.data.frame(compan_animal_only_grouped_trunc)
compan_animal_only_grouped_trunc = compan_animal_only_grouped_trunc[order(compan_animal_only_grouped_trunc$max_median,decreasing=TRUE),]

compan_animal_only_grouped_trunc = compan_animal_only_grouped_trunc %>% mutate(SGB_novelty = case_when(grepl("_SGB", compan_animal_only_grouped_trunc$variable) ~ "uSGB"))
compan_animal_only_grouped_trunc$SGB_novelty <- compan_animal_only_grouped_trunc$SGB_novelty %>% replace_na('kSGB')

compan_animal_only_grouped_trunc = compan_animal_only_grouped_trunc %>% group_by(SGB_novelty) %>% slice(1:8) %>% ungroup()
compan_animal_only_grouped_trunc = as.data.frame(compan_animal_only_grouped_trunc)
rownames(compan_animal_only_grouped_trunc) = compan_animal_only_grouped_trunc$variable
compan_animal_only_grouped_trunc = keep_SGB(compan_animal_only_grouped_trunc)
compan_animal_only_grouped_trunc$abbrev_tax = rownames(compan_animal_only_grouped_trunc)

compan_animal_only_grouped_trunc20 = compan_animal_only_grouped[compan_animal_only_grouped$variable %in% compan_animal_only_grouped_trunc$variable,] 
compan_animal_only_grouped_trunc20$abbrev_tax = compan_animal_only_grouped_trunc$abbrev_tax[match(compan_animal_only_grouped_trunc20$variable,compan_animal_only_grouped_trunc$variable)]
compan_animal_only_grouped_trunc20$SGB_novelty = compan_animal_only_grouped_trunc$SGB_novelty[match(compan_animal_only_grouped_trunc20$variable,compan_animal_only_grouped_trunc$variable)]
compan_animal_only_grouped_trunc20$abbrev_tax = gsub("_"," ",compan_animal_only_grouped_trunc20$abbrev_tax)

### unique species plots
unique_df$value[unique_df$value==0] = NA
unique_df = unique_df[!(unique_df$variable %in% dog_ & unique_df$host == "human"),]
unique_df = unique_df[!(unique_df$variable %in% dog_ & unique_df$host == "cat"),]

unique_df = unique_df[!(unique_df$variable %in% cat_ & unique_df$host == "human"),]
unique_df = unique_df[!(unique_df$variable %in% cat_ & unique_df$host == "dog"),]

unique_df = unique_df[!(unique_df$variable %in% human_ & unique_df$host == "dog"),]
unique_df = unique_df[!(unique_df$variable %in% human_ & unique_df$host == "cat"),]

unique_grouped = unique_df %>%
  dplyr::group_by(host, variable) %>%
  dplyr::summarise(median_abundance = median(value, na.rm=TRUE))
unique_grouped = as.data.frame(unique_grouped)

unique_grouped_trunc = unique_grouped %>%
  dplyr::group_by(host,variable) %>%
  dplyr::summarise(max_median = max(median_abundance))
unique_grouped_trunc = as.data.frame(unique_grouped_trunc)

unique_grouped_trunc = unique_grouped_trunc %>% mutate(SGB_novelty = case_when(grepl("_SGB", unique_grouped_trunc$variable) ~ "uSGB"))
unique_grouped_trunc$SGB_novelty <- unique_grouped_trunc$SGB_novelty %>% replace_na('kSGB')

unique_grouped_trunc = unique_grouped_trunc %>% dplyr::group_by(host) %>% dplyr::arrange(desc(max_median),.by_group = TRUE) %>% top_n(20,max_median)

unique_grouped_trunc = unique_grouped_trunc %>% group_by(host,SGB_novelty) %>% slice(1:2) %>% ungroup()
unique_grouped_trunc30 = unique_grouped[unique_grouped$variable %in% unique_grouped_trunc$variable,]
unique_grouped_trunc30 = as.data.frame(unique_grouped_trunc30)
rownames(unique_grouped_trunc30) = unique_grouped_trunc30$variable
unique_grouped_trunc30 = keep_SGB(unique_grouped_trunc30)
unique_grouped_trunc30$abbrev_tax = rownames(unique_grouped_trunc30)
unique_grouped_trunc30 = unique_grouped_trunc30 %>% mutate(SGB_novelty = case_when(grepl("_SGB", unique_grouped_trunc30$variable) ~ "uSGB"))
unique_grouped_trunc30$SGB_novelty <- unique_grouped_trunc30$SGB_novelty %>% replace_na('kSGB')
unique_grouped_trunc30$abbrev_tax = gsub("_"," ",unique_grouped_trunc30$abbrev_tax)

shared_plot = ggplot(shared_all_grouped_trunc20, aes(x=reorder(abbrev_tax,-median_abundance), y=log10(median_abundance),color=host)) + 
  geom_point(size=3,alpha=.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  scale_color_manual(values = c("cat" = "indianred",
                                "dog"="#E69F00",
                                "human"="#56B4E9")) +
  theme(legend.position="none") +
  facet_wrap(~ SGB_novelty,ncol=2,scales='free_x') + 
  scale_y_continuous(limits=c(-4.5,-1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none') +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

ca_plot = ggplot(compan_animal_only_grouped_trunc20, aes(x=reorder(abbrev_tax,-median_abundance), y=log10(median_abundance),color=host)) + 
  geom_point(size=3,alpha=.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  scale_color_manual(values = c("cat" = "indianred",
                                "dog"="#E69F00",
                                "human"="#56B4E9")) +
  #theme(legend.position="none") +
  facet_wrap(~ SGB_novelty,ncol=2,scales='free_x') + 
  scale_y_continuous(limits=c(-4.5,-1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

unique_plot = ggplot(unique_grouped_trunc30, aes(x=reorder(abbrev_tax,-median_abundance), y=log10(median_abundance),color=host)) + 
  geom_point(size=3,alpha=.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  scale_color_manual(values = c("cat" = "indianred",
                                "dog"="#E69F00",
                                "human"="#56B4E9")) +
  theme(legend.position="none") +
  facet_wrap(~ SGB_novelty,ncol=2,scales='free_x') + 
  scale_y_continuous(limits=c(-4.5,-1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none',
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  ylab("") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

cairo_pdf(file.path(output_dir,"nicheSpecificity_top_panel_Jun23tax.pdf"),width = 10,height = 5)
cowplot::plot_grid(shared_plot,ca_plot,unique_plot,nrow=1,ncol = 3,rel_widths = c(1,1,1),align='hv')
dev.off()