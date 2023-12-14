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
setwd("/Users/tobynbranck/Documents/pets_meta/")

# Import the companion animal metaphlan 4 taxonomic profiles relative abundance
metaphlan4 = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/metaphlan_v4_Oct22/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F))
names(metaphlan4) = gsub("_taxonomic_profile","",names(metaphlan4))
rownames(metaphlan4) = metaphlan4[,1]

# The keep_SGB function renames the taxonomy to the lowest known taxonomic level and preserves the SGB number (function written by Jacob Nearing)
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

# since the metaphlan4 table is stratified by taxonomic level, this next step greps for the SGBs only, which are noted by the "t__" prefix
# puts on a 0-1 scale (default is 0-100) and renormalizes (since we removed the unclassified row)
metaphlan4_t__ = metaphlan4 %>% filter(grepl("t__SGB", rownames(metaphlan4)))
metaphlan4_t__[,2:dim(metaphlan4_t__)[2]] = metaphlan4_t__[,2:dim(metaphlan4_t__)[2]]/100
metaphlan4_t__ = cbind(metaphlan4_t__$'# taxonomy',sweep(metaphlan4_t__[,2:dim(metaphlan4_t__)[2]],2,colSums(metaphlan4_t__[,2:dim(metaphlan4_t__)[2]]),'/'))

# Bring in metadata formatted using the metadata.R script
metadata_complete = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/all_metadata_preformatted.csv", header=T, check.names = F))
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
##########################################################

##substituting t__SGB for SP__ (as in "subspecies") so that I can grep for things that are unknown versus known SGBs
names(metaphlan4_t__)[1] = 'Taxonomy'
metaphlan4_t__$Taxonomy = gsub("t__SGB","SP__",metaphlan4_t__$Taxonomy)
#adding the unknown row (and putting the unclassified on 0-1 scale)
metaphlan4_t__[nrow(metaphlan4_t__) + 1,] = metaphlan4[4,]
metaphlan4_t__[nrow(metaphlan4_t__),][2:2549] = metaphlan4_t__[nrow(metaphlan4_t__),][2:2549]/100

###Bringing in MAGs for MAG quality plots
MAGs = as.data.frame(data.table::fread("/Users/tobynbranck/Downloads/checkm_qa_and_n50.tsv", header=T, check.names = F))
MAGs$sample = sub("\\..*", "", MAGs$bin_id)

#remove duplicated samples from the MAGs data; see explanation of duplicates below
duplicates = as.data.frame(data.table::fread("/Users/tobynbranck/Downloads/duplicated_samples.csv", header=T, check.names = F))
duplicates = as.list(duplicates$value)
MAGs <- MAGs[ ! MAGs$sample %in% duplicates, ]

#rename MAGs' corresponding samples to match metadata
MAGs$sample_match = MAGs$sample
MAGs$sample_match = gsub("re","",MAGs$sample_match)
MAGs$sample_match = gsub("_1","",MAGs$sample_match)
MAGs$study = metadata_complete$study_readable[match(MAGs$sample_match,metadata_complete$sample_id)]
MAGs = MAGs %>% mutate(sample = ifelse(!(study %in% c("DietInt_Study2", "DietInt_Study4")), 
                                   gsub("_L00[1-9]", "", sample), sample))

#keep only samples that exist in metadata 
MAGs = MAGs[MAGs$sample %in% metadata_complete$sample_id_metaphlan,]

#checks for and removes samples (columns) if 100% unknown or 0% relative abundance
all_unknown_list = names(metaphlan4_t__[,metaphlan4_t__[nrow(metaphlan4_t__),]>=1])
metaphlan4_t__ = metaphlan4_t__[, !colnames(metaphlan4_t__) %in% all_unknown_list]
MAGs <- MAGs[ ! MAGs$sample %in% all_unknown_list, ]

write.csv(MAGs,"/Users/tobynbranck/Documents/pets_meta/pets_MAGs_formatted_for_SGBtree.csv")
#write.csv(meta,"/Users/tobynbranck/Documents/pets_meta/pets_meta_formatted_for_SGBtree.csv")

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
write.csv(metaphlan4_t__[1:1483,],"/Users/tobynbranck/Documents/pets_meta/pets_metaphlan_v4_formatted_SGBtree.csv")

### getting raw datatable that has only the samples we want for the Supplement
metaphlan4_supp = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/metaphlan_v4_Oct22/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F))
names(metaphlan4_supp) = gsub("_taxonomic_profile","",names(metaphlan4_supp))
rownames(metaphlan4_supp) = metaphlan4_supp[,1]
metaphlan4_supp$'# taxonomy' = NULL
metaphlan4_supp = metaphlan4_supp[,colnames(metaphlan4_supp) %in% colnames(metaphlan4_t__)]
write.csv(metaphlan4_supp,"/Users/tobynbranck/Documents/pets_meta/rawPets_Metaphlan4_FilteredSamps.csv")

#get number MAGs
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
    legend.title = element_text(size=12.5, face="bold"),
    legend.text = element_text(size=12.5),
    legend.justification = c(1, 0),
    legend.position = c(1, 0.025),
    legend.background = element_rect(colour ="black", linetype = "solid", size=0.5)) +
  labs(x="Contamination (sqrt (%))", y="Completeness (%)", fill="MAG quality") +
  scale_x_sqrt(breaks=c(1,5,10,100)) +
  scale_y_continuous(breaks=c(50,90,100)) +
  scale_fill_manual(values = checkm_pal) +
  guides(fill = guide_legend(override.aes = list(size=2.5)))

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/figures/MAG_quality_scatter.pdf",width=6,height=5)
checkm_scatter
dev.off()

## number of low, medium, and high quality MAGs per host species ##
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

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/figures/MAG_quality_BAR.pdf",width=6,height=5)
checkm_bar
dev.off()

#distribution of percent unknown
perc_unknown = tail(metaphlan4_t__, n = 1)
pdf("/Users/tobynbranck/Documents/pets_meta/figures/distribution_percent_unknown.pdf")
u = hist(as.numeric(perc_unknown),
         main="distribution of percent unknown",
         xlab="percent unknown",
         xlim=c(0,100),breaks = 20)
text(u$mids,u$counts,labels=u$counts, adj=c(0.5, -0.5))
dev.off()

#pre-filtering
#how many unknowns vs. how many knowns?
sum(!grepl("_SGB",rownames(metaphlan4_t__)))
sum(grepl("_SGB", rownames(metaphlan4_t__)))

#### proportion of uSGBs per sample #####
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

pdf("/Users/tobynbranck/Documents/pets_meta/figures/Total_#SGBs_per_sample.pdf")
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

#### overlap scores
df_for_overlap = metaphlan4_t__
df_for_overlap = df_for_overlap[1:1483,1:2423]
df_for_overlap$SGB_novelty = NULL

df_for_overlap = as.data.frame(t(df_for_overlap))
df_for_overlap$host = metadata_complete$species[match(rownames(df_for_overlap),metadata_complete$sample)]
df_for_overlap$sample = rownames(df_for_overlap)

test = melt(df_for_overlap,id=c("sample","host"))
test_for_human = test

#### human metaphlan v4 output ####
human_metaphlan4 = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/hmp1_II/metaphlan_taxonomic_profiles.tsv", header=T, check.names = F))
madagascar_metaphlan4 = as.data.frame(data.table::fread("/Users/tobynbranck/Documents/pets_meta/CM_madagascar__metaphlan-4.0.4_vOct22_CHOCOPhlAnSGB_202212.tsv", header=T, check.names = F))
#removing Madagascar samples with unknown = 100
madagascar_metaphlan4 = madagascar_metaphlan4[ , -which(names(madagascar_metaphlan4) %in% c("A26_01_1FE","A37_02_1FE","A70_04_1FE"))]
names(madagascar_metaphlan4)[1] = "# taxonomy"
human_metaphlan4 = merge(human_metaphlan4, madagascar_metaphlan4, by = "# taxonomy", all = TRUE)
human_metaphlan4[is.na(human_metaphlan4)] = 0

names(human_metaphlan4) = gsub("_taxonomic_profile","",names(human_metaphlan4))
rownames(human_metaphlan4) = human_metaphlan4[,1]

human_metaphlan4_t__ = human_metaphlan4 %>% filter(grepl("t__SGB", rownames(human_metaphlan4)))
human_metaphlan4_t__$`# taxonomy` = gsub("t__SGB","SP__",human_metaphlan4_t__$`# taxonomy`)
rownames(human_metaphlan4_t__) = human_metaphlan4_t__$`# taxonomy`

human_metaphlan4_t__ = human_metaphlan4_t__[,names(human_metaphlan4_t__) %in% metadata_complete$sample]
human_metaphlan4_t__ = human_metaphlan4_t__/100
human_metaphlan4_t__ = sweep(human_metaphlan4_t__,2,colSums(human_metaphlan4_t__),'/')

write.csv(human_metaphlan4_t__,"/Users/tobynbranck/Documents/pets_meta/ALLhuman_metaphlan4_SGBs_formatted_SGBtree.csv")
#write.csv(human_meta,"/Users/tobynbranck/Documents/pets_meta/ALLhuman_metadata_formatted_SGBtree.csv")

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
for_uSGB_proportion_human$SGB = ifelse(grepl('_SGB',rownames(for_uSGB_proportion_human)),"uSGB","kSGB") #total 1992
for_uSGB_proportion_pets$SGB = ifelse(grepl('_SGB',rownames(for_uSGB_proportion_pets)),"uSGB","kSGB") #total 1051

####abundance-weighted version of proportion of uSGB per sample per host####
for_uSGB_proportion_human_weighted = for_uSGB_proportion_human
for_uSGB_proportion_pets_weighted = for_uSGB_proportion_pets[1:1483,]
for_uSGB_proportion_pets_weighted$SGB_novelty = NULL

for_uSGB_proportion_pets_weighted[ for_uSGB_proportion_pets_weighted<0.00001 ] <- 0
for_uSGB_proportion_human_weighted[ for_uSGB_proportion_human_weighted<0.00001 ] <- 0

pets_weighted_uSGB_prop = colSums(for_uSGB_proportion_pets_weighted[,1:2422][for_uSGB_proportion_pets_weighted$SGB=='uSGB',])
human_weighted_uSGB_prop = colSums(for_uSGB_proportion_human_weighted[,1:350][for_uSGB_proportion_human_weighted$SGB=='uSGB',])

pets_weighted_uSGB_prop = as.data.frame(pets_weighted_uSGB_prop)
human_weighted_uSGB_prop = as.data.frame(human_weighted_uSGB_prop)

names(human_weighted_uSGB_prop)[1] = "prop_uSGB"
names(pets_weighted_uSGB_prop)[1] = "prop_uSGB"

#setting up two columns, one for host species (cat, dog, human), and one for host species with human studies stratified
pets_weighted_uSGB_prop$host1 = metadata_complete$species[rownames(metadata_complete) %in% rownames(pets_weighted_uSGB_prop)]
pets_weighted_uSGB_prop$host2 = metadata_complete$species[rownames(metadata_complete) %in% rownames(pets_weighted_uSGB_prop)]

human_weighted_uSGB_prop$host1 = "human"
human_weighted_uSGB_prop$host2 = metadata_complete$study_readable[metadata_complete$sample %in% rownames(human_weighted_uSGB_prop)]

####Plot for proportions of uSGBs per sample, abundance weighted###################
proportion_combined = rbind(human_weighted_uSGB_prop,pets_weighted_uSGB_prop)
proportion_combined$host1 <- factor(proportion_combined$host1, levels=c("human","cat","dog"))
proportion_combined$host2 <- factor(proportion_combined$host2, levels=c("HMP1-II","Madagascar","cat","dog"))

prop_plot <- ggplot(proportion_combined, aes(x=host2, y=prop_uSGB,col=host2,alpha=host2)) + 
  geom_violin(linewidth=3) + coord_flip()

cairo_pdf("proportion_uSGBs_alldata_abund_weighted.pdf",width=5,height=5)
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
####################################################################################

human_overlap = human_metaphlan4_t__
human_overlap = as.data.frame(t(human_overlap))
human_overlap$host = "human"
human_overlap$sample = rownames(human_overlap)
human_test = melt(human_overlap,id=c("sample","host"))

combined_test = rbind(test_for_human,human_test)

#Campylobacter vignette
# campylobacter_df = combined_test %>% filter(str_detect(variable, "Campylobacter"))
# campylobacter_df = reshape(campylobacter_df, idvar = c("sample","host"), timevar = "variable", direction = "wide")
# campylobacter_df <- campylobacter_df[campylobacter_df$`value.k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_jejuni|t__19444`<=.02,]
# ggplot(campylobacter_df, aes(x = host, y = `value.k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_jejuni|t__19444`, color = host)) +
#   geom_boxplot()
# ggplot(campylobacter_df, aes(x=`value.k__Bacteria|p__Proteobacteria|c__Epsilonproteobacteria|o__Campylobacterales|f__Campylobacteraceae|g__Campylobacter|s__Campylobacter_jejuni|t__19444`, fill=host)) +
#   geom_density()


#counting how many samples (by host) an SGB is present at > 0.00001 relative abundance
combined_prev_abund_filtering = aggregate(value ~ host + variable, combined_test, function(x) sum(x > 0.00001, na.rm = TRUE))

names(combined_prev_abund_filtering)[3] = 'number_greater_0.00001'
combined_prev_abund_filtering$prevalence = ifelse(combined_prev_abund_filtering$host=="dog",combined_prev_abund_filtering$number_greater_0.00001/2056 , ifelse(combined_prev_abund_filtering$host=="cat",combined_prev_abund_filtering$number_greater_0.00001/367, 
                                                                                                                                                               combined_prev_abund_filtering$number_greater_0.00001/350))
combined_prev_abund_filtering_filtered = combined_prev_abund_filtering[combined_prev_abund_filtering$number_greater_0.00001 >=3, ] 

combined_prev_abund_filtering_filtered_dcast = dcast(combined_prev_abund_filtering_filtered, variable ~ host,value.var='number_greater_0.00001')
combined_prev_abund_filtering_filtered_dcast$variable = as.character(combined_prev_abund_filtering_filtered_dcast$variable)
write.csv(combined_prev_abund_filtering_filtered_dcast,"/Users/tobynbranck/Documents/pets_meta/SGB_Prevalence_By_Host.csv",row.names = FALSE)

########helicobacter table (Supplemental table)###########
helicob = combined_prev_abund_filtering_filtered_dcast[grep('Helicobacter',combined_prev_abund_filtering_filtered_dcast$variable),]
helicob$cat = helicob$cat/367
helicob$dog = helicob$dog/2055
helicob$human = helicob$human/350
helicob = helicob %>% mutate(across(where(is.numeric), ~ round(., 3)))
write.csv(helicob,"/Users/tobynbranck/Documents/pets_meta/helicobacter_df.csv",row.names = FALSE)
###########################################################

########campylobacter table (Supplemental table)###########
camp = combined_prev_abund_filtering_filtered_dcast[grep('Campylobacter_',combined_prev_abund_filtering_filtered_dcast$variable),]
camp$cat = camp$cat/367
camp$dog = camp$dog/2055
camp$human = camp$human/350
camp = camp %>% mutate(across(where(is.numeric), ~ round(., 3)))
write.csv(camp,"/Users/tobynbranck/Documents/pets_meta/campylobacter_df.csv",row.names = FALSE)
###########################################################

####venn diagram ##########################################
cat_list = combined_prev_abund_filtering_filtered_dcast$variable[!is.na(combined_prev_abund_filtering_filtered_dcast$cat)]
dog_list = combined_prev_abund_filtering_filtered_dcast$variable[!is.na(combined_prev_abund_filtering_filtered_dcast$dog)]
human_list = combined_prev_abund_filtering_filtered_dcast$variable[!is.na(combined_prev_abund_filtering_filtered_dcast$human)]

library(ggVennDiagram)
library(sf)
cat_dog = intersect(cat_list,dog_list)
cat_human = intersect(cat_list,human_list)
dog_human = intersect(dog_list,human_list)
cat_dog_human = intersect(cat_dog,human_list)

x = list(cat=cat_list,dog=dog_list,human=human_list)

venn <- ggVennDiagram::Venn(x)
d <- process_data(venn)
d2 <- process_data(venn)

d2@region <- st_polygonize(d@setEdge)

col <- c(cat = "indianred", dog = "#E69F00", human = "#56B4E9")

cairo_pdf("v4_VENN_0.00001_atLeast3samples.pdf",width = 5,height=5)
ggplot() +
  geom_sf(aes(fill = name), data = venn_region(d2)) +
  geom_sf(aes(color = name), data = venn_setedge(d)) +
  geom_sf_text(aes(label = name), data = venn_setlabel(d)) +
  geom_sf_text(aes(label = count), data = venn_region(d)) +
  scale_color_manual(values = ggplot2::alpha(col,.8)) +
  scale_fill_manual(values = ggplot2::alpha(col,.5)) +
  theme_void() +
  theme(legend.position = "none")
dev.off()
###########################################################

### Figure 3B and Supplemental plots for universally shared, companion animals only, host-specific; note - MEDIAN ABUNDANCE WHEN PRESENT
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

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/shared_all_hosts_top20SGBs.pdf",width=6,height=5)
ggplot(shared_all_grouped_trunc20, aes(x=reorder(abbrev_tax,-median_abundance), y=log10(median_abundance),color=host)) + 
  geom_point(size=3,alpha=.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  scale_color_manual(values = c("cat" = "indianred",
                                "dog"="#E69F00",
                                "human"="#56B4E9")) +
  theme(legend.position="none") +
  facet_wrap(~ SGB_novelty,ncol=2,scales='free_x') + 
  scale_y_continuous(limits=c(-4,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

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

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/CompAnim_shared_hosts_top20SGBs.pdf",width=6,height=5)
ggplot(compan_animal_only_grouped_trunc20, aes(x=reorder(abbrev_tax,-median_abundance), y=log10(median_abundance),color=host)) + 
  geom_point(size=3,alpha=.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  scale_color_manual(values = c("cat" = "indianred",
                                "dog"="#E69F00",
                                "human"="#56B4E9")) +
  theme(legend.position="none") +
  facet_wrap(~ SGB_novelty,ncol=2,scales='free_x') + 
  scale_y_continuous(limits=c(-4,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

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

cairo_pdf("/Users/tobynbranck/Documents/pets_meta/unique_top15SGBs.pdf",width=6,height=5)
ggplot(unique_grouped_trunc30, aes(x=reorder(abbrev_tax,-median_abundance), y=log10(median_abundance),color=host)) + 
  geom_point(size=3,alpha=.6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  xlab("") +
  scale_color_manual(values = c("cat" = "indianred",
                                "dog"="#E69F00",
                                "human"="#56B4E9")) +
  theme(legend.position="none") +
  facet_wrap(~ SGB_novelty,ncol=2,scales='free_x') + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

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

cairo_pdf("nicheSpecificity_top_panel.pdf",width = 10,height = 5)
cowplot::plot_grid(shared_plot,ca_plot,unique_plot,nrow=1,ncol = 3,rel_widths = c(1,1,1),align='hv')
dev.off()
#####################################################################################################

####### formatting metadata for running strainphlan on the cluster##############################
metadata_strainphlan = metadata_complete[,c('sample','species')]
metadata_strainphlan$sample = paste0(metadata_strainphlan$sample,"_bowtie2")

#write.table(metadata_strainphlan, 'metadata_strainphlan.txt', append = FALSE, sep = " ", dec = ".",
#            row.names = FALSE, col.names = TRUE)
metadata_strainphlan$sample = paste0(metadata_strainphlan$sample,".pkl")

fileConn<-file("/Users/tobynbranck/Documents/pets_meta/sampskeep.txt")
writeLines(metadata_strainphlan$sample, fileConn)
close(fileConn)
################################################################################################