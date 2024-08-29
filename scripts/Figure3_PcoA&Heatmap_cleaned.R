### Pet's pooled study: Figure 3 ###
### Last Edit: 08/19/24 ###

# rm(list = ls())

#######################################################################################################################################
## Packages
library(tidyverse)
library(vegan)
library(ggplot2)
library(ComplexHeatmap)
library(viridis)
library(gridExtra)
library(cowplot)
library(grid)
library(Maaslin2)
library(ggbreak) 
library(metagMisc)

#######################################################################################################################################
## Load the combined filtered MetaPhlAn4 table on SGB level
# input (total samples n = 2989)
animal = read_csv('rename/input/pets_metaphlan_v4_formatted_SGBtree_Jun23tax.csv') # n = 2639
human = read_csv('rename/input/ALLhuman_metaphlan4_SGBs_formatted_SGBtree_Jun23.csv') # n = 350

all_metadata = read_csv("Taxonomic profiling/metadata.csv") %>% # n = 2989
  mutate(weightpc = species,
         weightpc = ifelse(study_readable == "Madagascar", "Madagascar", weightpc),
         weightpc = ifelse(study_readable == "HMP1-II", "HMP1-II", weightpc)) %>%
  column_to_rownames("sample_id_metaphlan")

#######################################################################################################################################
## Combine all data together: data transformation
# animal matrix
animal_mat = animal %>%
  column_to_rownames("...1") %>%
  as.matrix() %>%
  t() 
# combine animal matrix with metadata
animal_meta = merge(all_metadata, animal_mat, by=0) %>%
  column_to_rownames("Row.names") # n = 2639

# madagascar human matrix
human_mat = human %>%
  column_to_rownames("...1") %>%
  as.matrix() %>%
  t()
# combine madagascar human matrix with metadata
human_meta = merge(all_metadata, human_mat, by=0) %>% 
  column_to_rownames("Row.names") # n = 350

# combine transformed data together
all_meta = full_join(human_meta, animal_meta) # n = 2989
all_meta[, 13:ncol(all_meta)][is.na(all_meta[, 13:ncol(all_meta)])] = 0 # First 12 columns for metadata
rownames(all_meta) = all_meta$sample_id

# save temporary matrix and metadata
pcoa_mat = all_meta[, 13:ncol(all_meta)]
pcoa_meta = all_meta[, 1:12]

# rename SGB to the closest known taxonomic level (function shared by Jacob Nearing)
keep_SGB <- function(bug4){
  temp <- bug4
  SGB_name <- names(bug4)
  #get suffix
  SGB_name_suffix <- gsub(".*SP__", "", SGB_name)
  #get prefix
  SGB_name_prefix <- sapply(SGB_name, function(x) str_extract(x, "^(.*?)(?=(\\|g__GGB|\\|o__OFGB|\\|c__CFGB|\\|p__PFGB|\\|SP__))"))
  #fix prefix to the last known assignment
  SGB_name_prefix <- gsub('.*[kpcofgst]__', "", SGB_name_prefix)
  #combined
  SGB_comb <- paste(SGB_name_prefix, "SGB", SGB_name_suffix, sep = "_")
  names(bug4) <- SGB_comb
  return(bug4)
}

pcoa_mat_rename = keep_SGB(pcoa_mat)
# ordered final matrix: pcoa_mat_rename
# ordered final metadata: pcoa_meta

#######################################################################################################################################
## Metadata visualization (Figure 1B and Supplemental Figure 1)
# sample (Figure 1B)
species.colors = c(cat = "indianred", dog = "#E69F00", human = "#56B4E0")

cairo_pdf("rename/figure/samples_ori.pdf", height = 6, width = 11.3, fallback_resolution = 300)

ggplot(all_metadata, aes(x=fct_infreq(study_readable), fill = species)) + 
  geom_bar(stat = "count") +
  geom_text(stat='count', aes(label=after_stat(count)), position = "stack", vjust= -0.5, size = 5) + 
  theme_classic() +
  theme(legend.text=element_text(size = 18), 
        legend.title=element_text(size = 18, face = "bold"), 
        legend.position = c(0.9, 0.5),
        plot.title = element_text(size = 18), 
        axis.text.x  = element_text(size = 10, angle = 45, hjust = 1), 
        axis.text.y  = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  labs(title = NULL, x=NULL, y= "Counts", fill = "Host Species") +
  scale_fill_manual(values = species.colors) +
  scale_y_break(c(330, 850)) + 
  scale_y_continuous(limits = c(0, 920), breaks = seq(0, 920, by = 100))

dev.off()

# age (Suppl. Fig 1A)
cairo_pdf("rename/figure/age.pdf", height = 6, width = 10, fallback_resolution = 300)

ggplot(all_metadata_CA, aes(x = species, y = age, color = species, fill = species)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 2, alpha = 0.7) +
  theme_classic() +
  theme(legend.text=element_text(size = 18), 
        legend.title=element_text(size = 18, face = "bold"), 
        plot.title = element_text(size = 18), 
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid = element_blank()) +
  labs(title = "Age", x=NULL, y= "Age (years)", color = "Host Species", fill = "Host Species") +
  scale_fill_manual(values = species.colors) +
  scale_color_manual(values = species.colors)

dev.off()

# weight (Suppl. Fig 1B)
cairo_pdf("rename/figure/weight.pdf", height = 6, width = 9, fallback_resolution = 300)

ggplot(all_metadata_CA, aes(x = species, y = weight, color = species, fill = species)) +
  geom_violin(alpha = 0.5) + 
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 2, alpha = 0.7) +
  theme_classic() +
  theme(legend.text=element_text(size = 18), 
        legend.title=element_text(size = 18, face = "bold"), 
        plot.title = element_text(size = 18), 
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid = element_blank()) +
  labs(title = "Weight", x = NULL, y = "Weight (kg)", color = "Host Species", fill = "Host Species") +
  scale_fill_manual(values = species.colors) +
  scale_color_manual(values = species.colors)

dev.off()

# sex (Suppl. Fig 1C)
all_metadata_sex = all_metadata %>% mutate(sex = ifelse(is.na(sex), "NA", as.character(sex)))

sex_df = all_metadata_sex[, 4:5] %>% as.data.frame()
sex_df_NA = sex_df[which(sex_df$sex == "NA"),]

cairo_pdf("rename/figure/sex.pdf", width = 10, height = 6, fallback_resolution = 300)

ggplot() + 
  geom_bar(sex_df, mapping = aes(x = species, fill = species, alpha = sex), stat = "count") +
  geom_bar(sex_df_NA, mapping = aes(x = species, fill = "gray")) + 
  theme_classic() +
  theme(legend.text=element_text(size = 18), legend.title=element_text(size = 18, face = "bold"), 
        plot.title = element_text(size = 18), 
        axis.text = element_text(size = 18), 
        axis.title =  element_text(size = 18),
        panel.grid = element_blank()) +
  labs(title = "Sex", x=NULL, y= "Counts", fill = "Host Species", alpha = "Sex") +
  scale_fill_manual(values = c(cat = "indianred", dog = "#E69F00", human = "#56B4E0")) +
  scale_alpha_manual(values = c(Female = 1, Male = 0.3, "NA" = 0.5))

dev.off()

# housing (Suppl. Fig 1D)
all_metadata_CA = all_metadata %>% filter(species != "human")

housing.colors= c(Facility = "#F2C143", Household = "#3A8A62", Stray = "#A4514F")

cairo_pdf("rename/figure/housing.pdf", height = 6, width = 9, fallback_resolution = 300)

ggplot(all_metadata_CA, aes(x = species, fill = `Facility/Household`)) + 
  geom_bar(stat = "count") +
  theme_classic() +
  theme(legend.text = element_text(size = 18), 
        legend.title = element_text(size = 18, face = "bold"), 
        legend.position = "right",
        plot.title = element_text(size = 18), 
        axis.text  = element_text(size = 18),
        axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        axis.text.y.right = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank()) +
  labs(title = "Housing", x = NULL, y = "Counts", fill = "Housing") +
  scale_fill_manual(values = housing.colors) 

dev.off()

#######################################################################################################################################
## MaAsLin2: Feature-wise testing for test differences in microbial abundances across host species
# save dataframe used for heatmap
all = cbind(pcoa_meta, pcoa_mat_rename)
# save(all, file = paste0("Taxonomic profiling/metaphlanv4/pcoa_all.Rdata"))

# format dataset by host species into long dataframe
cat = all[which(all$species=="cat"),]
tmp.cat = cat %>%
  select(-subject_id, -study_readable, -species, -sex, -age, -weight, -antibiotics, -medications, -`Hill's/Public`, -`Facility/Household`, -weightpc) %>% 
  arrange(sample_id) %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  select(species, sample_id, relab) %>%
  add_column(rep(cat$study_readable, each = 3682)) # number of species for formatting
names(tmp.cat) = c('species', 'sample_id', 'relab', 'study_readable')

dog = all[which(all$species=="dog"),]
tmp.dog = dog %>%
  select(-subject_id, -study_readable, -species, -sex, -age, -weight, -antibiotics, -medications, -`Hill's/Public`, -`Facility/Household`, -weightpc) %>% 
  arrange(sample_id) %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  select(species, sample_id, relab) %>%
  add_column(rep(dog$study_readable, each = 3682)) # number of species for formatting
names(tmp.dog) = c('species', 'sample_id', 'relab', 'study_readable')

human = all[which(all$species=="human"),]
tmp.human = human %>%
  select(-subject_id, -study_readable, -species, -sex, -age, -weight, -antibiotics, -medications, -`Hill's/Public`, -`Facility/Household`, -weightpc) %>% 
  arrange(sample_id) %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  select(species, sample_id, relab) %>%
  add_column(rep(human$study_readable, each = 3682)) # number of species for formatting
names(tmp.human) = c('species', 'sample_id', 'relab', 'study_readable')

# filter: relab >10^5 in at least 20 samples
species.list.filter = function(df){
  n_samples = nlevels(as.factor(df$sample_id))
  df_species = df %>% 
    filter(relab > 0.00001) %>%
    group_by(species) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 20) 
  dftop = df_species$species
}

filtered_cat_list = species.list.filter(tmp.cat)
filtered_dog_list = species.list.filter(tmp.dog)
filtered_human_list = species.list.filter(tmp.human)

top_all = c(filtered_cat_list, filtered_dog_list, filtered_human_list)
prevalent_species = unique(top_all)

df_numeric = all[,13:ncol(all)]
maaslin_mat = df_numeric[, names(df_numeric) %in% prevalent_species]

pcoa_meta_meta2 = pcoa_meta %>%   
  mutate(housing = `Facility/Household`,
         housing = ifelse(species == "human", "Household", housing))

pcoa_meta_mat_dog = data.frame(pcoa_meta_meta2, maaslin_mat) %>% filter(species == "dog") %>% select(-c(1:13))
pcoa_meta_meta_dog = data.frame(pcoa_meta_meta2, maaslin_mat) %>% filter(species == "dog") %>% select(c(1:13))

# run maaslin2 (2 random effects model)
Maaslin2(
  input_data = pcoa_meta_mat_dog, 
  input_metadata = pcoa_meta_meta_dog, 
  min_prevalence = 0,
  min_abundance = 0,
  output = "rename/double_re/maaslin_output_housing_dog", 
  fixed_effects = c("housing"),
  random_effects = c("study_readable", "subject_id"),
  reference = c("housing,Facility"))

Maaslin2(
  input_data = maaslin_mat,
  input_metadata = pcoa_meta_meta2,
  min_prevalence = 0,
  min_abundance = 0,
  output = "rename/double_re/maaslin_output_dog", 
  fixed_effects = c("species", "housing"),
  random_effects = c("study_readable", "subject_id"),
  reference = c("species,dog", "housing,Facility"))

Maaslin2(
  input_data = maaslin_mat,
  input_metadata = pcoa_meta_meta2,
  min_prevalence = 0,
  min_abundance = 0,
  output = "rename/double_re/maaslin_output_cat", 
  fixed_effects = c("species", "housing"),
  random_effects = c("study_readable", "subject_id"),
  reference = c("species,cat", "housing,Facility"))

Maaslin2(
  input_data = maaslin_mat,
  input_metadata = pcoa_meta_meta2,
  min_prevalence = 0,
  min_abundance = 0,
  output = "rename/double_re/maaslin_output_human", 
  fixed_effects = c("species", "housing"),
  random_effects = c("study_readable", "subject_id"),
  reference = c("species,human", "housing,Facility"))

# combine maaslin results for SGBs, called "SGB_maaslin_FDRcorrected_2re", random effect = study_readable, subject_id
dog_ref = read.csv('rename/double_re/maaslin_output_dog/all_results.tsv', head = TRUE, sep = "\t") %>% mutate(reference = "dog")
cat_ref = read.csv('rename/double_re/maaslin_output_cat/all_results.tsv', head = TRUE, sep = "\t") %>% mutate(reference = "cat")

maaslin_df = rbind(dog_ref, cat_ref) %>% 
  filter(!value == "Household") %>%
  filter(!value == "Stray") %>%
  group_by(feature) %>%
  mutate(ordered_val = pmin(value, reference),
         ordered_ref = pmax(value, reference)) %>%
  arrange(feature, ordered_val, ordered_ref) %>%
  filter(!duplicated(ordered_val) | !duplicated(ordered_ref)) %>%
  select(-ordered_val, -ordered_ref) %>%
  ungroup() %>%
  mutate(reference = paste0(reference, ",facility")) %>%
  mutate(fdr.pval = p.adjust(pval, method = "fdr")) %>%
  arrange(fdr.pval)

write.csv(maaslin_df,
          file = paste0("rename/double_re/SGB_maaslin_FDRcorrected_2re.csv"),
          row.names = F)  

#######################################################################################################################################
## PcoA: Frequency-corrected principal coordinates analysis (Figure 3A)
# calculate brat-curtis dissimilarity 
min_relab = min(pcoa_mat_rename[pcoa_mat_rename>0]) / 10 # set 0 to be 1/10 of the minimum relative abundance = 9.999994e-09
pcoa_mat_rename_pcoa = pcoa_mat_rename # save for aitchison distance 
pcoa_mat_rename_pcoa[pcoa_mat_rename_pcoa == 0] = min_relab # matrix

dist_pcoa_mat_rename = vegdist(pcoa_mat_rename, method = "bray") # or
dist_pcoa_mat_rename = vegdist(pcoa_mat_rename_pcoa, method = "aitchison", pseudocount = min_relab) 

weights = c(case_when(pcoa_meta$species == "human" ~ 1/350,
                      pcoa_meta$species == "cat" ~ 1/367,
                      pcoa_meta$species == "dog" ~ 1/2272)) # weighted by sample size

# calculate PC1 and PC2 scores
cmd_res_pcoa_mat_rename = wcmdscale(dist_pcoa_mat_rename, 
                                    k = 2,
                                    eig = TRUE,
                                    w = weights)

# explained variation percent
100 * cmd_res_pcoa_mat_rename$eig / sum(cmd_res_pcoa_mat_rename$eig)
percent_explained = metagMisc::eig_perc(cmd_res_pcoa_mat_rename$eig, positive = T, plot = F)
PC1 = percent_explained[1] %>% round(digits = 2) # 11.28 for bray; 19.2 for aitchison
PC2 = percent_explained[2] %>% round(digits = 2) # 8.07 for bray; 10.13 for aitchison

# select first two PC
pcoa_mat_2PC = tibble(PC1 = cmd_res_pcoa_mat_rename$points[,1], 
                      PC2 = cmd_res_pcoa_mat_rename$points[,2])

# calculate weighted average scores of bugs driving the cluster of PCoA
wascores = data.frame(wascores(pcoa_mat_2PC, pcoa_mat_rename))
rownames(wascores) = gsub("_", " ", rownames(wascores), fixed=TRUE)

# prepare two inputs for PcoA 
data_wa_sub = wascores[c('Bacteroidaceae SGB 1818', # <- Enterobacteriaceae SGB 1818
                         'Atopobiaceae SGB 14350', # <- Actinobacteria SGB 14350
                         'Phocaeicola vulgatus SGB 1814',
                         'Prevotella copri SGB 1626', #<- Prevotella copri clade A SGB 1626
                         'Bifidobacterium adolescentis SGB 17244'),] # <- # Bifidobacterium adolescentis SGB 17244 group
                         # selected by relative abundance clustering
rownames(data_wa_sub) = gsub("SGB ", "SGB", rownames(data_wa_sub))

pcoa_df = data.frame(pcoa_mat_2PC, pcoa_meta)

# color setting
species.colors = c(cat = "indianred", dog = "#E69F00", `HMP1-II` = "#56B4E0", Madagascar = "#0f3b50")

# main PCoA
cairo_pdf("rename/figure/weighted_pcoa_2xhuman_species_v4.pdf", width = 10, height = 7, fallback_resolution = 300)
cairo_pdf("rename/figure/weighted_pcoa_2xhuman_species_v4_aitchison_weighted.pdf", width = 10, height = 7, fallback_resolution = 300)
#cairo_pdf("rename/figure/weighted_pcoa_2xhuman_species_v4_aitchison_unweighted.pdf", width = 10, height = 7, fallback_resolution = 300)
#cairo_pdf("rename/figure/weighted_pcoa_2xhuman_species_v4_unweighted.pdf", width = 10, height = 7, fallback_resolution = 300)

ggplot(data = pcoa_df %>% arrange(species)) + 
  geom_point(size = 3, shape = 21, aes(x = PC1, y = PC2, fill = weightpc), alpha = 0.7) + 
  geom_text(data = data_wa_sub, aes(x = PC1, y = PC2, label = rownames(data_wa_sub)), 
            vjust = 0, hjust = 0.5, size = 6, fontface = "bold", color = "black") +
  theme_classic() + 
  theme(legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16, face = "bold"), 
        legend.position = "right",
        axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values = species.colors) + 
  labs(title = NULL, x = paste0("PCoA1 [", PC1,"%]"), y = paste0("PCoA2 [",PC2,"%]"), fill = "Host Species")

dev.off()

#######################################################################################################################################
## PERMANOVA with BETADISPER
# univariate
pcoa_meta_permanova = pcoa_meta %>% mutate(housing_dog = `Facility/Household`,
                                           housing_all = `Facility/Household`)
PERMANOVA.list = c("study_readable", "species", "housing_dog", "housing_all", "age")

PERMANOVA_results = NULL

for (variable in PERMANOVA.list) {
  if (variable == "housing_dog") {
    meta.tmp = pcoa_meta_permanova %>% select(all_of(c("sample_id", "species", variable))) %>% filter(species == "dog") %>% select(-species)
  } else {
    meta.tmp = pcoa_meta_permanova %>% select(all_of(c("sample_id", variable)))
  }
  
  meta.tmp = meta.tmp[complete.cases(meta.tmp),]
  
  common.sps = intersect(rownames(pcoa_mat_rename), meta.tmp$sample_id)
  
  mat.tmp = pcoa_mat_rename[match(common.sps, rownames(pcoa_mat_rename)),]
  meta.tmp = meta.tmp[match(common.sps, meta.tmp$sample_id),]
  
  if (variable == "age") {
    variable.tmp = meta.tmp[[variable]]
  } else {
    variable.tmp = factor(meta.tmp[[variable]])
  }
  
  dist_mat = vegdist(mat.tmp, method = "bray") 
  
  dispersion = betadisper(dist_mat, variable.tmp)
  
  anova.variable = anova(dispersion) # ANOVA's p-value is not significant meaning that group dispersions are homogenous ("Null hypothesis of no difference in dispersion between groups"
  tukey.variable = TukeyHSD(dispersion) # Tukey's test can be done to see if and which groups differ in relation to their variances 
  
  fml = paste0("mat.tmp ~ ", variable)
  Adonis = adonis2(as.formula(fml), data = meta.tmp, permutations = 999, method = "bray")
  
  res = c("variable" = variable,
          "dispersion.p" = anova.variable$`Pr(>F)`[1],
          "adonis.r2" = Adonis$R2[1],
          "adonis.p" = Adonis$`Pr(>F)`[1],
          "note" = ifelse(variable == "housing_dog", "dogs only", "all hosts"))
  
  PERMANOVA_results = bind_rows(PERMANOVA_results, res)
}

write.csv(PERMANOVA_results,
          file = paste0("rename/SGB_PERMANOVA_univariate.csv"),
          row.names = F)  

# multivariate
# housing_dog
meta.tmp = pcoa_meta_permanova %>% select(all_of(c("sample_id", "study_readable", "housing_dog", "age", "species"))) %>% filter(species == "dog") %>% select(-species)
meta.tmp = meta.tmp[complete.cases(meta.tmp),]
common.sps = intersect(rownames(pcoa_mat_rename), meta.tmp$sample_id)

mat.tmp = pcoa_mat_rename[match(common.sps, rownames(pcoa_mat_rename)),]
meta.tmp = meta.tmp[match(common.sps, meta.tmp$sample_id),] %>% select(-sample_id)

dist_mat = vegdist(mat.tmp, method = "bray") 

study_readable.tmp = factor(meta.tmp[["study_readable"]])
housing_dog.tmp = factor(meta.tmp[["housing_dog"]])
age.tmp = meta.tmp[["age"]]

study_readable.dispersion = betadisper(dist_mat, study_readable.tmp)
housing_dog.dispersion = betadisper(dist_mat, housing_dog.tmp)
age.dispersion = betadisper(dist_mat, age.tmp)

study_readable.variable = anova(study_readable.dispersion)$`Pr(>F)`[1]
housing_dog.variable = anova(housing_dog.dispersion)$`Pr(>F)`[1]
age.variable = anova(age.dispersion)$`Pr(>F)`[1]

Adonis_dog = adonis2(mat.tmp ~ study_readable + housing_dog + age, data = meta.tmp, by = "margin", permutations = 999, method = "bray")
res_dog = data.frame("variable" = c("study_readable", "housing_dog", "age"),
                     "dispersion.p" = c(study_readable.variable, housing_dog.variable, age.variable),
                     "adonis.r2" = Adonis_dog$R2[1:3],
                     "adonis.p" = Adonis_dog$`Pr(>F)`[1:3],
                     "note" = "only dogs")

write.csv(res_dog,
          file = paste0("rename/SGB_PERMANOVA_onlydogs_multivariate.csv"),
          row.names = F)  

# housing_all
meta.tmp = pcoa_meta_permanova %>% select(all_of(c("sample_id", "study_readable", "housing_all", "age", "species")))
meta.tmp = meta.tmp[complete.cases(meta.tmp),]
common.sps = intersect(rownames(pcoa_mat_rename), meta.tmp$sample_id)

mat.tmp = pcoa_mat_rename[match(common.sps, rownames(pcoa_mat_rename)),]
meta.tmp = meta.tmp[match(common.sps, meta.tmp$sample_id),] %>% select(-sample_id)

dist_mat = vegdist(mat.tmp, method = "bray") 
  
study_readable.tmp = factor(meta.tmp[["study_readable"]])
housing_all.tmp = factor(meta.tmp[["housing_all"]])
age.tmp = meta.tmp[["age"]]
species.tmp = factor(meta.tmp[["species"]])

study_readable.dispersion = betadisper(dist_mat, study_readable.tmp)
housing_all.dispersion = betadisper(dist_mat, housing_all.tmp)
age.dispersion = betadisper(dist_mat, age.tmp)
species.dispersion = betadisper(dist_mat, species.tmp)

study_readable.variable = anova(study_readable.dispersion)$`Pr(>F)`[1]
housing_all.variable = anova(housing_all.dispersion)$`Pr(>F)`[1]
age.variable = anova(age.dispersion)$`Pr(>F)`[1]
species.variable = anova(species.dispersion)$`Pr(>F)`[1]


Adonis_all = adonis2(mat.tmp ~ study_readable + housing_all + age + species, data = meta.tmp, by = "margin", permutations = 999, method = "bray")
    
res_all = data.frame("variable" = c("study_readable", "housing_all", "age", "species"),
                     "dispersion.p" = c(study_readable.variable, housing_all.variable, age.variable, species.variable),
                     "adonis.r2" = Adonis_all$R2[1:4],
                     "adonis.p" = Adonis_all$`Pr(>F)`[1:4],
                     "note" = "all hosts")

write.csv(res_all,
          file = paste0("rename/SGB_PERMANOVA_allhosts_multivariate.csv"),
          row.names = F)  

# for all samples
meta.tmp = pcoa_meta_permanova %>% select(all_of(c("sample_id", "study_readable", "species")))
meta.tmp = meta.tmp[complete.cases(meta.tmp),]
common.sps = intersect(rownames(pcoa_mat_rename), meta.tmp$sample_id)

mat.tmp = pcoa_mat_rename[match(common.sps, rownames(pcoa_mat_rename)),]
meta.tmp = meta.tmp[match(common.sps, meta.tmp$sample_id),] %>% select(-sample_id)

dist_mat = vegdist(mat.tmp, method = "bray") 

study_readable.tmp = factor(meta.tmp[["study_readable"]])
species.tmp = factor(meta.tmp[["species"]])

study_readable.dispersion = betadisper(dist_mat, study_readable.tmp)
species.dispersion = betadisper(dist_mat, species.tmp)

study_readable.variable = anova(study_readable.dispersion)$`Pr(>F)`[1]
species.variable = anova(species.dispersion)$`Pr(>F)`[1]

Adonis_species = adonis2(mat.tmp ~ study_readable + species, data = meta.tmp, by = "margin", permutations = 999, method = "bray")
res_species = data.frame("variable" = c("study_readable", "species"),
                         "dispersion.p" = c(study_readable.variable, species.variable),
                         "adonis.r2" = Adonis_species$R2[1:2],
                         "adonis.p" = Adonis_species$`Pr(>F)`[1:2],
                         "note" = "all hosts")

# only dog
meta.tmp = pcoa_meta_permanova %>% select(all_of(c("sample_id", "housing_dog", "age", "species"))) %>% filter(species == "dog") %>% select(-species)
meta.tmp = meta.tmp[complete.cases(meta.tmp),]
common.sps = intersect(rownames(pcoa_mat_rename), meta.tmp$sample_id)

mat.tmp = pcoa_mat_rename[match(common.sps, rownames(pcoa_mat_rename)),]
meta.tmp = meta.tmp[match(common.sps, meta.tmp$sample_id),] %>% select(-sample_id)

dist_mat = vegdist(mat.tmp, method = "bray") 

housing_dog.tmp = factor(meta.tmp[["housing_dog"]])
age.tmp = meta.tmp[["age"]]

housing_dog.dispersion = betadisper(dist_mat, housing_dog.tmp)
age.dispersion = betadisper(dist_mat, age.tmp)

housing_dog.variable = anova(housing_dog.dispersion)$`Pr(>F)`[1]
age.variable = anova(age.dispersion)$`Pr(>F)`[1]

Adonis_dog2 = adonis2(mat.tmp ~ age + housing_dog, data = meta.tmp, by = "margin", permutations = 999, method = "bray")
res_dog2 = data.frame("variable" = c("age", "housing_dog"),
                      "dispersion.p" = c(age.variable, housing_dog.variable),
                      "adonis.r2" = Adonis_dog$R2[1:2],
                      "adonis.p" = Adonis_dog$`Pr(>F)`[1:2],
                      "note" = "only dogs")

###########################################################################################################
## Heatmap (Figure 3C)
# filter criteria: relab > 10^5 in at least 10% of samples within the host species in at least 2 studies
species.list.filter.abundant = function(df){
  n_samples = nlevels(as.factor(df$sample_id))
  df_species = df %>% 
    filter(relab > 0.00001) %>%
    group_by(species) %>%
    tally() %>%
    ungroup() %>%
    filter(n > 0.1 * n_samples) 
  dftop = df_species$species
}

find.2studies = function(list){
  top_all = unlist(list)
  prevalent_species = top_all[duplicated(top_all)]
  top = unique(prevalent_species)
}

# cat top 10 list
cat = all[which(all$species == "cat"),]

tmp.cat = cat %>%
  select(-subject_id, -study_readable, -species, -sex, -age, -weight, -antibiotics, -medications, -`Hill's/Public`, -`Facility/Household`, -weightpc) %>% 
  arrange(sample_id) %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  select(species, sample_id, relab) %>%
  add_column(rep(cat$study_readable, each = 3682)) 
names(tmp.cat) = c('species', 'sample_id', 'relab', 'study_readable')

# apply for each study
cattop_10 = list()

for (i in unique(tmp.cat$study_readable)) {
  cat_i = tmp.cat[tmp.cat$study_readable == i,]
  cattop_10[[i]] = species.list.filter.abundant(cat_i)
}

cattop.10 = find.2studies(cattop_10)

cat_numeric = cat[, 13:ncol(cat)]
top_cat = cat_numeric[, names(cat_numeric) %in% cattop.10 ]

cattop_ori = names(sort(colSums(top_cat), decreasing = TRUE))[1:10] # find top 10 
#cattop_ori = names(sort(colSums(top_cat), decreasing = TRUE))[1:8] # find top 8

# dog top 10 list
dog = all[which(all$species=="dog"),]

tmp.dog = dog %>%
  select(-subject_id, -study_readable, -species, -sex, -age, -weight, -antibiotics, -medications, -`Hill's/Public`, -`Facility/Household`, -weightpc) %>% 
  arrange(sample_id) %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  select(species, sample_id, relab) %>%
  add_column(rep(dog$study_readable, each = 3682)) 
names(tmp.dog) = c('species', 'sample_id', 'relab', 'study_readable')

# apply for each study
dogtop_10 = list()

for (i in unique(tmp.dog$study_readable)) {
  dog_i = tmp.dog[tmp.dog$study_readable == i,]
  dogtop_10[[i]] = species.list.filter.abundant(dog_i)
}

dogtop.10 = find.2studies(dogtop_10)

dog_numeric = dog[, 13:ncol(dog)]
top_dog = dog_numeric[, names(dog_numeric) %in% dogtop.10] # dogtop.25 is the same as dogtop.10

dogtop_ori = names(sort(colSums(top_dog), decreasing = TRUE))[1:10] # find top 10 
#dogtop_ori = names(sort(colSums(top_dog), decreasing = TRUE))[1:8] # find top 8

# human top 10 list
human = all[which(all$species=="human"),]

tmp.human = human %>%
  select(-subject_id, -study_readable, -species, -sex, -age, -weight, -antibiotics, -medications, -`Hill's/Public`, -`Facility/Household`, -weightpc) %>% 
  arrange(sample_id) %>%
  pivot_longer(!sample_id, names_to = "species", values_to = "relab") %>%
  select(species, sample_id, relab) %>%
  add_column(rep(human$study_readable, each = 3682)) 
names(tmp.human) = c('species', 'sample_id', 'relab', 'study_readable')

# apply for each study
humantop_10 = list()

for (i in unique(tmp.human$study_readable)) {
  human_i = tmp.human[tmp.human$study_readable == i,]
  humantop_10[[i]] = species.list.filter.abundant(human_i)
}

humantop.10 = find.2studies(humantop_10)

human_numeric = human[, 13:ncol(human)]
top_human = human_numeric[,names(human_numeric) %in% humantop.10] 

humantop_ori = names(sort(colSums(top_human), decreasing = TRUE))[1:10] # find top 10 
#humantop_ori = names(sort(colSums(top_human), decreasing = TRUE))[1:8] # find top 8

# final lists across host
topspecies = c(cattop_ori, dogtop_ori, humantop_ori) %>% unique() # 27 (top 10) or 21 (top 8)
# save(topspecies, file = paste0("Taxonomic profiling/metaphlanv4/top27_list_heatmap.Rdata"))

# sort out columns based on final lists 
df_numeric = all[, 13:ncol(all)]
top = df_numeric[, names(df_numeric) %in% topspecies]

# transpose for heatmap
top_t = t(top) 
rownames(top_t) = gsub("_", " ", rownames(top_t), fixed=TRUE)
rownames(top_t) = gsub("SGB ", "SGB", rownames(top_t), fixed=TRUE)
min_relab = min(top_t[top_t>0]) / 10 # set 0 to be 1/10 of the minimum relative abundance = 2.000001e-08
logtop_t = top_t # save for log10 transformation
logtop_t[logtop_t == 0] = min_relab # matrix

metadata = pcoa_meta %>%
  select(`Study ID` = "study_readable", `Host Species` = "species", "Facility/Household", "Hill's/Public")
metadata[is.na(metadata)] = "NA" # metadata (also column annotation)

uSGB.meta = data.frame("old_bug_id" = rownames(top_t)) %>%
  mutate(SGB = case_when(grepl("SGB14350", `old_bug_id`) ~ "uSGB",
                         grepl("SGB53888", `old_bug_id`) ~ "uSGB",
                         grepl("SGB15313", `old_bug_id`) ~ "uSGB",
                         grepl("SGB1818", `old_bug_id`) ~ "uSGB",
                         grepl("SGB1481", `old_bug_id`) ~ "uSGB",
                         TRUE ~ "kSGB")) %>%
  column_to_rownames('old_bug_id') # row annotation

# plot heatmap
ann_colors = list(
  `Study ID` = c(`Madagascar` = "#0F3B50", `HMP1-II` = "#56B4E0", `Cross-sectional Study2` = "#B7AABA", `Cross-sectional Study3` = "#D44AB4", 
                 `Deusch et al. (2014)` = "#639CCF", `Coelho et al.` = "#4D643A", `Yarlagadda et al.` = "#B5E6D4",
                 `Allaway et al.` = "#31378C", `Deusch et al. (2015)` = "#6D81AE", `DietInt Study4` = "#798C7E", `DietInt Study5` = "#B05C63",
                 `Cross-sectional Study1` = "#FFFFB3", `DietInt Study1` = "#d69d68", `DietInt Study3` = "#7eda86", `DietInt Study2` = "#fbc086",
                 `Young et al.` = "#74D0D1", `Ateba et al.` = "#E0A8F5", `Liu et al.` = "#D2BED7", `Tanprasertsuk et al.` = "#C7A679",
                 `Ma et al.` = "#612141", `Xu et al.` = "#AE4109", `Alessandri et al.` = "#B435EE", `Wang et al.` = "#E9FDC3",
                 `Maldonado-Contreras et al.` = "#C7F4FC"),
  `Host Species` = c(dog = "#E69F00", cat = "indianred", human = "#56B4E0"),
  `Facility/Household` = c(Facility = "#F2C143", Household = "#3A8A62", Stray = "#A4514F", "NA" = "gray"),
  `Hill's/Public` = c(`Hill's` = "#A07855" , `Public` = "#D4B996", "NA" = "gray"),
  SGB = c(uSGB = "#FF0000", kSGB = "#000000"))

cairo_pdf("rename/figure/heatmap_filtered_top27_2Xhuman_v4_naturalcluster.pdf", width = 20, height = 9, fallback_resolution = 300)
#cairo_pdf("rename/figure/heatmap_filtered_top21_2Xhuman_v4_naturalcluster.pdf", width = 20, height = 9, fallback_resolution = 300)

pheatmap(mat = log10(logtop_t), 
         show_colnames = FALSE,
         color = inferno(10),
         border_color = NA, 
         annotation_col = metadata,
         annotation_row = uSGB.meta,
         annotation_colors = ann_colors,
         cellheight = 18, 
         cellwidth = 0.27,
         fontsize = 12,
         legend = TRUE,
         name = "log10 Scale Relative Abundance",
         annotation_names_row = FALSE)

dev.off()

#######################################################################################################################################
## Boxplot: microbial differences between dog housing conditions (Supplemental Figure 4) 
# format dataframe
df = data.frame(pcoa_meta_meta2, maaslin_mat) %>% 
  filter(species == "dog") %>%
  select(c(housing, Firmicutes_SGB_105987, Bifidobacterium_pseudolongum_SGB_17279, 
           Collinsella_SGB14744_SGB_14744,
           # <- Coriobacteriaceae_SGB14744,
           #Firmicutes_SGB_4052,
           # <- Moraxellaceae_SGB4052,
           Lactobacillaceae_SGB_7083,
           Tyzzerella_nexilis_SGB_4588, Bilophila_wadsworthia_SGB_15452, Succinivibrionaceae_SGB_3675,
           Slackia_SGB14780_SGB_14780, Turicibacter_sanguinis_SGB_6846, Helicobacter_bilis_SGB_19390))

# load maaslin result
dog_housing = read.csv('rename/double_re/maaslin_output_housing_dog/all_results.tsv', head = TRUE, sep = "\t")

# color
housing.colors= c(Facility = "#F2C143", Household = "#3A8A62", Stray = "#A4514F")

# boxplot function
plot_box = function(df, bug){
  df_bug = cbind(df['housing'], relab = df[bug]) %>% add_column(bug = names(df[bug]))
  names(df_bug) = c("Housing", "relab", "bug")
  
  maaslin_tmp = dog_housing[grepl(bug, dog_housing$feature),]
  coef1 = maaslin_tmp %>% filter(value == "Household") %>% pull(coef) %>% round(digits = 3)
  qval1 = maaslin_tmp %>% filter(value == "Household") %>% pull(qval) %>% round(digits = 3)
  
  coef2 = maaslin_tmp %>% filter(value == "Stray") %>% pull(coef) %>% round(digits = 3)
  qval2 = maaslin_tmp %>% filter(value == "Stray") %>% pull(qval) %>% round(digits = 3)
  
  ggplot(df_bug, aes(x = Housing, y = relab, fill = Housing)) +
    geom_boxplot(outlier.alpha = 0.0, na.rm = TRUE, alpha = .5, show.legend = FALSE) +
    geom_point(alpha = 0.5, size = 4, shape = 21, stroke = 0.15, color = 'black', position = ggplot2::position_jitterdodge(), show.legend = FALSE) +
    annotate("text", label = paste0("Household vs. Facility:", "\ncoef = ", coef1, "; FDR p-value = ", qval1,
                                    "\nStray vs. Facility:", "\ncoef = ", coef2, "; FDR p-value = ", qval2), 
             size = 5, hjust = 0, vjust = 1, fontface = "italic", x = "Facility", y = max(df_bug$relab)) +
    labs(x = NULL, y = paste0(bug)) +
    theme_classic() +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 18)) +
    scale_fill_manual(values = housing.colors)
}

p1 = plot_box(df, bug = 'Firmicutes_SGB_105987')
p2 = plot_box(df, bug = 'Bifidobacterium_pseudolongum_SGB_17279')
p3 = plot_box(df, bug = 'Collinsella_SGB14744_SGB_14744')
p4 = plot_box(df, bug = 'Lactobacillaceae_SGB_7083')
p5 = plot_box(df, bug = 'Tyzzerella_nexilis_SGB_4588')
p6 = plot_box(df, bug = 'Bilophila_wadsworthia_SGB_15452')
p7 = plot_box(df, bug = 'Succinivibrionaceae_SGB_3675')
p8 = plot_box(df, bug = 'Turicibacter_sanguinis_SGB_6846')
p9 = plot_box(df, bug = 'Slackia_SGB14780_SGB_14780')
p10 = plot_box(df, bug = 'Helicobacter_bilis_SGB_19390')

cairo_pdf("rename/figure/box_housing.pdf", width = 22.5, height = 8, fallback_resolution = 300)

main = grid.arrange(p6, p1, p3, p2, p7,
                    p9, p5, p10, p8, p4, ncol = 5)

dev.off()

#######################################################################################################################################
## Density plot: distribution patterns of microbial species across host species (Supplemental Figure 5)
# format dataframe and select species of interest
metadata_in = all[c(1,4)]
names(metadata_in) = c("sample_id", "Host Species")

df_numeric = all[,13:ncol(all)]
colnames(df_numeric) = gsub("_", " ", colnames(df_numeric), fixed=TRUE) # matrix

topspecies = gsub("_", " ", topspecies, fixed=TRUE) # top27 identified from heatmap section
interest = c("Phocaeicola vulgatus SGB 1814", "Collinsella intestinalis SGB 14741", 
             "Oscillospiraceae SGB 15313", "Succinivibrionaceae SGB 3677",
             "Bifidobacterium pseudolongum SGB 17279", "Megasphaera elsdenii SGB 5862", 
             "Atopobiaceae SGB 14350",
             # <- "Actinobacteria SGB 14350", 
             "Acidaminococcaceae SGB 5749",
             # <- "Firmicutes SGB 5749", 
             "Helicobacter canis SGB 21969", 
             "Prevotella copri SGB 1644",
             # <- "Prevotella copri clade C SGB 1644",
             "Firmicutes SGB 6260"
             # <- "Enterococcaceae SGB 6260"
) # bugs among the interest

all_species = unique(c(topspecies, interest)) # n = 32

top = df_numeric[, names(df_numeric) %in% all_species]

df = cbind(metadata_in, top)

human_df = df %>% filter(`Host Species` == "human")
for (i in 3:ncol(human_df)) {
  if (sum(human_df[[i]] > 0.00001) >= 3) {
    human_df[[i]][human_df[[i]] < 0.00001] <- NA
  } else {
    human_df[[i]] <- NA
  }
}

cat_df = df %>% filter(`Host Species` == "cat")
for (i in 3:ncol(cat_df)) {
  if (sum(cat_df[[i]] > 0.00001) >= 3) {
    cat_df[[i]][cat_df[[i]] < 0.00001] <- NA
  } else {
    cat_df[[i]] <- NA
  }
}

dog_df = df %>% filter(`Host Species` == "dog")
for (i in 3:ncol(dog_df)) {
  if (sum(dog_df[[i]] > 0.00001) >= 3) {
    dog_df[[i]][dog_df[[i]] < 0.00001] <- NA
  } else {
    dog_df[[i]] <- NA
  }
}

df = rbind(human_df, cat_df, dog_df)

df_set1 = df %>% 
  select(c(1:2, 3:10)) %>%
  pivot_longer(!c(1:2), names_to = "species", values_to = "relab")

df_set2 = df %>% 
  select(c(1:2, 11:18)) %>%
  pivot_longer(!c(1:2), names_to = "species", values_to = "relab")

df_set3 = df %>% 
  select(c(1:2, 19:26)) %>%
  pivot_longer(!c(1:2), names_to = "species", values_to = "relab")

df_set4 = df %>% 
  select(c(1:2, 27:34)) %>%
  pivot_longer(!c(1:2), names_to = "species", values_to = "relab")

# log10scale set1
species.colors = c(cat = "indianred", dog = "#E69F00", human = "#56B4E0")

df_set1$relab = log10(df_set1$relab)
df_set1$species = factor(df_set1$species, levels = unique(df_set1$species))

density_plot1 = ggplot(df_set1, aes(x = relab, fill = `Host Species`)) + 
  geom_density(aes(y=after_stat(density)), alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = species), size = 10, x = -5, y = 1.5, hjust = 0, vjust = 0.5, check_overlap = TRUE, fontface = "italic") +
  labs(title=NULL, x="log10 scale Relative abundance", y=NULL, fill = "Host Species") +
  theme_bw() +
  theme(axis.text=element_text(size= 30),
        legend.title=element_text(size = 30, face = "bold"), 
        legend.text= element_text(size = 30),
        axis.title = element_text(size = 35, face = "bold", color = "white"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = species.colors) +
  xlim(-5, 1) + 
  ylim(0, 1.6) + 
  facet_wrap(~ species, scale = "free_y", ncol= 1, strip.position="right")

# missing (proportion of zero) barplot
{
  missing_df_dog =  df %>% 
    select(c(1:2, 3:10)) %>% 
    filter(`Host Species` == "dog")
  missing_df_dog[-c(1:2)] = ifelse(is.na(missing_df_dog[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_dog[3:10])
  proportion =  (colSums(missing_df_dog[3:10])/nrow(missing_df_dog[3:10]))*100
  host = c("dog")
  
  dog = data.frame(species, proportion, host)
  
  missing_df_cat = df %>% 
    select(c(1:2, 3:10)) %>% 
    filter(`Host Species` == "cat")
  missing_df_cat[-c(1:2)] = ifelse(is.na(missing_df_cat[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_cat[3:10])
  proportion =  (colSums(missing_df_cat[3:10])/nrow(missing_df_cat[3:10]))*100
  host = c("cat")
  
  cat = data.frame(species, proportion, host)
  
  missing_df_human = df %>% 
    select(c(1:2, 3:10)) %>% 
    filter(`Host Species` == "human")
  missing_df_human[-c(1:2)] = ifelse(is.na(missing_df_human[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_human[3:10])
  proportion =  (colSums(missing_df_human[3:10])/nrow(missing_df_human[3:10]))*100
  host = c("human")
  
  human = data.frame(species, proportion, host)
  
  missing_df_p = rbind(cat, dog, human)
  missing_df_p$species = factor(missing_df_p$species, levels = unique(df_set1$species))
  
  missing_barplot1 = ggplot(missing_df_p, aes(x = host, y = proportion, fill = host)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 30, color = "white"),
      axis.text.y = element_text(size = 30),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 35, face = "bold", color = "white"),
      axis.title.x = element_text(size = 35, face = "bold", color = "white"),
      strip.text = element_blank(),
      strip.background = element_blank()) +
    scale_fill_manual(values = species.colors) +
    ylim(0, 100) +
    labs(y = "Proportion (%) of missing value", fill = NULL, x = "Host Species", title = NULL) +
    facet_wrap(~ species, ncol= 1, strip.position="right") +
    guides(fill ="none")
  }

plot1 = grid.arrange(missing_barplot1, density_plot1, ncol = 2, widths = c(0.5, 1.5))

# log10scale
df_set2$relab = log10(df_set2$relab)
df_set2$species = factor(df_set2$species, levels = unique(df_set2$species))

density_plot2 = ggplot(df_set2, aes(x = relab, fill = `Host Species`)) + 
  geom_density(aes(y=after_stat(density)), alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = species), size = 10, x = -5, y = 1.5, hjust = 0, vjust = 0.5, check_overlap = TRUE, fontface = "italic") +
  labs(title=NULL, x="log10 scale Relative abundance", y=NULL, fill = "Host Species") +
  theme_bw() +
  theme(axis.text=element_text(size= 30),
        legend.title=element_text(size = 30, face = "bold"), 
        legend.text= element_text(size = 30),
        axis.title = element_text(size = 35, face = "bold", color = "white"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = species.colors) +
  xlim(-5, 1) + 
  ylim(0, 1.6) + 
  facet_wrap(~ species, scale = "free_y", ncol= 1, strip.position="right")

# missing (proportion of zero) barplot
{
  missing_df_dog = df %>% 
    select(c(1:2, 11:18)) %>% 
    filter(`Host Species` == "dog")
  missing_df_dog[-c(1:2)] = ifelse(is.na(missing_df_dog[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_dog[3:10])
  proportion =  (colSums(missing_df_dog[3:10])/nrow(missing_df_dog[3:10]))*100
  host = c("dog")
  
  dog = data.frame(species, proportion, host)
  
  missing_df_cat = df %>% 
    select(c(1:2, 11:18)) %>% 
    filter(`Host Species` == "cat")
  missing_df_cat[-c(1:2)] = ifelse(is.na(missing_df_cat[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_cat[3:10])
  proportion =  (colSums(missing_df_cat[3:10])/nrow(missing_df_cat[3:10]))*100
  host = c("cat")
  
  cat = data.frame(species, proportion, host)
  
  missing_df_human = df %>% 
    select(c(1:2, 11:18)) %>% 
    filter(`Host Species` == "human")
  missing_df_human[-c(1:2)] = ifelse(is.na(missing_df_human[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_human[3:10])
  proportion =  (colSums(missing_df_human[3:10])/nrow(missing_df_human[3:10]))*100
  host = c("human")
  
  human = data.frame(species, proportion, host)
  
  missing_df_p = rbind(cat, dog, human)
  missing_df_p$species = factor(missing_df_p$species, levels = unique(df_set2$species))
  
  missing_barplot2 = ggplot(missing_df_p, aes(x = host, y = proportion, fill = host)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 30, color = "white"),
      axis.text.y = element_text(size = 30),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 35, face = "bold", color = "white"),
      axis.title.x = element_text(size = 35, face = "bold", color = "white"),
      strip.text = element_blank(),
      strip.background = element_blank()) +
    scale_fill_manual(values = species.colors) +
    ylim(0, 100) +
    labs(y = "Proportion (%) of missing value", fill = NULL, x = "Host Species", title = NULL) +
    facet_wrap(~ species, ncol= 1, strip.position="right") +
    guides(fill ="none")
  }

plot2 = grid.arrange(missing_barplot2, density_plot2, ncol = 2, widths = c(0.5, 1.5))

# log10scale
df_set3$relab = log10(df_set3$relab)
df_set3$species = factor(df_set3$species, levels = unique(df_set3$species))

density_plot3 = ggplot(df_set3, aes(x = relab, fill = `Host Species`)) + 
  geom_density(aes(y=after_stat(density)), alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = species), size = 10, x = -5, y = 1.5, hjust = 0, vjust = 0.5, check_overlap = TRUE, fontface = "italic") +
  labs(title=NULL, x="log10 scale Relative abundance", y=NULL, fill = "Host Species") +
  theme_bw() +
  theme(axis.text=element_text(size= 30),
        legend.title=element_text(size = 30, face = "bold"), 
        legend.text= element_text(size = 30),
        axis.title = element_text(size = 35, face = "bold", color = "white"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = species.colors) +
  xlim(-5, 1) + 
  ylim(0, 1.6) + 
  facet_wrap(~ species, scale = "free_y", ncol= 1, strip.position="right")

# missing (proportion of zero) barplot
{
  missing_df_dog = df %>% 
    select(c(1:2, 19:26)) %>% 
    filter(`Host Species` == "dog")
  missing_df_dog[-c(1:2)] = ifelse(is.na(missing_df_dog[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_dog[3:10])
  proportion =  (colSums(missing_df_dog[3:10])/nrow(missing_df_dog[3:10]))*100
  host = c("dog")
  
  dog = data.frame(species, proportion, host)
  
  missing_df_cat = df %>% 
    select(c(1:2, 19:26)) %>% 
    filter(`Host Species` == "cat")
  missing_df_cat[-c(1:2)] = ifelse(is.na(missing_df_cat[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_cat[3:10])
  proportion =  (colSums(missing_df_cat[3:10])/nrow(missing_df_cat[3:10]))*100
  host = c("cat")
  
  cat = data.frame(species, proportion, host)
  
  missing_df_human = df %>% 
    select(c(1:2, 19:26)) %>% 
    filter(`Host Species` == "human")
  missing_df_human[-c(1:2)] = ifelse(is.na(missing_df_human[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_human[3:10])
  proportion =  (colSums(missing_df_human[3:10])/nrow(missing_df_human[3:10]))*100
  host = c("human")
  
  human = data.frame(species, proportion, host)
  
  missing_df_p = rbind(cat, dog, human)
  missing_df_p$species = factor(missing_df_p$species, levels = unique(df_set3$species))
  
  missing_barplot3 = ggplot(missing_df_p, aes(x = host, y = proportion, fill = host)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 30, color = "white"),
      axis.text.y = element_text(size = 30),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 35, face = "bold", color = "white"),
      axis.title.x = element_text(size = 35, face = "bold", color = "white"),
      strip.text = element_blank(),
      strip.background = element_blank()) +
    scale_fill_manual(values = species.colors) +
    ylim(0, 100) +
    labs(y = "Proportion (%) of missing value", fill = NULL, x = "Host Species", title = NULL) +
    facet_wrap(~ species, ncol= 1, strip.position="right") +
    guides(fill ="none")
  }

plot3 = grid.arrange(missing_barplot3, density_plot3, ncol = 2, widths = c(0.5, 1.5))

# log10scale
df_set4$relab = log10(df_set4$relab)
df_set4$species = factor(df_set4$species, levels = unique(df_set4$species))

density_plot4 = ggplot(df_set4, aes(x = relab, fill = `Host Species`)) + 
  geom_density(aes(y=after_stat(density)), alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = species), size = 10, x = -5, y = 1.5, hjust = 0, vjust = 0.5, check_overlap = TRUE, fontface = "italic") +
  labs(title=NULL, x="log10 scale Relative abundance", y=NULL, fill = "Host Species") +
  theme_bw() +
  theme(axis.text=element_text(size= 30),
        legend.title=element_text(size = 30, face = "bold"), 
        legend.text= element_text(size = 30),
        axis.title = element_text(size = 35, face = "bold", color = "white"),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  scale_fill_manual(values = species.colors) +
  xlim(-5, 1) + 
  ylim(0, 1.6) + 
  facet_wrap(~ species, scale = "free_y", ncol= 1, strip.position="right")

# missing (proportion of zero) barplot
{
  missing_df_dog = df %>% 
    select(c(1:2, 27:34)) %>% 
    filter(`Host Species` == "dog")
  missing_df_dog[-c(1:2)] = ifelse(is.na(missing_df_dog[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_dog[3:10])
  proportion =  (colSums(missing_df_dog[3:10])/nrow(missing_df_dog[3:10]))*100
  host = c("dog")
  
  dog = data.frame(species, proportion, host)
  
  missing_df_cat = df %>% 
    select(c(1:2, 27:34)) %>% 
    filter(`Host Species` == "cat")
  missing_df_cat[-c(1:2)] = ifelse(is.na(missing_df_cat[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_cat[3:10])
  proportion =  (colSums(missing_df_cat[3:10])/nrow(missing_df_cat[3:10]))*100
  host = c("cat")
  
  cat = data.frame(species, proportion, host)
  
  missing_df_human = df %>% 
    select(c(1:2, 27:34)) %>% 
    filter(`Host Species` == "human")
  missing_df_human[-c(1:2)] = ifelse(is.na(missing_df_human[-c(1:2)]), 1, 0)
  
  species = colnames(missing_df_human[3:10])
  proportion =  (colSums(missing_df_human[3:10])/nrow(missing_df_human[3:10]))*100
  host = c("human")
  
  human = data.frame(species, proportion, host)
  
  missing_df_p = rbind(cat, dog, human)
  missing_df_p$species = factor(missing_df_p$species, levels = unique(df_set4$species))
  
  missing_barplot4 = ggplot(missing_df_p, aes(x = host, y = proportion, fill = host)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 30, color = "white"),
      axis.text.y = element_text(size = 30),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 35, face = "bold", color = "white"),
      axis.title.x = element_text(size = 35, face = "bold", color = "white"),
      strip.text = element_blank(),
      strip.background = element_blank()) +
    scale_fill_manual(values = species.colors) +
    ylim(0, 100) +
    labs(y = "Proportion (%) of missing value", fill = NULL, x = "Host Species", title = NULL) +
    facet_wrap(~ species, ncol= 1, strip.position="right") +
    guides(fill ="none")
  }

plot4 = grid.arrange(missing_barplot4, density_plot4, ncol = 2, widths = c(0.5, 1.5))

legend = get_legend(density_plot1) # set any density_plot1/2/3/4 show.legend = T to get the legend
bottom = textGrob("log10 Scale Relative Abundance", gp = gpar(fontsize = 60, fontface = "bold"))
yleft = textGrob("Proportion (%) of Missing Values", rot = 90, gp = gpar(fontsize = 60, fontface = "bold"))

plot = grid.arrange(plot1, plot2, plot3, plot4, legend, ncol = 5,
                    bottom = bottom,
                    left = yleft, widths = c(1, 1, 1, 1, 0.3))

ggsave(plot, file = paste0("rename/figure/density_all_supplement.pdf"), height = 30, width = 43, dpi = 300, limitsize = FALSE)

#######################################################################################################################################
#################################################################end###################################################################
#######################################################################################################################################