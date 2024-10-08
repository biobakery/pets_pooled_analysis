### Pet's pooled study: metadata ###
### Last Edit: 08/19/24 ###

# rm(list = ls())

#######################################################################################################################################
## packages
library(tidyverse)
library(dplyr)

## input
animal_public = read.csv('metadata/WGX_studylist.csv') 
# total 14 public studies with WGS samples only
# Allaway et al (PRJPRJEB34360) has total 64 samples deposited in ENA. Among all, 62 are WGS samples (the rest 2 are 16S).
   # 14 (*2) of them are duplicates of the same sample. We combine the dulipcates and rename the sample by the first of the two duplicates.
   # The final sample number of Allaway et al is 48.
   # Duplicates are: "ERR4183453|ERR4183458|ERR4183469|ERR4183472|ERR4183476|ERR4183479|ERR4183481|ERR4183483|ERR4183485|ERR4183487|ERR4183489|ERR4183492|ERR4183494|ERR4183496"
   # In addition, the study does not provde subject_id
# Coelho et al (PRJEB20308) has total 3096 samples deposited in ENA. All samples are WGS.
   # Each sample has 24 duplicates. We combine the duplicates and rename the sample by their unique accession (sample_accession) provided through ENA. 
   # The final sample number of Coelho et al is 129.
animal_newstudy1 = animal_public %>% filter(PMID == 29669589) %>% select(sample_accession_MGX)
animal_newstudy2 = animal_public %>% filter(PMID == 35506233) %>% select(sample_accession_MGX)

animal_newstudy_together = rbind(animal_newstudy1, animal_newstudy2)
animal_newstudy_together$sample_id_metaphlan = animal_newstudy_together$sample_accession_MGX
names(animal_newstudy_together) = c("sample_id_metaphlan", "sample_id")

# check actual sample_id used in metaphlanv4 
animal = read_csv('Taxonomic profiling/metaphlanv4/merged_taxonomic_SGB_animal_profile_normalized.csv') %>%
  names() %>% .[-1] %>% as_tibble() %>% setNames("sample_id_metaphlan") %>%
  mutate(sample_id = sample_id_metaphlan,
         sample_id = gsub("re", "", fixed = TRUE,
                          gsub("_1", "", fixed = TRUE, sample_id)))
animal = rbind(animal, animal_newstudy_together)

animal_id = read_csv('Taxonomic profiling/old public studies + Hills studies/Hills_samples_individual_animal_info.csv')
animal_all = read_csv('Taxonomic profiling/old public studies + Hills studies/formatted_metadata.csv')
# reasons of using formatted_metadata: 1) it contains only samples are actually used in taxonomic profiles (ie., subsets of Hill's samples)
#                                      2) it has the matched sample_id with the taxonomic profile (ie., some samples have "re")
# see how to format animal_all in the "skip" section
animal_Hills = read_csv("Taxonomic profiling/old public studies + Hills studies/new_Hills_Metadata.csv")
human_HMP = read_tsv('Taxonomic profiling/metaphlanv4/hmp1-ii_public_metadata_tidy.tsv')
human_hmp = read_csv('Taxonomic profiling/metaphlanv4/merged_taxonomic_SGB_human_profile_normalized.csv') %>%
  names() %>% .[-1] %>% as_tibble() %>% setNames("sample_id")
human_MAD = read_delim('Taxonomic profiling/metaphlanv4/filereport_read_run_PRJNA485056_tsv.txt')

## transform data
# animal_public
animal_public = animal_public %>%
  mutate(sex = case_when(sex == "male" ~ "Male",
                         sex == "female" ~ "Female"),
         `Hill's/Public` = c("Public"),
         `Facility/Household` = case_when(PMID == "32303546" ~ "Facility",
                                         PMID == "33676490" ~ "Facility",
                                         PMID == "32730134" ~ "Household",
                                         PMID == "25010839" ~ "Facility",
                                         PMID == "26659594" ~ "Facility",
                                         PMID == "27876765" ~ "Facility",
                                         PMID == "30680877" ~ "Household",
                                         PMID == "31021333" ~ "Facility",
                                         PMID == "31472697" ~ "Household",
                                         PMID == "33239463" ~ "Facility",
                                         PMID == "35467389" ~ "Facility",
                                         PMID == "34249503" ~ "Household",
                                         PMID == 29669589 ~ "Facility",
                                         PMID == 35506233 ~ case_when(batch == "India Dog" ~ case_when(database == "Street" ~ "Stray",
                                                                                                       database == "Shelter" ~ "Facility"),
                                                                      batch == "LD" ~ "Household",
                                                                      batch == "Exp1" ~ "Facility",
                                                                      batch == "Exp2" ~ "Facility",
                                                                      batch == "SA Dog" ~ "Household",
                                                                      TRUE ~ NA_character_)),
         antibiotics = case_when(antibiotics == "n" ~ "No",
                                 antibiotics == "y" ~ "Yes"),
         non_abx_med = case_when(non_abx_med == "n" ~ "No",
                                 non_abx_med == "y" ~ "Yes"),
         study_readable = case_when(PMID == 35467389 ~ "Ma et al.",
                                    PMID == 34249503 ~ "Tanprasertsuk et al.",
                                    PMID == 33676490 ~ "Liu et al.",
                                    PMID == 32730134 ~ "Maldonado-Contreras et al.",
                                    PMID == 25010839 ~ "Deusch et al. (2014)",
                                    PMID == 26659594 ~ "Deusch et al. (2015)",
                                    PMID == 27876765 ~ "Young et al.",
                                    PMID == 30680877 ~ "Alessandri et al.",
                                    PMID == 31021333 ~ "Xu et al.",
                                    PMID == 31472697 ~ "Wang et al.",
                                    PMID == 32303546 ~ "Allaway et al.",
                                    PMID == 33239463 ~ "Ateba et al.",
                                    PMID == 29669589 ~ "Coelho et al.",
                                    PMID == 35506233 ~ "Yarlagadda et al.")) %>%
  mutate(subject_id = case_when(PMID == 32303546 ~ NA_character_, # name from ENA report is for sample_id but not subject_id
                                TRUE ~ subject_accession),
         subject_id = ifelse(subject_id == "Fat louis", "Fat Louis", subject_id)) %>%
  select(c(sample_id = "sample_accession_MGX", "subject_id", 
           "study_readable", "species", "sex", 
           age = "age_at_sampling", "weight", 
          "antibiotics", "medications" = "non_abx_med",
          "Hill's/Public", "Facility/Household")) %>%
  filter(sample_id != "SRR14842322")
# %>% select(c(sample = "sample_accession_MGX", "Gender" = "gender", "Age @ collection" = "age_at_sampling", "Body weight (kg) @ collection" = "weight", "Antibiotics (>=4 weeks) (1 = Yes; 0= No)" = "antibiotics", "ANY Meds (>=4 weeks) (1 = yes; 0 = No)" = "non_abx_med", "Internal (I)/Household (E)" = "Internal/Household"))

# animal_Hills
animal_Hills =  animal_Hills %>%
  mutate(sample_id = 
           gsub("_R1_001.fastq.gz", "",
                gsub("_R2_001.fastq.gz", "", 
                     gsub("_R1_001b.fastq.gz", "",
                          gsub("_R2_001b.fastq.gz", "", 
                               gsub("_R2.fastq.gz", "",
                                    gsub("_R1.fastq.gz", "", `FASTQ Sequence Filename`)))))),
         sample_id = ifelse(!(Directory %in% c("Brent_FibersStudy", "Study_109490")), 
                            gsub("_L00[1-9]", "", sample_id), sample_id),
         study_readable = case_when(Directory == "107357_107351_107425" ~ "DietInt Study1",
                                    Directory == "Study_109490" ~ "DietInt Study2",
                                    Directory == "Derm_Keto_109432" ~ "DietInt Study3",
                                    Directory == "Brent_FibersStudy" ~ "DietInt Study4",
                                    Directory == "TOBI_Study" ~ "DietInt Study5",
                                    Directory == "SP_LifeStage" ~ "Cross-sectional Study1",
                                    Directory == "1000DogStudy" ~ "Cross-sectional Study2",
                                    Directory == "rawdata" ~ "Cross-sectional Study3",
                                    TRUE ~ as.character(Directory)),
         `Antibiotics (>=4 weeks) (1 = Yes; 0= No)` = case_when(`Antibiotics (>=4 weeks) (1 = Yes; 0= No)` == "0" ~ "No",
                                                                `Antibiotics (>=4 weeks) (1 = Yes; 0= No)` == "1" ~ "Yes"),
         `ANY Meds (>=4 weeks) (1 = yes; 0 = No)` = case_when(`ANY Meds (>=4 weeks) (1 = yes; 0 = No)` == "0" ~ "No",
                                                              `ANY Meds (>=4 weeks) (1 = yes; 0 = No)` == "1" ~ "Yes"),
         species = case_when(Species == "Cat" ~ "cat",
                             Species == "Dog" ~ "dog"),
         `Hill's/Public` = c("Hill's"),
         `Facility/Household` = ifelse(`Internal (I)/External (E)` == "I", "Facility", "Household")) %>%
  separate(Gender, into = c("gender", "neutered"), sep= " ") %>%
  full_join(animal_id, by = c("sample_id" = "Sample ID")) %>%
  select(c(sample_id = "sample_id", subject_id = "animal_id", "study_readable", 
           "species", "sex" = gender, 
           age = "Age @ collection", weight = "Body weight (kg) @ collection", 
           antibiotics = "Antibiotics (>=4 weeks) (1 = Yes; 0= No)", medications = "ANY Meds (>=4 weeks) (1 = yes; 0 = No)",
           "Hill's/Public", "Facility/Household")) %>%
  unique() 

#################################################### skip #################################################### 
# animal_all
names(animal_all)[names(animal_all) == '...1'] <- 'sample_id'
animal_all$sample_id = gsub("reSRR8861304.1" , "reSRR8861304", animal_all$sample_id)
animal_all = animal_all %>%
  mutate(sample = gsub("reSRR8861304.1" , "reSRR8861304", sample),
         sample = gsub("re", "", sample),
         sample = gsub("_1", "", sample),
         `Antibiotics (>=4 weeks) (1 = Yes; 0= No)` = case_when(`Antibiotics (>=4 weeks) (1 = Yes; 0= No)` == "n" ~ "No",
                                                                `Antibiotics (>=4 weeks) (1 = Yes; 0= No)` == "y" ~ "Yes"),
         `ANY Meds (>=4 weeks) (1 = yes; 0 = No)` = case_when(`ANY Meds (>=4 weeks) (1 = yes; 0 = No)` == "n" ~ "No",
                                                              `ANY Meds (>=4 weeks) (1 = yes; 0 = No)` == "y" ~ "Yes"),
         `Internal (I)/External (E)` = ifelse(`Internal (I)/External (E)` == "I", "Internal", "External"),
         study_readable = case_when(study_readable == "107357_107351_107425" ~ "DietInt Study1",
                                    study_readable == "Study_109490" ~ "DietInt Study2",
                                    study_readable == "Derm_Keto_109432" ~ "DietInt Study3",
                                    study_readable == "Brent_FibersStudy" ~ "DietInt Study4",
                                    study_readable == "TOBI_Study" ~ "DietInt Study5",
                                    study_readable == "SP_LifeStage" ~ "Cross-sectional Study1",
                                    study_readable == "1000DogStudy" ~ "Cross-sectional Study2",
                                    study_readable == "rawdata" ~ "Cross-sectional Study3",
                                    TRUE ~ as.character(study_readable))) %>%
  separate(Gender, into = c("Gender", "neutered"), sep= " ") %>%
  select(-neutered)
         
animal_metadata = left_join(animal_all, animal_public, by = c('sample')) %>%
  mutate(`Age @ collection` = case_when(is.na(`Age @ collection.x`) ~ `Age @ collection.y`, TRUE ~ `Age @ collection.x`),
         Gender = case_when(is.na(Gender.x) ~ Gender.y, TRUE ~ Gender.x),
         `Body weight (kg) @ collection` = case_when(is.na(`Body weight (kg) @ collection.x`) ~ `Body weight (kg) @ collection.y`, TRUE ~ `Body weight (kg) @ collection.x`),
         `Antibiotics (>=4 weeks) (1 = Yes; 0= No)` = case_when(is.na(`Antibiotics (>=4 weeks) (1 = Yes; 0= No).x`) ~ `Antibiotics (>=4 weeks) (1 = Yes; 0= No).y`, TRUE ~ `Antibiotics (>=4 weeks) (1 = Yes; 0= No).x`),
         `ANY Meds (>=4 weeks) (1 = yes; 0 = No)` = case_when(is.na(`ANY Meds (>=4 weeks) (1 = yes; 0 = No).x`) ~ `ANY Meds (>=4 weeks) (1 = yes; 0 = No).y`, TRUE ~ `ANY Meds (>=4 weeks) (1 = yes; 0 = No).x`),
         `Internal (I)/External (E)` = case_when(is.na(`Internal (I)/External (E).x`) ~ `Internal (I)/External (E).y`, TRUE ~ `Internal (I)/External (E).x`)) %>% 
  select(-ends_with('.x'), -ends_with('.y'))
########################################################################################################

# human_HMP
human_HMP = human_HMP %>% 
  filter(STArea == "Gut" & VISNO == 1 & SRS != "#N/A") %>%
  add_column(subject_id = NA_character_, study_readable = "HMP1-II", species = "human") %>%
  select(c(sample_id = 'SRS', 'subject_id' = 'SRS', 'study_readable', 'species', "sex" = 'Gender')) %>%
  right_join(human_hmp) %>%
  filter(!is.na(study_readable)) # 238
# human_HMP2 = human_HMP
# check setdiff(human_HMP2$sample_id, human_HMP$sample_id)
# n = 7 in metadata do not have metaphlan results 

# human_MAD
human_MAD = human_MAD %>%
  select(c(sample_id = "sample_alias")) %>%
  mutate(subject_id = sample_id, 
         study_readable = "Madagascar",
         species = "human") # 112

# combine all
animal_metadata = rbind(animal_public, animal_Hills) %>% right_join(animal) %>% filter(!is.na(study_readable))
# animal_metadata2 = rbind(animal_public, animal_Hills) %>% left_join(animal)
# setdiff(animal_metadata2, animal_metadata)
# n = 1 (ABCD_S70) in metaphlan results  does not have metadata 
# n = 126 (with 3M/6M/_b) in metadata do not have metaphlan results 

all_metadata = full_join(human_HMP, human_MAD) %>% full_join(animal_metadata) %>%
  mutate(sample_id_metaphlan = ifelse(is.na(sample_id_metaphlan), sample_id, sample_id_metaphlan)) %>%
  select(sample_id, sample_id_metaphlan, everything())

# all_metadata = all_metadata[,-c(5:12)]

# names(all_metadata) = c("sample_id", "study_readable", "species", "gender", "age", "weight", "antibiotics", "medications", "Internal/External")

write.csv(all_metadata, file = paste0("Taxonomic profiling/metadata.csv"),
          row.names = F)

# 2989

