# Organize TraitZOO data
# February 17 2026
# P. Pata
#
# Datasets: Perhirin et al. 2025 and Limoine et al. 2025

library(tidyverse)
library(here)
`%notin%` <- Negate(`%in%`)

# Additional GZTDv1 files
s.format <- read.csv(here("data_input/trait_dataset_standard_format_20230628.csv"))[-1,]
# Taxonomy table
taxonomy.2023 <- read.csv(here("data_input/taxonomy_table_20230628.csv"))



# A. Check Perhirin et al. 2025 data
# Margaux already organized it to follow the standardized format.
# Need to check verbatim names and units  + taxonomy
# TODO Follow up with Margaux about other trait metadata + references

fname <- "C:/Sync/Dissertation_Study_2/Functional traits lit and data/New_Trait_Data_For_Database/2025 data to add/Perhirin 2025 Fecal Pellets/1225_dataFaecalPelletZoo_Perhirin.csv"

df.per <- read.csv(fname) %>% 
  select(-c(X.1, X, catalogNumber, lifeStage))

# Margaux's standardization
AA <- distinct(df.per, traitName, traitUnit, verbatimTraitName, verbatimTraitUnit)

# Check taxonomy
BB <- distinct(df.per, scientificName, verbatimScientificName,
               kingdom, phylum, class, order, family, genus) %>% 
  arrange(scientificName) %>% 
  # Rename some scientificNames to match the taxonomy table
  mutate(scientificName = case_match(
    scientificName,
    "Acartia bifilosa" ~ "Acartia (Acanthacartia) bifilosa",
    "Acartia clausi" ~ "Acartia (Acartiura) clausi",
    "Acartia tonsa" ~ "Acartia (Acanthacartia) tonsa",
    "Eucalanus pileatus" ~ "Subeucalanus pileatus",
    "Iasis zonaria" ~ "Soestia zonaria",
    "Larvacea" ~ "Appendicularia",
    "Metridia lucens" ~ "Metridia lucens lucens",
    "Oikopleura dioica" ~ "Oikopleura (Vexillaria) dioica",
    "Oikopleura vanhoeffeni" ~ "Oikopleura (Vexillaria) vanhoeffeni",
    "Ostracods" ~ "Ostracoda",
    .default = scientificName
  ))

taxonomy.to.get <- BB %>% 
  filter(scientificName %notin% taxonomy.2023$scientificName)

# Manually assign some lifeStage information from the verbatimScientificName 
CC <- BB %>% 
  filter(grepl("C5|C6F|CIV|CV|CVI|F|M", verbatimScientificName) == T)

df.per <- df.per %>% 
  mutate(verbatimLifeStage = case_match(
    verbatimScientificName,
    "Acartia bifilosa F " ~ "F",
    "Acartia bifilosa M " ~ "M",
    "Calanus finmarchicus CV " ~ "CV",
    "Calanus finmarchicus F " ~ "F",
    "Calanus glacialis CIV " ~ "CIV",
    "Calanus glacialis CIV-CV " ~ "CIV-CV",
    "Calanus glacialis CV " ~ "CV",
    "Calanus glacialis CV-CVI " ~ "CV-CVI",
    "Calanus glacialis C5 " ~ "C5",
    "Calanus hyperboreus CV " ~ "CV",
    "Calanus hyperboreus C6F " ~ "C6F",
    "Calanus hyperboreus CIV " ~ "CIV",
    "Calanus hyperboreus F " ~ "F",
    "Centropages furcatus CV " ~ "CV",
    "Eucalanus bungii C6F " ~ "C6F",
    "Eucalanus pileatus CV " ~ "CV",
    "Gaetanus variabilis C5M/C6F " ~ "C5M/C6F",
    "Metridia longa C6F " ~ "C6F",
    "Metridia okhotensis C6F " ~ "C6F",
    "Neocalanus cristatus C5 " ~ "C5",
    "Neocalanus flemingeri C5 " ~ "C5",
    "Neocalanus plumchrus C5 " ~ "C5",
    "Oithona similis F " ~ "F",
    "Paracalanus sp. CV " ~ "CV",
    "Pseudocalanus elongatus F " ~ "F",
    "Undinula sp. CIV " ~ "CIV",
    .default = NA
  ))



# B. Check Limoine et al. 2025 data
fname <- "C:/Sync/Dissertation_Study_2/Functional traits lit and data/New_Trait_Data_For_Database/2025 data to add/Lemoine 2025/data_compilation_merged.xlsx"

df.lim <- readxl::read_xlsx(fname)


# C. Get new taxonomy from WoRMS


# D. Export results
save(df.per, taxonomy.per.lim)

