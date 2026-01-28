# Explore contents of Arctic Traits Dataset 
# January 20, 2026
# P. Pata
#
# Literature-based trait dataset of marine holoplankton from the European Arctic and sub-Arctic regions
# by T. Durazzano and J. Titocci
# Source: https://data.lifewatchitaly.eu/entities/dataset/269327ce-690f-493d-b129-03d7a6abb1ea

library(tidyverse)
`%notin%` <- Negate(`%in%`)

data.fname <- "C:/Users/patri/OneDrive - Université Laval/New_Trait_Data_For_Database/Durazzano_holoplankton_Arctic_traits/Dataset_marine_holoplankton_Arctic.csv"
refs.fname <- "C:/Users/patri/OneDrive - Université Laval/New_Trait_Data_For_Database/Durazzano_holoplankton_Arctic_traits/Dataset_marine_holoplankton_Arctic_reference_key.csv"

# Additional GZTDv1 files
s.format <- read.csv(here("data_input/trait_dataset_standard_format_20230628.csv"))[-1,]
# Taxonomy table
taxonomy.2023 <- read.csv(here("data_input/taxonomy_table_20230628.csv"))


dataset <- read.csv(data.fname, sep = ";")

# Edit the reference list
refs <- read.csv(refs.fname, sep = ";") %>% 
  select(-c(X:X.19)) %>% 
  filter(References != "") %>% 
  # 7 references have repeated rows and are filtered out, the full refences are similar but could be formatted slightly differently.
  distinct(DOI, .keep_all = T)

colnames(dataset)


# PP Notes: Data table is in a wide format with catalog number by species and lifestage/sex. 
# Columns 1-13 are taxa metadata. 
# The trait records generally follow, Trait | Citation | Occurrence Remarks format.
# Note that records can be multiple for a species-lifestage if from different references or with different comments.
# Data usually numeric or binary, except for geographicDistribution and environmentalPosition
# After elongation, rare instances of no citation (eg, Themisto compressa mean body length which is midpoint of min and max)


# --- Lengthen the data table ---- 

# This extraction only works until column 67. For reproductive traits (68 to 83), the references do not match all the columns and can include minimum and maximum values. Columns 84 to 89 and categorical traits.
# Simple numeric traits
AA1 <- dataset[,1:67] %>% 
  pivot_longer(
    cols = -(1:13),
    names_to = c("traitName", ".value"), 
    names_pattern = "^(.*?)(?:_(bibliographicCitation|occurrenceRemarks))?$",
    names_transform = list(.value = ~ ifelse(. == "", "value", .))
  ) %>% 
  filter(!is.na(value))

# Categorical traits
AA2 <- dataset[,c(1:13, 84: 89)] %>% 
  pivot_longer(
    cols = -(1:13),
    names_to = c("traitName", ".value"), 
    names_pattern = "^(.*?)(?:_(bibliographicCitation|occurrenceRemarks))?$",
    names_transform = list(.value = ~ ifelse(. == "", "value", .))
  ) %>% 
  filter(!is.na(value)) %>% 
  filter(value != "")

# Reproductive traits (binary, 1:true)
BB.reproMode <- dataset[,c(1:13,68:73)] %>% 
  pivot_longer(
    cols = -c(1:13, 18, 19),
    names_to = "traitName",
    values_to = "value"
  ) %>% 
  rename(bibliographicCitation = "asexual_bibliographicCitation",
         occurrenceRemarks = "asexual_occurrenceRemarks") %>% 
  filter(!is.na(value))

BB.eggSize <- dataset[,c(1:13,74:78)] %>% 
  pivot_longer(
    cols = c(meanEggSize, minimumEggSize, maximumEggSize),
    names_to = "traitName",
    values_to = "value"
  ) %>% 
  rename(bibliographicCitation = "eggSize_bibliographicCitation",
         occurrenceRemarks = "eggSize_occurrenceRemarks") %>% 
  filter(!is.na(value))

BB.clutchSize <- dataset[,c(1:13,79:83)] %>% 
  pivot_longer(
    cols = c(meanClutchSize, minimumClutchSize, maximumClutchSize),
    names_to = "traitName",
    values_to = "value"
  ) %>% 
  rename(bibliographicCitation = "clutchSize_bibliographicCitation",
         occurrenceRemarks = "clutchSize_occurrenceRemarks") %>% 
  filter(!is.na(value))

# Combine all data frames (set value to character)
all.data <- bind_rows(AA1, BB.clutchSize, BB.eggSize, BB.reproMode) %>% 
  mutate(valueType = "numeric") %>% 
  mutate(value = as.character(value)) %>% 
  mutate(secondaryReference = "Durazzano2026",
         secondaryReferenceDOI = "https://doi.org/10.48372/zzwt-r154") %>% 
  # Add the categorical data
  bind_rows(AA2 %>% 
              mutate(valueType = "categorical")) %>% 
  
  # remove records with no bibliographic citations
  filter(!is.na(bibliographicCitation)) %>% 
  
  # Set verbatim lifestage as combination of lifeStage and sex columns
  mutate(verbatimLifeStage = paste(lifeStage, sex)) %>% 
  mutate(verbatimLifeStage = str_remove(verbatimLifeStage, 
                                        " Undetermined sex")) %>% 
  # acceptedNameUsageID only retains the aphiaID numbers
  mutate(acceptedNameUsageID = gsub("\\D", "", acceptedNameUsageID)) %>% 
  # scientificName is acceptedNameUsage
  # TODO revise the records for 3 species that changed names
  rename(verbatimScientificName = acceptedNameUsage) %>% 
  mutate(scientificName = verbatimScientificName) %>% 
  # acceptedNameUsageID includes  scientificNameAuthorship
  mutate(scientificNameAuthorship = str_replace(scientificNameAuthorship, 
                                                "^(?![\\(]).*(?<![\\)])$", "(\\0)")) %>% 
  mutate(acceptedNameUsage = paste(scientificName, scientificNameAuthorship)) %>% 
  select(-c(scientificNameAuthorship, specificEpithet)) %>% 

# rename columns to verbatim values and comments
  rename(verbatimTraitName = traitName,
         verbatimTraitValue = value,
         verbatimNotes = occurrenceRemarks,
         primaryReferenceDOI = bibliographicCitation) %>% 

# Update primary reference based on DOI and also includes verbatimLocality
    # TODO figure out why not all DOIs are joined properly
  left_join(refs %>% 
              select(-higherGeographyID),
            join_by(primaryReferenceDOI == DOI)) %>% 
    
  # Remove the catalogNumber here
    select(-catalogNumber) %>% 
    
    # exclude records from the Pata & Hunt (2023) GZTDv1 
    filter(grepl("10.1002/lno.12478",primaryReferenceDOI) == F) 
    
# Reorganize columns to match the standard format
all.data.standardized <- s.format %>% 
  bind_rows(all.data)

# verify aphiaIDs and references are in the reference list
# TODO Append these new taxa to the overall taxonomy table
new_taxa <- all.data %>% 
  distinct(scientificName, acceptedNameUsageID, 
           kingdom, phylum, class, order, family, genus) %>% 
  filter(acceptedNameUsageID %notin% taxonomy.2023$aphiaID)

# Which columns in the standard GZTDv1 format are not in the long data table?
setdiff(colnames(s.format), colnames(all.data))

# TODO need verbatimTraitUnits
# TODO Columns to standardize: lifeStage, traitName, traitValue, traitUnit

# TODO update new references in the overall reference table, make sure full references are not duplicated.
# TODO Establish the reference codes according to FirstAuthorYear


rm(AA1, AA2, BB.clutchSize, BB.eggSize, BB.reproMode)