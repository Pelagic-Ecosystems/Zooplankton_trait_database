# Explore contents of Arctic Traits Dataset 
# January 20, 2026
# P. Pata
#
# Literature-based trait dataset of marine holoplankton from the European Arctic and sub-Arctic regions
# by T. Durazzano and J. Titocci
# Source: https://data.lifewatchitaly.eu/entities/dataset/269327ce-690f-493d-b129-03d7a6abb1ea

data.fname <- "C:/Users/patri/OneDrive - Université Laval/New_Trait_Data_For_Database/Durazzano_holoplankton_Arctic_traits/Dataset_marine_holoplankton_Arctic.csv"

dataset <- read.csv(data.fname, sep = ";")

colnames(dataset)


# PP Notes: Data table is in a wide format with catalog number by species and lifestage/sex. 
# Columns 1-13 are taxa metadata. 
# The trait records generally follow, Trait | Citation | Occurrence Remarks format.
# Note that records can be multiple for a species-lifestage if from different references or with different comments.
# Data usually numeric or binary, except for geographicDistribution and environmentalPosition
# After elongation, rare instances of no citation (eg, Themisto compressa mean body length which is midpoint of min and max)


# Try to lengthen the data table by trait and then removing blanks.

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
  mutate(secondaryReference = "Durazzano2026") %>% 
  bind_rows(AA2 %>% 
              mutate(valueType = "categorical")) %>% 
  
  # remove no bibliographic citations
  filter(!is.na(bibliographicCitation)) %>% 

# rename columns to verbatim values and comments
  rename(verbatimTraitName = traitName,
         verbatimTraitValue = value,
         verbatimNotes = occurrenceRemarks,
         primaryReferenceDOI = bibliographicCitation)

# TODO verify aphiaIDs and references are in the reference list

rm(AA1, AA2, BB.clutchSize, BB.eggSize, BB.reproMode)



# # Inspect taxa-stage with multiple records (for 5 taxa)
# view(all.data %>% 
#        group_by(acceptedNameUsage, lifeStage, sex, traitName) %>% 
#        filter(n() > 1))

