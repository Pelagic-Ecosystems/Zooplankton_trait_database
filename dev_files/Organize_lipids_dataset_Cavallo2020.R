# Lipids trait data standardization
# January 27, 2026
# P. Pata
#
# Data from Cavallo and Peck (2020)
# Cavallo, Alessandro, and Lloyd S Peck. “Lipid Storage Patterns in Marine Copepods: Environmental, Ecological, and Intrinsic Drivers.” ICES Journal of Marine Science 77, no. 5 (2020): 1589–601. https://doi.org/10.1093/icesjms/fsaa070.

# Data provided as supplementary excel spreadsheet. Raw Data in Supp Table 1.
# Supp Table 6 lists the lipid isolation and lipid class analysis names per reference.
# Separate Reference list marked by numbered IDs. This sheet was edited slightly to include primaryReference ID and the publication year of one of the papers.

# Table S1 Description:
# Supplementary Table 1. Raw lipid content data, as extracted from the scientific literature (individual sources corresponding to the Reference IDs listed here available in the "Reference IDs" tab). Collection latitude was often presented as a range, so so the mid-range values are presented in the Latitude column. Latitudinal zone categorisations were made as follows: Polar (66.5-90° N or S), Temperate (23.5-66.5° N or S), and Tropical (23.5° N – 23.5° S). Depth zones: epipelagic (E; 0-200 m), mesopelagic (M; 200-1000 m), bathypelagic (B, 1000-4000 m). Abbreviations: TL = total lipid, DW = dry weight, TAG = triacylglycerol, WE = wax ester, Y = dormancy present, N = dormancy absent.

library(tidyverse)
# source("toolkit.R")


getref <- function(refval = NA, refs.list = refs) {
  # If multiple vals seperated by comma, separate
  rfs <- as.numeric(unlist(str_split(refval, ",")))
  # Extract from the refs list
  ref.sub <- refs.list %>% 
    filter(Ref.code %in% rfs) 
  df <- paste0(ref.sub$ref_id, collapse = "; ")
  return(df)
}


data.fname <- "C:/Users/patri/OneDrive - Université Laval/New_Trait_Data_For_Database/Cavallo Peck Lipids/ICES_Supplementary_Tables_ed.xlsx"

dataset <- readxl::read_xlsx(data.fname,
                             sheet = "Supp. Table 1", skip = 1)

methods <- readxl::read_xlsx(data.fname,
                             sheet = "Supp. Table 6", skip = 1)

refs <- readxl::read_xlsx(data.fname,
                             sheet = "Reference IDs") %>% 
  rename(Ref.code = refID, ref_id = primaryReference)


# Add measurement methods in comments?
colnames(dataset)

Cavallo.long <- dataset %>% 
  rename(verbatimScientificName = Species,
         Ref.code = `Reference ID`,
         verbatimLifeStage = Stage,
         decimalLatitude = Latitude,
         verbatimLocality = `Collection area/coordinates`) %>% 
  # combine life stage and sex columns
  mutate(verbatimLifeStage = if_else(verbatimLifeStage == "Adult",
                                     Sex, verbatimLifeStage)) %>% 
  select(-Sex) %>% 
  # match references by refID
  group_by(Ref.code) %>% 
  mutate(primaryReference = getref(Ref.code)) %>% 
  ungroup() %>% 
  mutate(secondaryReference = "Cavallo2020",
         valueType = "numeric") %>%
  
  # For size assoc: dry weight from TL_ug and TL%DW
  mutate(sizeAssocValue = `TL (µg/individual)` / (`TL (% DW)`/100) / 1000) %>% 
  mutate(notes = if_else(!is.na(sizeAssocValue),
           "Associated size was back calculated or TL and TL%DW.", NA),
         sizeAssocValue = if_else(!is.na(sizeAssocValue),
                                  "dryWeight", NA), 
         sizeAssocUnit = if_else(!is.na(sizeAssocValue),
                                 "mg", NA),
         sizeAssocReference = if_else(!is.na(sizeAssocValue), 
                                      primaryReference, NA) ) %>% 
  
  # Convert WE (%TL) to numeric since one row is Trace
  mutate(`WE (% TL)` = as.numeric(`WE (% TL)`)) %>% 
  
  # Derive WE weight from WE%TL and TL_ug
  mutate(`WE (µg/individual)`= (`WE (% TL)`/100) * `TL (µg/individual)`, 
         `TAG (µg/individual)` = (`TAG (% TL)`/100) * `TL (µg/individual)`) %>% 
  
  # Organize some of the notes
  mutate(verbatimNotes = paste("Collection depth (m):", `Collection depth (m)`,
                               ", Latitudinal zone:", `Latitudinal zone`)) %>% 
  
  # remove feeding guild and dormancy, and some notes
  select(-c(`Feeding guild`:`Reference ID for dormancy categorisation`,
            `Latitudinal zone`, `Collection depth (m)`, `Depth zone`)) %>% 
  
  
  pivot_longer(cols = c(`TL (µg/individual)`:`WE (% TL)`,
                        `WE (µg/individual)`, `TAG (µg/individual)`),
               names_to = "verbatimTraitName",
               values_to = "verbatimTraitValue") %>% 
  # remove NAs
  filter(!is.na(verbatimTraitValue)) %>% 
  
  # Correct the trait names
  mutate(verbatimTraitName = str_replace(verbatimTraitName,
                                         "TL", "total lipid"),
         verbatimTraitName = str_replace(verbatimTraitName,
                                         "DW", "dry weight"),
         verbatimTraitName = str_replace(verbatimTraitName,
                                         "WE", "wax ester"),
         verbatimTraitName = str_replace(verbatimTraitName,
                                         "TAG", "triacylglycerol")) %>% 
  
  # Establish units
  mutate(verbatimTraitUnit = str_extract(verbatimTraitName, "(?<=\\()[^)]+(?=\\))")) 


# Dormant copepods
dormants <- dataset %>% 
  distinct(Species, # Stage,
           Dormancy, `Reference ID for dormancy categorisation`) %>% 
  filter(!is.na(Dormancy)) %>% 
  group_by(`Reference ID for dormancy categorisation`) %>% 
  mutate(primaryReference = getref(`Reference ID for dormancy categorisation`)) %>% 
  ungroup() %>% 
  arrange(Dormancy) 
  
# Temporary save before standardization
save(Cavallo.long, dormants,
     file = "C:/Users/patri/OneDrive - Université Laval/Laval_Postdoc/R_Worskapce/Laval_Arctic_Traits/data/Cavallo2020_lipid_data_temp.RData")



