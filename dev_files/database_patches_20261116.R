# Data patches
# Created: November 16, 2026
# By: P. Pata
#
#
# Corrections for reference list updates and AphiaID updates
require(tidyverse)
library(readxl)
library(worrms)
`%notin%` <- Negate(`%in%`)

# Load the level 1 data
lvl1 <- read.csv("data_input/Trait_dataset_level1/trait_dataset_level1-2023-08-15.csv") 
# TODO load level 2 data and apply similar curation


# The temporary updated files. Check the modify dates!
new.refs <- readxl::read_xlsx("data_input/reference_list_20261109.xlsx")

taxonomy <- readxl::read_excel("data_input/taxonomy_table_20261109.xlsx",
                              range = cell_cols(2:20)) 



# **** Updating the taxonomy table ****
new.taxonomy <- taxonomy %>% 
  filter(!is.na(isUpdated2026)) 

# Search WoRMS for the rest of the taxonomic details especially updated: acceptedName Usage, valid_authority, parentNameUsageID, acceptedNameUsageID
new.worms <- worrms::wm_record(id = unique(new.taxonomy$aphiaID)) %>% 
  distinct(AphiaID, valid_name, valid_authority, scientificname,
           taxonRankID, rank, parentNameUsageID, 
           kingdom, phylum, class, order, family, genus, lsid) %>% 
  mutate(acceptedNameUsage = paste(valid_name, valid_authority)) %>% 
  rename(aphiaID = AphiaID, species = scientificname,
         scientificName = valid_name,
         taxonRank = rank, acceptedNameUsageID = lsid) 

# Update the taxonomy file
updated.taxonomy.sub <- new.taxonomy %>% 
  select(taxonID, aphiaID, verbatimScientificName, majorgroup) %>% 
  left_join(new.worms, by = join_by(aphiaID))

updated.taxonomy <- taxonomy %>% 
  filter(is.na(isUpdated2026)) %>% 
  select(-isUpdated2026) %>% 
  bind_rows(updated.taxonomy.sub) %>% 
  # change duplicated taxonIDs
  mutate(taxonID = case_when(
    taxonID == 6185 ~ 6184,
    taxonID == 6193 ~ 6192,
    taxonID == 6195 ~ 6194,
    .default = taxonID)) %>% 
  arrange(taxonID)

# TODO export new taxonomy file... move updates working document to deprecated
# TODO deal with other changes in taxonomy


# **** Updating the taxonomy information in level 1 data ****

#  Update the relevant rows in the curated trait data tables
lvl1.taxonomy.to.change <- lvl1 %>% 
  filter(taxonID %in% updated.taxonomy.sub$taxonID) %>%   
  # TODO change taxonIDs that are actually duplicates, 6185 -> 6184; 6193 -> 6192; 6195 -> 6194
  mutate(taxonID = case_when(
    taxonID == 6185 ~ 6184,
    taxonID == 6193 ~ 6192,
    taxonID == 6195 ~ 6194,
    .default = taxonID)) %>% 
  # remove columns that are in the taxonomy dataframe, and replace these by taxonID
  select(-c(scientificName, acceptedNameUsageID, acceptedNameUsage,
            taxonRank, kingdom, phylum, class, order, family, genus)) %>% 
  left_join(updated.taxonomy.sub %>% 
              select(taxonID,
                     scientificName, acceptedNameUsageID, acceptedNameUsage,
                     taxonRank, kingdom, phylum, class, order, family, genus),
            by = join_by(taxonID))
  
# Change the rest of the database details to match the edits
lvl1.updated <- lvl1 %>% 
  filter(catalogNumber %notin% lvl1.taxonomy.to.change$catalogNumber) %>% 
  bind_rows(lvl1.taxonomy.to.change)

# rm(updated.taxonomy.sub, new.worms, new.taxonomy)
# rm(lvl1.taxonomy.to.change, lvl1, taxonomy)

# **** Updating the reference IDs (primaryReference, secondaryReference) and DOIs ****
new.refs <- new.refs %>% 
  # TODO include newly changed reference IDs
  filter(!is.na(`newReferenceDOI`))

# TODO apply any changes in referenceDOI first

# For dataframe lvl1.updated, get the rows with refIDs found in the new.refs column, then switch the referenceDOIs
AA.prim <- lvl1.updated %>% 
  mutate(referenceCode = primaryReference) %>% 
  left_join(new.refs %>% 
              select(referenceCode, newReferenceDOI), 
            by = join_by(referenceCode)) %>% 
  # # To simply check and verify
  # filter(primaryReference %in% new.refs$referenceCode) %>%
  # relocate(primaryReference, referenceCode, primaryReferenceDOI, newReferenceDOI) %>% 
  mutate(primaryReferenceDOI = if_else(!is.na(newReferenceDOI),
                                       newReferenceDOI, primaryReferenceDOI))

AA.seco <- lvl1.updated %>% 
  mutate(referenceCode = secondaryReference) %>% 
  left_join(new.refs %>% 
              select(referenceCode, newReferenceDOI), 
            by = join_by(referenceCode)) %>% 
  # # To simply check and verify
  # filter(secondaryReference %in% new.refs$referenceCode) %>%
  # relocate(secondaryReference, referenceCode, secondaryReferenceDOI, newReferenceDOI) %>% 
  mutate(secondaryReferenceDOI = if_else(!is.na(newReferenceDOI),
                                       newReferenceDOI, secondaryReferenceDOI))

