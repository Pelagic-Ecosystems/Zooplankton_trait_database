# Data patches
# Created: December 28, 2023
# By: P. Pata
#
#
# These are data correction patches applied to the trait database that involve curations done after September 2023. The records removed are based on literature but after analysis by P.Pata, were found to not be the best representation of trait values for the respective species. Justifications of these curations are provided.
#
# Note that these are temporary patches until the next version of the database.
#
require(tidyverse)
`%notin%` <- Negate(`%in%`)

# Load the level 1 data
lvl1 <- read.csv("data_input/Trait_dataset_level1/trait_dataset_level1-2023-08-15.csv") %>% 
  
  # ** Manual curation **
  
  # Remove the records of Gigantocypris muelleri because this is the only species with phosphate composition for Ostracods which is anomalously high in the analysis and the nitrogen composition is also anomalously high compared to the other ostracods- making estimates of stoichiometry and excretion rates for this species not representative of ostracoda in this analysis.
  filter(scientificName != "Gigantocypris muelleri") %>% 
  # Exclude the Clio pyramidata excretion rates which were observed for Clio sulcuta and measured at <0 degC in the antarctic resulting in really high values when standardized to 15degC. The total N and P values are also very high.
  filter(!(scientificName == "Clio pyramidata" & 
             traitName %in% c("excretionRateN_15C","excretionRateN_WSDW_15C",
                              "excretionRateP_15C","excretionRateP_WSDW_15C",
                              "nitrogenTotal", "phosphorusTotal"))) %>% 
  # Exclude the carbon weight of "Euphilomedes producta" because it is 4 orders of magnitude too small compared to the carbon weight of Euphilomedes interpuncta which has more trait data. This results in an erroneous calculated WSC respiration rate that is 3 orders of magnitude greater than any organism.
  filter(!(scientificName == "Euphilomedes producta" & traitName == "carbonWeight")) %>% 
  filter(!(scientificName == "Mysida" & traitName == "carbonWeight")) %>% # very low carbon weights
  filter(!(scientificName == "Clytia" & traitName == "nitrogenPDW")) %>% # value too high but most species have low values
  filter(!(scientificName == "Rathkea octopunctata" & traitName == "nitrogenPDW" 
           & primaryReference == "Matsakis1999")) %>% 
  filter(!(scientificName == "Acartia" & traitName == "phosphorusPDW")) %>% 
  # exclude some very high C:N records which are beyond the range for similar taxa
  filter(catalogNumber %notin% c("38-5223-4","38-2654-1",
                                 "38-630-14","38-630-31")) %>% 
  # exclude clearance rates for Tortanus (Boreotortanus) discaudatus which is a magnitude higher than other Tortanus values
  filter(!(scientificName == "Tortanus (Boreotortanus) discaudatus" &
             traitName %in% c("clearanceRate_15C","clearanceRate_WSC_15C"))) 

# # exclude the WSC respiration rate of Arietellus plumifer because it is an order of magnitude lower than the values of all taxa analyzed
# filter((scientificName == "Arietellus plumifer" & traitName %in%  
#           c("respirationRate_15C","respirationRate_WSC_15C"))) 

# Load the level 2 data
lvl2 <- read.csv("data_input/Trait_dataset_level2/trait_dataset_level2-2023-09-14.csv") %>% 
  # Exclude particularly high carbonPDW for Thysanoessa raschii which was calculated based on literature values of dry weight and carbon weight
  filter(catalogNumber %notin% c("23-3057-1"))
  
  