# Trait dataset paper tool kit
# Created by P. Pata
# Last modified April 13, 2023
#
# ***** Trait dataset tool kit *****
# 
# This file contains a set of functions that were utilized in developing and 
# curating the zooplankton trait dataset. A list of the functions with short
# descriptions are reported below.
# 
# ** Data curation functions **
#
# pooled.sd(SDarray, Narray): Calculates the pooled SD from two SDs. This is 
#   similar to the weighted mean. Note that this will only work if both arrays 
#   have N > 1.
# getSpeciesMeanSD(traitValue, dispersionSD, individualCount): Calculate the 
#   weighted means and SDs of traits when the dataset is a mix of individual
#   values and averaged records.
# geomean(x): Calculates the geometric mean of an array.
# scaleFun(x): Returns a number with two decimal places.
# convertRateT2WS(df, sizeAssoc, trtName1, trtName2, trtUnit):
#   Converts a by individual rate to a weight-specific rate.
# convertPW2Total(df, sizeAssoc, trtName1, trtName2, trtUnit = "mg"):
#   Converts percent composition value to total bulk composition value.
# convertTotal2PW(df, sizeAssoc, trtName1, trtName2, trtUnit = "percent"):
#   Converts total bulk composition to percent composition relative to a weight.
# mergeTraitDetails(df, trait.directory): Combine duplicated trait records of 
#   the same value but from different references.
# summarise.sizeAssoc(df): Calculates the average associated size of each trait
#   and taxon pair in the df. This returns a sizeAssoc data frame with the mean
#   N, and sd values which would be merged into the main df.
# updateIDs(data, trait.directory): Updates the trait and taxon-trait IDs.
# standardizeTable(file.list, s.format, taxonomy, lifestagelist, trait.directory):
#   Reorders the column and converts the names and formats of the traits.
# annotateLifeStage(data, lifestagelist): Modifies how the life stage was coded,
#   assigns the life stage ID, and records the original life stage.
# annotateTaxonomy(data, taxonomy): Assigns the taxonomic ID and the taxonomic 
#   ranks it belongs to. The recorded scientific name would be stored in the 
#   column verbatimScientificName.
# adjustClass(data, s.format): Changes the data class for numerical columns into
#   character variables to facilitate merging of data frames into a single file.
# cleanStrings(A): Cleans spaces, NAs, and blanks from strings for notes and 
#   references. Also removes duplicated strings separated by ;.
# cleanStringsList(A): Apply cleanStrings() to a list.
# standard_temp(rate0, t0, t, Q10): Calculate the value of a rate to a reference
#   temperature using a Q10 coefficient.
# cleanScientificName(name): Cleans the verbatim scientific name from trailing 
#   spaces and life stage information. Also removes the sp, aff, and cf texts.
# assignMajorGroup(taxonomy): Assigns the zooplankton major group in a taxonomy
#   data frame.
# getMaxObsNum(traitName, taxonID, df): Gets the maximum observation number for 
#    a trait and species.
#
#
# ** Data imputation functions **
# 
# getGroupLevelValue(taxon, trait, gen.level, trait.df, taxonomy.df)): Calculates
#   the group-level mean, standard deviation, and number of species-level 
#   observations for generalization of a given trait for a taxon. 
#
# conv.allom(W,a,b,base): General equation for allometric conversions. 
#   The default base is the natural log.
# getpval(x) : Calculates a p value from a regression model F statistic object.
#
# getRegressionModel(df, grp, X, Y, model): Derives the regression model between  
#   two trait variables for a zooplankton group. 
# calculateFromModel(df, model, excludeWithLit, applyToGeneralized, 
#   excludeCalculated): Calculates a trait using a regression model.
# plotAllometric(df, grp, X, Y): Plots the distribution of values between two 
#   traits and a regression line. [!!! Should be sensitive to if OLS or RMA regressions.]
# plotRegModel(model): Plots a simple regression model result with 95% confidence
#   intervals. 
# calculate.WSRates(trait.df, traits.calculated, trait.X, trait.Y): Calculates 
#   the weight specific rate for traits with calculate rate values derived
#   from regression equations.
#
# *** Functions for trait correlations ***
#
# calc_p_value(x, y): Calculates the p.value from either a t.test or the c
#    or.test functions. Need to manually switch 
# get_pairwise.N(x,y): Calculates the number of pair of two variables with data.
# widen_traits(trait.table, trait.list): Subsets the trait dataset to only 
#    the traits specified in the list. Only retains taxon information.
# correlate_traits(trait.table, trait.list): Calculates correlations between
#    variables in trait.list and returns a long data frame of the correlation,
#    N data pairs, and p-value.
# corphylo.tests(trait.num.phylo, trait.list, phylo.tree, zoop_vcv) :
#   This runs different correlation tests using the phyr::cor_phylo() function.
#   Note that this temporarily requires the phylo.tree and vcv objects to 
#   explore the effect of which phylo signal to use. This returns a long table
#   or correlation results.


require(dplyr)
require(lmodel2)
require(lubridate)
require(corrr)

# Helper functions ----

`%notin%` <- Negate(`%in%`)
scaleFUN <- function(x) sprintf("%.2f", x)
geomean <- function(x) {  exp(mean(log(x), na.rm = TRUE))  }

# Transformations for proportional or percent data. 
#  According to http://strata.uga.edu/8370/rtips/proportions.html. Logit is 
#  preferred for regressions because range is (-inf, inf), while arcsine is 
#  preferred for ordination and clustering because range is [0, 1] which 
#  works for multivariate methods. Both are linear at 0.3-0.7 range, logit 
#  gives a stronger transformation at the ends.

# arcsine transformation 
asinTransform <- function(x) { 
  # if percent data, change to proportion
  if (max(x) > 1 & min(x) >= 0 & max(x) <= 100) {
    x <- x/100
  } else if(min(x) < 0 | max(x) > 100){
    stop("Error: Data must either be proportion [0,1] or percentage [0, 100]")
  }
  asin(sqrt(x)) 
}

# logit transformation 
logitTransform <- function(x) { 
  # if percent data, change to proportion
  if (max(x) > 1 & min(x) >= 0 & max(x) <= 100) {
    x <- x/100
  } else if(min(x) <= 0 | max(x) > 100){
    stop("Error: Data must either be proportion [0,1] or percentage [0, 100]")
  }
  log(x/(1-x)) 
}


# pooled sd
pooled.sd <- function(SDarray, Narray){
  sqrt(sum((Narray-1)*(SDarray)^2) / (sum(Narray) - length(Narray)) )
  # # if inputs are: x,y,Nx,Ny
  # sqrt( ((Nx-1)*x^2 + (Ny-1)*y^2) / (Nx+Ny-2))
}

# Data curation functions ----

# Calculate the weighted means and SDs of traits when the dataset is a mix of 
#   individual values and averaged records.
getSpeciesMeanSD <- function(traitValue, dispersionSD, individualCount){
  # extract indiv values and grouped values
  val.indiv <- traitValue[which(individualCount == 1)]
  val.group.mean <- traitValue[which(individualCount > 1)]
  val.group.N <- individualCount[which(individualCount > 1)]
  val.group.sd <- dispersionSD[which(individualCount > 1)]
  # calculate mean of the individual values
  val.indiv.mean <- mean(val.indiv)
  val.indiv.N <- length(val.indiv)
  val.indiv.sd <- sd(val.indiv, na.rm = TRUE)
  
  # calculate N
  N <- sum(individualCount)
  # calculate weighted mean
  mean <- weighted.mean(c(val.indiv.mean, val.group.mean),
                        c(val.indiv.N, val.group.N))
  # calculate pooled sd
  # catch if val.indiv.N = 1, calculate the pooled sd from the groups vars only, 
  #  note that the mean and N are of all though. This is obviously not accurate 
  #  but would be a better estimate compared to calculating the unweighted sd, 
  #  especially when the group.N sums to a large sample.
  if (val.indiv.N == 1) {
    sd <- pooled.sd(c(val.group.sd),
                    c(val.group.N))
  } else {
    sd <- pooled.sd(c(val.indiv.sd, val.group.sd),
                    c(val.indiv.N, val.group.N))
  }
  # calculate standard error
  se <- sd / sqrt(N)
  
  
  return(list("N" = N, "mean" = mean, "sd" = sd, "se" = se))
}


# for converting individual rates to weight specific
convertRateInd2WS <- function(df, sizeAssoc, trtName1, trtName2, trtUnit) {
  D.1 <- df %>% 
    filter(traitName == trtName1)
  D.2 <- df %>% 
    filter(traitName == trtName2) %>% 
    filter(taxonID %notin% D.1$taxonID) %>% 
    # Derived columns are the original columns
    mutate(verbatimScientificName = scientificName,
           verbatimTraitName = traitName,
           verbatimTraitValue = as.character(traitValue),
           verbatimTraitUnit = traitUnit,
           verbatimLifeStage = lifeStage,
           verbatimNotes = notes,
           aggregateMeasure = FALSE,
           isDerived = TRUE,
           catalogSource = catalogNumber) %>% 
    select(-starts_with("sizeAssoc")) %>% 
    left_join(sizeAssoc, by = "taxonID") %>% 
    filter(!is.na(sizeAssocValue)) %>% 
    group_by(traitID, taxonID) %>% 
    mutate(traitName = trtName1,
           traitValue = traitValue/sizeAssocValue,
           traitUnit = trtUnit,
           basisOfRecord = "derived from related trait",
           dispersionSD = NaN,
           individualCount = 1,
           notes = paste0("Calculated from ",trtName1," and ",
                          trtName2,"; based on", catalogNumber),
           isDerived = TRUE,
           maxObsNum = 0) %>% 
    ungroup() %>% 
    standardizeID(trait.directory)
}


# for converting weight specific rates to individual-specific
convertRateWS2Ind <- function(df, sizeAssoc, trtName1, trtName2, trtUnit) {
  D.1 <- df %>% 
    filter(traitName == trtName1)
  D.2 <- df %>% 
    filter(traitName == trtName2) %>% 
    filter(taxonID %notin% D.1$taxonID) %>% 
    # Derived columns are the original columns
    mutate(verbatimScientificName = scientificName,
           verbatimTraitName = traitName,
           verbatimTraitValue = as.character(traitValue),
           verbatimTraitUnit = traitUnit,
           verbatimLifeStage = lifeStage,
           verbatimNotes = notes,
           aggregateMeasure = FALSE,
           isDerived = TRUE,
           catalogSource = catalogNumber) %>% 
    select(-starts_with("sizeAssoc")) %>% 
    left_join(sizeAssoc, by = "taxonID") %>% 
    filter(!is.na(sizeAssocValue)) %>% 
    group_by(traitID, taxonID) %>% 
    mutate(traitName = trtName1,
           traitValue = traitValue*sizeAssocValue,
           traitUnit = trtUnit,
           basisOfRecord = "derived from related trait",
           dispersionSD = NaN,
           individualCount = 1,
           notes = paste0("Calculated from ",trtName1," and ",
                          trtName2,"; based on", catalogNumber),
           isDerived = TRUE,
           maxObsNum = 0) %>% 
    ungroup() %>% 
    standardizeID(trait.directory)
}

# for converting between bulk and percent composition
convertPW2Total <- function(df, sizeAssoc, trtName1, trtName2, trtUnit = "mg"){
  D.1 <- df %>% 
    filter(traitName == trtName1)
  D.2 <- df %>% 
    filter(traitName  == trtName2) %>% 
    filter(taxonID %notin% D.1$taxonID) %>% 
    # Derived columns are the original columns
    mutate(verbatimScientificName = scientificName,
           verbatimTraitName = traitName,
           verbatimTraitValue = as.character(traitValue),
           verbatimTraitUnit = traitUnit,
           verbatimLifeStage = lifeStage,
           verbatimNotes = notes,
           aggregateMeasure = FALSE,
           isDerived = TRUE,
           catalogSource = catalogNumber) %>% 
    select(-starts_with("sizeAssoc")) %>% 
    left_join(sizeAssoc, by = "taxonID") %>% 
    filter(!is.na(sizeAssocValue)) %>% 
    group_by(traitID, taxonID) %>% 
    mutate(traitName = trtName1,
           traitValue = traitValue*sizeAssocValue / 100,
           traitUnit = trtUnit,
           basisOfRecord = "derived from related trait",
           dispersionSD = NaN,
           individualCount = 1,
           notes = paste0("Calculated from ",trtName1," and ",
                          trtName2,"; based on", catalogNumber),
           isDerived = TRUE,
           maxObsNum = 0) %>% 
    ungroup() %>% 
    standardizeID(trait.directory)
    # updateIDs(trait.directory)
}  


convertTotal2PW <- function(df, sizeAssoc, trtName1, trtName2, trtUnit = "percent") {
  D.1 <- df %>% 
    filter(traitName == trtName1)
  D.2 <- df %>% 
    filter(traitName  == trtName2) %>% 
    filter(taxonID %notin% D.1$taxonID) %>% 
    # Derived columns are the original columns
    mutate(verbatimScientificName = scientificName,
           verbatimTraitName = traitName,
           verbatimTraitValue = as.character(traitValue),
           verbatimTraitUnit = traitUnit,
           verbatimLifeStage = lifeStage,
           verbatimNotes = notes,
           aggregateMeasure = FALSE,
           isDerived = TRUE,
           catalogSource = catalogNumber) %>% 
    select(-starts_with("sizeAssoc")) %>% 
    left_join(sizeAssoc, by = "taxonID") %>% 
    filter(!is.na(sizeAssocValue)) %>% 
    group_by(traitID, taxonID) %>% 
    mutate(traitName = trtName1,
           traitValue = traitValue/sizeAssocValue * 100,
           traitUnit = trtUnit,
           basisOfRecord = "derived from related trait",
           dispersionSD = NaN,
           individualCount = 1, 
           isDerived = TRUE,
           notes = paste0("Calculated from ",trtName1," and ",
                          trtName2,"; based on", catalogNumber),
           maxObsNum = 0) %>% 
    ungroup() %>% 
    # updateIDs(trait.directory) %>% 
    standardizeID(trait.directory) %>% 
    filter(traitValue <= 100)
}

calculateRatio <- function(df, ratio, trtName1, trtName2){
  D.0 <- df %>% 
    filter(traitName == ratio)
  D.1 <- df %>% 
    filter(traitName == trtName1) %>% 
    filter(taxonID %notin% D.0$taxonID) %>% 
    select(catalogNumber, taxonID, traitValue1 = traitValue, scientificName, 
          taxonRank, acceptedNameUsageID, acceptedNameUsage, kingdom, phylum, 
           class, order, family, genus, majorgroup)
  D.2 <- df %>% 
    filter(traitName == trtName2) %>% 
    filter(taxonID %notin% D.0$taxonID) %>% 
    select(catalogNumber, taxonID, traitValue2 = traitValue)
  D.3 <- inner_join(D.1, D.2, by = "taxonID") %>% 
    mutate(traitValue = traitValue1/traitValue2, traitName = ratio, 
           traitUnit = NA, 
           valueType = "numeric",
           basisOfRecord = "derived from related trait",
           notes = paste0("Calculated from ",trtName1," and ",
                          trtName2),
           catalogSource = paste0(catalogNumber.x, "; ", catalogNumber.y),
           isDerived = TRUE,
           observationNumber = 1, maxObsNum = 1, traitID = NA) %>% 
    # Above is ratio by weight so convert to molar ratio
    mutate(traitValue = if_else(traitName == "ratioCN", 
                                traitValue * (14.0067/12.011), traitValue)) %>% 
    mutate(traitValue = if_else(traitName == "ratioCP", 
                                traitValue * (30.97376/12.011), traitValue)) %>% 
    mutate(traitValue = if_else(traitName == "ratioNP", 
                                traitValue * (30.97376/14.0067), traitValue)) %>% 
    # exclude values outside the range of observed values
    filter(traitValue >= min(D.0$traitValue) & traitValue <= max(D.0$traitValue)) %>% 
    select(-c(catalogNumber.x, catalogNumber.y, traitValue1, traitValue2)) %>% 
    standardizeID(trait.directory)
}

# For merging duplicated trait values or taxa-level averages. Note that this should be called for a grouped dataframe.
mergeTraitDetails <- function(df, trait.directory, rev.by = "P. Pata") {
  df <- df %>% 
    mutate(dup = n())
  
  df.same <- df %>% 
    filter(dup == 1) %>% 
    mutate(aggregateMeasure = FALSE) %>% 
    ungroup()
  
  df.updated <- df %>% 
    filter(dup > 1) %>% 
    mutate(aggregateMeasure = TRUE)  %>% 
    mutate(primaryReference = paste(primaryReference, collapse = "; "),
           primaryReferenceDOI = paste(primaryReferenceDOI, collapse = "; "),
           secondaryReference = paste(secondaryReference, collapse = "; "),
           secondaryReferenceDOI = paste(secondaryReferenceDOI, collapse = "; "),
           lifeStage = paste(lifeStage, collapse = "; "),
           sizeType = paste(sizeType, collapse = "; "), 
           # sizeAssocName = paste(sizeAssocName, collapse = "; "), 
           # sizeAssocUnit = paste(sizeAssocUnit, collapse = "; "), 
           # sizeAssocValue = paste(sizeAssocValue, collapse = "; "), 
           # sizeAssocReference = paste(sizeAssocReference, collapse = "; "), 
           verbatimLocality = paste(verbatimLocality, collapse = "; "), 
           # decimalLongitude = paste(decimalLongitude, collapse = "; "), 
           # decimalLatitude = paste(decimalLatitude, collapse = "; "), 
           notes = paste(notes, collapse = "; "),
           basisOfRecord = paste(basisOfRecord, collapse = "; "), 
           basisOfRecordDescription = paste(basisOfRecordDescription, collapse = "; "),
           # Updated March 22, 2023: Now the verbatim records are not retained
           #   when merging rows. Instead the data can be extracted from the
           #   catalogSource.
           verbatimScientificName = NA, verbatimTraitName = NA,
           verbatimTraitValue = NA, verbatimTraitUnit = NA,
           verbatimLifeStage = NA, verbatimTemperature = NA,
           verbatimNotes = NA, 
           # verbatimScientificName = paste(verbatimScientificName, collapse = "; "), 
           # verbatimTraitName = paste(verbatimTraitName, collapse = "; "), 
           # verbatimTraitValue = paste(verbatimTraitValue), 
           # verbatimTraitUnit = paste(verbatimTraitUnit, collapse = "; "), 
           # verbatimLifeStage = paste(verbatimLifeStage, collapse = "; "), 
           # verbatimNotes = paste(verbatimNotes, collapse = "; "),
           catalogSource = paste(catalogNumber, collapse = "; ") ) %>% 
    # clean strings
    mutate(primaryReference = cleanStrings(primaryReference),
           primaryReferenceDOI = cleanStrings(primaryReferenceDOI),
           secondaryReference = cleanStrings(secondaryReference),
           secondaryReferenceDOI = cleanStrings(secondaryReferenceDOI),
           lifeStage = cleanStrings(lifeStage),
           sizeType = cleanStrings(sizeType), 
           # sizeAssocName = cleanStrings(sizeAssocName), 
           # sizeAssocUnit = cleanStrings(sizeAssocUnit), 
           # sizeAssocValue = cleanStrings(sizeAssocValue), 
           # sizeAssocReference = cleanStrings(sizeAssocReference), 
           verbatimLocality = cleanStrings(verbatimLocality), 
           # decimalLongitude = cleanStrings(decimalLongitude), 
           # decimalLatitude = cleanStrings(decimalLatitude), 
           notes = cleanStrings(notes),
           basisOfRecord = cleanStrings(basisOfRecord), 
           basisOfRecordDescription = cleanStrings(basisOfRecordDescription), 
           # verbatimScientificName = cleanStrings(verbatimScientificName), 
           # verbatimTraitName = cleanStrings(verbatimTraitName), 
           # verbatimTraitValue = cleanStrings(verbatimTraitValue), 
           # verbatimTraitUnit = cleanStrings(verbatimTraitUnit), 
           # verbatimLifeStage = cleanStrings(verbatimLifeStage), 
           # verbatimNotes = cleanStrings(verbatimNotes),
           catalogSource = cleanStrings(catalogNumber) )%>% 
    distinct(catalogSource, .keep_all = TRUE) %>% 
    # mutate(observationNumber = NA) %>% 
    # updateIDs() %>% 
    mutate(notes = cleanStrings(paste("Trait value merged from multiple records.; ", notes))) %>% 
    # update upload date
    mutate(uploadBy = rev.by, uploadDate = as.character(ymd(Sys.Date()))) %>% 
    ungroup() 
  
  df <- bind_rows(df.same, df.updated) %>% 
    select(-dup)
}

# This function removes the values in the verbatim columns when processing
#   level 2 or 3 data such as merging multiple records or deriving the mean.
#   As a fail-safe, this only deletes the verbatim records when the value for
#   catalogSource is not empty, i.e., the record was based on something else.
removeVerbatimRecords <- function(df) {
  # Preserve the order of rows
  df <- df %>% 
    rownames_to_column(var = "rowIndex")
  
  # Separate into if merged/derived or not
  df.A <- df %>% 
    filter(is.na(catalogSource))
  
  df.B <- df %>% 
  filter(!is.na(catalogSource)) %>% 
    mutate(verbatimScientificName = NA, verbatimTraitName = NA,
           verbatimTraitValue = NA, verbatimTraitUnit = NA,
           verbatimLifeStage = NA, verbatimTemperature = NA,
           verbatimNotes = NA)
  
  # merge back and rearrange rows
  df <- bind_rows(df.A, df.B) %>% 
    arrange(rowIndex) %>% 
    select(-rowIndex)
  
  return(df)
}

# Function for selecting the average associated size for a trait value. If there
#   are multiple size types for a trait, the one with the highest sample size 
#   will be selected. References will also be summarized. Note that this may be 
#   a useful reference but it does not correspond to the species-average size 
#   value. For analyses requiring comparing the trait with size, it might be 
#   better to use the species-averaged size instead of the averaged sizeAssocValue.
summarise.sizeAssoc <- function(df) {
  A <- df %>% 
    select(traitID, taxonID,
           sizeAssocName, sizeAssocValue, sizeAssocUnit, sizeAssocReference) %>% 
    filter(!is.na(sizeAssocValue)) %>% 
    group_by(traitID, taxonID, sizeAssocName, sizeAssocUnit) %>% 
    # mutate(N = n(),
    #        Nmax = max(N)) %>% 
    # ungroup() %>% 
    # filter(N == Nmax) %>% 
    # group_by(traitID, taxonID, traitUnit) %>% 
    mutate(sizeAssocN = n(),
           sizeAssocSD = sd(sizeAssocValue, na.rm = TRUE),
           sizeAssocValue = mean(sizeAssocValue, na.rm = TRUE),
           sizeAssocReference = cleanStrings(paste(sizeAssocReference, 
                                                   collapse = "; "))) %>% 
    ungroup() %>% 
    distinct(traitID, taxonID, 
             sizeAssocName, sizeAssocValue, sizeAssocUnit, sizeAssocReference,
             sizeAssocN, sizeAssocSD) %>% 
    # arrange and get first instance
    group_by(traitID, taxonID) %>% 
    arrange(sizeAssocN) %>% 
    filter(row_number() == 1)
}

# Update the trait ID and unit when standardizing raw data
standardizeID <- function(data, trait.directory){
  col.order <- colnames(data)
  data <- data %>% 
    select(-c(traitID)) %>% 
    left_join(distinct(trait.directory, traitID, traitName), by = "traitName")  %>% 
    relocate(col.order)
}
standardizeUnit <- function(data, trait.directory){
  col.order <- colnames(data)
  data <- data %>% 
    select(-c(traitUnit)) %>% 
    left_join(distinct(trait.directory, traitName, traitUnit), by = "traitName")  %>% 
    relocate(col.order)
}

# For updating trait IDs and taxon-trait-IDs 
# master.ID.List is a dataframe that keep tracks of the maximum observation 
#   of each traitID and taxonID. Need to regenerate this df after every time 
#   updateIDs is called
updateIDs <- function(data) {
  # Create and update a masterIDList file which is loaded and updated by this function.
  #  This would be a safer approach when records are curated so a trait and species 
  #  may completely be excluded and thus the maxObsNum resets to 0.
  load("data_input/master_id_list.RData")
  
  col.order <- colnames(data)
  data <- data %>% 
    select(-maxObsNum) %>% 
    left_join(master.ID.List, by = c("traitID", "taxonID")) %>% 
    group_by(traitID,taxonID) %>% 
    # if there is no maxObsNum set to zero
    mutate(maxObsNum = if_else(is.na(maxObsNum), 0, maxObsNum)) %>% 
    # if there is no observation number, assign one
    mutate(observationNumber = if_else(!is.na(observationNumber),
                                       observationNumber, 
                                       maxObsNum + row_number())) %>% 
    # update the maxObsNum to either be the original value or the max updated
    #   observation number.
    mutate(maxObsNum = max(c(observationNumber, maxObsNum))) %>% 
    # if there is no catalogNumber, assign one
    mutate(catalogNumber = if_else(is.na(catalogNumber), 
                                  paste0(traitID,"-",taxonID,"-",observationNumber),
                                  catalogNumber)) %>% 
    ungroup() %>% 
    relocate(col.order)
  
  # Update the list of maxObsNum
  # master.ID.List <- master.ID.List %>% 
  master.ID.List <- master.ID.List %>% 
    bind_rows(distinct(data, traitID, taxonID, maxObsNum)) %>% 
    group_by(traitID, taxonID) %>% 
    summarise(maxObsNum = max(maxObsNum), .groups = "drop")
  
  save(master.ID.List, file = "data_input/master_id_list.RData")
  
  # print(paste0("Updated master ID list. N rows = ", nrow(master.ID.List),
  #              " N maxObsNum = ", length(unique(master.ID.List$maxObsNum))))
  
  return(data)
}

# Functions for standardizing trait files 
standardizeTable <- function(file.list, s.format, taxonomy, lifestagelist,
                             trait.directory, rev.by = "P. Pata") {
  data <- s.format
    # # These columns are set during database curation
    # dplyr::select(-c(catalogSource, maxObsNum, aggregateMeasure))
  
  for(i in file.list) {
    d <- openxlsx::read.xlsx(paste0(infol,i)) %>% 
      mutate(secondaryReference = as.character(secondaryReference))
    if("traitUnit" %notin% colnames(d)){
      d <- d %>% 
        mutate(traitUnit = NA)
    }
    if("verbatimTraitUnit" %notin% colnames(d)){
      d <- d %>% 
        mutate(verbatimTraitUnit = NA)
    }
    d <- d %>%  
      mutate(traitValue = as.character(traitValue),
             verbatimTraitValue = as.character(verbatimTraitValue)) %>% 
      # set traitID, traitName based on trait directory
      select(-c(traitName, traitUnit)) %>%
      left_join(select(trait.directory, traitID, traitName,
                       verbatimTraitName, verbatimTraitUnit),
                by = c("verbatimTraitName","verbatimTraitUnit"))
    
    if("individualCount" %in% colnames(d)) {
      d$individualCount <- as.numeric(str_trim(d$individualCount))
    }
    if("assocTemperature" %in% colnames(d)) {
      d$assocTemperature <- as.character(d$assocTemperature)
      d$verbatimTemperature <- as.character(d$assocTemperature)
    }
    if("lifeStage" %in% colnames(d)) {
      d$lifeStage <- as.character(d$lifeStage)
    }
    if("verbatimLocality" %in% colnames(d)) {
      d$verbatimLocality <- as.character(d$verbatimLocality)
    }
    if("notes" %in% colnames(d)) {
      d$notes <- as.character(d$notes)
    }
    # Update some fields to match DWC (March 16, 2023)
    if("longitude" %in% colnames(d)) {
      d <- d %>% 
        rename(decimalLongitude = longitude)
    }
    if("latitude" %in% colnames(d)) {
      d <- d %>% 
        rename(decimalLatitude = latitude)
    }
    if("rank" %in% colnames(d)){
      d <- d %>% 
        rename(taxonRank = rank)     
    }
    
    data <- data %>% 
      bind_rows(d)
  }
  # Filter duplicates (Updated Feb 11 2023)
  data <- data %>% 
    distinct(verbatimScientificName, verbatimTraitName, verbatimTraitUnit, 
             verbatimTraitValue, lifeStage,
             primaryReference, secondaryReference,
             assocTemperature, .keep_all = TRUE)
    
  # Annotate taxonomy information
  data <- annotateTaxonomy(data, taxonomy) %>% 
    # Set stage ID
    annotateLifeStage(lifestagelist) %>% 
    # Record file generation time stamp
    mutate(uploadBy = rev.by, uploadDate = ymd(Sys.Date())) %>% 
    # Remove not in mesozooplankton major groups of interest
    filter(!is.na(majorgroup)) %>% 
    # assign catalogNumber
    group_by(traitID,taxonID) %>% 
    mutate(observationNumber = row_number()) %>% 
    mutate(catalogNumber = paste0(traitID,"-",taxonID,"-",observationNumber)) %>% 
    mutate(maxObsNum = max(observationNumber)) %>% 
    ungroup() %>% 
    # trait unit = verbatim trait unit if empty
    mutate(traitUnit = if_else(is.na(traitUnit), verbatimTraitUnit, traitUnit)) %>% 
    # assign trait value source
    mutate(basisOfRecord = "literature") %>% 
    relocate(colnames(s.format))
}

annotateLifeStage <- function(data, lifestagelist) {
  data <- data %>% 
    mutate(verbatimLifeStage = lifeStage) %>% 
    select(-c(lifeStage, stageID)) %>% 
    left_join(lifestagelist, by = "verbatimLifeStage")
}

annotateTaxonomy <- function(data, taxonomy) {
  data <- data %>% 
    select(-c(taxonID, scientificName, acceptedNameUsageID, acceptedNameUsage,
              taxonRank, kingdom, phylum, class, order, family, genus,
               majorgroup)) %>% 
    mutate(vsn_ed = cleanScientificName(verbatimScientificName)) %>% 
    left_join(select(taxonomy, verbatimScientificName, scientificName, 
                     taxonID, acceptedNameUsageID, 
                     acceptedNameUsage,
                     taxonRank, kingdom, phylum, class, order, family, genus, 
                     majorgroup), 
              by = c("vsn_ed" = "verbatimScientificName")) %>%
    select(-vsn_ed)
}

# Function for changing class types to merge separate level 1 data files.
# Need to improve this code to include all columns
adjustClass <- function(data, s.format){
  data <- data %>% 
    mutate(traitUnit = as.character(traitUnit),
           lifeStage = as.character(lifeStage),
           assocTemperature = as.character(assocTemperature),
           primaryReferenceDOI = as.character(primaryReferenceDOI),
           secondaryReferenceDOI = as.character(secondaryReferenceDOI),
           uploadDate = as.character(uploadDate),
           sizeType = as.character(sizeType),
           sizeAssocName = as.character(sizeAssocName),
           sizeAssocUnit = as.character(sizeAssocUnit),
           sizeAssocReference = as.character(sizeAssocReference),
           verbatimLocality = as.character(verbatimLocality),
           # species = as.character(species),
           basisOfRecordDescription = as.character(basisOfRecordDescription),
           verbatimTraitUnit = as.character(verbatimTraitUnit),
           verbatimLifeStage = as.character(verbatimLifeStage),
           verbatimTemperature = as.character(verbatimTemperature),
           verbatimNotes = as.character(verbatimNotes),
           catalogSource = as.character(catalogSource),
           aggregateMeasure = as.logical(aggregateMeasure),
           isDerived = as.logical(isDerived),
           errorRisk = as.numeric(errorRisk),
           errorRiskRank = as.character(errorRiskRank),
           notes = as.character(notes))
  
  # for (i in c(1:ncol(data))) {
  #   if(colnames(data[i]) == colnames(s.format[i])){
  #     if (class(s.format[,i]) == "character") {
  #       data[,i] <- as.character(data[,i])
  #     }
  #     if (class(s.format[,i]) == "numeric") {
  #       data[,i] <- as.numeric(data[,i])
  #     }
  #   }
  # }
}

# function for cleaning concatenated strings
cleanStrings <- function(A){
  A <- str_split(A, pattern =";", simplify = TRUE) 
  A <- str_trim(A)
  A <- sort(A)
  A <- str_replace(A, "NA", replacement = "")
  A <- A[A != '']
  A <- unique(as.list(A))
  A <- paste(A, collapse = "; ")
  A <- as.character(A)
}
cleanStringsList <- function(A){
  A <- lapply(A, cleanStrings)
}

standard_temp <- function(rate0, t0, t = 15, Q10 = 2.8) {
  10^(log10(rate0) + log10(Q10)*((t-t0)/10) )
}

cleanScientificName <- function(name){
  # remove trailing *sp., ?-. symbols, and anything inside parenthesis
  # remove trailing life stage information
  name <- gsub(" aff\\.","", name)
  name <- gsub(" cf\\.","", name)
  name <- gsub(" C4| C5| CI| CIV| CV| III| IV| VI| NI| NII| NIII| NIV| NV","",name)
  name <- str_replace(name, "\\s*\\([^\\)]+\\)", "")
  name <- str_replace_all(name, "[[:punct:]]", "")
  name <- str_replace_all(name, "[:digit:]", "")
  name <- gsub(" sp$","", name)
  name <- gsub(" spp$","", name)
  name <- gsub(" V$","", name)
  name <- str_replace(name, " female","")
  name <- str_replace(name, " male","")
  name <- str_replace(name, "Agg","")
  name <- str_replace(name, "aggregate","")
  name <- str_replace(name, "solitary","")
  name <- str_squish(str_trim(name))
  name
}

# Assign major group based on a row/s of taxonomy information
assignMajorGroup <- function(taxonomy){
  taxonomy <- taxonomy %>% 
    mutate(majorgroup = "") %>% 
    mutate(majorgroup = if_else(class %in% c("Polychaeta"),"Polychaete",majorgroup)) %>% 
    mutate(majorgroup = if_else(phylum %in% c("Chaetognatha"),"Chaetognath",majorgroup)) %>% 
    mutate(majorgroup = if_else(class %in% c("Branchiopoda"),"Cladoceran",majorgroup)) %>% 
    mutate(majorgroup = if_else(phylum %in% c("Ctenophora"),"Ctenophore",majorgroup)) %>% 
    mutate(majorgroup = if_else(class %in% c("Scyphozoa"),"Scyphomedusae",majorgroup)) %>%  
    mutate(majorgroup = if_else(order %in% c("Siphonophorae"),"Siphonophore",majorgroup)) %>%  
    mutate(majorgroup = if_else(order %in% c("Narcomedusae","Leptothecata",
                                             "Trachymedusae","Limnomedusae",
                                             "Anthoathecata"),
                                "Hydromedusae",majorgroup)) %>% 
    mutate(majorgroup = if_else(order %in% c("Pteropoda"),"Pteropod",majorgroup)) %>% 
    mutate(majorgroup = if_else(class %in% c("Thaliacea"), "Thaliacean",majorgroup)) %>% 
    mutate(majorgroup = if_else(class %in% c("Appendicularia"),
                                "Appendicularian",majorgroup)) %>% 
    mutate(majorgroup = if_else(class %in% c("Ostracoda"), "Ostracod",majorgroup)) %>% 
    mutate(majorgroup = if_else(class %in% c("Copepoda"),"Non-calanoid",majorgroup)) %>% 
    mutate(majorgroup = if_else(order %in% c("Calanoida"),"Calanoid",majorgroup)) %>% 
    mutate(majorgroup = if_else(order %in% c("Amphipoda"),"Amphipod",majorgroup)) %>% 
    mutate(majorgroup = if_else(order %in% c("Decapoda"),"Decapod",majorgroup)) %>% 
    mutate(majorgroup = if_else(order %in% c("Euphausiacea"),"Euphausiid",majorgroup)) %>%
    mutate(majorgroup = if_else(order %in% c("Mysida","Mysidacea"),
                                "Mysid",majorgroup)) 
  return(taxonomy$majorgroup)
}

getMaxObsNum <- function(traitName, taxonID, df) {
  df <- df %>% 
    filter(traitName == traitName & taxonID == taxonID)
  if (nrow(df) > 0) {
    i <- max(df$observationNumber)
  } else {
    i <- 0
  }
  i
}


# Data imputation functions ----

# This function calculates the mean, standard deviation, and number of 
#   observations for generalization of a given trait for a taxon. This will not 
#   generate an updated catalogNumber and this should created if including the 
#   generalized trait observation in the overall trait table.
getGroupLevelValue <- function(taxon, trait, gen.level = "genus", trait.df, taxonomy.df){
  # check if trait is numeric or binary
  trait.type <- trait.df %>% 
    filter(traitName == trait) %>% 
    filter(row_number()==1) %>% 
    select(valueType)
  
  if (trait.type$valueType %notin% c("numeric","binary")) {
    stop("Error: Please select a trait that is either numeric or binary.")
  }
  
  
  # get taxonomy details of a taxon
  taxon.details <- taxonomy.df %>% 
    filter(scientificName == taxon) %>% 
    distinct(taxonID, .keep_all = TRUE)
  
  if (nrow(taxon.details) == 0) {
    stop("Error: Please select a scientific name that is found in the taxonomy data frame.")
  }
  
  
  # get a list of species that are in the focus taxon's group
  if (gen.level == "genus"){
    taxon.group <- taxonomy.df %>% 
      filter(genus == taxon.details$genus)
  } else if (gen.level == "family") {
    taxon.group <- taxonomy.df %>% 
      filter(family == taxon.details$family)
  } else if (gen.level == "order") {
    taxon.group <- taxonomy.df %>% 
      filter(order == taxon.details$order)
  } else {
    stop("Error: Please select a generalization level that is either genus, family, or order.")
  }
  
  taxon.group <- taxon.group %>%  
    filter(taxonRank %in% c("Species", "Subspecies")) %>% 
    distinct(taxonID, scientificName)
  
  # filter the trait dataframe to select the trait and the species of interest
  trait.sub <- trait.df %>% 
    filter(traitName == trait) %>% 
    filter(taxonID %in% taxon.group$taxonID)
  
  # If there are no other taxa found for this species, return with a blank and provide a warning that there was no generalization.
  if(nrow(trait.sub) == 0){
    warning("Warning: There are no records of the selected trait at the selected group level.")
    return(NA)
  }
  
  # Calculate the mean and standard deviation across all members in a group and return a row of generalized trait value following the standardized format. 
  generalized.trait <- trait.sub %>% 
    mutate(traitValue = as.numeric(traitValue)) %>% 
    mutate(dispersionSD = sd(traitValue, na.rm = TRUE), individualCount = n(),
           traitValue = mean(traitValue), 
           basisOfRecord = "generalized", 
           notes = paste0("Trait value generalized from the ",gen.level," level average.")) %>% 
    distinct(traitID, traitName, traitValue, traitUnit, valueType,
             assocTemperature, dispersionSD, individualCount, 
             notes, basisOfRecord) %>% 
    bind_cols(select(taxon.details, -c(verbatimScientificName)))
  
  return(generalized.trait)
}


# general equation for allometric conversion, base defaults to natural log
conv.allom <- function(W,a,b,base = exp(1)) {
  base^(a + (b*log(W,base)))
}

getpval <- function(x){ 
  p <- pf(x[1],x[2],x[3],lower.tail = F) 
  attributes(p) <- NULL
  return(p)
}


# This is a wrapper for an OLS regression and does not calculate confidence intervals.
getRegressionModel <- function(df, grp = "All", X, Y, base = "10") {
  trait.sub <- df %>% 
    filter(traitName %in% c(all_of(X), all_of(Y)))
  if (grp != "All") {
    trait.sub <- trait.sub %>% 
      filter(str_detect(group, grp)) 
  }
  trait.sub <- trait.sub %>% 
    filter(taxonRank %in% c("Subspecies","Species")) %>% 
    select(taxonID, scientificName, traitName, traitValue, majorgroup) %>% 
    pivot_wider(names_from = traitName, values_from = traitValue) %>% 
    filter(!is.na(get(X)) & !is.na(get(Y))) %>% 
    relocate(X, Y)
  
  if (nrow(trait.sub) > 3) {
    # Calculate OLS
    if(base == "10") {
      reg <- lm(log10(get(Y)) ~ log10(get(X)), data = trait.sub)
    } else if (base == "e") {
      reg <- lm(log(get(Y)) ~ log(get(X)), data = trait.sub)
    } else {
      stop("Error: Please select either base 10 or e.")
    }
    
    reg.results.size <- data.frame(grp = grp, X = X,  Y = Y, 
                                   a = reg$coefficients[1],
                                   b = reg$coefficients[2],
                                   n = nrow(trait.sub),
                                   R2 = summary(reg)$adj.r.squared,
                                   RSE = summary(reg)$sigma,
                                   pval = getpval(summary(reg)$fstatistic), 
                                   model = "OLS", base = as.character(base),
                                   minX = min(trait.sub[,1]), maxX = max(trait.sub[,1]))
  } else {
    reg.results.size <- data.frame(grp = grp, X = X,  Y = Y, 
                                   a = NA, b = NA, n = nrow(trait.sub),
                                   R2 = NA, RSE = NA, pval = NA, 
                                   model = "OLS", base = as.character(base),
                                   minX = NA, maxX = NA)
  }
  
}

# Calculates regression model using the lmodel2() function that allows major axis regressions
# df is a standardized trait dataframe
getRegressionModel2 <- function(df, grp = "All", X, Y, 
                               model = "OLS", base = "10"){
  trait.sub <- df %>% 
    filter(traitName %in% c(all_of(X), all_of(Y)))
  if (grp != "All") {
    trait.sub <- trait.sub %>% 
      filter(str_detect(group, grp)) 
  }
  trait.sub <- trait.sub %>% 
    filter(taxonRank %in% c("Subspecies","Species")) %>% 
    select(taxonID, scientificName, traitName, traitValue, majorgroup) %>% 
    pivot_wider(names_from = traitName, values_from = traitValue) %>% 
    filter(!is.na(get(X)) & !is.na(get(Y))) %>% 
    relocate(X, Y)
  
  # This calculates both the OLS and RMA
  if(base == "10") {
    reg <- lmodel2(log10(get(Y)) ~ log10(get(X)), data = trait.sub)
  } else if (base == "e") {
    reg <- lmodel2(log(get(Y)) ~ log(get(X)), data = trait.sub)
  } else {
    stop("Error: Please select either base 10 or e.")
  }

  if (model == "OLS"){
    ii <- 1
  } else if (model =="RMA") {
    ii <- 3
  } else {
    stop("Error: Please select either OLS or RMA regression model.")
  }
  
  # return(reg)
  
  reg.results.size <- data.frame(grp = grp, X = X,  Y = Y, 
                                 a = reg$regression.results$Intercept[ii], 
                                 b = reg$regression.results$Slope[ii], 
                                 a.ci.2.5 = reg$confidence.intervals$`2.5%-Intercept`[ii],
                                 a.ci.97.5 = reg$confidence.intervals$`97.5%-Intercept`[ii],
                                 b.ci.2.5 = reg$confidence.intervals$`2.5%-Slope`[ii],
                                 b.ci.97.5 = reg$confidence.intervals$`97.5%-Slope`[ii],
                                 n = nrow(trait.sub), 
                                 R2 = reg$rsquare,
                                 # RSE = reg$sigma,
                                 pval = reg$P.param, 
                                 model = model, base = as.character(base),
                                 minX = min(trait.sub[,1]), maxX = max(trait.sub[,1]))
}

# calculate values based on an allometric conversion model
calculateFromModel <- function(df, model, trait.directory, excludeWithLit = TRUE,
                               applyToGeneralized = FALSE,
                               excludeCalculated = TRUE) {
  # make sure this is just one model
  stopifnot(nrow(model) == 1)
  
  # error if the base of the logarithm is not 10 or e (natural log).
  if (!(model$base == "10" | model$base == "e" | model$base == "ln")) {
    stop("Error: Please select model with either base 10 or e.")
  }
  
  df.withdata <- df %>% 
    # filter(str_detect(basisOfRecord, "literature|calculated taxon average| midpoint")) %>% 
    filter(traitName %in% model$Y)
  
  df.calc <- df %>% 
    # Assign the verbatim trait information as the calculated values
    mutate(verbatimTraitName = traitName, verbatimTraitUnit = traitUnit,
           verbatimTraitValue = as.character(traitValue), verbatimNotes = notes) %>% 
    select(-c(primaryReferenceDOI, secondaryReferenceDOI, catalogNumber, 
              sizeAssocName, sizeAssocUnit, sizeAssocValue,
              sizeAssocN, sizeAssocSD, assocTemperature,
              sizeAssocReference, verbatimLocality, 
              decimalLongitude, decimalLatitude, notes,
              aggregateMeasure, isDerived, catalogSource, basisOfRecordDescription)) %>%  
    filter(str_detect(group, model$grp),
           traitName %in% model$X) %>% 
    mutate(traitValue = as.numeric(traitValue)) %>% 
    filter(traitValue >= model$minX,
           traitValue <= model$maxX) 
  
  # If the predictor is size trait, assign it as sizeAssoc
  if (model$X %in% c("bodyLengthMax", "carbonWeight", "dryWeight", "wetWeight")){
    df.calc <- df.calc %>%
      mutate(sizeAssocName = traitName,
             sizeAssocUnit = traitUnit, sizeAssocValue = traitValue,
             sizeAssocSD = dispersionSD, sizeAssocN = individualCount,
             sizeAssocReference = if_else(!is.na(secondaryReference),
                                          secondaryReference, primaryReference))
  }
  # If y variable is a rate, assign temperature to default of 15C, if this is
  #   not the case, need to revise externally.
  if (grepl("_15C",model$Y)){
    df.calc <- df.calc %>% 
      mutate(assocTemperature = 15)
  }
  
  if(excludeCalculated == TRUE) {
    df.calc <- df.calc %>% 
      filter(!(basisOfRecord == "calculated from regression"))
  }
  df.calc <- df.calc %>% 
    dplyr::select(-c(primaryReference, secondaryReference, dispersionSD, individualCount))
  
  if(applyToGeneralized == TRUE) {
    df.calc <- df.calc 
      # apply only genus level generalization if for level 3
      # filter(!(basisOfRecord == "generalized" & valueRank != "Genus"))
  } else {
    df.calc <- df.calc %>% 
      filter(!(basisOfRecord == "generalized"))
  }
  
  # calculate for traits - note that traits here are log-transformed to a 
  #  specific base. None of the models involve percent or ratio traits, and 
  #  those would need a separate transformation if ever.
  if (model$base == "10") {
    df.calc <- df.calc %>% 
      mutate(traitName = model$Y,
             traitValue =  conv.allom(traitValue, model$a, model$b, base = 10))
  } else if (model$base == "e" | model$base == "ln") {
    mutate(traitName = model$Y,
           traitValue =  conv.allom(traitValue, model$a, model$b))
  }
  
  df.calc <- df.calc %>% 
    mutate(basisOfRecord = "calculated from regression",
           notes = paste0("Value calculated from ",model$X, " using the equation: ",
                          "y = ",model$base,"^(",sprintf("%.3f",model$a),
                          "+(",sprintf("%.3f",model$b),
                          "*log(x,",model$base,")))."),
           uploadBy = "P.Pata", 
           uploadDate = as.character(ymd(Sys.Date()))) %>% 
    # TODO revise the next three lines below
    # group_by(taxonID) %>% 
    # mutate(maxObsNum = getMaxObsNum(traitName, taxonID, df)) %>% 
    # ungroup() %>% 
    # If rate, note that assocTemperature was standardized to 15C
    mutate(assocTemperature = if_else(str_detect(traitName,"Rate"),
                                           15, -999),
           assocTemperature = na_if(assocTemperature, -999),
           assocTemperature = as.character(assocTemperature),
           verbatimTemperature = "NA")
  
  # exclude taxa which already have literature/derived trait information
  if(excludeWithLit == TRUE){
    df.calc <- df.calc %>% 
      filter(taxonID %notin% df.withdata$taxonID)
  }
  
  df.calc <- df.calc %>% 
    # update units and traitTaxonIDs
    standardizeID(trait.directory) %>% 
    standardizeUnit(trait.directory) %>% 
    group_by(traitID,taxonID) %>% 
    # TODO revise the next three lines below
    mutate(observationNumber = maxObsNum + row_number()) %>% 
    mutate(maxObsNum = max(observationNumber)) %>% 
    mutate(catalogNumber = paste0(traitID,"-",taxonID,"-",observationNumber)) %>% 
    mutate(isDerived = TRUE) %>% 
    ungroup()
}


# function for plot
plotAllometric <- function(df, grp, X, Y, base = "10") {
  trait.sub <- df %>% 
    filter(traitName %in% c(X,Y)) %>% 
    filter(str_detect(group, grp)) %>% 
    filter(taxonRank %in% c("Subspecies","Species")) %>% 
    select(taxonID, scientificName, traitName, traitValue, majorgroup) %>% 
    pivot_wider(names_from = traitName, values_from = traitValue) %>% 
    filter(!is.na(get(X)) & !is.na(get(Y))) 
  
  if (base == "10") {
    base.trans <-  "log10"
  } else if(base == "e") {
    base.trans <-  "log"
  } else {
    stop("Error: Please select the axis scales to either be log10 or ln scaled.")
  }
  
  ggplot(trait.sub, aes(x = get(X), y = get(Y))) +
    geom_point(aes(color = majorgroup, text = scientificName)) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_x_continuous(labels = scaleFUN, trans = base.trans) +
    scale_y_continuous(labels = scaleFUN, trans = base.trans) +
    xlab(X) + ylab(Y) +
    theme_bw() + 
    stat_regline_equation(label.y = 1) +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             label.y = 1.5) 
}

# Plot regression line with confidence intervals - this only plots the function 
#  but not the data points.
plotRegModel <- function(model){
  model <- model[1,]
  
  if (model$base == "10") {
    base.trans <-  "log10"
    base.allom <- 10
  } else if(model$base == "e") {
    base.trans <-  "log"
    base.allom <- exp(1)
  } else {
    stop("Error: Please select the axis scales to either be log10 or ln scaled.")
  }
  
  df <- data.frame( x = seq(model$minX, model$maxX, len = 100) ) %>% 
    mutate(y = conv.allom(x, model$a, model$b, base = base.allom),
           y.lower = conv.allom(x, model$a.ci.2.5, model$b.ci.2.5, base = base.allom),
           y.upper = conv.allom(x, model$a.ci.97.5, model$b.ci.97.5, base = base.allom))
  
  g <- ggplot(df, aes(x, y)) +
    geom_line() + 
    geom_line(aes(y = y.lower), lty = 2, color = "gray") +
    geom_line(aes(y = y.upper), lty = 2, color = "gray") +
    theme_bw() +
    xlab(model$X) + ylab(model$Y) +
    scale_x_continuous(labels = scaleFUN, trans = base.trans) +
    scale_y_continuous(labels = scaleFUN, trans = base.trans)
  
  return(g)
}

# Updated February 13, 2023 - uses the associated size value
calculate.WSRates <- function(traits.calculated, trait.X, trait.Y) {
  # Prepare for updates in traitName and traitUnit
  if (trait.X == "carbonWeight") {
    name.suffix <- "_WSC_15C"
    unit.suffix <- "mg C^-1 h"
  } else if (trait.X == "dryWeight") {
    name.suffix <- "_WSDW_15C"
    unit.suffix <- "mg^-1 h"
  } else if (trait.X == "wetWeight") {
    name.suffix <- "_WSWW_15C"
    unit.suffix <- "mg^-1 h"
  } else {
    stop("Error: Please select either carbonWeight, dryWeight, or wetWeight.")
  }
  
  traits.calculated <- traits.calculated %>% 
    filter(traitName %in% trait.Y) %>% 
    # Assign the verbatim trait information as the calculated values
    mutate(verbatimTraitName = traitName, verbatimTraitUnit = traitUnit,
           verbatimTraitValue = traitValue, verbatimNotes = notes,
           traitValue = as.numeric(traitValue),
           sizeAssocValue = as.numeric(sizeAssocValue)) %>% 
    group_by(catalogNumber) %>% 
    filter(!is.na(traitValue) & !is.na(sizeAssocValue)) %>% 
    mutate(traitValue = traitValue / sizeAssocValue,
           traitName = str_replace(traitName,"_15C", name.suffix),
           traitUnit = str_replace(traitUnit,"h", unit.suffix),
           notes = "Calculated from per individual rate and associated size value.") %>% 
    standardizeID(trait.directory) %>% 
    mutate(traitValue = as.character(traitValue))
}

# Generalize categorical traits
# This function generalizes the binary version of a categorical trait for missing
#  species based on records of categorical traits. The results are proportions 
#  between 0 to 1.0 representing how likely a category is present in a taxonomic 
#  group. The number of species this was based on is stored in traitValues.
generalizeCategoricalTrait <- function(df, binaryPrefix, taxonomy) {
  cat.traits <- df %>% 
    filter(valueType == "binary") %>% 
    mutate(traitValue = as.numeric(traitValue))
  
  # Per trait, per taxa (e.g. genus), get proportions of a categorical trait
  genus.binary <- cat.traits %>% 
    filter(grepl(binaryPrefix,traitName)) %>% 
    # group at a taxonomic level
    group_by(genus, traitName) %>% 
    summarise(traitValue = sum(traitValue), nspecies = n(), .groups = "drop") %>% 
    group_by(genus) %>% 
    mutate(total.obs = sum(traitValue)) %>% 
    ungroup() %>% 
    # convert to proportions
    mutate(traitValue = traitValue / total.obs)
  
  # generalize to missing species levels?
  missing.cat <- taxonomy %>% 
    filter(taxonRank %in% c("Species","Subspecies")) %>% 
    # Only calculate for species which there observed trait values
    filter(genus %in% filter(df, grepl(binaryPrefix,traitName))$genus) %>% 
    # Exclude species which we already have data for
    filter(taxonID %notin% filter(df, grepl(binaryPrefix,traitName))$taxonID) 
  
  # Generalize the categorical traits to missing species
  trait.generalized <- missing.cat %>% 
    distinct(taxonID, scientificName, genus) %>% 
    inner_join(genus.binary, by = "genus", multiple = "all") %>% 
    mutate(notes = paste0("generalized at genus level based on ", nspecies, 
                          " species and ", total.obs, " records"),
           individualCount = total.obs) %>% 
    select(-c(nspecies, total.obs)) %>% 
    left_join(select(taxonomy, taxonRank, acceptedNameUsageID, 
                     acceptedNameUsage, kingdom, phylum, 
                     class, order, family, genus, species, majorgroup, taxonID), 
              by = c("taxonID", "genus")) %>% 
    mutate(traitID = NA,traitUnit = NA,
           aggregateMeasure = FALSE, isDerived = TRUE, valueType = "binary",
           basisOfRecord = "generalized") %>% 
    standardizeID(trait.directory) %>% 
    mutate(traitValue = as.character(traitValue))
}


# Functions for trait correlations ----

# Functions for colpair_map()
calc_p_value <- function(vec_a, vec_b) {
  t.test(vec_a, vec_b)$p.value # The t.test here is not the same as cor.test
  # cor.test(vec_a, vec_b)$p.value
}

get_pairwise.N <- function(vec_a, vec_b) {
  sum(!is.na(vec_a) & !is.na(vec_b))
}

# Widens a table for listed traits. Column names do not change and SD is excluded.
widen_traits <- function(trait.table, trait.list) {
  trait.table %>% 
    filter(traitName %in% trait.list) %>% 
    dplyr::select(taxonID, scientificName, majorgroup, acceptedNameUsageID, 
                  traitName, traitValue) %>% 
    pivot_wider(names_from = traitName, values_from = traitValue) %>% 
    relocate(all_of(trait.list)) 
}

# This function returns a wide table of 2 selected traits for analysis. The SE
#   is calculated here which will be useful for analysis that need the
#   measurement error as input. This loses information on trait ID and units.
widen_traits_SE <- function(trait.table, trait.list) {
  
  val <- trait.table %>%
    filter(traitName %in% trait.list) %>%
    dplyr::select(taxonID, scientificName, majorgroup, acceptedNameUsageID,
                  traitName, traitValue) %>%
    pivot_wider(names_from = traitName, values_from = traitValue) %>% 
    relocate(all_of(trait.list)) %>% 
    # rename to trait.x and trait.y
    rename(trait.x = all_of(trait.list[1]), trait.y = all_of(trait.list[2]))
  
  sd <- trait.table %>%
    filter(traitName %in% trait.list) %>%
    dplyr::select(taxonID, scientificName, majorgroup, acceptedNameUsageID,
                  traitName, dispersionSD) %>%
    pivot_wider(names_from = traitName, values_from = dispersionSD) %>% 
    relocate(all_of(trait.list)) %>% 
    # rename to sd.x and sd.y
    rename(sd.x = all_of(trait.list[1]), sd.y = all_of(trait.list[2]))
  
  
  n <- trait.table %>%
    filter(traitName %in% trait.list) %>%
    dplyr::select(taxonID, scientificName, majorgroup, acceptedNameUsageID,
                  traitName, individualCount) %>%
    pivot_wider(names_from = traitName, values_from = individualCount) %>% 
    relocate(all_of(trait.list)) %>% 
    # rename to n.x and n.y
    rename(n.x = all_of(trait.list[1]), n.y = all_of(trait.list[2])) %>% 
    # if N is na, N = 1
    mutate(n.x = if_else(is.na(n.x), 1, n.x),
           n.y = if_else(is.na(n.y), 1, n.y))
  
  wide.table <- left_join(val, sd, by = c("taxonID", "scientificName", 
                                          "majorgroup", "acceptedNameUsageID")) %>% 
    left_join(n, by = c("taxonID", "scientificName", 
                        "majorgroup", "acceptedNameUsageID")) %>% 
    # This section assigns the SE for instances when SE is NaN. Here, we use 
    #   the mean of available SEs rather than assuming no error.
    # Consider SD = 0 as erroneous, which can happen when N=2 or values are the same
    mutate(sd.x = if_else(sd.x == 0, NaN, sd.x),
           sd.y = if_else(sd.y == 0, NaN, sd.y)) %>% 
    mutate(se.x = sd.x / sqrt(n.x),
           se.y = sd.y / sqrt(n.y)) %>% 
    # Calculate mean SE across all species for a particular trait and assign 
    #  this to rows with missing SE
    mutate(mean.se.x = mean(se.x, na.rm = TRUE),
           mean.se.y = mean(se.y, na.rm = TRUE),
           se.x = if_else(is.na(se.x), mean.se.x, se.x),
           se.y = if_else(is.na(se.y), mean.se.y, se.y)) %>% 
    dplyr::select(-c(mean.se.x, mean.se.y))
  
  return(wide.table)
}


correlate_traits <- function(trait.table, trait.list,
                             logtransX = FALSE, logtransY = FALSE) {
  A <- widen_traits(trait.table, trait.list) %>% 
    # filter(!is.na(get(trait.list[1])) & !is.na(get(trait.list[2]))) %>% 
    dplyr::select(all_of(trait.list)) 
  
  # If log10 transform the X or Y variable
  if (logtransX == TRUE) {
    A[,1] = log10(A[,1])
  }
  if (logtransY == TRUE) {
    A[,2] = log10(A[,2])
  }
  
  # with the corrr::correlate() function, multiple traits can be analyzed with missing points. This returns a correlation matrix which can be stretched into a long data frame.
  # Correlate, get p-value, and N
  B <- corrr::correlate(A, use = "pairwise.complete.obs", 
                        method = "pearson", quiet = TRUE) %>% 
    corrr::stretch(na.rm = TRUE, remove.dups = TRUE) 
  
  C <- corrr::colpair_map(A, get_pairwise.N) %>% 
    corrr::stretch(na.rm = TRUE, remove.dups = TRUE) %>% 
    rename(N = r)
  
  D <- corrr::colpair_map(A, calc_p_value) %>%  
    corrr::stretch(na.rm = TRUE, remove.dups = TRUE) %>% 
    rename(pval = r)
  
  E <- left_join(B,C, by = c("x","y")) %>% 
    left_join(D, by = c("x","y")) %>% 
    arrange(-r)
  
  return(E)
}

# This runs the cor.phylo test with or without a phylogenetic signal and returns
#   a long table of correlation results.
corphylo.tests <- function(trait.num.phylo, trait.list, phylo.tree, 
                           logtransX = FALSE, logtransY = FALSE) {
  # *** Prepare trait and phylo tree subsets ***
  # Filter and organize the data set
  trait.sub <- widen_traits_SE(trait.num.phylo, trait.list)
  # *** Log transform the trait values ****
  # This is often applied to the trait dataset since most do not have a normal
  #   distribution. For the current trait dataset, the exception are
  #   traitName == "energyDensityVolume", ratios, maybe composition 
  if (logtransX == TRUE) {
    trait.sub$trait.x <- log10(trait.sub$trait.x)
  }
  if (logtransY == TRUE) {
    trait.sub$trait.y <- log10(trait.sub$trait.y)
  }
  
  trait.sub <- trait.sub %>% 
    # Exclude pairs with missing information. Note that because the exclusion 
    #   happens after calling the widen_traits_SE() function, the average SE 
    #   used to pad the NaN values are based on the entire available data for a trait.
    filter(!is.na(trait.x) & !is.na(trait.y)) %>% 
    # add ott_id
    left_join(dplyr::select(phylo.table, ott_id, taxonID, color), by = "taxonID") %>% 
    mutate(ott_id = as.character(ott_id))
  
  # if sample size is too small, exit
  if(nrow(trait.sub) < 10) {
    results <- NA
    warning(paste("Not calculated for",trait.list[1],"and",trait.list[2],"because sample size <10."))
    return(results)
  }
  
  # Prune the phylo tree to only include species in the trait subset
  phylo.sub <- phylo.tree %>% 
    ape::keep.tip(trait.sub$ott_id)
  # recalculate distances
  phylo.sub <- ape::compute.brlen(phylo.sub, method = 'Grafen')
  
  # Create a no phylogenetic effect tree
  star <- ape::stree(length(phylo.sub$tip.label))
  star$edge.length <- array(1, dim = c(length(phylo.sub$tip.label), 1))
  star$tip.label <- phylo.sub$tip.label
  
  # Reorder trait.sub to match phylo.sub tips. Useful for plotting
  trait.sub <- data.frame(ott_id = phylo.sub$tip.label) %>% 
    left_join(trait.sub, by = "ott_id")
  
  # *** Run correlations ***
  # 1. cor_phylo without phylo and ME
  cor.nophylo.nome <- cor_phylo(variates = ~ trait.x + trait.y,
                                species = ~ ott_id,
                                phy = star,
                                data = trait.sub) # boot = 100
  
  # 3. cor_phylo with phylo but no ME
  cor.phylo.nome <- cor_phylo(variates = ~ trait.x + trait.y,
                              species = ~ ott_id,
                              phy = phylo.sub,
                              data = trait.sub)
  
  # *** Summarize results into table ***
  model.type <- c("no.phylo.no.me","phylo.no.me")
  # correlation coefficients
  res.r <- data.frame(trait1 = trait.list[1], trait2 = trait.list[2],
                      N.pairs = nrow(trait.sub))
  res.r <- bind_cols(res.r,
                     data.frame(model.type = model.type),
                     data.frame(cor.phylo = rbind(cor.nophylo.nome$corrs[1,2],
                                                  cor.phylo.nome$corrs[1,2])))
  
  # phylogenetic effect from OU (Ornstein-Uhlenbeck) process
  res.d <- cbind(model.type,t(cbind(cor.nophylo.nome$d, cor.phylo.nome$d)) ) %>% 
    as.data.frame() %>% 
    pivot_longer(cols = c("trait.x","trait.y"), 
                 names_to = "trait", values_to = "d")
  # regression coefficients
  res.B <- cbind(model.type,t(cbind(cor.nophylo.nome$B[,1], cor.phylo.nome$B[,1])) ) %>% 
    as.data.frame() %>% 
    rename("trait.x" = "trait.x_0", "trait.y" = "trait.y_0" )%>% 
    pivot_longer(cols = c("trait.x","trait.y"), 
                 names_to = "trait", values_to = "B") 
  res.pval <- cbind(model.type,t(cbind(cor.nophylo.nome$B[,4], cor.phylo.nome$B[,4])) ) %>% 
    as.data.frame() %>% 
    rename("trait.x" = "trait.x_0", "trait.y" = "trait.y_0" )%>% 
    pivot_longer(cols = c("trait.x","trait.y"), 
                 names_to = "trait", values_to = "pval")
  
  # Get metrics to compare models: loglike (higher better), AIC and BIC (lower better).
  #  note that adding more predictor variables would increase the log-likelihood
  #  value even if the predictor is not statistically significant so caution on 
  #  comparing models. Complex models are also preferred in AIC. BIC penalizes 
  #  complex models, so this prefers simpler ones.
  res.loglik <- data.frame(log.lik = 
                             rbind(cor.nophylo.nome$logLik, cor.phylo.nome$logLik) )
  res.AIC <- data.frame(AIC = 
                          rbind(cor.nophylo.nome$AIC, cor.phylo.nome$AIC))
  res.BIC <- data.frame(BIC = 
                          rbind(cor.nophylo.nome$BIC, cor.phylo.nome$BIC))
  
  # Add general taxonomic groupings (list of groups and % calanoid species)
  mg.list <- trait.sub %>% 
    group_by(majorgroup) %>% 
    summarise(Nspecies = n()) %>% 
    mutate(Ntotal = sum(Nspecies) )
  perc.calanoid <- filter(mg.list, majorgroup == "Calanoid") %>% 
    summarise(perc.cal = Nspecies/Ntotal*100)
  if (nrow(perc.calanoid) == 0) {
    perc.calanoid <- data.frame(perc.cal = 0)
  }
  mg.list <- mg.list %>% 
    arrange(-Nspecies) %>% 
    mutate(list = paste0(majorgroup," (",Nspecies,")")) %>% 
    summarise(list = paste(list, collapse = ", "))
  
  # Wrap the results into a wide table that could have iterations every row
  results.r <- bind_cols(res.r, res.loglik, res.AIC, res.BIC,
                         perc.calanoid, mg.list) %>% 
    mutate(analysis = paste0(trait.list[1],"-",trait.list[2])) %>% 
    relocate(analysis)
  results.phylo = left_join(res.d, res.B, by = c("model.type","trait")) %>% 
    left_join(res.pval, by = c("model.type","trait")) %>% 
    mutate(trait = str_replace(trait,"trait.x",trait.list[1]),
           trait = str_replace(trait,"trait.y",trait.list[2])) %>% 
    mutate(analysis = paste0(trait.list[1],"-",trait.list[2])) %>% 
    relocate(analysis)
  
  results <- left_join(results.r, results.phylo, 
                       by = c("analysis","model.type")) %>% 
    rename(trait.coef = trait)
  
  rm(results.r, results.phylo, res.r, res.d, res.B, res.pval,
     res.loglik, res.AIC, res.BIC, mg.list, perc.calanoid)
  return(results)
}


# This runs different correlation tests using the phyr::cor_phylo() function.
#   Note that this temporarily requires the phylo.tree and vcv objects to 
#   explore the effect of which phylo signal to use. This returns a long table
#   or correlation results. 
# This version tries to explore ME and a subset phylo vcv rather than a pruned tree.
# Last updated: Jan 23, 2023.
corphylo.tests.extra <- function(trait.num.phylo, trait.list, phylo.tree, 
                                 zoop_vcv, 
                                 logtransX = FALSE, logtransY = FALSE) {
  # *** Prepare trait and phylo tree subsets ***
  # Filter and organize the data set
  trait.sub <- widen_traits_SE(trait.num.phylo, trait.list)
  # *** Log transform the trait values and MEs ****
  # This is often applied to the trait dataset since most do not have a normal
  #   distribution. For the current trait dataset, the exception are
  #   traitName == "energyDensityVolume", ratios, maybe composition (?)
  # If log10 transform the X or Y variable, also transform the sd. 
  #   *Note that is is erroneous when sd<0 because a log transformation would 
  #    be negative and the magnitude of the error is larger.
  if (logtransX == TRUE) {
    trait.sub$trait.x <- log10(trait.sub$trait.x)
    # recalculate the sd (this is not the best estimator)
    # TODO improve estimate for sd of log-transformed means
    trait.sub$sd.x <- log10(trait.sub$sd.x)
  }
  if (logtransY == TRUE) {
    trait.sub$trait.y <- log10(trait.sub$trait.y)
    # recalculate the sd
    trait.sub$sd.y <- log10(trait.sub$sd.y)
  }
  
  trait.sub <- trait.sub %>% 
    # Recalculate the SE . If the SE is NaN/0, use the mean SE of the trait
    mutate(se.x = sd.x / sqrt(n.x),
           se.y = sd.y / sqrt(n.y)) %>% 
    mutate(mean.se.x = mean(se.x, na.rm = TRUE),
           mean.se.y = mean(se.y, na.rm = TRUE),
           se.x = if_else(is.na(se.x), mean.se.x, se.x),
           se.y = if_else(is.na(se.y), mean.se.y, se.y)) %>% 
    dplyr::select(-c(mean.se.x, mean.se.y)) %>% 
    # Exclude pairs with missing information. Note that because the exclusion 
    #   happens after calling the widen_traits_SE() function, the average SE 
    #   used to pad the NaN values are based on the entire available data for a trait.
    filter(!is.na(trait.x) & !is.na(trait.y)) %>% 
    
    # add ott_id
    left_join(dplyr::select(phylo.table, ott_id, taxonID, color), by = "taxonID") %>% 
    mutate(ott_id = as.character(ott_id))
  
  # Prune the phylo tree to only include species in the trait subset
  phylo.sub <- phylo.tree %>% 
    ape::keep.tip(trait.sub$ott_id)
  # Create a no phylogenetic effect tree
  star <- ape::stree(length(phylo.sub$tip.label))
  star$edge.length <- array(1, dim = c(length(phylo.sub$tip.label), 1))
  star$tip.label <- phylo.sub$tip.label
  
  # Reorder trait.sub to match phylo.sub tips. Useful for plotting
  trait.sub <- data.frame(ott_id = phylo.sub$tip.label) %>% 
    left_join(trait.sub, by = "ott_id")
  
  # VCV.sub
  ii <- which(row.names(zoop_vcv) %in% phylo.sub$tip.label)
  zoop.vcv.sub <- zoop_vcv[ii,ii]
  
  # *** Run correlations ***
  # 1. Regular correlation
  cor.reg <- correlate_traits(trait.num.phylo, trait.list,
                              logtransX, logtransY)
  
  # 2. cor_phylo without phylo and ME
  cor.nophylo.nome <- cor_phylo(variates = ~ trait.x + trait.y,
                                species = ~ ott_id,
                                phy = star,
                                data = trait.sub) # boot = 100
  
  # 3. cor_phylo with phylo but no ME
  cor.phylo.nome <- cor_phylo(variates = ~ trait.x + trait.y,
                              species = ~ ott_id,
                              phy = phylo.sub,
                              data = trait.sub) 
  
  # 4. cor_phylo with phylo and ME
  cor.phylo.me <- cor_phylo(variates = ~ trait.x + trait.y,
                            species = ~ ott_id,
                            phy = phylo.sub,
                            meas_errors = list(trait.x ~ se.x, 
                                               trait.y ~ se.y),
                            data = trait.sub) 
  # 5. cor_phylo with phylo and no ME but phylo is subset of the vcv
  # Seems like very small difference but test anyway
  cor.phylo.nome.vcv <- cor_phylo(variates = ~ trait.x + trait.y,
                                  species = ~ ott_id,
                                  phy = zoop.vcv.sub,
                                  data = trait.sub) 
  
  # *** Summarize results into table ***
  model.type <- c("no.phylo.no.me","phylo.no.me","phylo.me","phylo.no.me.vcv")
  # correlation coefficients
  res.r <- data.frame(trait1 = trait.list[1], trait2 = trait.list[2],
                      N.pairs = nrow(trait.sub), reg.cor = cor.reg$r)
  res.r <- bind_cols(res.r,
                     data.frame(model.type = model.type),
                     data.frame(cor.phylo = rbind(cor.nophylo.nome$corrs[1,2],
                                                  cor.phylo.nome$corrs[1,2],
                                                  cor.phylo.me$corrs[1,2],
                                                  cor.phylo.nome.vcv$corrs[1,2])))
  
  # phylogenetic effect from OU (Ornstein-Uhlenbeck) process
  res.d <- cbind(model.type,t(cbind(cor.nophylo.nome$d, cor.phylo.nome$d, 
                                    cor.phylo.me$d, cor.phylo.nome.vcv$d)) ) %>% 
    as.data.frame() %>% 
    pivot_longer(cols = c("trait.x","trait.y"), 
                 names_to = "trait", values_to = "d")
  # regression coefficients
  res.B <- cbind(model.type,t(cbind(cor.nophylo.nome$B[,1], cor.phylo.nome$B[,1], 
                                    cor.phylo.me$B[,1], cor.phylo.nome.vcv$B[,1])) ) %>% 
    as.data.frame() %>% 
    rename("trait.x" = "trait.x_0", "trait.y" = "trait.y_0" )%>% 
    pivot_longer(cols = c("trait.x","trait.y"), 
                 names_to = "trait", values_to = "B") 
  res.pval <- cbind(model.type,t(cbind(cor.nophylo.nome$B[,4], cor.phylo.nome$B[,4], 
                                       cor.phylo.me$B[,4], cor.phylo.nome.vcv$B[,4])) ) %>% 
    as.data.frame() %>% 
    rename("trait.x" = "trait.x_0", "trait.y" = "trait.y_0" )%>% 
    pivot_longer(cols = c("trait.x","trait.y"), 
                 names_to = "trait", values_to = "pval")
  
  # Get metrics to compare models: loglike (higher better), AIC and BIC (lower better).
  #  note that adding more predictor variables would increase the log-likelihood
  #  value even if the predictor is not statistically significant so caution on 
  #  comparing models. Complex models are also preferred in AIC. BIC penalizes 
  #  complex models, so this prefers simpler ones.
  res.loglik <- data.frame(log.lik = 
                             rbind(cor.nophylo.nome$logLik, cor.phylo.nome$logLik,
                                   cor.phylo.me$logLik, cor.phylo.nome.vcv$logLik) )
  res.AIC <- data.frame(AIC = 
                          rbind(cor.nophylo.nome$AIC, cor.phylo.nome$AIC,
                                cor.phylo.me$AIC, cor.phylo.nome.vcv$AIC))
  res.BIC <- data.frame(BIC = 
                          rbind(cor.nophylo.nome$BIC, cor.phylo.nome$BIC,
                                cor.phylo.me$BIC, cor.phylo.nome.vcv$BIC))
  
  # Add general taxonomic groupings (list of groups and % calanoid species)
  mg.list <- trait.sub %>% 
    group_by(majorgroup) %>% 
    summarise(Nspecies = n()) %>% 
    mutate(Ntotal = sum(Nspecies) )
  perc.calanoid <- filter(mg.list, majorgroup == "Calanoid") %>% 
    summarise(perc.cal = Nspecies/Ntotal*100)
  mg.list <- mg.list %>% 
    arrange(-Nspecies) %>% 
    mutate(list = paste0(majorgroup," (",Nspecies,")")) %>% 
    summarise(list = paste(list, collapse = ", "))
  
  # Wrap the results into a wide table that could have iterations every row
  results.r <- bind_cols(res.r, res.loglik, res.AIC, res.BIC,
                         perc.calanoid, mg.list) %>% 
    mutate(analysis = paste0(trait.list[1],"-",trait.list[2])) %>% 
    relocate(analysis)
  results.phylo = left_join(res.d, res.B, by = c("model.type","trait")) %>% 
    left_join(res.pval, by = c("model.type","trait")) %>% 
    mutate(trait = str_replace(trait,"trait.x",trait.list[1]),
           trait = str_replace(trait,"trait.y",trait.list[2])) %>% 
    mutate(analysis = paste0(trait.list[1],"-",trait.list[2])) %>% 
    relocate(analysis)
  
  results <- left_join(results.r, results.phylo, 
                       by = c("analysis","model.type")) %>% 
    rename(trait.coef = trait)
  
  rm(results.r, results.phylo, res.r, res.d, res.B, res.pval,
     res.loglik, res.AIC, res.BIC, mg.list, perc.calanoid)
  return(results)
  
  # # If extracting p-vals wit B coef
  # cor.type <- c("no.phylo.no.me", "no.phylo.no.me","phylo.no.me",
  #               "phylo.no.me","phylo.me","phylo.me","phylo.nome.vcv")
  # res.B <- cbind(cor.type, trait = trait.list,
  #                rbind(cor.nophylo.nome$B, cor.phylo.nome$B, cor.phylo.me$B,
  #                      cor.phylo.nome.vcv$B))
  # rownames(res.B) <- c(1:8)
  
}


# Wrapper function for phylolm using the OU model. Calculates the regression
#  with and without the phylogenetic signal for the specified pair of variables.
#  The first variable in trait.list is the predictor variables.
phylolm.test <- function(data, trait.list, phylo.tree, base = "raw"){
  # log transform the data when needed
  if(base == "10") {
    data <- data %>% 
      mutate(traitValue = log10(traitValue))
  } else if(base == "e" | base == "ln"){
    data <- data %>% 
      mutate(traitValue = log(traitValue))
  }
  
  # Subset the trait data
  data <- widen_traits_SE(data, trait.list) %>% 
    filter(!is.na(trait.x) & !is.na(trait.y)) %>%
    dplyr::select(trait.x, trait.y, taxonID, scientificName, majorgroup) %>% 
    # add ott_id
    left_join(dplyr::select(phylo.table, ott_id, taxonID, color), by = "taxonID") %>% 
    mutate(ott_id = as.character(ott_id))
  
  # Sample size must at least be 3 species
  if(nrow(data) >= 3) {
    
    # Add general taxonomic groupings (list of groups and % calanoid species)
    mg.list <- data %>% 
      group_by(majorgroup) %>% 
      summarise(Nspecies = n()) %>% 
      mutate(Ntotal = sum(Nspecies) )
    perc.calanoid <- filter(mg.list, majorgroup == "Calanoid") %>% 
      summarise(perc.cal = Nspecies/Ntotal*100)
    if(nrow(perc.calanoid) == 0) { # especially for gelatinous only groups
      perc.calanoid <- data.frame(perc.cal = 0)
    }
    mg.list <- mg.list %>% 
      arrange(-Nspecies) %>% 
      mutate(list = paste0(majorgroup," (",Nspecies,")")) %>% 
      summarise(list = paste(list, collapse = ", "))
    
    # Prune tree and recalculate branch distances
    # species composition of the subset
    phylo.sub <- phylo.tree %>% 
      ape::keep.tip(unique(data$ott_id))
    phylo.sub <- ape::compute.brlen(phylo.sub, method = 'Grafen')
    phylo.sub.table <- data.frame(ott_id = as.numeric(phylo.sub$tip.label)) %>% 
      left_join(dplyr::select(phylo.table, ott_id, scientificName, taxonID, color), 
                by = "ott_id")
    # Create a no phylogenetic effect tree
    star <- ape::stree(length(phylo.sub$tip.label))
    star$edge.length <- array(1, dim = c(length(phylo.sub$tip.label), 1))
    star$tip.label <- phylo.sub$tip.label
    # Reorder data to match phylo.sub tips. Useful for plotting
    data <- data.frame(ott_id = phylo.sub$tip.label) %>% 
      left_join(data, by = "ott_id") %>% 
      column_to_rownames("ott_id")
    
    # Calculate regressions
    # 1. phylolm with phylo effects
    reg1 <- phylolm(trait.y ~ trait.x, phy=phylo.sub, 
                    model = "OUrandomRoot", data=data, 
                    lower.bound = 0, upper.bound = 10000,
                    boot = 100)
    # 2. phylolm with no phylo effects = a regular OLS regression
    reg2 <- phylolm(trait.y ~ trait.x, phy=star, 
                    model = "OUfixedRoot", data=data,
                    boot = 100)
    
    # Organize results into a table
    model.type <- c("reg.phylo","reg.no.phylo")
    results <- data.frame(model.type = model.type, # add test.num
                          trait1 = trait.list[1], trait2 = trait.list[2],
                          N.pairs = nrow(data),
                          perc.calanoid = perc.calanoid,
                          r.squared = c(reg1$r.squared, reg2$r.squared),
                          adj.r.squared = c(reg1$adj.r.squared, reg2$adj.r.squared),
                          aic = c(reg1$aic, reg2$aic),
                          loglik = c(reg1$logLik, reg2$logLik),
                          sigma = c(reg1$sigma2, reg2$sigma2),
                          alpha = c(reg1$optpar, reg2$optpar),
                          intercept = c(reg1$coefficients[1], reg2$coefficients[1]),
                          slope = c(reg1$coefficients[2], reg2$coefficients[2]),
                          intercept.se = c(summary(reg1)$coefficients[1,2],
                                           summary(reg2)$coefficients[1,2]),
                          slope.se = c(summary(reg1)$coefficients[2,2],
                                       summary(reg2)$coefficients[2,2]),
                          
                          
                          intercept.ci.lo = c(summary(reg1)$coefficients[1,4],
                                              summary(reg2)$coefficients[1,4]),
                          slope.ci.lo = c(summary(reg1)$coefficients[2,4],
                                          summary(reg2)$coefficients[2,4]),
                          intercept.ci.up = c(summary(reg1)$coefficients[1,5],
                                              summary(reg2)$coefficients[1,5]),
                          slope.ci.up = c(summary(reg1)$coefficients[2,5],
                                          summary(reg2)$coefficients[2,5]),
                          
                          intercept.pval = c(summary(reg1)$coefficients[1,6],
                                             summary(reg2)$coefficients[1,6]),
                          slope.pval = c(summary(reg1)$coefficients[2,6],
                                         summary(reg2)$coefficients[2,6]),
                          # Add min and max of predictors
                          minX = 10^min(data$trait.x), maxX = 10^max(data$trait.x),
                          base = "10",
                          mg.list = mg.list)
  } else {
    results <- data.frame(model.type = c("reg.phylo","reg.no.phylo"), # add test.num
                          trait1 = trait.list[1], trait2 = trait.list[2],
                          N.pairs = nrow(data), adj.r.squared = NA,
                          minX = NA, maxX = NA, slope = NA, intercept = NA, base = "10")
  }
  return(results)
}