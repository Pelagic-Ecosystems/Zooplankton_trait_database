# Zooplankton_trait_database
The global zooplankton trait database developed by Pata and Hunt. This repository also contains scripts for generating figures summarizing and analyzing the contents of the database, and a template for subsetting trait information from the database and estimating missing trait values.

This zooplankton trait database is currently stored in the data_input/ folder. The database is stored in two levels: level 1 contain the standardized and curated individual-level trait records and level 2 contain the species-level trait records. Necessary supporting data files to curate and analyze the database contents are provided in the data_input/ folder. The data_output/ folder contains the stored results of the analysis of the level 2 data for phylogenetic trait-trait correlations, regression models of allometric scaling relationships between traits, and a data file of the raw results of the trait value estimation tests.

The R script files provide the scripts for generating most of the figures in the Pata and Hunt (for submission) manuscript and the scripts used in processing and analyzing the trait database. The latter are not executed by default as these may take up considerable computation time. In displaying some of the figures, the previously calculated results in the data_output/ folder are loaded. The toolkit.R file contains useful functions for data curation and analysis.

The figures in the manuscript can be generated from the following files:
- Fig 3 - 2-Database summary figures and tables.Rmd
- Fig 4 - 3-Database summary phylo dendrogram.Rmd (partly, annotated externally)
- Fig 5 - 2-Database summary figures and tables.Rmd
- Fig 6 - 4-Association between traits.Rmd
- Fig 7 - 6-Filling in taxonomic gaps tests and figures.Rmd
- Fig 8 - 6-Filling in taxonomic gaps tests and figures.Rmd
- Supp Fig 1 - 3-Database summary phylo dendrogram.Rmd
- Supp Fig 2 - 2-Database summary figures and tables.Rmd
- Supp Fig 3 - 6-Filling in taxonomic gaps tests and figures.Rmd

Please note that the zooplankton trait database is continuously being developed and will eventually be hosted in a different repository. We welcome any feedback and support to improve the contents and structure of the database. We are also open to collaborations in applying and analyzing the database for your study region or research interests and/or to assisting you in extracting trait information from the database. Please send us an email at p.pata@oceans.ubc.ca.
