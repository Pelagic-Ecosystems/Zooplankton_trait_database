# Trait database curation shiny app
# Created by: P. Pata
# Last updated: February 27, 2024
# 
# This app explores the contents of the global zooplankton trait database. It provides some basic options to filter the contents of the dataset. The number of data points (by species) is displayed in a figure and a direct download button is provided.


# # To download package
# install.packages("remotes")
# remotes::install_github("rstudio/shinyuieditor")
# # To edit and run shiny app
# shinyuieditor::launch_editor(app_loc = here::here("ZTD_app/"))
# shiny::runApp("ZTD_app")



library(shiny)
library(plotly)
library(gridlayout)
library(bslib)
library(tidyverse)

# UI Layout --------------------------------------------------

ui <- grid_page(
  layout = c(
    "header  header",
    "sidebar plot  "
  ),
  row_sizes = c(
    "55px",
    "1fr"
  ),
  col_sizes = c(
    "500px",
    "1fr"
  ),
  gap_size = "1rem",
  grid_card(
    area = "sidebar",
    card_header("Filter the database"),
    card_body(
      radioButtons(
        selected = "2",
        inputId = "myTraitLevel",
        label = "Select the level of trait data you want to work with:",
        choices = list(
          "Level 1: Individual level" = "1",
          "Level 2: Species level" = "2"
        ),
        width = "100%"
      ),
      checkboxGroupInput(
        selected = c(
          "Amphipod",
          "Euphausiid",
          "Pteropod"
        ),
        inputId = "myMajorGroupBox",
        label = "Select the major groups:",
        choices = list(
          "Calanoid copepods" = "Calanoid",
          "Non-calanoid copepods" = "Non-calanoid",
          "Amphipods" = "Amphipod",
          "Euphausiids" = "Euphausiid",
          "Pteropods" = "Pteropod"
        ),
        width = "100%"
      ),
      checkboxGroupInput(
        selected = c(
          "dryWeight",
          "respirationRate_15C"
        ),
        inputId = "myTraitBox",
        label = "Select the traits:",
        choices = list(
          "Body length" = "bodyLengthMax",
          "Dry weight" = "dryWeight",
          "Respiration rate (individual)" = "respirationRate_15C",
          "Respiration rate (weight specific)" = "respirationRate_WSC_15C",
          "Feeding mode" = "feedingMode",
          "Trophic group" = "trophicGroup"
        )
      ),
      checkboxGroupInput(
        selected = c(
          "scientificName","traitName","traitValue","traitUnit",
          "valueType","primaryReference","secondaryReference"
        ),
        inputId = "myFieldBox",
        label = "Select the database fields to download:",
        choices = list(
          "Scientific name" = "scientificName",
          "Trait" = "traitName",
          "Trait value" = "traitValue",
          "Trait unit" = "traitUnit",
          "Value type" = "valueType",
          "Primary reference" = "primaryReference",
          "Secondary reference" = "secondaryReference",
          "Life stage" = "lifeStage",
          "Associated temperature" = "assocTemperature"
        )
      )
    )
  ),
  grid_card_text(
    area = "header",
    content = "Global Zooplankton Trait Database",
    alignment = "center",
    is_title = FALSE
  ),
  grid_card(
    area = "plot",
    card_header("Contents of data filtered"),
    card_body(
      plotlyOutput(
        outputId = "plot",
        width = "100%",
        height = "100%"
      ),
      downloadButton(
        outputId = "myDownloadButton",
        label = "Download to csv file",
        class = NULL,
        icon = shiny::icon("download")
      )
    )
  )
)

# Functions --------------------------------------------------

## Function to open and filter the appropriate trait data source
filterTraitData <- function(traitLevel, majorGroupList, traitList) {
  # Open the trait database for whatever level
  if (traitLevel == 1){
    traitDB <- read.csv("../data_input/Trait_dataset_level1/trait_dataset_level1-2023-08-15.csv")
  }
  else {
    traitDB <- read.csv("../data_input/Trait_dataset_level2/trait_dataset_level2-2023-09-14.csv")
  }
  
  traitDB <- traitDB %>% 
    filter(majorgroup %in% majorGroupList) %>% 
    filter(traitName %in% traitList)
  
}

plotDataHist <- function(traitDB){
  trait.directory <- read.csv("../data_input/trait_directory_20230628.csv") %>% 
    distinct(traitID, .keep_all = TRUE)
    
  majorgroup.colors <- data.frame(
    majorgroup = c("Calanoid","Non-calanoid","Mysid","Amphipod","Decapod","Euphausiid",
                   "Ostracod","Polychaete","Pteropod",
                   "Chaetognath","Appendicularian","Thaliacean","Cladoceran",
                   "Hydromedusae","Siphonophore","Scyphomedusae","Ctenophore"),
    color = c("blue","navy","slateblue2","cornflowerblue","cyan2","forestgreen",
                    "seagreen3","tan1","green",
                    "darkorchid4","violetred","hotpink3","saddlebrown",
                    "tomato1","goldenrod3","yellow2","red4"))
                    
  # To slightly change the labels for the plot
  majorgroup.colors.ed.labels <-  c("Calanoids","Non-calanoids","Mysids",
                                    "Amphipods","Decapods","Euphausiids",
                                    "Ostracods","Polychaetes","Pteropods",
                                    "Chaetognaths","Appendicularians","Thaliaceans",
                                    "Cladocerans","Hydromedusae",
                                    "Siphonophores","Scyphozoans","Ctenophores")
  
  trait.sum <- traitDB %>% 
    filter(valueType %in% c("numeric","categorical")) %>% 
    dplyr::select(traitID, traitName, taxonID, scientificName, majorgroup, taxonRank,
                  primaryReference, secondaryReference) %>% 
    left_join(distinct(trait.directory,
                       traitBucket, traitNum, trait, traitID),
              by = "traitID") %>% 
    # only consider species and subspecies level when presenting results
    filter(taxonRank %in% c("Species","Subspecies")) 
  
  # Number of species with at least one trait info
  length(unique(trait.sum$taxonID))
  
  # count number of unique species-trait records (note that this is not traitName thus nitrogenPDW and nitrogenTotal are the same trait:nitrogen content)
  trait.sum <- trait.sum %>% 
    distinct(traitBucket, traitNum, traitName, majorgroup, taxonID) %>% 
    group_by(traitBucket, traitNum, traitName, majorgroup) %>% 
    # Number of species-trait per major group
    summarise(Nsp.mg = n()) %>% 
    group_by(traitBucket, traitNum, traitName) %>% 
    # Number of species per trait
    mutate(Nspecies = sum(Nsp.mg)) %>% 
    arrange(traitNum) %>% 
    ungroup() %>% 
    # arrange by trait bucket and decreasing total value
    mutate(majorgroup = factor(majorgroup, levels = majorgroup.colors$majorgroup),
           traitBucket = factor(traitBucket, 
                                levels = c("morphological","biochemical composition",
                                           "physiological","behavioral",
                                           "life history"))) %>% 
    mutate(traitName = fct_reorder(traitName, Nspecies))
  
  
  bplot <- ggplot(trait.sum, aes(y = traitName, x = Nsp.mg, 
                                 fill = majorgroup)) +
    geom_bar(stat = "identity") +
    
    # Add number at end of bar
    geom_text(aes(label = after_stat(x), group = traitName),
              stat = "summary", fun = sum, 
              size = 5) +
    scale_fill_manual(values = majorgroup.colors$color,
                      limits = majorgroup.colors$majorgroup,
                      labels = majorgroup.colors.ed.labels,
                      name = "Major group") +
    # coord_flip() +
    ylab("Trait") + xlab("Number of species") +
    # facets proportioned by number of traits
    facet_grid(traitBucket ~., scales = "free_y", space = "free", switch = "x")  +
    theme_bw() +
    theme(strip.placement = "outside", text = element_text(size = 12),
          legend.position = "right",
          legend.background = element_rect(colour = NA)) 
}


# Server ------------------------------------------------------------

server <- function(input, output) {

   # Display an overview of the data subset contents
  output$plot <- renderPlotly({
    traitDB <- filterTraitData(input$myTraitLevel,
                               input$myMajorGroupBox,
                               input$myTraitBox)
    
    plotDataHist(traitDB)
  })
  
  # Download the data to a system folder
  output$myDownloadButton <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("trait_database_subset_",Sys.Date(),".csv")
    },
    content = function(file) {
      # Select columns and write dataset to file
      traitDB <- filterTraitData(input$myTraitLevel,
                                 input$myMajorGroupBox,
                                 input$myTraitBox) %>% 
        dplyr::select(any_of(input$myFieldBox))
      
      write.csv(traitDB, file)
    }
  )
}

shinyApp(ui, server)
  

