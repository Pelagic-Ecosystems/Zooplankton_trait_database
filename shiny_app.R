# Create Shiny App
# Created by P. Pata

# install.packages("remotes")
# remotes::install_github("rstudio/shinyuieditor")

shinyuieditor::launch_editor(app_loc = here::here("ZTD_app/"))

shiny::runApp("ZTD_app")
