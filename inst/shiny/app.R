require(DT)
require(shiny)
require(shinyjs)
require(shinythemes)
require(ggplot2)
require(plotly)
require(data.table)
require(dplyr)
require(tidyr)
require(universalmotif)
require(motifStack)
require(JASPAR2020)
require(biomaRt)
require(GenomicFeatures)
require(ChIPseeker)
require(AnnotationHub)
require(AnnotationDbi)
require(magrittr)

ui <- navbarPage(
  title = "EpiScope",
  id="EpiScope",
  fluid=TRUE,
  theme = shinytheme("yeti"),

  source(file.path("ui", "ui_main_page.R"),  local = TRUE)$value,
  source(file.path("ui", "ui_motif_enrichment.R"),  local = TRUE)$value


)

server <- function(input, output, session) {
  source(file.path("./server/", "server.R"),  local = TRUE)$value
}

shinyApp(ui = ui, server = server)





