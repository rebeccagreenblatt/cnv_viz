library(shiny)
library(dplyr)
library(plotly)

choices <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  
  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px; padding-top: 10px} .selectize-dropdown { font-size: 12px; line-height: 12px;}"),
  
  sidebarLayout(
    sidebarPanel(width = 2,
      selectInput(inputId = "chr",
                  label = "",
                  choices = choices,
                  selected = "all")
    ),
    
    mainPanel(width = 10,
      h2("Sample"),
      plotlyOutput("chr_plot", width = "100%"),
      br(),
      plotlyOutput("gene_plot", width = "100%")
    )
  )
    
)


