library(shiny)
library(dplyr)
library(plotly)


karyotype_filename <- "CPDV183182_20179_D1_20007_B_20184_200206_S4.final-diagram.pdf"
sample_name <- gsub(".final-diagram.pdf", "", karyotype_filename)

choices <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  
  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px; padding-top: 10px} .selectize-dropdown { font-size: 12px; line-height: 12px;}"),
  
  sidebarLayout(
    sidebarPanel(width = 2,
      selectInput(inputId = "chr",
                  label = "",
                  choices = choices,
                  selected = "all"),
      a("karyotype",target="_blank",href=karyotype_filename)
    ),
    
    mainPanel(width = 10,
      h2(sample_name),
      plotlyOutput("chr_plot", width = "100%"),
      br(),
      plotlyOutput("selected_plot")
    )
  )
    
)


