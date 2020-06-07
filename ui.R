library(shiny)
library(dplyr)
library(plotly)

karyotype_filename <- "CPDV183182_20179_D1_20007_B_20184_200206_S4.final-diagram.pdf"
sample_name <- gsub(".final-diagram.pdf", "", karyotype_filename)

path <-"../../cnv_viz_demo/rg_viz"
test_files <- list.files(path)
cnr_filename <- paste0(path, "/", test_files[158])
cnr <- read.table(file = cnr_filename, sep = '\t', header = TRUE)
cnr_target <- filter(cnr, !gene %in% c("Antitarget", ".")) %>% select(chromosome, gene)
genes <- c("", sort(unique(cnr_target$gene)))

choices <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  
  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px; padding-top: 10px} .selectize-dropdown { font-size: 12px; line-height: 12px;}"),
  
  sidebarLayout(
    sidebarPanel(width = 2,
      selectInput(inputId = "chr",
                  label = "chromosome",
                  choices = choices,
                  selected = "all"),
      selectizeInput(inputId = "gene",
                     label = "gene",
                     choices = genes,
                     selected = ""),
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


