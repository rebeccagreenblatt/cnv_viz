library(shiny)
library(dplyr)
library(plotly)

path <-"../../cnv_viz_demo/rg_viz"
test_files <- list.files(path)

karyotype_filename <- test_files[157]
karyotype_path <- paste0(path, "/", karyotype_filename)
file.copy(karyotype_path, "www")
sample_name <- gsub(".final-diagram.pdf", "", karyotype_filename)

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
      h3(sample_name),
      plotlyOutput("chr_plot", width = "100%"),
      br(),
      conditionalPanel("input.chr != 'all'", plotlyOutput("selected_plot"))
    )
  )
    
)


