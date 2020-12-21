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

cbio_studies <- read.csv("../cbio_study_names.csv")

choices <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  
  tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px; padding-top: 10px} .selectize-dropdown { font-size: 12px; line-height: 12px;}"),
  
  sidebarLayout(
    sidebarPanel(width = 3,
      selectInput(inputId = "chr",
                  label = "chromosome",
                  choices = choices,
                  selected = "all"),
      selectizeInput(inputId = "gene",
                     label = "gene",
                     choices = genes,
                     selected = ""),
      a("Karyotype",target="_blank",href=karyotype_filename),
      br(), br(), 
      a("Launch GISTIC", target = "_blank", href = 'http://portals.broadinstitute.org/tcga/gistic/browseGisticByGene#'),
      br(), br(), 
      img(src = "cnviz_legend.png", width = "100%"),
      br(), br()
    ),
    
    
    
    mainPanel(width = 9,
              tabsetPanel(type = "tabs",
              tabPanel("Patient Data",
      h3(sample_name),
      plotlyOutput("chr_plot", width = "100%"),
      br(),
      conditionalPanel("input.chr != 'all'", plotlyOutput("selected_plot")),
     ),
      tabPanel("TCGA Pan-Cancer Atlas 2018 Data",  selectizeInput(inputId = "cancer",
                                            label = "cancer",
                                            choices = c("", cbio_studies$Cancer),
                                            selected = ""),
               dataTableOutput("cbioOutput")))
    )
  )
    
)


