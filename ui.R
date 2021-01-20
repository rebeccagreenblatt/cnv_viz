library(shiny)
library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(plotly)

path <- "../../1218_rg_cnv_viz/"
karyotype_filename <- "CPDV193638.cnv.pdf"
karyotype_path <- paste0(path, karyotype_filename)
file.copy(karyotype_path, "www")

cnr_name <- "../../CPDC181710/CPDC181710.final.cnr"
cnr_file <- read.table(cnr_name, sep = '\t', header = TRUE)
cnr_target <- filter(cnr_file, !gene %in% c("Antitarget", ".")) %>% select(chromosome, gene)
genes <- c("", sort(unique(cnr_target$gene)))

sample_name <- tail(strsplit(cnr_name, "\\.|\\/")[[1]], 3)[1]

cbio_studies <- read.csv("../cbio_study_names.csv")

choices <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  
  #tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px; padding-top: 10px} .selectize-dropdown { font-size: 12px; line-height: 12px;}"),
  
  navbarPage("CNViz",
             tabPanel("Patient Data (adjusted)", fluid=TRUE,
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          selectInput(inputId = "adj_chr",
                                      label = "chromosome",
                                      choices = choices,
                                      selected = "all"),
                          selectizeInput(inputId = "adj_gene",
                                         label = "gene",
                                         choices = genes,
                                         selected = ""),
                          br(), 
                          img(src = "../adjusted_legend.jpg", width = "100%")
                        ),
                        mainPanel(
                          h3(sample_name),
                          tableOutput("meta"),
                          textOutput("comment"),
                          br(),
                          plotlyOutput("adj_chr_plot", width = "100%"),
                          textOutput("adj_copies"),
                          dataTableOutput("mutations"),
                          br(),
                          conditionalPanel("input.adj_chr != 'all'", plotlyOutput("adj_selected_plot"))
                        )
                      )),
             tabPanel("Patient Data (unadjusted)", fluid=TRUE, 
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
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
                          #a("Launch GISTIC", target = "_blank", href = 'http://portals.broadinstitute.org/tcga/gistic/browseGisticByGene#'),
                          #br(), br(), 
                          img(src = "../unadjusted_legend.jpg", width = "100%"),
                          br(), br()
                        ),
                        mainPanel(
                          width = 9,
                          h3(sample_name),
                          plotlyOutput("chr_plot", width = "100%"),
                          br(),
                          conditionalPanel("input.chr != 'all'", plotlyOutput("selected_plot"))
                        )
                      )),
             tabPanel("TCGA Pan-Cancer Atlas 2018 Data", fluid=TRUE,
                      selectizeInput(inputId = "cancer", label = "cancer",
                                     choices = c("", cbio_studies$Cancer),
                                     selected = ""),
                      dataTableOutput("cbioOutput")))
)
