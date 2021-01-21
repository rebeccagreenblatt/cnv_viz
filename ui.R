library(shiny)
library(dplyr)
library(plotly)
library(cBioPortalData)
library(DT)
library(stringr)
options(dplyr.summarise.inform = FALSE)

sample <- "CPDC160895"
sample_path <- "../../CPDC160895/"

karyotype_filename <- paste0(sample, ".cnv.pdf")
karyotype_file <- paste0(sample_path, karyotype_filename)
file.copy(karyotype_file, "www")

cnr_file <- read.table(paste0(sample_path, sample, ".final.cnr"), sep = '\t', header = TRUE)
cnr_target <- filter(cnr_file, !gene %in% c("Antitarget", ".")) %>% select(chromosome, gene)
genes <- c("", sort(unique(cnr_target$gene)))

adjusted_gene_file <- read.csv(paste0(sample_path, sample, "_genes.csv"))
adj_genes <- c("", sort(unique(adjusted_gene_file$gene.symbol)))

cbio <- cBioPortal()
cbio_studies <- read.csv("cbio_study_names.csv")

chrs <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  
  #tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 12px; padding-top: 10px} .selectize-dropdown { font-size: 12px; line-height: 12px;}"),
  
  navbarPage("CNViz",
             tabPanel("Patient Data (adjusted)", fluid=TRUE,
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          selectInput(inputId = "adj_chr",
                                      label = "chromosome",
                                      choices = chrs,
                                      selected = "all"),
                          selectizeInput(inputId = "adj_gene",
                                         label = "gene",
                                         choices = adj_genes,
                                         selected = ""),
                          br(), 
                          img(src = "../adjusted_legend.jpg", width = "100%")
                        ),
                        mainPanel(
                          h3(sample),
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
                                      choices = chrs,
                                      selected = "all"),
                          selectizeInput(inputId = "gene",
                                         label = "gene",
                                         choices = genes,
                                         selected = ""),
                          a("Karyotype",target="_blank",href=karyotype_filename),
                          br(), br(), 
                          img(src = "../unadjusted_legend.jpg", width = "100%"),
                          br(), br()
                        ),
                        mainPanel(
                          width = 9,
                          h3(sample),
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
