library(shiny)
library(dplyr)
library(plotly)
#library(cBioPortalData)
library(DT)
options(dplyr.summarise.inform = FALSE)

# sample <- "CPDC160895"
# sample_path <- paste0("../../CPDC160895/")
#sample <- "CPDC181710"
#sample_path <- paste0("../../CPDC181710/")
sample <- "CPDV193638"
sample_path <- paste0("../../CPDV193638/")
#sample <- "CPDC192023_xGen_20002_CapC_G8_PSL_20070_SEQ_200082_S36"
#sample_path <- paste0("../../EGFRvIII_CDNK2A/",sample, "/")

cnr_file <- read.table(paste0(sample_path, sample, ".final.cnr"), sep = '\t', header = TRUE)
cnr_target <- filter(cnr_file, !gene %in% c("Antitarget", ".")) %>% select(chromosome, gene)
genes <- c("", sort(unique(cnr_target$gene)))

adjusted_gene_file <- read.csv(paste0(sample_path, sample, "_genes.csv")) %>% filter(!gene.symbol %in% c("Antitarget", "."))
adj_genes <- c("", sort(unique(adjusted_gene_file$gene.symbol)))

cbio <- cBioPortal()
cbio_studies <- read.csv("cbio_study_names.csv")

chrs <- c("all", paste0("chr", "1":"22"), "chrX", "chrY")

fluidPage(
  navbarPage("CNViz",
             #textInput("sample_id", "Sample ID", value = "", placeholder = "e.g. CPDC181710"),
             tabPanel("Patient Data (adjusted)", fluid=TRUE,
                      conditionalPanel("output.purecn_passed == 'true'", sidebarLayout(
                          sidebarPanel(
                          width = 3,
                          selectInput(inputId = "adj_chr",
                                      label = "chromosome",
                                      choices = chrs,
                                      selected = "all"),
                          selectInput(inputId = "adj_gene_panel",
                                      label = "gene panel",
                                      choices = c("all genes", "CPD Hematologic Malignancies Panel", "CPD Solid Tumor Panel"),
                                      selected = "all_genes"),
                          selectizeInput(inputId = "adj_gene",
                                         label = "gene",
                                         choices = adj_genes,
                                         selected = ""),
                          downloadButton("karyotype", "karyotype"),
                          br(), br(),
                          img(src = "../adjusted_legend.jpg", width = "100%")
                        ),
                        mainPanel(
                          h3(sample),
                          tableOutput("meta"),
                          textOutput("comment"),
                          br(),
                          column(12, plotlyOutput("adj_chr_plot", width = "100%")),
                          conditionalPanel("input.adj_chr != 'all'",
                                           column(11, offset = 1, 
                                                  strong(textOutput("adj_copies")),
                                                  dataTableOutput("mutations"),
                                                  plotlyOutput("adj_selected_plot"),
                                                  style = 'padding:20px'))
                        ))),
                      conditionalPanel("output.purecn_passed == 'false'", p("PureCN failed for this sample. See Information page."))
                      ),
             tabPanel("Patient Data (unadjusted)", fluid=TRUE, 
                      sidebarLayout(
                        sidebarPanel(
                          width = 3,
                          selectInput(inputId = "chr",
                                      label = "chromosome",
                                      choices = chrs,
                                      selected = "all"),
                          selectInput(inputId = "gene_panel",
                                      label = "gene panel",
                                      choices = c("all genes", "CPD Hematologic Malignancies Panel", "CPD Solid Tumor Panel"),
                                      selected = "all_genes"),
                          selectizeInput(inputId = "gene",
                                         label = "gene",
                                         choices = genes,
                                         selected = ""),
                          br(), br(), 
                          img(src = "../unadjusted_legend.jpg", width = "100%"),
                          br(), br()
                        ),
                        mainPanel(
                          width = 9,
                          h3(sample),
                          column(12, plotlyOutput("chr_plot", width = "100%")),
                          conditionalPanel("input.chr != 'all'",
                                           column(11, offset = 1, 
                                                  plotlyOutput("selected_plot"),
                                                  style = 'padding:20px'))
                        )
                      )),
             tabPanel("TCGA Pan-Cancer Atlas 2018 Data", fluid=TRUE,
                      selectizeInput(inputId = "cancer", label = "cancer",
                                     choices = c("", cbio_studies$Cancer),
                                     selected = ""),
                      dataTableOutput("cbioOutput")),
             tabPanel(icon("info-circle", class = NULL, lib = "font-awesome"), fluid=TRUE,
                      h5("Patient Data (Unadjusted)"),
                      p("Unadjusted data is obtained from CNVkit. CNVkit takes hybrid capture DNA sequencing data and infers copy number. CNViz rounds copy number estimates to the nearest integer. Unadjusted refers to the fact that this data is not adjusted for the sample's purity or the tumor's ploidy. Instead, it uses the raw sequencing data. This tab includes probe-level data for each gene."),
                      p(em("Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014). CNVkit: Genome-wide copy number detection and visualization from targeted sequencing. PLOS Computational Biology 12(4):e1004873.")),
                      div(style="display: inline-block;", a("CNVkit paper,", target ="_blank", href = "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004873")),
                      div(style="display: inline-block;", a("CNVkit documentation", target ="_blank", href = "https://cnvkit.readthedocs.io/en/stable/")),
                      br(), br(),
                      h5("Patient Data (Adjusted)"),
                      p("Adjusted data is obtained from PureCN, using data generated by CNVkit. PureCN estimates sample purity, tumor ploidy, loss of heterozygosity (LOH) and integer copy number. Somatic mutations are also inferred. Of note, germline mutations are not displayed. Purity is PureCN's estimate for the percentage of the sample containing tumor cells. Ploidy is estimated by PureCN, then rounded to the nearest whole number. Only segments where copy number differs from tumor ploidy or those with LOH are highlighted. LOH is determined using SNPs, so if a gene or segment is displayed as having 0 or 1 copies, it still may not be marked as LOH if there were no SNPs in that gene or segment. If a sample is marked as Flagged, this means PureCN has deemed the results as unreliable. If a sample is marked as Failed, this means PureCN was unable to estimate purity and ploidy. Adjusted allelic fraction refers to the allelic fraction we would expect to see in the sample if these are the correct purity, ploidy and copy number estimates. Lastly, a copy number greater than 64 will have log2 value > 5, and thus not fit on the plot - it is artifically brought into view, but the hover information will be correct. Copy number of 0 will be displayed at -2.5, as it would otherwise be negative infinity."),
                      p(strong("If this panel does not appear, this means PureCN failed for this sample. See the unadjusted panel, but know the results are not adjusted for purity or ploidy.")),
                      p(em("Riester, M., Singh, A.P., Brannon, A.R. et al. (2016). PureCN: copy number calling and SNV classification using targeted short read sequencing. Source Code Biol Med 11(1), 13.")),
                      div(style="display: inline-block;", a("PureCN paper,", target ="_blank", href = "https://scfbm.biomedcentral.com/articles/10.1186/s13029-016-0060-z")),
                      div(style="display: inline-block;", a("PureCN documentation", target ="_blank", href = "https://rdrr.io/bioc/PureCN/f/inst/doc/PureCN.pdf")),
                      br(), br(),
                      h5("TCGA Pan-Cancer Atlas Data"),
                      p("TCGA Pan-Cancer Atlas Data was obtained from cBioPortal's R package cBioPortalData. This same information can be found on cBioPortal's website. Similar information from many additional studies can also be found on their website. The copy number data displayed was generated by the GISTIC algorithm."),
                      p(" - Deep Deletion indicates a deep loss, possibly a homozygous deletion"),
                      p(" - Shallow Deletion indicates a shallow loss, possibley a heterozygous deletion"),
                      p(" - Gain indicates a low-level gain (a few additional copies, often broad)"),
                      p(" - Amplification indicate a high-level amplification (more copies, often focal)"),
                      div(style="display: inline-block;", a("cBioPortal website,", target ="_blank", href = "https://www.cbioportal.org/")),
                      div(style="display: inline-block;", a("cBioPortal FAQ,", target ="_blank", href = "https://docs.cbioportal.org/1.-general/faq")),
                      div(style="display: inline-block;", a("cbioportaldata R Package,", target ="_blank", href = "https://bioconductor.org/packages/release/bioc/html/cBioPortalData.html")),
                      div(style="display: inline-block;", a("GISTIC paper", target ="_blank", href = "https://pubmed.ncbi.nlm.nih.gov/18077431/")),
                      br(), br(), br(), br(), br(), br(), br(), br(),
                      p(em("Questions? Feedback? Please email rebecca.greenblatt@pennmedicine.upenn.edu")),
                      div(style="display: inline-block;", a("CNViz GitHub", target ="_blank", href = "https://github.com/rebeccagreenblatt/cnv_viz"))
                      
                      ))
)