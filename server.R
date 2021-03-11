library(karyoploteR)
library(CopyNumberPlots)
library(GenomicRanges)

# sample <- "CPDC160895"
# sample_path <- paste0("../../CPDC160895/")
#sample <- "CPDC181710"
#sample_path <- paste0("../../CPDC181710/")
sample <- "CPDV193638"
sample_path <- paste0("../../CPDV193638/")

source("plot_function.R")

cnr_file <- read.table(paste0(sample_path, sample, ".final.cnr"), sep = '\t', header = TRUE)
cns_file <- read.table(paste0(sample_path, sample, ".final.cns"), sep='\t', header = TRUE)
adjusted_gene_file <- read.csv(paste0(sample_path, sample, "_genes.csv"))
loh_file <- read.csv(paste0(sample_path, sample, "_loh.csv")) #not using dna seg file, only using loh seg file for adjusted segments
variants_file <- read.csv(paste0(sample_path, sample, "_variants.csv"))
metadata_file <- read.csv(paste0(sample_path, sample, ".csv"))

green <- "#009E73" # normal
blue <- "#0072B2" # out of range genes
pink <- "#CC79A7" #mutation
orange <- "#D55E00" #segment
black <- "#000000" #LOH
white <- "#FFFFFF"

cbio <- cBioPortal()
cbio_studies <- read.csv("cbio_study_names.csv")

heme <- read.csv("hematologic_malignancies_panel.csv")
solid <- read.csv("solid_tumor_panel.csv")

function(input, output, session) {
  
  cnr <- cnr_file
  cnr_target <- filter(cnr, !gene %in% c("Antitarget", ".")) %>% 
    mutate(m_probe = (start + end)/2, 
           cn_est = 2*2^log2)
  
  by_gene <- cnr_target %>% group_by(chromosome, gene) %>%
    summarise(s = min(start), e = max(end),
              mean_cn_est = weighted.mean(cn_est, weight),
              total_weight = sum(weight))%>%
    mutate(mean_log2 = log(mean_cn_est/2,2)) %>%
    mutate(m = (s+e)/2, 
           cn = round((2^mean_log2)*2)) %>%
    mutate(blue = as.numeric(mean_log2 < -0.41 | mean_log2 > 0.32),
           copies = ifelse(cn == 1, " copy", " copies")) 
  
  by_gene$mean_log2 <- ifelse(by_gene$mean_log2 < -2.5, -2.5, by_gene$mean_log2) #keep all points in range of y axis
  by_gene$mean_log2 <- ifelse(by_gene$mean_log2 > 5, 5, by_gene$mean_log2)
  
  cns <- filter(cns_file, gene != "-")
  gl <- filter(cns, log2 < -0.41 | log2 > 0.32) #segments with copy # < 1.5 or > 2.5
  gl$log2 <- ifelse(gl$log2 < -2.5, -2.5, gl$log2) #keep all points in range of y axis
  
  by_gene$marker_color <- ifelse(by_gene$cn == 2, white, blue)
  by_gene$outline_color <- ifelse(by_gene$marker_color == white, green, blue)
  
  gene_list <- reactive({ switch(input$gene_panel,
                                 "all genes" = by_gene$gene,
                                 "CPD Hematologic Malignancies Panel" = heme$gene,
                                 "CPD Solid Tumor Panel" = solid$gene) 
  })
  
  output$chr_plot <- renderPlotly({ getPlots(by_gene, gl, cytoband_data, ploidy = 2, gene_list(), input$chr, adjusted = FALSE) })
  
  observe({
    updateSelectInput(session, "gene", choices = switch(input$gene_panel,
                                                            "all genes" = c("", sort(by_gene$gene)),
                                                            "CPD Hematologic Malignancies Panel" = c("", sort(heme$gene)),
                                                            "CPD Solid Tumor Panel" = c("", sort(solid$gene)))
    )})
  
  observe({
    updateSelectInput(session, "chr", selected = as.character(by_gene[by_gene$gene==input$gene,]$chromosome[1]))
  })
  
  d <- reactive ({ event_data(event="plotly_click", source = "a")[[3]] })
  
  observeEvent(event_data(event = "plotly_click", source = "a"), {
    updateSelectizeInput(session, "gene", selected = by_gene[by_gene$m == d(),]$gene[1])
  })
  
  observeEvent(input$chr, {
    updateSelectizeInput(session, "gene", selected = 
                           ifelse(by_gene[by_gene$gene == input$gene,"chromosome"][1] == input$chr, input$gene, ""))
  })
  
  gene_data <- eventReactive(input$gene,{
    filter(cnr_target, gene == input$gene)
  })
  
  test1 <- reactive({ nrow(gene_data()) })
  
  cn <- reactive({ by_gene[by_gene$gene == input$gene,]$cn[1] })
  copy <- reactive({ ifelse(cn() == 1, "copy", "copies")})
  
  seg <- reactive({ filter(cns, chromosome == input$chr, start == d() | end == d()) })
  seg_genes <- reactive({ strsplit(seg()$gene[1], ",")[[1]] })
  seg_dat <- reactive({ data.frame(gene = unique(seg_genes())) })
  seg_gene_dat <- reactive({ inner_join(seg_dat(), cnr_target, by =c("gene")) })
  
  test2 <- reactive({ nrow(seg_gene_dat()) })
  
  plot_check <- reactive({ max(input$gene != "", d()) >= 1 })
  
  output$selected_plot <- renderPlotly({
    req(plot_check())
    if(test1() > 0 & input$chr != "all"){
      plot_ly(height = 250, type = 'scatter', mode = 'markers') %>%
        add_trace(x = gene_data()$m_probe, 
                  y = gene_data()$log2, 
                  marker = list(color='yellow', size = 15*gene_data()$weight),
                  showlegend = F) %>%
        add_segments(x = min(gene_data()$start), xend = max(gene_data()$end), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        add_segments(x = min(gene_data()$start), xend = max(gene_data()$end), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        layout(title = paste0("probe data: ", gene_data()$gene[1], " (", cn(), " ", copy(),")"), xaxis = list(tickfont = list(size = 6), range = min(gene_data()$s), max(gene_data()$e)), yaxis=list(tickfont = list(size = 6), title = "log(2) copy number ratio", titlefont = list(size = 8), range = c(-3, 6)), paper_bgcolor='#fafafa', plot_bgcolor='#fafafa',margin = list(t = 80))
    } else if(test2() > 0 & input$chr != "all"){
      plot_ly(height = 250, type = 'scatter', mode = 'markers') %>%
        add_trace(x = seg_gene_dat()$m_probe, 
                  y = seg_gene_dat()$log2, 
                  hoverinfo = 'text',
                  text = seg_gene_dat()$gene,
                  marker = list(color='red', size = seg_gene_dat()$weight*15),
                  showlegend = F) %>%
        add_segments(x = min(seg_gene_dat()$start), xend = max(seg_gene_dat()$end), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        add_segments(x = min(seg_gene_dat()$start), xend = max(seg_gene_dat()$end), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        layout(title = paste0("probe data: segment (", min(seg_gene_dat()$start) ,"-", max(seg_gene_dat()$end),")"), xaxis = list(tickfont = list(size = 6), range = min(seg_gene_dat()$s), max(seg_gene_dat()$e)), yaxis=list(tickfont = list(size = 6), title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)), paper_bgcolor='#fafafa', plot_bgcolor='#fafafa',margin = list(t = 80))
    } 
  })
  
  ##cbio table
  if(exists("cbio")){
  cbio_studyId <- reactive({cbio_studies$studyId[cbio_studies$Cancer == input$cancer]})
  cbio_table <- reactive({ 
    getDataByGenePanel(api = cbio, 
                       studyId = cbio_studyId(), 
                       genePanelId = "IMPACT468",
                       molecularProfileId = paste0(cbio_studyId(), "_gistic"), 
                       sampleListId = paste0(cbio_studyId(), "_cna"))
  })
  cbio_dat <- reactive({ data.frame(cbio_table()[[1]], stringsAsFactors = FALSE) })
  output$cbioOutput <- DT::renderDataTable({
    req(nchar(input$cancer)>0)
    datatable(cbio_dat() %>% group_by(hugoGeneSymbol) %>% 
                summarise(Gain = sum(value ==1)/n(), 
                          Amplification = sum(value == 2)/n(), 
                          ShallowDeletion = sum(value == -1)/n(), 
                          DeepDeletion = sum(value == -2)/n()),
              rownames = FALSE) %>% 
      formatPercentage(c("Gain", "Amplification", "ShallowDeletion", "DeepDeletion"), 2)
  })
}
  
  ##adjusted data
  adj_bygene <- adjusted_gene_file %>% filter(!(gene.symbol %in% c("Antitarget", ".")), C.flagged %in% c(NA, FALSE))
  colnames(adj_bygene)[colnames(adj_bygene) == "gene.symbol"] <- "gene"
  colnames(adj_bygene)[colnames(adj_bygene) == "chr"] <- "chromosome"
  colnames(adj_bygene)[colnames(adj_bygene) == "number.targets"] <- "total_weight"
  adj_bygene$mean_log2 <- round(log(adj_bygene$C/2, 2),2)
  adj_bygene$mean_log2 <- ifelse(adj_bygene$C < 0.4, -2.5, adj_bygene$mean_log2) #ensure all points are within range of y axis
  adj_bygene$mean_log2 <- ifelse(adj_bygene$mean_log2 > 5, 5, adj_bygene$mean_log2) #ensure all points are within range of y axis
  adj_bygene$cn <- round(adj_bygene$C)
  adj_bygene$copies <- ifelse(adj_bygene$cn == 1, " copy", " copies")
  adj_bygene$loh <- ifelse(is.na(adj_bygene$loh), FALSE, adj_bygene$loh)
  
  adj_bygene$all_genes <- rep(TRUE, nrow(adj_bygene))
  adj_bygene <- left_join(adj_bygene, heme, by=c("gene"))
  adj_bygene <- left_join(adj_bygene, solid, by=c("gene"))
  
  adj_cns <- loh_file %>% select(chr, start, end, type, C) %>% mutate(cn = C, chromosome = chr)
  adj_cns$log2 <- log(adj_cns$cn/2, 2)
  adj_cns$log2 <- ifelse(adj_cns$cn < 0.4, -2.5, adj_cns$log2)
  adj_cns$loh <- adj_cns$type %in% c("COPY-NEUTRAL LOH", "LOH", "WHOLE ARM COPY-NEUTRAL LOH", "WHOLE ARM LOH")
  
  variants <- variants_file
  variants$AR <- round(variants$AR, 2)
  variants$ML.AR <- round(variants$ML.AR, 2)
  variants$ML.AR <- ifelse(variants$ML.C == 0, "N/A", variants$ML.AR) #shouldn't have an allelic fraction if copy number is truly 0
  mutations <- filter(variants, ML.SOMATIC == TRUE & FLAGGED == FALSE & !(gene.symbol %in% c(NA, "<NA>", "Antitarget"))) %>% 
    select(gene.symbol, ID, depth, AR, ML.AR)
  colnames(mutations) <- c("Gene", "ID", "Depth", "Raw Allelic Fraction", "Adjusted Allelic Fraction (adjusted for purity, ploidy and copy # estimates)") 
  
  varsome_links <- c()
  for(i in 1:nrow(mutations)){
    t <-  strsplit(mutations$ID[i], "chr|\\:|\\_|\\/")[[1]]
    link <- paste0("https://varsome.com/variant/hg38/", t[2], "%3A", t[3], "%3A", t[4], "%3A", t[5])
    varsome_links <- c(varsome_links, link)
  }
  mutations$ID <- paste0("<a href='", varsome_links, "' target='_blank'>", mutations$ID, "</a>" )
  
  genes_wmutations <- data.frame("gene" = unique(mutations$Gene), "mutation" = TRUE)
  
  adj_bygene <- left_join(adj_bygene, genes_wmutations, by = c("gene"))
  adj_bygene$mutation <- ifelse(is.na(adj_bygene$mutation), FALSE, adj_bygene$mutation)
  
  metadata <- metadata_file
  metadata$Sex = ifelse(metadata$Sex == TRUE, "MALE", "FEMALE")
  metadata$Ploidy = round(metadata$Ploidy)
  output$meta <- renderTable(metadata %>% select(Purity, Ploidy, Sex, Contamination, Flagged, Failed), bordered = TRUE)
  output$comment <- renderText(ifelse(!is.na(metadata$Comment[1]), paste0("Comment: ", metadata$Comment[1]), ""))
  ploidy <- metadata$Ploidy
  
  adj_by_gene <- adj_bygene %>% 
    mutate(m = (start+end)/2)#, cn = C)
  
  adj_gl <- filter(adj_cns, log2 < log((ploidy-0.5)/2,2) | log2 > log((ploidy+0.5)/2,2) | loh == TRUE) #only show segments +/- 0.5 of plidy AND any LOH segments
  
  adj_colors <- ifelse(adj_by_gene$cn == 2, white, blue)
  adj_colors <- ifelse(adj_by_gene$loh == TRUE, black, adj_colors)
  adj_colors <- ifelse(adj_by_gene$mutation == TRUE, pink, adj_colors)
  adj_colors2 <- ifelse(adj_colors == white, green, adj_colors)
  adj_colors2 <- ifelse(adj_colors == pink & adj_by_gene$loh == TRUE, black, adj_colors2)
  adj_by_gene$marker_color <- adj_colors
  adj_by_gene$outline_color <- adj_colors2
  
  cytoband_data <- read.csv("cytoband_data.csv")
  
  adj_gene_list <- reactive({ switch(input$adj_gene_panel,
                                  "all genes" = adj_by_gene$gene,
                                  "CPD Hematologic Malignancies Panel" = heme$gene,
                                  "CPD Solid Tumor Panel" = solid$gene) 
  })
  
  output$adj_chr_plot <- renderPlotly({ getPlots(adj_by_gene, adj_gl, cytoband_data, ploidy, adj_gene_list(), input$adj_chr, adjusted = TRUE) })
  
  
  observe({
    updateSelectInput(session, "adj_gene", choices = switch(input$adj_gene_panel,
                                                  "all genes" = c("", sort(adj_by_gene$gene)),
                                                  "CPD Hematologic Malignancies Panel" = c("", sort(heme$gene)),
                                                  "CPD Solid Tumor Panel" = c("", sort(solid$gene)))
    )})
  
  
  observe({
    updateSelectInput(session, "adj_chr", selected = as.character(adj_by_gene[adj_by_gene$gene==input$adj_gene,]$chromosome[1]))
  })
  
  adj_d <- reactive ({ event_data(event="plotly_click", source = "b")[[3]] }) 
  
  chromosome_names <- c(paste0("chr", "1":"22"), "chrX", "chrY")
  
  observeEvent(event_data(event = "plotly_click", source = "b"), {
    updateSelectizeInput(session, "adj_gene", selected = adj_by_gene[adj_by_gene$m == adj_d(),]$gene[1])
    updateSelectizeInput(session, "adj_chr", selected = chromosome_names[filter(adj_gl, start == adj_d() | end == adj_d())$chromosome[1]])
  })
  
  observeEvent(input$adj_chr, {
    updateSelectizeInput(session, "adj_gene", selected = 
                           ifelse(adj_by_gene[adj_by_gene$gene == input$adj_gene,"chromosome"][1] == input$adj_chr, input$adj_gene, ""))
  })
  
  adjusted_cn <- eventReactive(input$adj_gene,{
    if(nchar(input$adj_gene) > 1){
      paste0(input$adj_gene, " (", round(filter(adj_by_gene, gene == input$adj_gene)$start/1e6,2), "M", "-", round(filter(adj_by_gene, gene == input$adj_gene)$end/1e6,2), "M", ")", ": ", 
             round(filter(adj_by_gene, gene == input$adj_gene)$cn,2), 
             ifelse(filter(adj_by_gene, gene == input$adj_gene)$cn == 1, " copy", " copies"),
             ifelse(filter(adj_by_gene, gene == input$adj_gene)$loh == TRUE, ", suspected LOH", ""))
    } else paste0("")
  })
  
  output$adj_copies <- renderText(adjusted_cn())
  
  gene_mutations <- eventReactive(input$adj_gene,{
    filter(mutations, Gene == input$adj_gene)
  })
  
  output$mutations <- DT::renderDataTable(if(nrow(gene_mutations()) > 0){
    return(datatable(gene_mutations(),escape=FALSE, options = list(dom = 't', columnDefs = list(list(className = 'dt-center', targets = c(1:5)))))) 
  } else return(data.frame()))
  
  adj_seg <- reactive({ filter(adj_gl, chromosome == gsub("adj_", "", input$adj_chr), start == adj_d() | end == adj_d()) })
  seg_start <- reactive({ adj_seg()$start[1] })
  seg_end <- reactive({ adj_seg()$end[1] })
  adj_seg_gene_dat <- reactive({ 
    if(nrow(adj_seg()) > 0){
      return(filter(adj_by_gene, chromosome == input$adj_chr & start >= seg_start() & end <= seg_end())) 
    } else return(data.frame())
  }) 
  
  seg_cytoband_dat <- reactive({
    filter(cytoband_data, chrom == input$adj_chr & seg_start() < chromEnd & seg_end() > chromStart)
  })

  adj_test2 <- reactive({ nrow(adj_seg_gene_dat()) > 0 })
  
  adj_plot_check <- reactive({ adj_d() >= 1 })
  
  output$adj_selected_plot <- renderPlotly({
    req(adj_plot_check() & input$adj_gene == "")
    if(adj_test2() & input$adj_chr != "all"){
      plot_ly(height = 250, type = 'scatter', mode = 'markers') %>%
        add_trace(x = adj_seg_gene_dat()$m, 
                  y = adj_seg_gene_dat()$mean_log2, 
                  hoverinfo = 'text',
                  text = paste0(adj_seg_gene_dat()$gene, " (", adj_seg_gene_dat()$cn, " copies)"),
                  marker = list(color= adj_seg_gene_dat()$marker_color,
                                size = 3*log(adj_seg_gene_dat()$total_weight+1),
                                line = list(color = adj_seg_gene_dat()$outline_color)),
                  showlegend = F) %>%
        add_segments(x = adj_seg()$start[1], xend = adj_seg()$end[1], y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        add_segments(x = adj_seg()$start[1], xend = adj_seg()$end[1], y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        layout(title = paste0("segment: ", adj_seg()$start[1], "-", adj_seg()$end[1]), 
               xaxis = list(range = c(min(seg_cytoband_dat()$chromStart), max(adj_seg_gene_dat()$end)),
                                    ticktext = as.list(seg_cytoband_dat()$name), tickvals = as.list(seg_cytoband_dat()$chromStart),
                                    tickfont = list(size = 8), tickmode = "array", tickangle = 270, side = "top"), 
               margin = list(t = 80),
               yaxis=list(tickfont = list(size = 6), title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
    } 
  })
  
  ##karyotype file
  output$karyotype <- downloadHandler(
    filename = paste0(sample, "_karyotype.pdf"),
    content = function(file) {
      karyo_data <- adj_cns %>% select(chr, start, end, cn, loh)
      karyo_data$cn <- ifelse(karyo_data$cn > 6, 6, karyo_data$cn)
      granges <- GenomicRanges::makeGRangesFromDataFrame(karyo_data, keep.extra.columns = TRUE, ignore.strand = TRUE)
      pdf(file)
      kp <- karyoploteR::plotKaryotype("hg38", plot.type = 2)
      CopyNumberPlots::plotCopyNumberCalls(kp, cn.calls = granges, labels = "", label.cex = 0, cn.colors = "red_blue", loh.color = "green")
      graphics::legend(x = "bottomright", fill = c("#EE0000", "#FFC1C1", "#E0E0E0", "#B2DFEE", "#87CEFA", "#1E90FF", "#0000FF", "green"),
                       legend = c("0", "1", "2", "3", "4", "5", "6+", "loh"), title = "copies", bty = "n")
      dev.off()
    })
  
}
