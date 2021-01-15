#path <-"../../cnv_viz_demo/rg_viz"
#test_files <- list.files(path)
#cnr_filename <- paste0(path, "/", test_files[158])
#cns_filename <- paste0(path, "/", test_files[159])
#cnr_filename <- "CPDV183182_20179_D1_20007_B_20184_200206_S4.final.cnr"
#cns_filename <- "CPDV183182_20179_D1_20007_B_20184_200206_S4.final.cns"
cnr_filename <- "../../1218_rg_cnv_viz/CPDV193638.final.cnr"
cns_filename <- "../../1218_rg_cnv_viz/CPDV193638.final.cns"

green <- "#009E73" # normal
blue <- "#0072B2" # out of range genes
pink <- "#CC79A7" #mutation
orange <- "#D55E00" #segment
black <- "#000000" #LOH

library(cBioPortalData)
cbio <- cBioPortal()
cbio_studies <- read.csv("../cbio_study_names.csv")
library(DT)
library(stringr)
#source("../adjusted_panel.R")

function(input, output, session) {
  
  values <- reactiveValues()
  
  output$karyotype <- renderText({
    return(paste('<iframe style="height:600px; width:100%" src="', karyotype_filename, '"></iframe>', sep = ""))
  })
  
  cnr <- read.table(file = cnr_filename, sep = '\t', header = TRUE)
  cnr_target <- filter(cnr, !gene %in% c("Antitarget", "."))
  cnr_target <- mutate(cnr_target, m_probe = (start + end)/2)
  #new
  cnr_target$cn_est <- 2^cnr_target$log2*2
  
  probes <- filter(cnr_target, log2 < -1 | log2 > 0.585) %>% ## log(1/2,2), log(3/2,2)
    group_by(gene) %>% summarise(probe_oor = 1)
  
  by_gene <- cnr_target %>% group_by(chromosome, gene) %>%
    summarise(s = min(start), e = max(end),
              #mean_log2 = weighted.mean(log2, weight),
              ##
              mean_cn_est = weighted.mean(cn_est, weight),
              sd_log2 = sd(log2),
              ##
              total_weight = sum(weight))%>%
    mutate(mean_log2 = log(mean_cn_est/2,2)) %>%
    mutate(m = (s+e)/2, cn = round((2^mean_log2)*2)) %>%
    mutate(blue = as.numeric(mean_log2 < -0.41 | mean_log2 > 0.32)) %>%
    left_join(probes, by = c("gene")) %>% #instead of probe should change to cnr_target? or do we even need to join by probes at all?? 
    #mutate(pink = as.numeric(blue == 0 & probe_oor == 1))
    mutate(pink = ifelse(blue == 0 & sd_log2 > 0.3, 1, 0))
  
  by_gene$mean_log2 <- ifelse(by_gene$mean_log2 < -2, -2, by_gene$mean_log2)
  
  cns <- read.table(cns_filename, header = TRUE)
  gene_ranges <- c()
  for(i in 1:nrow(cns)) { 
    if(unlist(strsplit(cns$gene[i], ","))[1] != "-") {
      gene_ranges <- c(gene_ranges, paste0(unlist(strsplit(cns$gene[i], ","))[1], " - ", tail(unlist(strsplit(cns$gene[i], ",")),1))) }
    else gene_ranges <- c(gene_ranges, "-")
  }
  cns$gene_ranges <- gene_ranges
  gl <- filter(cns, log2 < -0.41 | log2 > 0.32)
  gl$log2 <- ifelse(gl$log2 < -2, -2, gl$log2)
  
  chromosomes <- c(paste0("chr", "1":"22"), "chrX", "chrY")
  for(i in c(1:length(chromosomes))){
    
    chr <- filter(by_gene, chromosome == chromosomes[i])
    assign(chromosomes[i], chr)
    
    chr_gl <- filter(gl, chromosome == chromosomes[i])
    assign(paste0(chromosomes[i], "_gl"), chr_gl)
    
    colors <- ifelse(get(chromosomes[i])$mean_log2 > 0.32 | get(chromosomes[i])$mean_log2 < -0.41, blue, 
                     ifelse(get(chromosomes[i])$pink == 1, "hotpink", "white"))
    colors2 <- ifelse(colors == "white", green, colors) #what?!!!
    
    plot <- plot_ly(source = "a", type = 'scatter', mode = 'markers') %>%
      add_trace(x = get(chromosomes[i])$m, 
                y = get(chromosomes[i])$mean_log2, 
                text = paste0(get(chromosomes[i])$gene),
                hoverinfo = 'text',
                marker = list(color = colors, 
                              line = list(color = green), 
                              #line = list(color = colors2),
                              size = 3*log(get(chromosomes[i])$total_weight+1)),
                showlegend = F) %>%
      add_segments(x = 0, xend = max(get(chromosomes[i])$m), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      add_segments(x = 0, xend = max(get(chromosomes[i])$m), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      layout(annotations = list(x = 40e6 , y = 6, text = chromosomes[i], showarrow= F), yaxis=list(title = "log(2) copy number ratio", titlefont = list(size = 8), range = c(-3, 6)), xaxis= list(range = c(0,250e6)))
    
    if(nrow(get(paste0(chromosomes[i], "_gl")))>0){
      for(j in 1:nrow(get(paste0(chromosomes[i], "_gl")))){
        plot <- plot %>% add_segments(x = get(paste0(chromosomes[i], "_gl"))$start[j], xend = get(paste0(chromosomes[i], "_gl"))$end[j], y = get(paste0(chromosomes[i], "_gl"))$log2[j], yend = get(paste0(chromosomes[i], "_gl"))$log2[j], 
                                      line = list(color = orange, width = 3), showlegend = F,
                                      name = get(paste0(chromosomes[i], "_gl"))$gene_ranges[j])
        #name = get(paste0(chromosomes[i], "_gl"))$gene)
      }
    }
    
    assign(paste0(chromosomes[i],"_plot"), plot)
    
  }
  
  all_plot <- subplot(chr1_plot, chr2_plot, chr3_plot,
                      chr4_plot, chr5_plot, chr6_plot,
                      chr7_plot, chr8_plot, chr9_plot,
                      chr10_plot, chr11_plot, chr12_plot,
                      chr13_plot, chr14_plot, chr15_plot,
                      chr16_plot, chr17_plot, chr18_plot,
                      chr19_plot, chr20_plot, chr21_plot,
                      chr22_plot, chrX_plot, chrY_plot,
                      nrows=9, shareY = TRUE, shareX = TRUE) %>%
    layout(autosize = F, height = 1200)
  
  plot_todisplay <- reactive({ get(paste0(input$chr, "_plot")) })
  
  output$chr_plot <- renderPlotly(plot_todisplay())
  
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
  
  seg <- reactive({ filter(cns, start == d() | end == d()) })
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
                  #hoverinfo = 'text',
                  marker = list(color='yellow', size = 15*gene_data()$weight),
                  showlegend = F) %>%
        add_segments(x = min(gene_data()$start), xend = max(gene_data()$end), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        add_segments(x = min(gene_data()$start), xend = max(gene_data()$end), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        layout(title = paste0(gene_data()$gene[1], " (", cn(), " ", copy(),")"), xaxis = list(tickfont = list(size = 6), range = min(gene_data()$s), max(gene_data()$e)), yaxis=list(tickfont = list(size = 6), title = "log(2) copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
    } else if(test2() > 0 & input$chr != "all"){
    #else if(input$chr != "all"){
      plot_ly(height = 250, type = 'scatter', mode = 'markers') %>%
        add_trace(x = seg_gene_dat()$m_probe, 
                  y = seg_gene_dat()$log2, 
                  hoverinfo = 'text',
                  text = seg_gene_dat()$gene,
                  marker = list(color='red', size = seg_gene_dat()$weight*15),
                  showlegend = F) %>%
        add_segments(x = min(seg_gene_dat()$start), xend = max(seg_gene_dat()$end), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        add_segments(x = min(seg_gene_dat()$start), xend = max(seg_gene_dat()$end), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        layout(title = paste0("segment"), xaxis = list(tickfont = list(size = 6), range = min(seg_gene_dat()$s), max(seg_gene_dat()$e)), yaxis=list(tickfont = list(size = 6), title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
    } 
  })
  
  ######add cbio table
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
  
  ## adjusted
  # get rid of flagged genes?? 
  adj_bygene <- read.csv("../../1218_rg_cnv_viz/CPDV193638_genes.csv") %>% filter(!(gene.symbol %in% c("Antitarget", ".")))
  colnames(adj_bygene)[colnames(adj_bygene) == "gene.symbol"] <- "gene"
  colnames(adj_bygene)[colnames(adj_bygene) == "chr"] <- "chromosome"
  colnames(adj_bygene)[colnames(adj_bygene) == "number.targets"] <- "total_weight"
  adj_bygene$mean_log2 <- round(log(adj_bygene$C/2, 2),2)
  adj_bygene$mean_log2 <- ifelse(adj_bygene$C == 0, -2, adj_bygene$mean_log2)
  adj_bygene$mean_log2 <- ifelse(adj_bygene$mean_log2 > 5, 5, adj_bygene$mean_log2)
  adj_bygene$C <- round(adj_bygene$C)
  
  adj_cns <- read.table("../../1218_rg_cnv_viz/CPDV193638_dnacopy.seg", sep= "\t", header=TRUE)
  colnames(adj_cns)[colnames(adj_cns) == "chrom"] <- "chromosome"
  colnames(adj_cns)[colnames(adj_cns) == "loc.start"] <- "start"
  colnames(adj_cns)[colnames(adj_cns) == "loc.end"] <- "end"
  adj_cns$log2 <- log(adj_cns$C/2, 2)
  adj_cns$log2 <- ifelse(adj_cns$C == 0, -2, adj_cns$log2) 
  adj_cns$C <- round(adj_cns$C)
  
  loh <- read.csv("../../1218_rg_cnv_viz/CPDV193638_loh.csv") %>% select(chr, start, end, type)
  adj_cns <- left_join(adj_cns, loh, by = c("start", "end"))
  adj_cns$loh <- adj_cns$type %in% c("COPY-NEUTRAL LOH", "LOH", "WHOLE ARM COPY-NEUTRAL LOH", "WHOLE ARM LOH")
  
  variants <- read.csv("../../1218_rg_cnv_viz/CPDV193638_variants.csv")
  variants$AR <- round(variants$AR, 2)
  variants$ML.AR <- round(variants$ML.AR, 2)
  variants$ML.AR <- ifelse(variants$ML.C == 0, "N/A", variants$ML.AR) #shouldn't have an allelic fraction if copy number is truly 0...
  mutations <- filter(variants, ML.SOMATIC == TRUE & FLAGGED == FALSE & !(gene.symbol %in% c(NA, "<NA>", "Antitarget"))) %>% 
    #select(gene.symbol, ID, CN.SUBCLONAL, CELLFRACTION, AR.ADJUSTED, ML.LOH)
    select(gene.symbol, ID, depth, AR, ML.AR)
  #colnames(mutations) <- c("Gene", "ID", "Subclonal", "Cell Fraction", "Allelic Fraction", "LOH")
  colnames(mutations) <- c("Gene", "ID", "Depth", "Raw Allelic Fraction", "Allelic Fraction Adjusted for Copy # Change") 
  
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
  
  metadata <- read.csv("../../1218_rg_cnv_viz/CPDV193638.csv")
  metadata$Sex = ifelse(metadata$Sex == TRUE, "MALE", "FEMALE")
  output$meta <- renderTable(metadata %>% select(Purity, Ploidy, Sex, Contamination, Flagged, Failed))
  output$comment <- renderText(ifelse(!is.na(metadata$Comment[1]), paste0("Comment: ", metadata$Comment[1]), ""))
  
  adj_by_gene <- adj_bygene %>% 
    #mutate(blue = as.numeric(mean_log2 < -0.41 | mean_log2 > 0.32)) %>%
    mutate(m = (start+end)/2, cn = C)
  
  max_xranges <- adj_by_gene %>% group_by(chromosome) %>% summarise(max = max(end))
  
  #adj_cns <- adj_cns %>% mutate(log2 = log(C/2,2))
  adj_gl <- filter(adj_cns, log2 < -0.41 | log2 > 0.32 | loh == TRUE)
  #adj_gl <- adj_cns
  
  adj_colors <- ifelse(adj_by_gene$C == 2, "white", blue)
  adj_colors <- ifelse(adj_by_gene$loh == TRUE, black, adj_colors)
  adj_colors <- ifelse(adj_by_gene$mutation == TRUE, pink, adj_colors)
  adj_colors2 <- ifelse(adj_colors == "white", green, adj_colors)
  adj_colors2 <- ifelse(adj_colors == pink & adj_by_gene$loh == TRUE, black, adj_colors2)
  adj_by_gene$marker_color <- adj_colors
  adj_by_gene$outline_color <- adj_colors2
  
  adj_chromosomes <- c(paste0("adj_chr", "1":"22"), "adj_chrX", "adj_chrY")
  for(i in c(1:length(adj_chromosomes))){
    
    adj_chr <- filter(adj_by_gene, chromosome == gsub("adj_", "", adj_chromosomes[i]))
    assign(paste0(adj_chromosomes[i]), adj_chr)
    
    max_xrange <- max_xranges$max[i]
    #can you replace below with just adj_chr isntead of get(... etc)??
    
    adj_chr_gl <- filter(adj_gl, chromosome == gsub("adj_chr", "", adj_chromosomes[i]))
    assign(paste0(adj_chromosomes[i], "_gl"), adj_chr_gl)
    
    seg_colors <- ifelse(get(paste0(adj_chromosomes[i],"_gl"))$loh == TRUE, black, orange)
    
    adj_plot <- plot_ly(source = "b", type = 'scatter', mode = 'markers') %>%
      add_trace(x = get(adj_chromosomes[i])$m, 
                y = get(adj_chromosomes[i])$mean_log2, 
                text = paste0(get(adj_chromosomes[i])$gene, " (", round(get(adj_chromosomes[i])$C,2), " copies)"),
                hoverinfo = 'text',
                marker = list(color = #adj_colors,
                                get(adj_chromosomes[i])$marker_color,
                              #line = list(color = green), 
                              line = list(color = 
                                            #adj_colors2),
                                            get(adj_chromosomes[i])$outline_color),
                              size = 2*log(get(adj_chromosomes[i])$total_weight+1)),
                showlegend = F) %>%
      add_segments(x = 0, xend = max(c(get(adj_chromosomes[i])$m,0)), 
                   y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      add_segments(x = 0, xend = max(c(get(adj_chromosomes[i])$m,0)), 
                   y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      layout(yaxis=list(title = "log(2) copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
    
    if(nrow(get(paste0(adj_chromosomes[i], "_gl")))>0){
      for(j in 1:nrow(get(paste0(adj_chromosomes[i], "_gl")))){
        adj_plot <- adj_plot %>% add_segments(x = get(paste0(adj_chromosomes[i], "_gl"))$start[j], xend = get(paste0(adj_chromosomes[i], "_gl"))$end[j], y = get(paste0(adj_chromosomes[i], "_gl"))$log2[j], yend = get(paste0(adj_chromosomes[i], "_gl"))$log2[j], color = I(seg_colors[j]), 
                                      line = list(
                                        #color = seg_colors[j], 
                                        width = 3), showlegend = F)
        #,name = get(paste0(adj_chromosomes[i], "_gl"))$gene_ranges[j])
      }
    }
    
    cytoband_data <- read.csv("../cytoband_data.csv")
    cytoband_chrom <- filter(cytoband_data, chrom == gsub("adj_", "", adj_chromosomes[i]))
    
    subplot <- adj_plot %>% layout(
      annotations = list(x = 40e6 , y = 6, text = gsub("adj_", "", adj_chromosomes[i]), showarrow= F), 
      xaxis = list(range = c(0, 250e6), dtick = 100e6), yaxis = list(range(-3,6)))
    
    adj_plot <- adj_plot %>% 
      # add_lines(x=c(0, max_xrange),
      #                                y = c(-3, -3),
      #           xaxis = "x2",
      #                               showlegend = FALSE,
      #                              line=list(color = "black"),
      #                             hoverinfo = 'skip') %>%
      layout(
        #xaxis = list(showgrid = FALSE, range = c(0, max_xrange)),
        xaxis = list(scaleanchor = "x", ticktext = as.list(cytoband_chrom$name), tickvals = as.list(cytoband_chrom$chromStart),
                      tickfont = list(size = 8), tickmode = "array", tickangle = 270, side = "top"),
      #domain = c(0, max(get(adj_chromosomes[i])$end))),
      #yaxis=list(title = "log(2) copy number ratio", titlefont = list(size = 8), range = c(-3, 6)),
      margin = list(t = 80))
    
    
    assign(paste0(adj_chromosomes[i],"_plot"), adj_plot)
    assign(paste0(adj_chromosomes[i], "_subplot"), subplot)
    
  }
  
  adj_all_plot <- subplot(adj_chr1_subplot, adj_chr2_subplot, adj_chr3_subplot,
                          adj_chr4_subplot, adj_chr5_subplot, adj_chr6_subplot,
                          adj_chr7_subplot, adj_chr8_subplot, adj_chr9_subplot,
                          adj_chr10_subplot, adj_chr11_subplot, adj_chr12_subplot,
                          adj_chr13_subplot, adj_chr14_subplot, adj_chr15_subplot,
                          adj_chr16_subplot, adj_chr17_subplot, adj_chr18_subplot,
                          adj_chr19_subplot, adj_chr20_subplot, adj_chr21_subplot,
                          adj_chr22_subplot, adj_chrX_subplot, adj_chrY_subplot,
                          nrows=9, shareY = TRUE, 
                          shareX = TRUE) %>%
    layout(autosize = F, height = 1200)
  
  adj_plot_todisplay <- reactive({ get(paste0("adj_", input$adj_chr, "_plot")) })
  
  output$adj_chr_plot <- renderPlotly({ adj_plot_todisplay() })
  
  observe({
    updateSelectInput(session, "adj_chr", selected = as.character(adj_by_gene[adj_by_gene$gene==input$adj_gene,]$chromosome[1]))
  })
  
  adj_d <- reactive ({ event_data(event="plotly_click", source = "b")[[3]] }) 
  
  observeEvent(event_data(event = "plotly_click", source = "b"), {
    updateSelectizeInput(session, "adj_gene", selected = adj_by_gene[adj_by_gene$m == adj_d(),]$gene[1])
  })
  
  observeEvent(input$adj_chr, {
    updateSelectizeInput(session, "adj_gene", selected = 
                           ifelse(adj_by_gene[adj_by_gene$gene == input$adj_gene,"chromosome"][1] == input$adj_chr, input$adj_gene, ""))
  })
  
  #adj_cn <- reactive({ adj_by_gene[adj_by_gene$gene == input$adj_gene,]$cn[1] })
  #adj_copy <- reactive({ ifelse(adj_cn() == 1, "copy", "copies")})
  
  adjusted_cn <- eventReactive(input$adj_gene,{
    paste0(round(filter(adj_by_gene, gene == input$adj_gene)$C,2), ifelse(filter(adj_by_gene, gene == input$adj_gene)$C == 1, " copy", " copies"))
  })
  
  output$adj_copies <- renderText(adjusted_cn())
  
  gene_mutations <- eventReactive(input$adj_gene,{
    filter(mutations, Gene == input$adj_gene)
  })
  
  output$mutations <- DT::renderDataTable(if(nrow(gene_mutations()) > 0){
    return(datatable(gene_mutations(),escape=FALSE,options = list(dom = 't', columnDefs = list(list(className = 'dt-center', targets = c(1:5)))))) 
  } else return(data.frame()))
  
  adj_seg <- reactive({ filter(adj_gl, start == adj_d() | end == adj_d()) })
  seg_start <- reactive({ adj_seg()$start[1] })
  seg_end <- reactive({ adj_seg()$end[1] })
  adj_seg_gene_dat <- reactive({ 
    if(nrow(adj_seg()) > 0){
      return(filter(adj_by_gene, chromosome == input$adj_chr & start >= seg_start() & end <= seg_end())) 
    } else return(data.frame())
  }) 
  #%>%mutate(m = (loc.start+loc.end)/2, log2 = 2*2^C)})
  #seg_dat <- reactive({ data.frame(gene = unique(seg_genes())) })
  #seg_gene_dat <- reactive({ inner_join(seg_dat(), cnr_target, by =c("gene")) })

  adj_test2 <- reactive({ nrow(adj_seg_gene_dat()) > 0 })
  #adj_test2 <- reactive({ nrow(adj_seg()) })
  
  adj_plot_check <- reactive({ adj_d() >= 1 })
  
  output$adj_selected_plot <- renderPlotly({
    req(adj_plot_check())
    if(adj_test2() & input$adj_chr != "all"){
      plot_ly(height = 250, type = 'scatter', mode = 'markers') %>%
        add_trace(x = adj_seg_gene_dat()$m, 
                  y = adj_seg_gene_dat()$mean_log2, 
                  hoverinfo = 'text',
                  text = adj_seg_gene_dat()$gene,
                  marker = list(color=#adj_seg_gene_colors(),
                                  adj_seg_gene_dat()$marker_color,
                                size = 3*log(adj_seg_gene_dat()$total_weight+1),
                                line = list(color = 
                                              #adj_seg_gene_colors2())),
                                              adj_seg_gene_dat()$outline_color)),
                  showlegend = F) %>%
        add_segments(x = min(adj_seg_gene_dat()$start), xend = max(adj_seg_gene_dat()$end), y = -0.41, yend = -0.41, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        add_segments(x = min(adj_seg_gene_dat()$start), xend = max(adj_seg_gene_dat()$end), y = 0.32, yend = 0.32, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
        layout(title = paste0("segment"), xaxis = list(tickfont = list(size = 6), range = min(adj_seg_gene_dat()$start), max(adj_seg_gene_dat()$end)), yaxis=list(tickfont = list(size = 6), title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
    } 
  })
  
}

# get rid of "." gene
# show copy number for unadjusted data on label
