path <-"../../cnv_viz_demo/rg_viz"
test_files <- list.files(path)
cnr_filename <- paste0(path, "/", test_files[158])
cns_filename <- paste0(path, "/", test_files[159])
#cnr_filename <- "CPDV183182_20179_D1_20007_B_20184_200206_S4.final.cnr"
#cns_filename <- "CPDV183182_20179_D1_20007_B_20184_200206_S4.final.cns"

function(input, output, session) {
  
  output$karyotype <- renderText({
    return(paste('<iframe style="height:600px; width:100%" src="', karyotype_filename, '"></iframe>', sep = ""))
  })
  
  cnr <- read.table(file = cnr_filename, sep = '\t', header = TRUE)
  cnr_target <- filter(cnr, !gene %in% c("Antitarget", "."))
  cnr_target <- mutate(cnr_target, m_probe = (start + end)/2)

  by_gene <- cnr_target %>% group_by(chromosome, gene) %>%
    summarise(s = min(start), e = max(end),
              mean_log2 = weighted.mean(log2, weight),
              total_weight = sum(weight))%>%
    mutate(m = (s+e)/2, cn = round((2^mean_log2)*2))
  
  cns <- read.table(cns_filename, header = TRUE)
  gene_ranges <- c()
  for(i in 1:nrow(cns)) { 
    if(unlist(strsplit(cns$gene[i], ","))[1] != "-") {
      gene_ranges <- c(gene_ranges, paste0(unlist(strsplit(cns$gene[i], ","))[1], " - ", tail(unlist(strsplit(cns$gene[i], ",")),1))) }
    else gene_ranges <- c(gene_ranges, "-")
  }
  cns$gene_ranges <- gene_ranges
  gl <- filter(cns, log2 < -1 | log2 > 1)
  
  chromosomes <- c(paste0("chr", "1":"22"), "chrX", "chrY")
  for(i in c(1:length(chromosomes))){
    
    chr <- filter(by_gene, chromosome == chromosomes[i])
    assign(chromosomes[i], chr)
    
    chr_gl <- filter(gl, chromosome == chromosomes[i])
    assign(paste0(chromosomes[i], "_gl"), chr_gl)
    
    colors <- ifelse(get(chromosomes[i])$mean_log2 > 1 | get(chromosomes[i])$mean_log2 < -1, "blue", "white")
    
    plot <- plot_ly(type = 'scatter', mode = 'markers') %>%
      add_trace(x = get(chromosomes[i])$m, 
                y = get(chromosomes[i])$mean_log2, 
                text = paste0(get(chromosomes[i])$gene),
                hoverinfo = 'text',
                marker = list(line = list(color = "green"), color=colors, size = 3*log(get(chromosomes[i])$total_weight+1)),
                showlegend = F) %>%
      add_segments(x = 0, xend = max(get(chromosomes[i])$m), y = -1, yend = -1, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      add_segments(x = 0, xend = max(get(chromosomes[i])$m), y = 1, yend = 1, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      layout(annotations = list(x = 40e6 , y = 6, text = chromosomes[i], showarrow= F), yaxis=list(title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)), xaxis= list(range = c(0,250e6)))
    
    if(nrow(get(paste0(chromosomes[i], "_gl")))>0){
      for(j in 1:nrow(get(paste0(chromosomes[i], "_gl")))){
        plot <- plot %>% add_segments(x = get(paste0(chromosomes[i], "_gl"))$start[j], xend = get(paste0(chromosomes[i], "_gl"))$end[j], y = get(paste0(chromosomes[i], "_gl"))$log2[j], yend = get(paste0(chromosomes[i], "_gl"))$log2[j], 
                                      line = list(color = "red", width = 6), showlegend = F,
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
  
  #trying to add rectangle highlighter
  #first get input gene start and stop location
  # gene_start <- eventReactive(input$gene, { by_gene[by_gene$gene == input$gene,"s"] })
  # gene_stop <- eventReactive(input$gene, { by_gene[by_gene$gene == input$gene, "e"] })
  # final_plot <- eventReactive(input$gene, { plot_todisplay() %>% add_segments(y = -3, yend = 6, x = gene_start(), yend = gene_end(), line = list(opacity = 0.2, color = "purple")) })
  # output$chr_plot <- renderPlotly({ final_plot() })
  
  output$chr_plot <- renderPlotly(plot_todisplay())
  
  observe({
    updateSelectInput(session, "chr", selected = as.character(by_gene[by_gene$gene==input$gene,]$chromosome[1]))
  })
  
  d <- reactive ({ event_data("plotly_click")[[3]] })
  
  observeEvent(event_data("plotly_click"), {
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
                    hoverinfo = 'text',
                    marker = list(color='yellow', size = gene_data()$weight*15),
                    showlegend = F) %>%
          add_segments(x = min(gene_data()$start), xend = max(gene_data()$end), y = -1, yend = -1, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
          add_segments(x = min(gene_data()$start), xend = max(gene_data()$end), y = 1, yend = 1, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
          layout(title = paste0(gene_data()$gene[1], " (", cn(), " ", copy(),")"), xaxis = list(tickfont = list(size = 6), range = min(gene_data()$s), max(gene_data()$e)), yaxis=list(tickfont = list(size = 6), title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
      } else if(test2() > 0 & input$chr != "all"){
          plot_ly(height = 250, type = 'scatter', mode = 'markers') %>%
            add_trace(x = seg_gene_dat()$m_probe, 
                      y = seg_gene_dat()$log2, 
                      hoverinfo = 'text',
                      text = seg_gene_dat()$gene,
                      marker = list(color='red', size = seg_gene_dat()$weight*15),
                      showlegend = F) %>%
            add_segments(x = min(seg_gene_dat()$start), xend = max(seg_gene_dat()$end), y = -1, yend = -1, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
            add_segments(x = min(seg_gene_dat()$start), xend = max(seg_gene_dat()$end), y = 1, yend = 1, line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
            layout(title = paste0("segment"), xaxis = list(tickfont = list(size = 6), range = min(seg_gene_dat()$s), max(seg_gene_dat()$e)), yaxis=list(tickfont = list(size = 6), title = "Copy number ratio", titlefont = list(size = 8), range = c(-3, 6)))
      } 
  }) 
  
}
