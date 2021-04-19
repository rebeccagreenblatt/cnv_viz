getPlots <- function(by_gene, gl, cytoband_data, ploidy, gene_list, chr_choice, adjusted){
  green <- "#009E73" # normal
  blue <- "#0072B2" # out of range genes
  pink <- "#CC79A7" #mutation
  orange <- "#D55E00" #segment
  black <- "#000000" #LOH

  if(length(gene_list) < 200){
    marker_size <- 3.25
  } else marker_size <- 2
  
  chromosomes <- c(paste0("chr", "1":"22"), "chrX", "chrY")
  for(i in c(1:length(chromosomes))){

    chr <- filter(by_gene, chromosome == gsub("", "", chromosomes[i]), gene %in% gene_list)
    assign(chromosomes[i], chr)

    chr_gl <- filter(gl, chromosome == gsub("", "", chromosomes[i]))
    assign(paste0(chromosomes[i], "_gl"), chr_gl)

    if(adjusted == TRUE){
      seg_colors <- ifelse(get(paste0(chromosomes[i],"_gl"))$loh == TRUE, black, orange)
      source = "b"
      } 
    else {
      seg_colors <- orange
      source = "a"
      }

    out_of_range <- ifelse(get(chromosomes[i])$cn > 64, " - log-2 value outside range of y axis", "")

    plot <- plot_ly(source = source, type = 'scatter', mode = 'markers') %>%#,
      add_trace(x = get(chromosomes[i])$m,
                y = get(chromosomes[i])$mean_log2,
                text = paste0(get(chromosomes[i])$gene, " (", round(get(chromosomes[i])$cn,2), " copies)", out_of_range),
                hoverinfo = 'text',
                marker = list(color = get(chromosomes[i])$marker_color,
                              line = list(color = get(chromosomes[i])$outline_color),
                              size = marker_size*log(get(chromosomes[i])$total_weight+1)),
                showlegend = F) %>%
      add_segments(x = 0, xend = max(c(get(chromosomes[i])$m,0)),
                   y = log((ploidy-0.5)/2, 2), yend = log((ploidy-0.5)/2, 2), line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      add_segments(x = 0, xend = max(c(get(chromosomes[i])$m,0)),
                   y = log((ploidy+0.5)/2, 2), yend = log((ploidy+0.5)/2, 2), line = list(color = "gray", width = 1, dash = "dot"), showlegend = F) %>%
      layout(yaxis=list(
        title = "log(2) copy number ratio", titlefont = list(size = 8), range = c(-3.5, 6)))

    if(nrow(get(paste0(chromosomes[i], "_gl")))>0){
      for(j in 1:nrow(get(paste0(chromosomes[i], "_gl")))){
        plot <- plot %>% add_segments(
          x = get(paste0(chromosomes[i], "_gl"))$start[j], xend = get(paste0(chromosomes[i], "_gl"))$end[j],
          y = get(paste0(chromosomes[i], "_gl"))$log2[j], yend = get(paste0(chromosomes[i], "_gl"))$log2[j],
          color = I(seg_colors[j]),
          text = paste0("segment: ", get(paste0(chromosomes[i], "_gl"))$start[j], "-", get(paste0(chromosomes[i], "_gl"))$end[j]),
          hoverinfo = 'text',
          line = list(width = 3),
          showlegend = F)
      }
    }

    cytoband_chrom <- filter(cytoband_data, chrom == gsub("", "", chromosomes[i]))

    subplot <- plot %>% layout(
      annotations = list(x = 40e6 , y = 6, text = gsub("", "", chromosomes[i]), showarrow= F),
      xaxis = list(range = c(0, 250e6), dtick = 100e6), yaxis = list(range(-3,6)))    #%>%

    plot <- plot %>%
      add_trace(x = cytoband_chrom$chromStart, y = 6, xaxis = 'x2', showlegend = F, marker = list(size = 0.1), hoverinfo = 'skip') %>%
      layout(title = gsub("", "", chromosomes[i]),
             xaxis = list(range = c(0, max(cytoband_chrom$chromEnd)), zeroline = TRUE, showline = TRUE),
             xaxis2 = list(range = c(0, max(cytoband_chrom$chromEnd)),
                           ticktext = as.list(cytoband_chrom$name), tickvals = as.list(cytoband_chrom$chromStart),
                           tickfont = list(size = 8), tickmode = "array", tickangle = 270, side = "top",
                           overlaying = 'x', zeroline = TRUE, autorange = FALSE, matches = 'x'),
             margin = list(t = 80))

    assign(paste0(chromosomes[i],"_plot"), plot)
    assign(paste0(chromosomes[i], "_subplot"), subplot)

  }

  all_plot <- subplot(
    chr1_subplot, chr2_subplot, chr3_subplot,
    chr4_subplot, chr5_subplot, chr6_subplot,
    chr7_subplot, chr8_subplot, chr9_subplot,
    chr10_subplot, chr11_subplot, chr12_subplot,
    chr13_subplot, chr14_subplot, chr15_subplot,
    chr16_subplot, chr17_subplot, chr18_subplot,
    chr19_subplot, chr20_subplot, chr21_subplot,
    chr22_subplot, chrX_subplot, chrY_subplot,
    nrows=8, shareY = TRUE, shareX = TRUE) %>% 
    layout(autosize = F, height = 1200)

  plot_todisplay <- get(paste0(gsub("_adj", "", chr_choice), "_plot"))
  return(plot_todisplay)

}
