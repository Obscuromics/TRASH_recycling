
# Read TRASH summary file
read_in_summary_file <- function(file, chr_lengths, plot.min = 30, plot.max = 300){
  repeats_tbl <- read.csv(file) # read file
  repeats_tbl$most.freq.value.N <- as.integer(repeats_tbl$most.freq.value.N)
  repeats_tbl <- repeats_tbl[repeats_tbl$name %in% names(chr_lengths),]
  repeats_l <- repeats_tbl[repeats_tbl$most.freq.value.N >= plot.min,]
  repeats_l <- repeats_l[repeats_l$most.freq.value.N <= plot.max,]
  return(repeats_l)
}

# generate summary tables (function copied from TRASH)
generate_summary_tbl <- function(repeats_l, chr_lengths, prefix, 
                                 plot.min = 50, plot.max = 300, out_dir = "summary_files",
                                 top = NULL, chr_fraction = NULL){
  # Identify main peaks
  repeat_no <- aggregate(repeats.identified ~ most.freq.value.N, repeats_l, sum)
  total_length <- repeat_no$most.freq.value.N * repeat_no$repeats.identified
  percent_fraction <- total_length*100 / sum(chr_lengths)
  repeat_sum <- cbind(repeat_no, total_length, percent_fraction)
  colnames(repeat_sum)[1] <- "repeat_length"
  colnames(repeat_sum)[2] <- "no_repeats"
  
  peaks <- repeat_sum[repeat_sum$percent_fraction >= min(length(chr_lengths) * 0.01, 0.2),]
  png(file = paste0(out_dir, "/", paste(prefix, "peaks.png", sep = "_")),
      width = 800, height = 600)
  plot(repeat_sum$repeat_length, repeat_sum$no_repeats, type = "h", 
       xlim = c(50, min(max(peaks$repeat_length)*2, plot.max)),
       xlab = "monomer length", ylab = "no of repeats")
  dev.off()
  
  peaks_merge <- peaks[order(peaks$no_repeats, decreasing = TRUE),]
  
  i = 1
  
  tryCatch({
    while(peaks_merge$repeat_length[i] > 0){
      main_peak <- peaks_merge$repeat_length[i]
      j = 1
      for(peak in repeat_sum$repeat_length){
        if(peak == main_peak - 2 | peak == main_peak - 1 | peak == main_peak + 1 | peak == main_peak + 2){
          peaks_merge[i,2] <- peaks_merge[i,2] + repeat_sum[j,2]
          peaks_merge[i,3] <- peaks_merge[i,3] + repeat_sum[j,3]
          peaks_merge[i,4] <- peaks_merge[i,4] + repeat_sum[j,4]
        }
        j = j + 1
      }
      peaks_merge <- subset(peaks_merge, repeat_length != main_peak - 2 & repeat_length != main_peak - 1 &
                              repeat_length != main_peak + 1 & repeat_length != main_peak + 2)
      i = i + 1
    }
  }, error=function(e){})
  
  repeat_merge <- repeat_sum
  for(i in 1:nrow(peaks_merge)){
    repeat_merge <- subset(repeat_merge, repeat_length != peaks_merge[i,1] 
                           & repeat_length != peaks_merge[i,1] - 2 & repeat_length != peaks_merge[i,1] - 1
                           & repeat_length != peaks_merge[i,1] + 1 & repeat_length != peaks_merge[i,1] + 2)
  }
  
  repeat_merge <- rbind(repeat_merge, peaks_merge)
  assemblyName <- prefix
  
  name_col <- rep(assemblyName, nrow(peaks_merge))
  export_peaks <- cbind(name_col, peaks_merge)
  colnames(export_peaks)[1] <- "genome"
  export_peaks <- export_peaks[order(export_peaks$percent_fraction, decreasing = TRUE),]
  rownames(export_peaks) <- NULL
  
  write.csv(export_peaks, file = paste0(out_dir, "/", paste(prefix, "peaks.csv", sep = "_")))
  png(file = paste0(out_dir, "/", paste(prefix, "peaks_m.png", sep = "_")),
      width = 800, height = 600)
  plot(repeat_merge$repeat_length, repeat_merge$no_repeats, type = "h", 
       xlim = c(50, min(max(peaks$repeat_length)*2, 2000)),
       xlab = "monomer length", ylab = "no of repeats")
  dev.off()
  
  # TAKE ONLY TOP REPEATS PER SET OF CHROMOSOMES (IGNORE IF NOT PROVIDED -> TAKE ALL)
  if(!is.null(top)){
    if(length(export_peaks$repeat_length) > top){
      export_peaks <- head(export_peaks, top)
    }
  }
  
  # FILTER PEAKS BASED ON CHROMOSOME FRACTION
  if(!is.null(chr_fraction)){
    export_peaks <- export_peaks %>% filter(percent_fraction >= chr_fraction)
  }
  
  return(export_peaks)
}

# Extract all consensus sequences for satellite repeats of a certain length 
extract_consensus_copies <- function(repeats_l, repeats_list, prefix){
  df <- repeats_l[repeats_l$most.freq.value.N %in% repeats_list,]
  satellites_all = NULL
  for(rep in unique(df$most.freq.value.N)){
    df_f <- df[df$most.freq.value.N == rep,]
    for(i in 1:nrow(df_f)){
      sequence <- df_f[i,]$consensus.primary
      chrom <- df_f[i,]$name
      start <- df_f[i,]$start
      end <- df_f[i,]$end
      rep_name <- paste(prefix, chrom, start, end, paste0("len", rep), sep = "_")
      
      if(!sequence == "none_identified"){
        satellites_all <- rbind(satellites_all, data.frame(
          rep_id = rep_name, rep_length = paste0("len", rep), sequence = sequence))
      }
    }
  }
  return(satellites_all)
}

# Generate list of kmers for a satellite repeat
generate_kmers <- function(sequence, k = 3){
  count = 1
  kmers <- list()
  while((count + k - 1) <= nchar(sequence)){
    kmer <- substr(sequence, count, count + k - 1)
    kmers <- append(kmers, kmer)
    count <- count + 1
  }
  return(unlist(kmers))
}

# Generate jaccard similarity index for a pair of sequences
jaccard_similarity <- function(seq_a, seq_b, kmers){
  intersection <- length(intersect(kmers[[seq_a]], kmers[[seq_b]]))
  union <- length(union(kmers[[seq_a]], kmers[[seq_b]]))
  similarity = intersection / union
  return(similarity)
}

# filter outliers based on median jaccard similarity
filter_outliers <- function(m){
  median_rep_dist <- NULL
  for(rep_id in rownames(m)){
    median_rep_dist <- rbind(median_rep_dist, 
                             data.frame(rep_id = rep_id, median_sim = median(m[rep_id,])))
  }
  boxplot(median_rep_dist$median_sim)
  out <- boxplot.stats(median_rep_dist$median_sim)$out
  if(!length(out) == 0){
    out_ind <- which(median_rep_dist$median_sim %in% c(out))
    #to_filter <- median_rep_dist[out_ind,] %>%
    # filter(median_sim < 0.5)
    to_filter <- median_rep_dist[out_ind,]
    m_index <- match(to_filter$rep_id, rownames(m))
    m_out <- m[-m_index, -m_index]
  }else{
    m_out <- m
  }
  return(m_out)
}

# predict an optimal number of clusters based on the height drop
calculate_n_clusters <- function(m){
  hc <- hclust(dist(m, method = "euclidean"), method = "complete")
  last_heights <- rev(tail(hc$height, 10))
  h_dist_all <- vector()
  h_dist_all <- append(h_dist_all, 0)
  for(i in 1:(length(last_heights)-1)){
    h_dist <- (last_heights[i] - last_heights[i+1])
    h_dist_all <- append(h_dist_all, h_dist)
  }
  i <- match(max(h_dist_all), h_dist_all)
  h <- ((last_heights[i-1] - last_heights[i]) / 2) + last_heights[i]
  clusters <- cutree(hc, h = h)
  return(clusters)
}

# generate a heatmap
generate_heatmap <- function(m, fontsize = 3){
  color <- bluered(100)
  breaks <- seq(0, 1, length.out = 100)
  legend_breaks <- seq(0, 1, length.out = 11)
  legend_labels <- as.character(legend_breaks)
  plt <- pheatmap::pheatmap(m, color = color, breaks = breaks, 
                                    legend_breaks = legend_breaks, legend_labels = legend_labels,
                                    fontsize = fontsize, method = "complete")
  return(plt)
}

