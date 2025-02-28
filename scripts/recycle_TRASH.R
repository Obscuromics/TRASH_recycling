
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("dplyr"))
library("pheatmap")
library(dendextend)
library("ggplot2")
library("gplots")
library("RColorBrewer")
library("glue")
library("stringr")

source('scripts/helper_functions.R') # import functions
#source('/Users/ab66/Documents/sanger_work/Tools/TRASH_recycling/scripts/helper_functions.R') # import functions

parser <- ArgumentParser()
parser$add_argument("-i", "--input", help = "TRASH output summary file")
parser$add_argument("-c", "--chrom", help = "Chromosome file")
parser$add_argument("--min", type = "integer", help = "Minimum repeat length")
parser$add_argument("--max", type = "integer", help = "Maximum repeat length")
parser$add_argument("-p", "--prefix", default = "sp", help = "Species prefix for the output files")
parser$add_argument("-o", "--out_dir", default = "summary_files", 
                    help = "Output directory for summary files")
args <- parser$parse_args()

# read arguments
trash_summary <- args$input
#trash_summary <- "/Users/ab66/Documents/sanger_work/Tools/TRASH_recycling/data/Summary.of.repetitive.regions.bcop.fa.csv"
chrom_file <- args$chrom
#chrom_file <- "/Users/ab66/Documents/sanger_work/Tools/TRASH_recycling/data/idBraCopr2.core.chroms.tsv"
plot.min <- args$min
#plot.min = 30
plot.max <- args$max
#plot.max = 300
prefix <- args$prefix
#prefix = "Bcop"

# output directories
#NOTE: create directories for output files
#summary_dir <- "/Users/ab66/Documents/sanger_work/Tools/TRASH_recycling/summary_files"
summary_dir <- args$out_dir

# process chromosomes lengths
chr_lengths_f <- read.csv(chrom_file, sep = "\t", header = TRUE)
chr_lengths <- list()
for(grp in unique(chr_lengths_f$group)){
  chr_grp <- chr_lengths_f[chr_lengths_f$group == grp,]
  chr_lengths [[grp]] <- c(chr_grp$length)
  names(chr_lengths [[grp]]) <- chr_grp$chrom
}

# read TRASH summary file and extract main peaks
chr_repeats_l <- list()
chr_peaks <- list()
all_repeats <- c()
for(grp in names(chr_lengths)){
  chr_repeats_l [[grp]] <- read_in_summary_file(file = trash_summary, chr_lengths = chr_lengths[[grp]])
  chr_peaks [[grp]] <- generate_summary_tbl(
    repeats_l = chr_repeats_l[[grp]], chr_lengths = chr_lengths[[grp]], 
    plot.min = plot.min, plot.max = plot.max,
    prefix = paste0(prefix, "_", grp), out_dir = summary_dir, chr_fraction = 0.1)
  all_repeats <- append(all_repeats, chr_peaks [[grp]]$repeat_length)
}

print("Processing satellite repeats with lengths...")
for(grp in names(chr_peaks)){
  print(paste(grp, paste(chr_peaks[[grp]]$repeat_length, collapse = ", "), sep = ": "))
}

# Extract consensus sequences for repeats of all lengths from all chromosomes
# need to check if repeats of the same lengths but from different chromosomes have the same consensus
all_repeats <- unique(all_repeats)
consensus_sequences <- NULL
for(grp in names(chr_repeats_l)){
  grp_consensuses <- extract_consensus_copies(
    repeats_l = chr_repeats_l[[grp]], repeats_list = all_repeats, prefix = prefix)
  consensus_sequences <- rbind(consensus_sequences, grp_consensuses)
}

# Compare repeat sequences of the same length with each other using jaccard similarity index
consensus_list <- consensus_sequences$sequence
names(consensus_list) <- consensus_sequences$rep_id
kmers <- sapply(consensus_list, generate_kmers, k = 6)
jaccard_m_all <- list()
for(rep in unique(consensus_sequences$rep_length)){
  consensus_flt <- consensus_sequences[consensus_sequences$rep_length == rep,]$rep_id
  m_size <- length(consensus_flt)
  jaccard_m <- matrix(nrow = m_size, ncol = m_size)
  colnames(jaccard_m) <- consensus_flt
  rownames(jaccard_m) <- consensus_flt
  
  for(seq_a in colnames(jaccard_m)){
    for(seq_b in rownames(jaccard_m)){
      a_index <- match(seq_a, colnames(jaccard_m))
      b_index <- match(seq_b, rownames(jaccard_m))
      ab_dist <- jaccard_similarity(seq_a, seq_b, kmers = kmers)
      jaccard_m[a_index, b_index] <- ab_dist
    }
  }
  jaccard_m_all[[rep]] <- jaccard_m
}

# generate heatmaps for each repeat group individuallys
for(rep in names(jaccard_m_all)){
  plt <- generate_heatmap(jaccard_m_all[[rep]])
  ggsave(plot = plt, 
         filename = paste(summary_dir, paste(prefix, rep, "jaccard_similarity.png", sep = "_"), 
                          sep = "/"), device = "png")
}

# extract good clusters (median jaccard similarity across the whole cluster is higher than 0.65)
median_jaccard_similarities <- sapply(jaccard_m_all, median)
clusters_flt <- names(median_jaccard_similarities[median_jaccard_similarities >= 0.65])
good_cluster_names <- names(jaccard_m_all[names(jaccard_m_all) %in% clusters_flt])

# check if these clusters have outliers and filter them out
jaccard_m_all_good <- list()
jaccard_m_all_good_rep_ids <- vector()
for(rep in good_cluster_names){
  jaccard_m_all_good[[rep]] <- filter_outliers(m = jaccard_m_all[[rep]])
  jaccard_m_all_good_rep_ids <- append(jaccard_m_all_good_rep_ids, rownames(jaccard_m_all_good[[rep]]))
  # generate a new heatmap
  plt <- generate_heatmap(jaccard_m_all_good[[rep]])
  ggsave(plot = plt, 
         filename = paste(summary_dir, paste(prefix, rep, "jaccard_similarity.f.png", sep = "_"), 
                          sep = "/"), device = "png")
}

# predict an optimal number of clusters for the rest based on the biggest height drop
other_clusters_to_process <- names(jaccard_m_all)[!names(jaccard_m_all) %in% good_cluster_names]

for(rep in other_clusters_to_process){
  m <- jaccard_m_all[[rep]]
  clusters <- calculate_n_clusters(m)
  for(j in unique(clusters)){
    rep_ids <- names(clusters[clusters == j])
    cluster_id <- paste0(rep, "_", j)
    
    # take only the clusters that have a median jaccard similarity more than 0.5
    # and discard the rest
    if(length(rep_ids) >= 5 & median(m[rep_ids, rep_ids]) >= 0.5){
      m_f <- filter_outliers(m = m[rep_ids, rep_ids])
      jaccard_m_all_good[[cluster_id]] <- m_f
      jaccard_m_all_good_rep_ids <- append(jaccard_m_all_good_rep_ids, 
                                           rownames(jaccard_m_all_good[[cluster_id]]))
      # generate a new heatmap
      plt <- generate_heatmap(jaccard_m_all_good[[cluster_id]])
      ggsave(plot = plt, 
             filename = paste(summary_dir, paste(prefix, cluster_id, "jaccard_similarity.f.png", sep = "_"), 
                              sep = "/"), device = "png")
      
    }else if(length(rep_ids) > 1 & median(m[rep_ids, rep_ids]) < 0.5){
      # check again if there are any subclusters?
      # generate a new heatmap
      plt <- generate_heatmap(m[rep_ids, rep_ids])
      ggsave(plot = plt, 
             filename = paste(summary_dir, paste(prefix, cluster_id, "jaccard_similarity.d.png", sep = "_"), 
                              sep = "/"), device = "png")
    }
  }
}

print(paste0("Number of sequences left: ", length(jaccard_m_all_good_rep_ids)))

# calculate jaccard similarity for all repeats of all lengths together
m_size <- length(jaccard_m_all_good_rep_ids)
jaccard_m_filtered <- matrix(nrow = m_size, ncol = m_size)
colnames(jaccard_m_filtered) <- jaccard_m_all_good_rep_ids
rownames(jaccard_m_filtered) <- jaccard_m_all_good_rep_ids
  
for(seq_a in colnames(jaccard_m_filtered)){
  for(seq_b in rownames(jaccard_m_filtered)){
    a_index <- match(seq_a, colnames(jaccard_m_filtered))
    b_index <- match(seq_b, rownames(jaccard_m_filtered))
    ab_dist <- jaccard_similarity(seq_a, seq_b, kmers = kmers)
    jaccard_m_filtered[a_index, b_index] <- ab_dist
  }
}

plt <- generate_heatmap(jaccard_m_filtered, fontsize = 1)
ggsave(plot = plt, 
       filename = paste(summary_dir, paste(prefix, "jaccard_similarity.all.png", sep = "_"), 
                        sep = "/"), device = "png", units = "cm", width = 40, height = 40, dpi = 900)

## plotting jaccard and kmer based distances (PCOA plot)
jaccard_dist <- 1 - jaccard_m_filtered
jaccard_pcoa <- cmdscale(jaccard_dist, k=2, eig=TRUE, add=TRUE)
jaccard_positions <- jaccard_pcoa$points
colnames(jaccard_positions) <- c("pcoa1", "pcoa2")
rownames(jaccard_positions) <- colnames(jaccard_dist)

# plot
percent_explained <- 100 * jaccard_pcoa$eig / sum(jaccard_pcoa$eig)
pretty_pe <- format(round(percent_explained[1:2], digits=2), nsmall=1, trim=TRUE)

labs <- c(glue("PCOA1 ({pretty_pe[1]}%)"), glue("PCOA2 ({pretty_pe[2]}%)"))

jaccard_positions <- jaccard_positions %>% as_tibble(rownames = "repeats")

jaccard_positions["genome"] <- 
  paste(str_split_fixed(jaccard_positions$repeats, "_", 4)[,2],
        str_split_fixed(jaccard_positions$repeats, "_", 4)[,3], sep = "_")

jaccard_positions$repeat_length <- sapply(str_split(jaccard_positions$repeats, '_'), tail, 1)
jaccard_positions$annot <- paste(jaccard_positions$genome, jaccard_positions$repeat_length, sep = "_")

jaccard_dist_plt <- jaccard_positions %>%
  ggplot(aes(x = pcoa1, y = pcoa2, col = repeat_length, shape = genome)) +
  geom_point(size = 3, alpha = 0.5) +
  labs(x=labs[1], y=labs[2]) +
  geom_text(data = jaccard_positions %>% filter(repeat_length %in% c(
    "len185", "len94", "len155", "len156", "len118", "len129")),
    aes(label = annot), size = 2,
    nudge_x = 0.1, nudge_y = 0.1, 
    check_overlap = T)

jaccard_dist_plt

ggsave(plot = jaccard_dist_plt, 
       filename = paste(summary_dir, paste(prefix, "jaccard_dist.all.png", sep = "_"), 
                        sep = "/"), device = "png", units = "cm", width = 20, height = 25)
