# This script processes chromosome data to calculate various statistics such as length, GC content,
# and gene-related metrics for different species. The script is divided into two main parts:
#
# 1. Calculating Chromosome Length and GC Content:
#    - The script reads chromosome sequences from .fna files located in the specified directory.
#    - For each species (excluding "Squalus_acanthias"), it calculates the total length and GC content
#      of each chromosome.
#    - The results are saved in a tab-separated file for each species.
#
# 2. Calculating Gene-related Metrics for Chromosomes:
#    - The script reads pre-calculated chromosome statistics from a file.
#    - Chromosomes are classified into micro, middle, and macro categories based on their size.
#    - For each species, the script reads gene annotations from a .gff file and calculates the median
#      gene length and gene number for each chromosome.
#    - The gene density (number of genes per Mb) is also calculated.
#    - The updated statistics are saved back to the file.
library(data.table)
library(dplyr)
library(Biostrings)
library(rtracklayer)
library(parallel)

# Calculate chromosome length and GC content
chrom_files <- list.files("/home/vonui/odp_shark/repre_sp_data/onlychrom", pattern = ".fna", full.names = TRUE)

for (i in chrom_files) {
  sp <- basename(i) %>% gsub(".fna", "", .)
  if (sp == "Squalus_acanthias") {
    next
  }
  sequences <- readDNAStringSet(i)
  seq_names <- names(sequences) %>% gsub(",", "", .)
  corre_dt <- data.table(seq_id = gsub(" .*", "", seq_names), chrom_id = gsub(".+(chr[a-z]* [^\\ ]+).+", "\\1", seq_names))
  mclapply(unique(corre_dt$chrom_id), function(x) {
    seq_length <- sum(width(sequences[grepl(paste0("\\b", x, "\\b"), seq_names)]))
    seq_gc <- letterFrequency(sequences[grepl(paste0("\\b", x, "\\b"), seq_names)], letters = c("G", "C")) %>% sum() / sum(width(sequences[grepl(paste0("\\b", x, "\\b"), seq_names)]))
    data.table(chrom_id = x, seq_size = seq_length / 1e6, seq_gc = seq_gc, species = sp)
  }, mc.cores = 10) %>% rbindlist() -> stats_dt_tmp
  fwrite(stats_dt_tmp, paste0("/home/vonui/odp_shark/repre_sp_data/chrom_stats/", sp, "_stats.tsv"), sep = "\t", quote = F, col.names = T)
}



# Calculate gene-related metrics for chromosomes
# identify micro-, middle-, macro-chromosomes based on size boundaries of 30 Mb and 70 Mb
# ref https://genome.cshlp.org/content/33/9/1527.full
library(ggplot2)
stats_dt <- fread("/home/vonui/odp_shark/repre_sp_data/represent_8sp_chromStats.tsv")
stats_dt$chrom_type <- cut(stats_dt$seq_size, breaks = c(0, 30, 70, Inf), labels = c("micro", "middle", "macro"))
stats_dt$median_gene_length <- 0
stats_dt$gene_number <- 0

# seq_size Mb, gene_length Kb
for (i in unique(stats_dt$species)) {
  #   i <- unique(stats_dt$species)[8]
  sp <- gsub("_chrom", "", i)
  gff <- import(paste0("/home/vonui/odp_shark/repre_sp_data/", sp, ".gff")) %>% as.data.table()
  fna <- readDNAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/onlychrom/", sp, "_chrom.fna"))
  names(fna)

  chromo_id <- names(fna) %>% gsub(" .*", "", .)
  chromo_names <- names(fna) %>%
    gsub(".+(chromosome [^ ]+) .+", "\\1", .) %>%
    gsub(",", "", .)
  chromo_namesDT <- data.table(id = chromo_id, name = chromo_names)
  gff_chrom <- gff[seqnames %in% chromo_id & type == "gene", c("seqnames", "width", "type")]
  gff_chrom$chrom_name <- chromo_namesDT[match(gff_chrom$seqnames, chromo_namesDT$id), "name"]
  gff_chrom %>%
    group_by(chrom_name) %>%
    summarise(median_gene_length = median(width) / 1e3, gene_number = n()) %>%
    ungroup() -> gff_chrom_stats
  stats_dt[species == i]$median_gene_length <- gff_chrom_stats[match(stats_dt[species == i]$chrom_id, gff_chrom_stats$chrom_name), "median_gene_length"]
  stats_dt[species == i]$gene_number <- gff_chrom_stats[match(stats_dt[species == i]$chrom_id, gff_chrom_stats$chrom_name), "gene_number"]
  stats_dt[species == i]
}

stats_dt$gene_density <- stats_dt$gene_number / stats_dt$seq_size # per Mb
fwrite(stats_dt, "/home/vonui/odp_shark/repre_sp_data/represent_8sp_chromStats.tsv", sep = "\t", quote = F, col.names = T)
