# This script prepares input files for the PAML yn00 program.
# It performs the following steps:
# 1. Loads necessary libraries: dplyr, data.table, Biostrings, and parallel.
# 2. Lists all .fasta files in the specified directory.
# 3. Extracts species names from the CDS_SingleCopy directory and formats them.
# 4. Defines a function `process_fun` to process each .fasta file:
#    a. Extracts the orthologous group (OG) name from the file name.
#    b. Creates a directory for the OG.
#    c. Reads the DNA sequences from the .fasta file.
#    d. Renames the sequences using the species names.
#    e. Writes the sequences to a .pml file in the OG directory.
#    f. Creates a control (.ctl) file for the PAML yn00 program with the necessary parameters.
# 5. Uses `mclapply` to apply the `process_fun` function to all .fasta files in parallel using 20 cores.
library(dplyr)
library(data.table)
library(Biostrings)
library(parallel)

OG_files <- list.files("/home/vonui/odp_shark/ka_ks/trimal_output_larger_100bp", pattern = ".fasta", full.names = T)
sp_names_order <- list.files("/home/vonui/odp_shark/ka_ks/CDS_SingleCopy", full.names = F) %>%
  gsub("._cds.fasta", "", .) %>%
  gsub("(\\w).+_(\\w\\w).+", "\\1\\2", .)
process_fun <- function(i) {
  OG_name <- gsub("_final_trimmed_cds.fasta", "", basename(i))
  system(paste0("mkdir /home/vonui/odp_shark/ka_ks/PAML/all_InAndOut/", OG_name))
  cds <- readDNAStringSet(i)
  names(cds) <- sp_names_order
  seq_names <- names(cds)
  seqs <- as.character(cds)
  nseq <- length(seqs)
  seqlen <- unique(width(cds))
  header <- sprintf("%d %d", nseq, seqlen)
  seq_lines <- sprintf("%s  %s", seq_names, seqs)
  writeLines(c(header, seq_lines), con = paste0("/home/vonui/odp_shark/ka_ks/PAML/all_InAndOut/", OG_name, "/", OG_name, ".pml"))
  # ctl文件
  seqfile <- paste0("seqfile = /home/vonui/odp_shark/ka_ks/PAML/all_InAndOut/", OG_name, "/", OG_name, ".pml")
  outfile <- paste0("outfile = /home/vonui/odp_shark/ka_ks/PAML/all_InAndOut/", OG_name, "/", OG_name, ".yn00")
  ctl <- c(
    seqfile,
    outfile,
    "verbose = 0",
    "icode = 0",
    "weighting = 0",
    "commonf3x4 = 0"
  )
  writeLines(ctl, con = paste0("/home/vonui/odp_shark/ka_ks/PAML/all_InAndOut/", OG_name, "/", OG_name, ".ctl"))
}

mclapply(OG_files, process_fun, mc.cores = 20)
