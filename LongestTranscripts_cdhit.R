# The main functions of this script are:
# 1. Extract the longest transcripts on chromosomes.
# 2. Remove duplicate protein sequences using the cd-hit tool.
#
# The main steps include:
# 1. Load necessary R packages: data.table, dplyr, Biostrings, rtracklayer.
# 2. Get the filenames of all species and remove the "_chrom.fna" part from the filenames.
# 3. For each species, perform the following operations:
#    - Read the DNA sequence file on the chromosomes.
#    - Read the GFF file and filter out the CDS (coding sequences) on the chromosomes.
#    - Filter the corresponding protein sequences based on the longest CDS.
#    - Write the filtered DNA sequences and protein sequences to new files.
#    - Write the filtered GFF information to new files.
# 4. Use the cd-hit tool to remove duplicate protein sequences.
#
# Note:
# - The script will skip processing for the species "Squalus_acanthias".
# - When processing protein sequences, if the number of filtered protein sequences does not match the number of unique protein sequence IDs, a warning message will be printed.
# - The path and parameters of the cd-hit tool need to be adjusted according to the actual situation.
library(data.table)
library(dplyr)
library(Biostrings)
library(rtracklayer)

####### extract longest transcripts
sp_names <- list.files("/home/vonui/odp_shark/repre_sp_data/onlychrom", pattern = ".fna") %>% gsub("_chrom.fna", "", .)

for (i in sp_names) {
  if (i == "Squalus_acanthias") {
    next
  }
  chrom_fna <- readDNAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/onlychrom/", i, "_chrom.fna"))
  # chrom_fna <- chrom_fna[grepl("^NC", names(chrom_fna))]
  chrom_names <- names(chrom_fna) %>% gsub(" .*", "", .)
  gff <- import(paste0("/home/vonui/odp_shark/repre_sp_data/", i, ".gff")) %>% as.data.table()
  gff <- gff[seqnames %in% chrom_names & type == "CDS", c("ID", "seqnames", "strand", "start", "end", "width")]
  gff %>%
    group_by(ID) %>%
    filter(width == max(width)) %>%
    ungroup() %>%
    select(-width) -> gff
  gff$ID <- gsub("cds-", "", gff$ID)
  protein_seq_names <- gff$ID
  protein <- readAAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/", i, ".pep"))
  names(protein) <- names(protein) %>% gsub(" .*", "", .)
  # names(protein) <- names(protein) %>% gsub("rna-", "", .)
  protein_filter <- protein[names(protein) %in% protein_seq_names]
  if (length(protein_filter) != length(unique(protein_seq_names))) {
    print(paste(i, length(protein_filter), length(unique(protein_seq_names))))
  }
  writeXStringSet(chrom_fna, paste0("/home/vonui/odp_shark/repre_sp_data/onlychrom/chrom_nounlocal_4odp/", i, "_chrom.fna"))
  writeXStringSet(protein_filter, paste0("/home/vonui/odp_shark/repre_sp_data/inchrom_protein/", i, "_chrom.pep"))
  fwrite(gff, paste0("/home/vonui/odp_shark/repre_sp_data/chrom_files_new/", i, ".chrom"), sep = "\t", quote = F, col.names = F)
}

####### cdhit
pep_names <- list.files("/home/vonui/odp_shark/repre_sp_data/inchrom_nounlocal_protein_longest", pattern = ".pep", full.names = TRUE)
pep_names <- pep_names[pep_names != "/home/vonui/odp_shark/repre_sp_data/Deania_profundorum.pep"]

for (i in pep_names) {
  cdhit <- paste0("~/anaconda3/envs/odp/bin/cd-hit", " -i ", i, " -o ", "/home/vonui/odp_shark/repre_sp_data_cdhit99/", basename(i), " -c 0.99 -n 5 -M 0 -T 40")
  system(cdhit)
}
