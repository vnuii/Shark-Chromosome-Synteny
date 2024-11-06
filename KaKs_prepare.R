# This script calculates Ks values according to the literature: https://genome.cshlp.org/content/33/9/1527.full
# It performs the following steps:
# 1. Extracts proteins located on chromosomes using SonicParanoid in a conda environment.
# 2. Aligns protein sequences using MAFFT in a BUSCO environment.
# 3. Extracts corresponding CDS sequences (nucleotides) for reverse translation using gffread.
# 4. Processes single-copy orthologous groups to prepare CDS sequences for trimming.
# 5. Trims CDS sequences based on the length of the original AA sequences.
# 6. Uses TrimAl to trim aligned protein sequences and backtranslate to CDS sequences.
# 7. Filters out sequences with internal stop codons and those shorter than 100bp.
library(data.table)
library(dplyr)
library(Biostrings)
library(rtracklayer)
library(parallel)

# extract proteins located on chromosomes, including unlocalized proteins
fna_files <- list.files("/home/vonui/odp_shark/repre_sp_data/onlychrom", pattern = ".fna", full.names = TRUE)
for (i in fna_files) {
  i <- fna_files[9]
  if (i == "/home/vonui/odp_shark/repre_sp_data/Callorhinchus_milii.fna") {
    next
  }
  sp <- basename(i) %>% gsub("_chrom.fna", "", .)
  fna <- readDNAStringSet(i)
  gff <- import(paste0("/home/vonui/odp_shark/repre_sp_data/", sp, ".gff")) %>% as.data.table()
  pep <- readAAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/", sp, ".pep"))
  names(pep) <- gsub(" .+", "", names(pep))
  fna_names <- names(fna) %>% gsub(" .+", "", .)
  gff_need <- gff[seqnames %in% fna_names & type == "CDS", c("seqnames", "type", "ID")]
  gff_need

  gff_need$ID <- gsub("cds-", "", gff_need$ID)
  # names(pep) <- gsub("rna-", "", names(pep))
  pep_need <- pep[names(pep) %in% gff_need$ID]
  pep
  pep_need

  writeXStringSet(pep_need, paste0("/home/vonui/odp_shark/repre_sp_data/inchrom_unlocal_protein/", sp, "_chrom.pep"))
}

# mafft
## mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)
singlecopy_files <- list.files("/home/vonui/odp_shark/ka_ks/orthofinder/Single_Copy_Orthologue_Sequences", pattern = ".fa", full.names = TRUE)
mclapply(singlecopy_files, function(i) {
  sp <- basename(i) %>% gsub(".fa", "", .)
  system(paste0("/home/vonui/anaconda3/envs/busco/bin/linsi --anysymbol ", i, " > /home/vonui/odp_shark/ka_ks/mafft/", sp, ".aln"))
}, mc.cores = 20)

# extract CDS sequences for corresponding genes (nucleotides) for reverse translation
## gffread -w /home/vonui/odp_shark/repre_sp_data/CDS/all_cds.fasta -g /home/vonui/odp_shark/repre_sp_data/Callorhinchus_milii.fna /home/vonui/odp_shark/repre_sp_data/Callorhinchus_milii.gff
rp_sps <- list.files("/home/vonui/odp_shark/repre_sp_data", pattern = "pep", full.names = TRUE)
rp_sps <- rp_sps[!rp_sps == "/home/vonui/odp_shark/repre_sp_data/Callorhinchus_milii.pep"]
rp_sps <- rp_sps[!rp_sps == "/home/vonui/odp_shark/repre_sp_data/Deania_profundorum_edited.pep"]
for (i in rp_sps) {
  sp <- basename(i) %>% gsub(".pep", "", .)
  system(paste0("/home/vonui/anaconda3/envs/busco/bin/gffread -x /home/vonui/odp_shark/repre_sp_data/CDS/", sp, "_cds.fasta -g /home/vonui/odp_shark/repre_sp_data/", sp, ".fna /home/vonui/odp_shark/repre_sp_data/", sp, ".gff"))
}

# prepare CDS sequences for trimming
OG_dt <- fread("/home/vonui/odp_shark/ka_ks/orthofinder/Orthogroups/Orthogroups.tsv")
singlecopy_OGs <- fread("/home/vonui/odp_shark/ka_ks/orthofinder/Orthogroups/Orthogroups_SingleCopyOrthologues.txt", header = FALSE) %>% pull(V1)
OG_dt <- OG_dt[Orthogroup %in% singlecopy_OGs]
spnames <- list.files("/home/vonui/odp_shark/repre_sp_data/CDS", pattern = ".fasta", full.names = FALSE) %>% gsub("_cds.fasta", "", .)
OG_files <- list.files("/home/vonui/odp_shark/ka_ks/mafft", pattern = ".aln", full.names = TRUE)
for (i in spnames) {
  i <- spnames[8]
  sp_singlecopy_genes <- OG_dt[[i]]
  cds <- readDNAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/CDS/", i, "_cds.fasta"))
  gff <- import(paste0("/home/vonui/odp_shark/repre_sp_data/", i, ".gff")) %>% as.data.table()
  protein <- readAAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/inchrom_unlocal_protein/", i, "_chrom.pep"))
  # names(cds) <- gsub("rna-", "", names(cds))
  gff_need <- gff[Parent %in% names(cds) & type == "CDS", c("seqnames", "type", "ID", "Parent")]
  gff_need$Parent <- unlist(gff_need$Parent)
  gff_need <- unique(gff_need)
  gff_need$ID <- gsub("cds-", "", gff_need$ID)
  # 二选一
  sum(!sp_singlecopy_genes %in% gsub("rna-", "", names(cds)))
  names(cds) <- gsub("rna-", "", names(cds))
  cds <- cds[names(cds) %in% sp_singlecopy_genes]
  cds
  writeXStringSet(cds, paste0("/home/vonui/odp_shark/ka_ks/CDS_SingleCopy/", i, "_cds.fasta"))
  # 二选一
  sum(!sp_singlecopy_genes %in% gff_need$ID)
  names(cds) <- gff_need$ID[match(names(cds), gff_need$Parent)]

  cds <- cds[names(cds) %in% sp_singlecopy_genes]
  cds
  writeXStringSet(cds, paste0("/home/vonui/odp_shark/ka_ks/CDS_SingleCopy/", i, "_cds.fasta"))
}

# delete the stop codon at the end
sp_names <- list.files("/home/vonui/odp_shark/ka_ks/CDS_SingleCopy_copy", pattern = ".fasta", full.names = TRUE)

for (i in sp_names) {
  sp <- basename(i) %>% gsub("_cds.fasta", "", .)
  cds <- readAAStringSet(i)
  pep <- readAAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/inchrom_unlocal_protein/", sp, "_chrom.pep"))
  pep <- pep[names(pep) %in% names(cds)]
  pep <- pep[match(names(cds), names(pep))]
  need_check <- names(cds)[width(cds) > width(pep) * 3]
  cds[need_check] <- narrow(cds[need_check], start = 1, end = width(cds[need_check]) - 3)
  sum(!width(cds) == width(pep) * 3)
  cds <- narrow(cds, start = 1, end = width(cds) - 3)
  writeXStringSet(cds, i)
}

sp_names <- list.files("/home/vonui/odp_shark/ka_ks/CDS_SingleCopy_copy", pattern = ".fasta", full.names = TRUE)
i <- sp_names[5]
sp <- basename(i) %>% gsub("_cds.fasta", "", .)
cds <- readAAStringSet(i)
# pep <- readAAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/inchrom_unlocal_protein/", sp, "_chrom.pep"))
pep <- readAAStringSet("/home/vonui/odp_shark/ka_ks/all_alignments.aln")
pep <- pep[names(pep) %in% names(cds)]
pep <- pep[match(names(cds), names(pep))]
# pep_2 <- readAAStringSet(paste0("/home/vonui/odp_shark/repre_sp_data/inchrom_unlocal_protein/", sp, "_chrom.pep"))
# pep_2 <- pep_2[names(pep_2) %in% names(cds)]
# pep_2 <- pep_2[match(names(cds), names(pep_2))]
need_check <- names(cds)[width(cds) > width(pep) * 3]

cds <- narrow(cds, start = 1, end = width(cds) - 3)
names(cds)[!names(cds) %in% names(pep)]
writeXStringSet(cds, i)

system("cat /home/vonui/odp_shark/ka_ks/CDS_SingleCopy_copy/*.fasta >> /home/vonui/odp_shark/ka_ks/CDS_SingleCopy_copy/all_cds.fasta")
trimmed_sequences <- readDNAStringSet("/home/vonui/odp_shark/ka_ks/CDS_SingleCopy_copy/all_cds.fasta")
process_file <- function(i) {
  OG_name <- basename(i) %>% gsub(".aln", "", .)
  single_protein <- readAAStringSet(i)
  cds_need <- trimmed_sequences[names(trimmed_sequences) %in% names(single_protein)]

  output_path <- paste0("/home/vonui/odp_shark/ka_ks/CDS4trimal/", OG_name, "_cds.fasta")
  writeXStringSet(cds_need, output_path)
}
mclapply(OG_files, process_file, mc.cores = 20)

## trimal
## trimal -in aligned_proteins.fasta -out trimmed_proteins.fasta -automated1 -backtrans cds_sequences.fasta
## trimal -in trimmed_cds.fasta -out final_trimmed_cds.fasta -nogaps
mclapply(OG_files, function(i) {
  OG_name <- basename(i) %>% gsub(".aln", "", .)
  trimal_commad1 <- paste0("/home/vonui/anaconda3/envs/busco/bin/trimal -in ", i, " -out /home/vonui/odp_shark/ka_ks/trimal_output/", OG_name, "_trimmed_cds.fasta -automated1 -backtrans /home/vonui/odp_shark/ka_ks/CDS4trimal/", OG_name, "_cds.fasta")
  trimal_commad2 <- paste0("/home/vonui/anaconda3/envs/busco/bin/trimal -in /home/vonui/odp_shark/ka_ks/trimal_output/", OG_name, "_trimmed_cds.fasta -out /home/vonui/odp_shark/ka_ks/trimal_output/", OG_name, "_final_trimmed_cds.fasta -nogaps")
  system(trimal_commad1)
  system(trimal_commad2)
}, mc.cores = 20)

# system(paste0("rm /home/vonui/odp_shark/ka_ks/trimal_output/*"))

# delete sequences with internal stop codons and those shorter than 100bp
completed_files <- list.files("/home/vonui/odp_shark/ka_ks/trimal_output", pattern = "final_trimmed_cds.fasta", full.names = F) %>% gsub("_final_trimmed_cds.fasta", "", .)
OG_names <- basename(OG_files) %>% gsub(".aln", "", .)
OG_rest <- OG_names[!OG_names %in% completed_files] # 有些文件没有完成即中间有终止子的
for (i in completed_files) {
  cds <- readDNAStringSet(paste0("/home/vonui/odp_shark/ka_ks/trimal_output/", i, "_final_trimmed_cds.fasta"))
  if (any(width(cds) > 100)) {
    system(paste0("cp /home/vonui/odp_shark/ka_ks/trimal_output/", i, "_final_trimmed_cds.fasta /home/vonui/odp_shark/ka_ks/trimal_output_larger_100bp/", i, "_final_trimmed_cds.fasta"))
  }
}
mclapply(completed_files, function(i) {
  cds <- readDNAStringSet(paste0("/home/vonui/odp_shark/ka_ks/trimal_output/", i, "_final_trimmed_cds.fasta"))
  if (any(width(cds) > 100)) {
    system(paste0("cp /home/vonui/odp_shark/ka_ks/trimal_output/", i, "_final_trimmed_cds.fasta /home/vonui/odp_shark/ka_ks/trimal_output_larger_100bp/", i, "_final_trimmed_cds.fasta"))
  }
}, mc.cores = 20)
