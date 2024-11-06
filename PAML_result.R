# This script processes PAML results to calculate Ka/Ks ratios for different species comparisons.
#
# Steps:
# 1. List all files in the specified directory with the pattern "yn00".
# 2. Filter out files with size 0.
# 3. Create a data table with species IDs and names.
# 4. Initialize an empty data table for storing Ka/Ks results.
# 5. Loop through each filtered file:
#    - Extract the OG name from the file name.
#    - Read the file content.
#    - Extract and clean the header information.
#    - Extract and clean the result data.
#    - Append the results to the KaKs_dt data table.
# 6. Define shallow groups for comparison.
# 7. Loop through each shallow group to:
#    - Filter and process data for deep sea comparisons.
#    - Filter and process data for middle sea comparisons.
# 8. Calculate medians of omega values for deep sea comparisons.
# 9. Write the final processed data to a TSV file.
library(data.table)
library(dplyr)

OG_files <- list.files("/home/vonui/odp_shark/ka_ks/PAML/all_InAndOut", pattern = "yn00", full.names = T, recursive = T)
files_info <- file.info(OG_files)
OG_files_need <- rownames(files_info[files_info$size > 0, ]) # 内部包含终止子的文件，剩下1783 --> 1754
id_dt <- data.table(id = c(1, 2, 3, 4, 5, 6, 7, 8), spname = c("Cca", "Dpr", "Hpe", "Hfr", "Ler", "Sac", "Sti", "Tob"))

KaKs_dt <- data.table(matrix(nrow = 0, ncol = 12))
setnames(KaKs_dt, c("sp1", "sp2", "S", "N", "t", "kappa", "omega", "dN", "dN_SE", "dS", "dS_SE", "OG"))
for (i in OG_files_need) {
  OG_name <- basename(i) %>% gsub(".yn00", "", .)
  yn00 <- readLines(i)
  header <- yn00[146] %>%
    strsplit(" ") %>%
    unlist() %>%
    .[. != "" & . != "+-"]
  header[1:2] <- c("sp1", "sp2")
  header[c(9, 11)] <- c("dN_SE", "dS_SE")
  result <- yn00[148:175] %>% strsplit(" ")
  result <- lapply(result, function(x) x[x != "" & x != "+-"] %>% unname())
  tmp_dt <- data.table(matrix(nrow = 0, ncol = length(header)))
  setnames(tmp_dt, header)
  for (j in result) {
    tmp_dt <- rbind(tmp_dt, as.list(as.numeric(j)))
  }
  # tmp_dt_subset <- tmp_dt[sp1 %in% id_dt$id & sp2 %in% id_dt$id]
  tmp_dt[, OG := OG_name]
  KaKs_dt <- rbind(KaKs_dt, tmp_dt)
}
# fwrite(KaKs_dt, "/home/vonui/odp_shark/ka_ks/PAML/KaKs_dt.tsv", sep = "\t")
shallow_groups <- list(c(1, 3), c(4, 5), c(7, 8))
plot_dt <- data.table()
for (i in shallow_groups) {
  i <- unlist(i)
  dmp <- KaKs_dt[sp1 %in% c(2, i) & sp2 %in% c(2, i)]
  dmp[, sp1 := id_dt[sp1, spname]]
  dmp[, sp2 := id_dt[sp2, spname]]
  dmp$group <- paste0(dmp$sp2, "_", dmp$sp1)
  dmp$type <- "deep_sea"
  dmp$color <- ifelse(dmp$sp1 == "Dpr" | dmp$sp2 == "Dpr", "deep", "shallow")
  dmp$comparison <- paste(id_dt$spname[i], collapse = "_")
  plot_dt <- rbind(plot_dt, dmp)
}
for (i in shallow_groups) {
  i <- unlist(i)
  dmp <- KaKs_dt[sp1 %in% c(6, i) & sp2 %in% c(6, i)]
  dmp[, sp1 := id_dt[sp1, spname]]
  dmp[, sp2 := id_dt[sp2, spname]]
  dmp$group <- paste0(dmp$sp2, "_", dmp$sp1)
  dmp$type <- "middle_sea"
  dmp$color <- ifelse(dmp$sp1 == "Sac" | dmp$sp2 == "Sac", "middle", "shallow")
  dmp$comparison <- paste(id_dt$spname[i], collapse = "_")
  plot_dt <- rbind(plot_dt, dmp)
}

medians <- plot_dt[type == "deep_sea", .(median_omega = median(omega)), by = comparison]
fwrite(plot_dt, "/home/vonui/odp_shark/ka_ks/PAML/plot_dt.tsv", sep = "\t")
