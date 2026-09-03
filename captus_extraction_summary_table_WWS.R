#JJ 15/1/26
#WW 25/6/26
# Script to produce summary table of base pair and reference length recovery

#set working directory
setwd("your/working/directory")

setwd("C:/Users/arist/OneDrive/Desktop/Tesis2/captus_final")
getwd()

#install packages
install.packages("tidyverse")
install.packages("dplyr")

library(tidyverse)
library(dplyr)
library(readr)

#remember to delete first two lines in the tsv file
stats <- read_tsv(file = 'captus-extract_stats.tsv')


summary_df <- stats %>%
  group_by(sample_name) %>%
  summarise(
    n_genes_recovered = n_distinct(locus),
    total_base_recovered = sum(ref_len_matched, na.rm = TRUE) * 3 ,
    .groups = "drop"
  )

write_delim(summary_df, file = 'summary_samples_captus_extract.tsv', delim = '\t')
