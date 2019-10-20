# Install packages
pacman::p_load("tidyverse")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the PATRIC wastewater isolates
patric <- read_csv("data/PATRIC_genome_wastewater.csv") %>%
  janitor::clean_names() %>%
 
head(patric)
write_csv(patric, "data/PATRIC_complete_genomes_wastewater.csv")
dim(patric)


patric <- read_csv("data/PATRIC_genome_wastewater.csv") %>%
  janitor::clean_names() %>%
  dplyr::filter(genome_status == "Complete") %>%
  dplyr::filter(!is.na(ref_seq_accessions))
write_csv(patric, "data/PATRIC_complete_RefSeq_genomes_wastewater.csv")
