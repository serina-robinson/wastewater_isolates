# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the Eawag BBD data
bbd <- read_excel("data/EawagBBD_reactions.xlsx", col_names = T) %>%
  janitor::clean_names() %>%
  dplyr::mutate(bbd_class = word(ec_code, sep = "\\.", 1)) %>%
  dplyr::filter(bbd_class != "EC Code") 
head(bbd)

# Limit to monooxygenases?
oxygenases <- bbd[grep("oxygenase", bbd$enzyme_name),]
oxygenases$enzyme_name

# Check unique E.C. numbers
table(oxygenases$ec_code)
oxygenases[grepl("-", oxygenases$ec_code),]

no_dashes <- oxygenases %>%
  dplyr::filter(!grepl("-", ec_code)) %>%
  dplyr::filter(!grepl("unspecific", enzyme_name))


# Parse the uniprot
uni <- read_delim("~/Documents/EAWAG/github/oxyphen/DATA/ec_uniprot.tsv", delim = "\t")
no_dashes$ec_code %in% uni$EC # all present! woo!

# Find all associated uniprot IDs
oxy_search <- uni[uni$EC %in% no_dashes$ec_code,] %>%
  dplyr::filter(!is.na(uniprot)) %>%
  dplyr::filter(nchar(uniprot) > 8) %>%
  dplyr::filter(!grepl("Transferred entry", description))
dim(oxy_search)
head(oxy_search)
write_csv(oxy_search, "data/BBD_oxygenases_in_uniprot.csv")
