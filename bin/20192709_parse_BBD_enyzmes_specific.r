# Install packages
pacman::p_load("tidyverse", "readxl")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the Eawag BBD data
bbd <- read_excel("data/EawagBBD_reactions.xlsx", col_names = T) %>%
  janitor::clean_names() %>%
  dplyr::mutate(bbd_class = word(ec_code, sep = "\\.", 1)) %>%
  dplyr::filter(bbd_class != "EC Code")  %>%
  dplyr::filter(bbd_class == 1)
head(bbd)

# Limit to oxygenases?
# oxygenases <- bbd
oxygenases <- bbd[grep("oxygenase", bbd$enzyme_name),]
# oxygenases$enzyme_name

# Check unique E.C. numbers
table(oxygenases$ec_code)
oxygenases[grepl("-", oxygenases$ec_code),]

no_dashes <- oxygenases %>%
  dplyr::filter(!grepl("\\.-\\.-", ec_code)) %>%
  dplyr::filter(!grepl("unspecific", enzyme_name)) %>%
  dplyr::mutate(unspecific_ec = paste0(word(ec_code, 1, sep = "\\."), ".", word(ec_code, 2, sep = "\\."), ".", word(ec_code, 3, sep = "\\.")))
head(no_dashes)

# Parse the uniprot
uni <- read_delim("~/Documents/EAWAG/github/oxyphen/DATA/ec_uniprot.tsv", delim = "\t") %>%
  dplyr::mutate(unspecific_ec = paste0(word(EC, 1, sep = "\\."), ".", word(EC, 2, sep = "\\."), ".", word(EC, 3, sep = "\\.")))

table(no_dashes$unspecific_ec %in% uni$unspecific_ec) # 5 not present


# Find all associated uniprot IDs
oxy_search <- uni[uni$unspecific_ec %in% no_dashes$unspecific_ec,] %>%
  dplyr::filter(!is.na(uniprot)) %>%
  dplyr::filter(nchar(uniprot) > 8) %>%
  dplyr::filter(!grepl("Transferred", description))


dim(oxy_search)
head(oxy_search)
write_csv(oxy_search, "data/20192909_BBD_oxygenases_only_in_uniprot.csv")


# Pick favorites

