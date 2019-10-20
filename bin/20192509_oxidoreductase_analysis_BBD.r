# Install packages
pacman::p_load("tidyverse")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the Eawag BBD data
bbd <- read_excel("data/EawagBBD_reactions.xlsx", col_names = T) %>%
  janitor::clean_names() %>%
  dplyr::mutate(bbd_class = word(ec_code, sep = "\\.", 1)) %>%
  dplyr::filter(bbd_class != "EC Code")


# Find all oxidoreductases
prop.table(table(bbd$bbd_class))

# Find all oxygenases
length(grep("oxygenase", bbd$enzyme_name))
240/(table(bbd$bbd_class)[1])
oxygenases <- bbd[grep("oxygenase", bbd$enzyme_name),]
write_csv(oxygenases, "output/EawagBBD_oxygenases.csv")
