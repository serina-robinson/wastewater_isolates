# Install packages
pacman::p_load("tidyverse", "readxl", "Biostrings")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the 2-ketoglutarate
KG <- readAAStringSet("~/Downloads/PF05118_seed.txt")
uniprot <- word(names(KG), sep = "_", 1)
uniprot6 <- uniprot[nchar(uniprot) == 6]
uniprot6
writeLines(uniprot6, "output/PF05118_alphaKG.txt")

# AMO
amo <- readAAStringSet("~/Downloads/PF02461_full.txt")
uniprot <- word(names(amo), sep = "_", 1)
amo
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
uniprot6
writeLines(uniprot6, "output/PF02461_AMO.txt")

# Rieske
rieske <- readAAStringSet("~/Downloads/PF00848_seed.txt")
uniprot <- word(names(rieske), sep = "_", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
uniprot6
writeLines(uniprot6, "output/PF00848_Rieske.txt")

# LPMO
lpmo <- readAAStringSet("~/Downloads/PF03443_seed.txt")
uniprot <- word(names(lpmo), sep = "_", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
uniprot6
writeLines(uniprot6, "output/PF03443_LPMO.txt")

# CYP p450
cyp450 <- readAAStringSet("~/Downloads/PF00067_rp15.txt")
head(cyp450)
names(cyp450)
uniprot <- word(names(cyp450), sep = "\\/", 1)
uniprot6 <- data.frame(uniprot[nchar(uniprot) == 6]) 
colnames(uniprot6) <- "uniprot6_1"
unip_fin <- uniprot6 %>%
  dplyr::mutate(uniprot3 = paste0(substr(uniprot6_1, 1, 2))) %>%
  distinct(uniprot3, .keep_all = T)
unip_fin
dim(unip_fin)
writeLines(as.character(unip_fin$uniprot6_1), "output/PF00067_CYP450.txt")

######### Max limit submitted to MetaHMM on 09/29/2019
# http://pitgroup.org/metahmm/

# soluble METHANE MONOOXYGENASE
sMMO <- readAAStringSet("~/Downloads/PF02964_uniprot.txt")
uniprot <- word(names(sMMO), sep = "\\/", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
writeLines(uniprot6, "output/PF02964_sMMO.txt")

# PF0866? sMMO not worth submitting
# FAD_binding_3 (PF01494)
# FMO-like (PF00743)


pf <- readAAStringSet("~/Downloads/PF03241_rp15.txt")
head(pf)
uniprot <- word(names(pf), sep = "\\/", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
length(uniprot6)
writeLines(uniprot6, "output/PF03241_accs.txt")

pf <- readAAStringSet("~/Downloads/PF00903_seed (1).txt")
uniprot <- word(names(pf), sep = "_", 1)
uniprot6 <- uniprot[nchar(uniprot) == 6]
length(uniprot6)
writeLines(uniprot6, "output/PF00903_accs.txt")

pf <- readAAStringSet("~/Downloads/PF01494_seed (1).txt")
uniprot <- word(names(pf), sep = "_", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
length(uniprot6)
writeLines(uniprot6, "output/PF01494_accs.txt")

pf <- readAAStringSet("~/Downloads/PF05138_seed.txt")
uniprot <- word(names(pf), sep = "_", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
length(uniprot6)
writeLines(uniprot6, "output/PF05138_accs.txt")

pf <- readAAStringSet("~/Downloads/PF02332_seed.txt")
uniprot <- word(names(pf), sep = "_", 1)
uniprot
uniprot6 <- uniprot[nchar(uniprot) == 6]
length(uniprot6)
writeLines(uniprot6, "output/PF02332_accs.txt")




