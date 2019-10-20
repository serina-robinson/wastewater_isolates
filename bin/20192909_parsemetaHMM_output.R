# Install packages
pacman::p_load("tidyverse", "readxl", "Biostrings", "DECIPHER")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the metaHMM jobs
job452 <- readAAStringSet("data/MetaHMM_jobs/job_452_a3aff3d2c5a3efbb9012461ea0f1dc48/output_unaligned.fasta")
names(job452) <- paste0(names(job452), "_PF05118_sludge")
length(job452) # 14 alpha keto-glutarate sequences

job453 <- readAAStringSet("data/MetaHMM_jobs/job_453_e7b98007034a7a223d031e012346f14c/output_unaligned.fasta")
head(job453)
names(job453) <- paste0(names(job453), "_PF02461_sludge")
length(job453) # no AMOs/PMOs

job454 <- readAAStringSet("data/MetaHMM_jobs/job_454_8e261066b8257aced38efb583a9acd71/output_unaligned.fasta")
names(job454) <- paste0(names(job454), "_PF00848_sludge")
length(job454) # 52 ring hydroxyl A

job455 <- readAAStringSet("data/MetaHMM_jobs/job_455_762ce71357d905baad49fbbbe8e2c907/output_unaligned.fasta")
names(job455) <- paste0(names(job455), "_PF03443_sludge")
length(job455) # no LPMOs present

job456 <- readAAStringSet("data/MetaHMM_jobs/job_456_8304afe4a3d5f45b2b5519ff0513adfc/output_unaligned.fasta")
names(job456) <- paste0(names(job456), "_PF00667_sludge")
length(job456) # 90 cytochrome p450s

job457 <- readAAStringSet("data/MetaHMM_jobs/job_457_2cd968913af4cfa9963af1c8e0f511de/output_unaligned.fasta")
names(job457) <- paste0(names(job457), "_PF02332_sludge")

job458 <- readAAStringSet("data/MetaHMM_jobs/job_458_b826a1252d141eb3e78f8920e86ded88/output_unaligned.fasta")
names(job458) <- paste0(names(job458), "_PF05138_sludge")

job459 <- readAAStringSet("data/MetaHMM_jobs/job_459_7a60f028f5f6e857a1872a0c4afa523a/output_unaligned.fasta")
names(job459) <- paste0(names(job459), "_PF01494_sludge")

job461 <- readAAStringSet("data/MetaHMM_jobs/job_461_cf5c27a9176f62ecfc4c0995e2d5a38f/output_unaligned.fasta")
names(job461) <- paste0(names(job461), "_PF03241_sludge")

# Combine everything
comb <- c(job452, job453, job454, job455, job456, job457, job458, job459, job461)
length(comb) # 156 sequences
comb_long <- comb[width(comb) > 80]
length(comb_long)


# Read in the 2-ketoglutarate
kg <- readAAStringSet("~/Downloads/PF05118_seed.txt")
names(kg) <- paste0(names(kg), "_PF05118")
length(kg)

# AMO
amo <- readAAStringSet("~/Downloads/PF02461_full.txt")
names(amo) <- paste0(names(amo), "_PF02461")
length(amo)

# Rieske
rieske <- readAAStringSet("~/Downloads/PF00848_seed.txt")
names(rieske) <- paste0(names(rieske), "_PF00848")
length(rieske)

# LPMO
lpmo <- readAAStringSet("~/Downloads/PF03443_seed.txt")
names(lpmo) <- paste0(names(lpmo), "_PF03443")
uniprot <- word(names(lpmo), sep = "_", 1)
uniprot
uniprot6 <- data.frame(uniprot[nchar(uniprot) == 6]) 
colnames(uniprot6) <- "uniprot6_1"
unip_fin <- uniprot6 %>%
  dplyr::mutate(uniprot3 = paste0(substr(uniprot6_1, 1, 2))) %>%
  distinct(uniprot3, .keep_all = T) 
nrow(unip_fin)
lpmo <- lpmo[grepl(paste0(unip_fin$uniprot6_1, collapse = "|"),  names(lpmo))]
length(lpmo)

# CYP p450
cyp450 <- readAAStringSet("~/Downloads/PF00067_seed.txt")
names(cyp450) <- paste0(names(cyp450), "_PF00067")
length(cyp450)

hpa <- readAAStringSet("~/Downloads/PF03241_seed.txt") 
names(hpa) <- paste0(names(hpa), "_PF03241")
length(hpa) # 147

gly <- readAAStringSet("~/Downloads/PF00903_seed (1).txt")
names(gly) <- paste0(names(gly), "_PF00903")
length(gly)

paag <- readAAStringSet("~/Downloads/PF05138_seed.txt")
names(paag) <- paste0(names(paag), "_PF05138")
length(paag)

fad <- readAAStringSet("~/Downloads/PF01494_seed (1).txt")
names(fad) <- paste0(names(fad), "_PF01494")
length(fad)

phenol <- readAAStringSet("~/Downloads/PF02332_seed.txt")
names(phenol) <- paste0(names(phenol), "_PF02332")
length(phenol)

pfams <- c(kg, amo, rieske, lpmo, cyp450, hpa, gly, paag, fad, phenol)

comb_everything <- c(comb_long, pfams)
comb_dedup <- comb_everything[!duplicated(comb_everything)]
length(comb_dedup)
table(word(names(comb_dedup), sep = "_", -1))

names(comb_dedup) <- gsub("\\/", "_", names(comb_dedup))
names(comb_dedup) <- gsub("-", "_", names(comb_dedup))
tail(names(comb_dedup))
summary(width(comb_dedup))
numseqs <- length(comb_dedup)
numseqs
writeXStringSet(comb_dedup, "output/1143_oxygenases_for_ssn.fasta")

######### Max limit submitted to MetaHMM on 09/29/2019
# http://pitgroup.org/metahmm/




