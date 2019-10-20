# fastafile, hmmfile, outfile, size_threshold
# Install packages
pacman::p_load("data.table", "tidyverse", "Biostrings","DECIPHER","igraph","RColorBrewer", "officer", "rvg")

# Set working directory
setwd("~/Documents/EAWAG/github/wastewater_isolates/")

# Read in the all-vs-all BLAST
allvall <- fread("data/1143_oxygenases_for_ssn.tsv", stringsAsFactors = F, data.table = F)
head(allvall)


namvec <- c(levels(as.factor(word(allvall$V1, sep = "_", -1))))
length(namvec)

pal <- c("#92d050", "#E41A1C", "goldenrod", "#984EA3", "#FF7F00", "blue1", "#8B4513", "#F781BF", 'cornflowerblue', 'seagreen4')
length(pal)
pal2 <- c(pal[1:10], "gray") 

length(namvec)
length(pal2)

source("lib/make_cluster_diagram.r")
shrt <- allvall[,c(1,2,11)]
head(shrt)

colnames(shrt) <- c("prot1", "prot2", "eval")

# gr <- make_cluster(shrt, thresh1=36, thresh2=36, bynum=1, namvec, pal2)
noprs <- shrt[!shrt$prot1 == shrt$prot2,]
noprs <- noprs[order(noprs$eval), ]

source("lib/color_mibig_ANL.r")
meta <- color_mibig(noprs, namvec, pal2)
metu <- unique(meta)

thresh <- 1

if (thresh >= 1) {
  net <- noprs[noprs$eval < (as.numeric(paste0("1.00e-", thresh))),]
}

if (thresh < 1) {
  net <- noprs[noprs$eval < as.numeric(thresh),]
}

g <- simplify(graph.data.frame(net, vertices = metu, directed = FALSE))
# g <- delete.vertices((g), degree(g)<1) # delete singletons

#Append metadata
V(g)$color <- V(g)$color
V(g)$size <- V(g)$size

l <- layout_components(g)


pdf(file = paste0("output/20191013_", nrow(net),"_oxygenases_SSN_e-",thresh,".pdf"), width = 20, height = 20)
# jpeg(file = paste0("output/",nrow(net),"_OleC_only_network_e-",thresh,"_1000px.jpg"),width = 1000, height = 1000)
# tiff(file = paste0("output/",nrow(net),"_AMP_binding_network_e-",thresh,"_300dpi_MIBiG_outliers_try2.tiff"), units = "in", width = 10, height = 10, res = 300)
par(mar=c(1,1,1,1))
pl<-plot(g, #vertex.label=ifelse(grepl("_OleC",V(g)$name), V(g)$name, NA),
         vertex.label = NA, vertex.label.size = 0.1,
         #vertex.label = ifelse(grepl("NRPS", V(g)$name), V(g)$name, NA),
         layout=l, edge.color="gray40", edge.width=0.1)
# title(main = paste0("BLAST e-value cut-off: 1e-",thresh))
dev.off() 

leg <- read_excel("data/pfam_key.xlsx")
# proteins <- unique(meta$fam)
proteins <- as.character(paste0(leg$name, ' (', leg$pfam, ')'))
proteins[1] <- gsub ("\\(sludge)", "", proteins[1])
colors<-unique(meta$color)
colors

pdf(file=paste0("output/",nrow(net),"_legend_sludge_metagenome.pdf"), width = 20, height = 5)
plot.new()
legend("bottomright",legend=proteins,fill = colors, bty="n")
dev.off()

