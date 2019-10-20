make_cluster<-function(shrt, thresh1, thresh2, bynum, namvec, pal2) { # input is an all vs. all table 
  # #Pull out columns to make a sequence similarity network

  #First remove all comparisons between identical proteins (closed loop nodes)
  
  noprs<-shrt[!shrt$prot1==shrt$prot2,]
  noprs<-noprs[order(noprs$eval),]
  
  source("src/color_mibig_ANL.r")
  meta<-color_mibig(noprs, namvec, pal2)
  
  #Make a unique data frame 
  metu<-unique(meta)
  
  for(i in seq(from = thresh1, to = thresh2, by = bynum)) {
    
  # Set a similarity threshold 
  thresh<-i
  net <- noprs[noprs$eval<(as.numeric(paste0("1.00e-",thresh))),]
  
  metu2 <- metu[order(metu$size),]
  # write_csv(net, "output/graph_net_test.csv")
  #Simplify the graph
  #gr<-simplify(graph.data.frame(net,vertices=metu2,directed = FALSE))
  g <- simplify(graph.data.frame(net, vertices=metu2, directed = FALSE))
  g <- delete.vertices((g), degree(g)<1)
  
  #Append metadata
  V(g)$color <- V(g)$color
  V(g)$size <- V(g)$size
  
  l <- layout_components(g)
  }
  
  ##Make the graph

  # tiff(file = paste0("output/",nrow(net),"_OleC_only_network_e-",thresh,"_300dpi.tiff"), units = "in", width = 10, height = 10, res = 300)
  # jpeg(file = paste0("output/", nrow(net),"_antismashv2_CreM_predictions_network_try3_e-",thresh,"_1000px_labeled.jpg"), width = 2000, height = 2000)
  # # pdf(file = paste0("output/",nrow(net),"_OleC_only_network_e-",thresh,".pdf"),width = 20, height = 20)
  # par(mar=c(1,1,1,1))
  #   pl<-plot(g, vertex.label = NA, 
  #            # vertex.label=ifelse(grepl("_MACS",V(g)$name), V(g)$name, NA), 
  #            layout=l, edge.color="gray40", edge.width=0.3)
  #   title(main = paste0("BLAST e-value cut-off: 1e-",thresh))
  # dev.off()
  # }
  
  #data<-toVisNetworkData(g)
  
  #visNetwork(nodes=data$nodes,edges=data$edges,height="500px") %>%
  #visEdges(smooth = FALSE) %>%
  #visPhysics(stabilization = FALSE) %>%
  #visSave(file = paste0("Biuret_hydrolase_network_e-",thresh,".html"), background = "white")
  
  # Make a legend
  colors<-unique(metu$color)
  proteins<-unique(metu$fam)
  clegend <- data.frame(colors, proteins)
  
  # Plot the legend
  pdf(file=paste0("output/", nrow(net), "_legend_functional_class_final.pdf"), width = 5, height = 10)
  plot.new()
  legend("bottomright",legend=proteins, fill = colors, bty="n") #horiz = T, nrow = 2)
  dev.off() 
  
  return(list(net=net,g=g, legend = clegend))
}
