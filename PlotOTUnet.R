#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Wrapper function to create a network with the n most abundant OTUs per sample in a sample by otu table. 
# Abundant taxa will be selected per sample, so that the final number of taxa in the plot may be larger than n.
# abudnance of OTUs per sampl will determine shape of network
# size of vertices determined by total relative abundance of OTUs
#
# input: 
# data - sample (colums) by OTUs (rows) table in sequence counts, row names are identical to rownames(tax)
# tax - OTU by taxon table, taxonomic levels are column names
# file - name of the image file to write the plot to
# abund - number of abundant taxa to plot
# tax.level - taxonomic resolution for legend
# transparent - should the transparency of the vertices be determined by total relative abundance, 
#  default TRUE, if FALSE vertices are not transparent
# clr - should clr-transformed abundance be used to determine the shape of the network
# colorPalette - range of colors for taxonomy legend
# sample.names - names of the sample vertices (default: colnames(data))
# net.layout - default: fruchterman.reingold with edge weights determined by OTU abundance per sample
# keep.layout - save layout table (coordinates of vertices) as R object and as file
#
# output:
# stacked barplot and corresponding legend in a separate window
#
# dependencies:
# The function was designed for usage in windows
# R packages: vegan, reshape, igraph, scales
#=============================================================================
#documentation end

PlotOTUnet <- function (data, tax, abund, 
                        file = "network.png",
                        tax.level = "class", transparent = T, clr = F,
                        colorPalette = c("violet","purple","blue","white","darkgreen","yellow","red","darkred"),
                        sample.names = colnames(data),
                        net.layout = layout.fruchterman.reingold(net.graph, weights = E(net.graph)$weights),
                        keep.layout = T) {
  
  # loading required packages
  if (!"vegan" %in% installed.packages()) {
    install.packages("vegan")
  }
  require(vegan)
  
  if (!"reshape" %in% installed.packages()) {
    install.packages("reshape")
  }
  require(reshape)
  
  if (!"igraph" %in% installed.packages()) {
    install.packages("igraph")
  }
  require(igraph)
  
  if (!"scales" %in% installed.packages()) {
    install.packages("scales")
  }
  require(scales)
  
  # input data
  Data0 <- data
  Tax0 <- tax
  Data.rel0 <- prop.table(as.matrix(Data0), 2) * 100
  tot.rel.abund0 <- rowSums(Data0) / sum(rowSums(Data0)) * 100
  
  # selecting names of abundant OTUs
  abund.names <- c()
  for (i in 1:ncol(Data.rel0)) {
    abund.names <- unique(c(abund.names, rownames(Data.rel0)[order(Data.rel0[, i], decreasing = T)][1:abund]))
  }
  
  # subsetting and sorting data
  # ascending order of total relative OTU abundance to ensure that larger vertices are plotted on top of smaller vertices
  tot.rel.abund <- tot.rel.abund0[abund.names]
  tot.rel.abund <- tot.rel.abund[order(tot.rel.abund)]
  
  Data <- Data0[names(tot.rel.abund), ]
  Tax <- Tax0[names(tot.rel.abund), ]
  
  Data.01 <- decostand(Data, "pa")
  Data.01.otu <- data.frame(otu = rownames(Data.01), Data.01)
  Data.melt.0 <- melt(Data.01.otu)
  Data.melt.1 <- Data.melt.0[Data.melt.0$value > 0, ]
  
  # creating graph
  net.graph <- graph.data.frame(Data.melt.1, directed = F)
  
  # edge weights
  if (clr == F) {
    Data.rel <- Data.rel0[names(tot.rel.abund), ]
    Data.rel.otu <- data.frame(otu = rownames(Data.rel), Data.rel)
    Data.melt.rel0 <- melt(Data.rel.otu)
    Data.melt.rel1 <- Data.melt.rel0[Data.melt.0$value > 0, ]
    E(net.graph)$weights <- rank(Data.melt.rel1$value)
  } else {
    Data.tmp <- Data + 0.5
    Data.clr <- Data.tmp
    for (i in 1:ncol(Data.tmp)) {
      log2gm <- mean(log2(Data.tmp[, i]))
      Data.clr[, i] <- log2(Data.tmp[, i]) - log2gm
    }
    Data.clr.otu <- data.frame(otu = rownames(Data.clr), Data.clr)
    Data.melt.clr0 <- melt(Data.clr.otu)
    Data.melt.clr1 <- Data.melt.clr0[Data.melt.0$value > 0, ]
    E(net.graph)$weights <- rank(Data.melt.clr1$value)
  }
  
  #setting vertex attributes
  Data.vertex <- data.frame(otu = rownames(Data),
                            tot.rel.abund = tot.rel.abund,
                            tax = droplevels(Tax[, colnames(Tax) == tax.level])
                            )
  Data.vertex$color <- Data.vertex$tax
  levels(Data.vertex$color) <- colorRampPalette(colorPalette)(length(levels(Data.vertex$tax)))
  if (transparent == T) {
    Data.vertex$color.alpha <- c()
    for (i in 1:nrow(Data.vertex)) {
      Data.vertex$color.alpha[i] <- adjustcolor(as.character(Data.vertex$color)[i],
                                                alpha.f = rescale(Data.vertex$tot.rel.abund, to = c(0, 1))[i]
                                                )
    }
  } else {
    Data.vertex$color.alpha <- as.character(Data.vertex$color)
  }
  
  # order Data.vertex to match graph object
  Data.vertex <- Data.vertex[match(V(net.graph)$name[-c((nrow(Data) + 1):length(V(net.graph)$name))],
                                   as.character(Data.vertex$otu)), ]
  
  V(net.graph)$size <- c(rescale(Data.vertex$tot.rel.abund, to = c(1, 25)), nchar(sample.names) * 2.4)
  V(net.graph)$size2 <- c(rescale(Data.vertex$tot.rel.abund, to = c(1, 25)), rep(4.5, ncol(Data)))
  V(net.graph)$color <- c(Data.vertex$color.alpha, rep("white", ncol(Data)))
  V(net.graph)$label <- c(rep(NA,nrow(Data)), sample.names)
  V(net.graph)$shape <- c(rep("circle", nrow(Data)), rep("rectangle", ncol(Data)))
  
  # plotting
  png(filename = file, width = 12, height = 9, units = "in", res = 600, type = "windows")
  layout(matrix(c(1, 2), 1, 2), widths = c(3, 1))
  par(mar = c(1, 1, 1, 1))
  plot.igraph(net.graph,
              vertex.size = V(net.graph)$size,
              vertex.size2 = V(net.graph)$size2,
              vertex.label = V(net.graph)$label,
              vertex.label.cex = 0.7,
              vertex.label.color = "black",
              vertex.label.family =  "sans",
              vertex.shape = V(net.graph)$shape,
              edge.color = "lightgrey",
              edge.width = 0.5,
              layout = net.layout
              )
  plot.new()
  par(mar=c(1,1,1,4))
  legend("left",
         pch = 22,
         col = "black",
         pt.bg = colorRampPalette(colorPalette)(length(levels(Data.vertex$tax))),
         legend = levels(Data.vertex$tax),
         pt.cex = 1,
         cex = 0.7,
         bty = "n"
         )
  dev.off()
  
  # save layout table
  if (keep.layout == T) {
    return(net.layout)
    write.table(net.layout, "layout.txt",sep = "\t",row.names = F, col.names = F)
  }
}