#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Function to create a dotplot of e.g. OTU abundances
#
# input: 
# Data - sample (colums) by taxa (rows) table (data.frame or matrix), or similar. Values will be rescaled to cex between 0.1 and max.size
# taxonomy - row labels for plot (character vector), will be used to order rows, default rownames(Data)
# condition - column labels for plot (character vector), default colnames(Data)
# color - color matrix for dotplot (character matrix), same dimensions as Data
# max.size - maximum cex for dots (default: 4)
# grid - should gridlines be plotted in the background (default: TRUE)
# 
# output:
# dotplot showing abundance (values) of each OTU per sample
#
# dependencies:
# scales (if not available, will be installed by running the function)
#=============================================================================
#documentation end


dotPlot <- function(Data, 
                    taxonomy = rownames(Data), 
                    condition = colnames(Data), 
                    color = matrix("black", nrow(Data), ncol(Data)),
                    max.size = 4, 
                    grid = T) {
  
  if (!"scales" %in% installed.packages()) {
    install.packages("scales")
  }
  require(scales)
  
  scaleData <- rescale(as.matrix(Data), to = c(0.1, max.size)) # rescale values to cex
  scaleData2 <- scaleData[order(taxonomy, decreasing = T),] # order rows by taxon name
  taxonomy2 <- taxonomy[order(taxonomy, decreasing = T)] # order row labels by taxon name
  color2 <- color[order(taxonomy, decreasing = T),]
  
  plot(0,
       0,
       type = "n",
       axes = F,
       ylab = "",
       xlab = "",
       ylim = c(1, nrow(scaleData2)),
       xlim = c(1, ncol(scaleData2)))
  axis(1, 
       las = 2,
       at = 1:ncol(scaleData2),
       labels = condition,
       cex.axis = 0.8)
  axis(2,
       las = 2,
       at = 1:nrow(scaleData2),
       labels = taxonomy2,
       cex.axis = 0.8)
  if(grid == T) {
    abline(h = 1:nrow(scaleData2), 
           v = 1:ncol(scaleData2), 
           col = "lightgrey", 
           lty = 3)
  }
  for(i in 1:nrow(scaleData2)) {
    points(1:ncol(scaleData2),
           rep(i, ncol(scaleData2)),
           bg = color2[i, ],
           cex = scaleData2[i, ],
           pch = 21,
           col = "white")
  }
}
