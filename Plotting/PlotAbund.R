#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenr�ck
# acknowledgements: 
# primary authority: Christiane Hassenr�ck
# other authorities: 
#=============================================================================
# File contents
# Wrapper function to create a stacked barplot of the n most abundant taxa per sample in a sample by taxa table. 
# Abundant taxa will be selected per sample, so that the final number of taxa in the plot may be larger than n.
#
# input: 
# relData - sample (colums) by taxa (rows) table in relative abundances
# abund - number of abundant taxa to plot
# 
# output:
# stacked barplot and corresponding legend in a separate window
#
# dependencies:
# The function was designed for usage in windows. The calling of a new graphics device on unix is different
#=============================================================================
#documentation end

PlotAbund <- function(relData, abund,
                      colorPalette = c("violet", "purple", "blue", "white", "darkgreen", "yellow", "red", "darkred")) {
  abund_names <- c()
  for (i in 1:ncol(relData)) {
    abund_names <- unique(c(abund_names, rownames(relData)[order(relData[, i], decreasing = T)][1:abund]))
  }
  abund_rel0 <- relData[abund_names, ]
  abund_rel <- rbind(abund_rel0, 100 - colSums(abund_rel0))
  rownames(abund_rel)[nrow(abund_rel)] <- "other"
  windows(width = 20, height = 10)
  par(mfrow = c(1, 2))
  barplot(abund_rel, 
          col = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
          ylim = c(0, 100), 
          las = 2
          )
  plot.new()
  legend("center", 
         pch = 22, 
         col = "black", 
         pt.bg = rev(c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey")),
         legend = rev(rownames(abund_rel)),
         pt.cex = 1.5, 
         cex = 0.8
         ) 
}