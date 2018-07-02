#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Wrapper function to create a stacked barplot of the n most abundant taxa per sample in a sample by taxa table. 
# Abundant taxa will be selected per sample, so that the final number of taxa in the plot may be larger than n.
#
# input: 
# relData - sample (colums) by taxa (rows) table in relative abundances
# abund - number of abundant taxa to plot
# margin - adjust plot margin for long labels
# 
# output:
# stacked barplot and corresponding legend in a separate window
#
# dependencies:
# The function was designed for usage in windows. The calling of a new graphics device on unix is different
#=============================================================================
#documentation end

PlotAbund <- function(relData, abund, 
                      margin = par()$mar,
                      colorPalette = c("violet", "purple", "blue", "white", "darkgreen", "yellow", "red", "darkred"),
                      plot.ratio = c(3, 1),
                      method = c("nmost", "precentage"),
                      open.window = T,
                      save.OTU = F,
                      sample.names = colnames(relData),
                      sort.taxa = F,
                      ...) {
  if(method == "nmost") {
    abund_names <- c()
    for (i in 1:ncol(relData)) {
      abund_names <- unique(c(abund_names, rownames(relData)[order(relData[, i], decreasing = T)][1:abund]))
    }
    abund_rel0 <- relData[abund_names, ]
  }
  if(method == "percentage") {
    abund_rel0 <- relData[apply(relData, 1, function(x) { max(x) >= abund }), ]
  }
  if (sort.taxa == T) {
    abund_rel0 <- abund_rel0[order(rownames(abund_rel0)), ]
  }
  abund_rel <- rbind(abund_rel0, 100 - colSums(abund_rel0))
  rownames(abund_rel)[nrow(abund_rel)] <- "other"
  abund_rel <- as.matrix(abund_rel)
  if (open.window == T) {
    windows(width = 20, height = 10)
  }
  # par(mfrow = c(1, 2), mar = margin, xpd = NA)
  layout(mat = matrix(c(1, 2), 1, 2, byrow = T), widths = plot.ratio)
  par(mar = margin)
  barplot(abund_rel,
          col = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
          ylim = c(0, 100), 
          las = 2,
          ylab = "Relative sequence abundance [%]",
          names.arg = sample.names,
          cex.names = 0.8,
          ...)
  par(mar = c(1,1,1,1))
  plot.new()
  legend("center", 
         pch = 22, 
         col = "black", 
         pt.bg = rev(c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey")),
         legend = rev(rownames(abund_rel)),
         pt.cex = 1.5, 
         cex = 0.8
         ) 
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  par(mfrow = c(1,1))
  
  if(save.OTU == T) {
    return(
      data.frame(
        abund_rel,
        color = c(colorRampPalette(colorPalette)(nrow(abund_rel) - 1), "darkgrey"),
        stringsAsFactors = F
      )
    )
  }
}
