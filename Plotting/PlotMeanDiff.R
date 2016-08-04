#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Wrapper function for a 2-panel plot (1) to show mean sequence proportion between 2 treatments as barplot, 
# and (2) to show the mean difference with 95% confidence interval
#
# input: 
# mean.diff - matrix with 5 columns: 
#   mean sequence proportion for treatment A
#   mean sequence proportion for treatment B
#   mean difference between treatment A and treatment B
#   lower confidence interval of mean difference
#   upper confidence interval of mean difference
# labels - labels of taxa, default: rownames(mean.diff)
# col1 - color of bars for condition A in barplot (character)
# col2 - color of bars for condition A in barplot (character)
# mar1 - margins of barplot, adjust for long labels
# xlim2 - range of axis for plotting mean difference (second plot)
# xlab1 - x axis label of barplot
# xlab2 - x axis label of second plot
# 
# output:
# 2-panel plot in separate window
#
# dependencies:
# The function was designed for usage in windows. The calling of a new graphics device on unix is different
#=============================================================================
#documentation end

PlotMeanDiff <- function(mean.diff, 
                         labels = rownames(mean.diff), 
                         col1 = "darkgrey", 
                         col2 = "black", 
                         mar1 = c(4,15,1,1),
                         xlim2 = if(min(mean.diff[, 3:5]) < 0 & max(mean.diff[, 3:5]) > 0) {
                             range(mean.diff[, 3:5])
                           } else {
                             if (min(mean.diff[, 3:5]) > 0) {
                               c(-0.1, max(mean.diff[, 3:5]))
                             }
                             if (max(mean.diff[, 3:5]) < 0) {
                               c(min(mean.diff[, 3:5]), 0.1)
                             }
                           }, 
                         xlab1 = "Sequence proportion [%]", 
                         xlab2 = "Log2-fold change") {
  windows(width = 20, height = (nrow(mean.diff) + 1))
  par(mfrow = c(1, 2), mar = mar1)
    
  mean.plot <- c()
  for(i in 1:nrow(mean.diff)) {
    mean.plot <- c(mean.plot, mean.diff[i, 1:2], 0)
  }
  bp <- barplot(rev(mean.plot), 
                horiz = T, 
                col = rev(rep(c(col1, col2, NA), nrow(mean.diff))), 
                names.arg = "", 
                border = NA,
                xlab = xlab1,
                yaxs = "i")
  y.coord <- apply(matrix(bp[-c(3 * (1:nrow(mean.diff)) - 2), 1], 
                          nrow(mean.diff), 
                          2, 
                          byrow = T),
                   1,
                   mean)
  mtext(rev(as.character(labels)), 
        side = 2, 
        at = y.coord, 
        las = 1,
        cex = 0.8)
  par(mar = c(4,2,1,1))
  plot(0,0,
       type = "n",
       ylim = c(min(bp) - 0.5, max(bp) + 0.5),
       xlim = xlim2,
       axes = F,
       ylab = "",
       xlab = xlab2,
       yaxs = "i")
  axis(1)
  abline(v = 0)
  segments(rev(mean.diff[, 4]),
           y.coord,
           rev(mean.diff[, 5]),
           y.coord)
  segments(rev(mean.diff[, 4]),
           y.coord - 0.3,
           rev(mean.diff[, 4]),
           y.coord + 0.3)
  segments(rev(mean.diff[, 5]),
           y.coord - 0.3,
           rev(mean.diff[, 5]),
           y.coord + 0.3)
  points(rev(mean.diff[, 3]),
         y.coord,
         cex = 2,
         pch = 1)
}
