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
# The following code can help you create the input for the plotting function
# Data is a data.frame with samples in columns and the taxa of interest in rows
# Treatment is a vector (factor with 2 levels) with the treatment for each sample
#   mean.diff <- matrix(NA, nrow(Data), 5)
#   rownames(mean.diff) <- rownames(Data)
#   colnames(mean.diff) <- c("mean.uninhabited", "mean.inhabited", "mean.diff", "CI.lower", "CI.upper")
#   mean.diff[, 1:2] <- t(apply(Data, 1, function(x) {by(x, Treatment, mean)}))
#   mean.diff[, 3] <- mean.diff[, 2] - mean.diff[, 1]
#   mean.diff[, 4:5] <- t(apply(Data, 1, function(x) {t.test(x ~ Treatment)$conf.int}))
# further arguments
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
                         open.window = T,
                         higher.tax = NULL,
                         set.layout = T,
                         legend = F,
                         cond = NULL,
                         panel.no = NULL,
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
  
  if(open.window == T) {
    windows(width = 20, height = (0.8 * nrow(mean.diff) + 1))
  }
  
  if(set.layout == T) {
    par(mfrow = c(1, 2))
  }
  
  par(mar = mar1)
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
                yaxs = "i",
                cex.axis = 0.8,
                lwd = 0.5)
  
  y.coord <- apply(matrix(bp[-c(3 * (1:nrow(mean.diff)) - 2), 1], 
                          nrow(mean.diff), 
                          2, 
                          byrow = T),
                   1,
                   mean)
  
  if(is.null(higher.tax)) {
    mtext(rev(as.character(labels)), 
          side = 2, 
          at = y.coord, 
          las = 1,
          cex = 0.5)
  } else {
    mtext(rev(as.character(labels)), 
          side = 2, 
          at = bp[seq(2, nrow(bp), 3), 1], 
          las = 1,
          cex = 0.5)
    mtext(rev(higher.tax), 
          side = 2, 
          at = bp[seq(3, nrow(bp), 3), 1], 
          las = 1,
          cex = 0.5)
  }
  
  if(!is.null(panel.no)) {
    par(xpd = NA)
    text(
      par("usr")[1] - 0.5 * par("usr")[2],
      par("usr")[4] + 2,
      labels = panel.no,
      font = 2,
      cex = 1.5
    )
    par(xpd = F)
  }
  
  if(legend == T) {
    par(xpd = NA)
    points(
      rep(par("usr")[1] - 0.5 * par("usr")[2], 2),
      c(par("usr")[3] - par("usr")[4] * 0.1, par("usr")[3] - par("usr")[4] * 0.25),
      pch = 15,
      cex = 2,
      col = c(col1, col2)
    )
    text(
      rep(par("usr")[1] - 0.5 * par("usr")[2], 2),
      c(par("usr")[3] - par("usr")[4] * 0.1, par("usr")[3] - par("usr")[4] * 0.25),
      labels = cond,
      cex = 1,
      pos = 4
    )
    par(xpd = F) 
  }
  
  par(mar = c(mar1[1],2,mar1[3],1))
  plot(0,0,
       type = "n",
       ylim = c(min(bp) - 0.5, max(bp) + 0.5),
       xlim = xlim2,
       axes = F,
       ylab = "",
       xlab = xlab2,
       yaxs = "i",
       cex.axis = 0.8)
  axis(1, cex.axis = 0.8, lwd = 0.5)
  abline(v = 0, lwd = 0.5)
  segments(rev(mean.diff[, 4]),
           y.coord,
           rev(mean.diff[, 5]),
           y.coord,
           lwd = 0.5)
  segments(rev(mean.diff[, 4]),
           y.coord - 0.3,
           rev(mean.diff[, 4]),
           y.coord + 0.3,
           lwd = 0.5)
  segments(rev(mean.diff[, 5]),
           y.coord - 0.3,
           rev(mean.diff[, 5]),
           y.coord + 0.3,
           lwd = 0.5)
  points(rev(mean.diff[, 3]),
         y.coord,
         cex = 1.5,
         pch = 16)
}
