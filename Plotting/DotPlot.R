# dotplot for OTU table (or similar)

dotPlot <- function(Data, taxonomy, condition = colnames(Data), color = rep("black", nrow(Data)), max.size = 4, grid = T) {
  if (!"scales" %in% installed.packages()) {
    install.packages("scales")
  }
  require(scales)
  
  scaleData <- rescale(as.matrix(Data), to = c(0.1, max.size))
  scaleData2 <- scaleData[order(taxonomy, decreasing = T),]
  taxonomy2 <- taxonomy[order(taxonomy, decreasing = T)]
  
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
           bg = color[i],
           cex = scaleData2[i, ],
           pch = 21,
           col = "white")
  }
}
