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
# legend.value - abundance values for legend
# taxonomy - row labels for plot (character vector), will be used to order rows, default rownames(Data)
# condition - column labels for plot (character vector), default colnames(Data)
# color - color matrix for dotplot (character matrix), same dimensions as Data
# min.size - minimum cex for dots (default: 0.1)
# max.size - maximum cex for dots (default: 4)
# grid - should gridlines be plotted in the background (default: TRUE)
# sort - should the rows of the input be ordered accoring to taxonomy (default TRUE)
# 
# output:
# dotplot showing abundance (values) of each OTU per sample
#
# dependencies:
# scales (if not available, will be installed by running the function)
#=============================================================================
#documentation end


dotPlot <- function(Data, legend.value,
                    taxonomy = rownames(Data), 
                    condition = colnames(Data), 
                    color = matrix("black", nrow(Data), ncol(Data)),
                    min.size = 0.1,
                    max.size = 4, 
                    grid = T,
                    sort = T,
                    filled = T) {
  
  if (!"scales" %in% installed.packages()) {
    install.packages("scales")
  }
  require(scales)
  
  scaleData <- rescale(as.matrix(Data), to = c(min.size, max.size)) # rescale values to cex
  if (sort == T) {
    scaleData2 <- scaleData[order(taxonomy, decreasing = T),] # order rows by taxon name
    taxonomy2 <- taxonomy[order(taxonomy, decreasing = T)] # order row labels by taxon name
    color2 <- color[order(taxonomy, decreasing = T),]
  } else {
    scaleData2 <- scaleData 
    taxonomy2 <- taxonomy
    color2 <- color
  }
  
  LM <- lm(as.vector(scaleData) ~ as.vector(as.matrix(Data)), na.action = na.omit)
  legend.size <- (LM$coefficients[2] * legend.value + LM$coefficients[1])

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
       cex.axis = 0.6)
  axis(4,
       las = 2,
       at = 1:nrow(scaleData2),
       labels = taxonomy2,
       cex.axis = 0.6)
  if(grid == T) {
    abline(h = 1:nrow(scaleData2), 
           v = 1:ncol(scaleData2), 
           col = "lightgrey", 
           lty = 1)
  }
  for(i in 1:nrow(scaleData2)) {
    if(filled == T) {
      points(1:ncol(scaleData2),
             rep(i, ncol(scaleData2)),
             bg = color2[i, ],
             col = "white",
             cex = scaleData2[i, ],
             pch = 21
      )
    } else {
      points(1:ncol(scaleData2),
             rep(i, ncol(scaleData2)),
             col = color2[i, ],
             cex = scaleData2[i, ],
             pch = 21
      )
    }
  }
  legend.x <- seq(par()$usr[2], 
                  par()$usr[2] + ncol(scaleData2)*2/10,
                  length.out = length(legend.size)
                  )
  legend.y <- rep(par()$usr[3],
                  length(legend.size)
                  )
  par(xpd = NA)
  points(legend.x, 
         legend.y, 
         pch = 16, 
         cex = legend.size
         )
  text(legend.x, legend.y, labels = legend.value, pos = 1, cex = 0.7)
  text(mean(legend.x), par()$usr[3] - nrow(scaleData2)/20, labels = "clr-transformed", cex = 0.8)
  text(mean(legend.x), par()$usr[3] - nrow(scaleData2)/20, labels = "sequence counts", cex = 0.8, pos = 1)
  par(xpd = F)
}
