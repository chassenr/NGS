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
                    filled = T,
                    show.legend = T,
                    scale.from = c(min(Data, na.rm = T), max(Data, na.rm = T))) {
  
  if (!"scales" %in% installed.packages()) {
    install.packages("scales")
  }
  require(scales)
  
  scaleData <- rescale(as.matrix(Data), from = scale.from, to = c(min.size, max.size)) # rescale values to cex
  if (sort == T) {
    scaleData2 <- scaleData[order(taxonomy, decreasing = T),] # order rows by taxon name
    taxonomy2 <- taxonomy[order(taxonomy, decreasing = T)] # order row labels by taxon name
    color2 <- color[order(taxonomy, decreasing = T),]
  } else {
    scaleData2 <- scaleData 
    taxonomy2 <- taxonomy
    color2 <- color
  }
  
  # LM <- lm(as.vector(scaleData) ~ as.vector(as.matrix(Data)), na.action = na.omit)
  legend.size <- rescale(legend.value, from = scale.from, to = c(min.size, max.size))

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
       cex.axis = 0.6,
       lwd = 0.5)
  axis(4,
       las = 2,
       at = 1:nrow(scaleData2),
       labels = taxonomy2,
       cex.axis = 0.6,
       lwd = 0.5)
  if(grid == T) {
    abline(h = 1:nrow(scaleData2), 
           v = 1:ncol(scaleData2), 
           col = "lightgrey", 
           lty = 1,
           lwd = 0.5)
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
  
  # https://stackoverflow.com/questions/29125019/get-margin-line-locations-mgp-in-user-coordinates
  line2user <- function(line, side) {
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(0:1, 'inches', 'user'))
    y_off <- diff(grconvertY(0:1, 'inches', 'user'))
    switch(side,
           `1` = par('usr')[3] - line * y_off * lh,
           `2` = par('usr')[1] - line * x_off * lh,
           `3` = par('usr')[4] + line * y_off * lh,
           `4` = par('usr')[2] + line * x_off * lh,
           stop("side must be 1, 2, 3, or 4", call.=FALSE))
  }
  
  legend.x <- seq(line2user(par("mar")[4] * 0.2, 4), 
                  line2user(par("mar")[4] * 0.8, 4),
                  length.out = length(legend.size)
                  )
  legend.y <- rep(line2user(par("mar")[1] * 0.2, 1),
                  length(legend.size)
                  )
  if(show.legend == T) {
    par(xpd = NA)
    points(legend.x, 
           legend.y, 
           pch = 16, 
           cex = legend.size
    )
    text(legend.x, line2user(par("mar")[1] * 0.4, 1), labels = legend.value, cex = 0.7)
    text(mean(legend.x), line2user(par("mar")[1] * 0.6, 1), labels = "clr-transformed", cex = 0.8)
    text(mean(legend.x), line2user(par("mar")[1] * 0.8, 1), labels = "sequence counts", cex = 0.8)
    par(xpd = F)
  }
}
