#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# alternative plotting function for iNEXT output that does not rely on ggplot2
# only 'rarefaction' plots supported
#
# input: 
# iNEXT() output
# hill: which Hill numbers to plot (default: 0, 1, 2)
# color: character vector giving color for each sample (default based on rainbow color palette)
# labels: customized sample labels (default: name names of input)
# selection: which samples to plot (default: plot all samples in the order they occur in the input)
#  selection can either be:
#  an integer vector with the column number in the iNEXT input (subsetting and ordering at the same time)
#  a character vector giving the sample names (subsetting and ordering at the same time)
#  a logical vector of the length == sample number
#
# output:
# rarefaction plots for specified Hill numbers with confidence interval
# plot created in a separate window
#
# dependencies:
# The function was designed for usage in windows. The calling of a new graphics device on unix is different
#=============================================================================
#documentation end


iNEXTplot <- function(
  input,
  hill = c(0, 1, 2),
  color = rainbow(length(input$iNextEst)),
  labels = names(input$iNextEst),
  selection = 1:length(input$iNextEst),
  window.new = T,
  ylim = c(0, max(as.matrix(input2$qD.UCL))),
  dynamic.ylim = F,
  xlim = c(0, max(as.matrix(input2$m)))
) {
  
  # format input
  input1 <- vector("list", length = 6)
  names(input1) <- c("m", "method", "order", "qD", "qD.LCL", "qD.UCL")
  # create subset for hill number estimates
  input1$qD <- sapply(input$iNextEst, function(x) { x$qD })
  # create subset for upper CI
  input1$qD.UCL <- sapply(input$iNextEst, function(x) { x$qD.UCL })
  # create subset for lower CI
  input1$qD.LCL <- sapply(input$iNextEst, function(x) { x$qD.LCL })
  # create subset for sequencing depth
  input1$m <- sapply(input$iNextEst, function(x) { x$m })
  # create subset for sequencing depth
  input1$method <- sapply(input$iNextEst, function(x) { x$method })
  # create subset for order of hill numbers
  input1$order <- sapply(input$iNextEst, function(x) { x$order })

  # apply selection
  input2 <- lapply(input1, function(x) { x[, selection] })
  color2 <- color[selection]
  labels2 <- labels[selection]
  
  # open new graphics device
  if (window.new == T) {
    windows(width = (length(hill) + 1) * 5, height = 10)
  }
  par(mfrow = c(1, length(hill) + 1))

  for (i in 1:length(hill)) {
    # adjust ylim
    if(dynamic.ylim == T) {
      ylim <- c(0, max(as.matrix(input2$qD.UCL[input2$order[, 1] == hill[i], ]))) 
    } 
    
    # create empty plot
    plot(
      0, 0,
      type = "n",
      xlim = xlim,
      ylim = ylim,
      xlab = "Number of sequences",
      ylab = paste("q = ", hill[i], sep = ""),
      axes = F
    )

    # add axes
    axis(1, cex.axis = 0.8)
    axis(2, cex.axis = 0.8, las = 2)

    # subset data to hill number of interest
    input.sub <- lapply(input2, function(x) { x[input2$order[, 1] == hill[i], ] })

    # add data
    for (j in 1:ncol(input.sub$qD)) {
      polygon(
        x = c(input.sub$m[, j], rev(input.sub$m[, j])),
        y = c(input.sub$qD.UCL[, j], rev(input.sub$qD.LCL[, j])),
        col = adjustcolor(color2[j], alpha.f = 0.3),
        border = NA
      )
      lines(
        input.sub$m[input.sub$method[, j] != "extrapolated", j],
        input.sub$qD[input.sub$method[, j] != "extrapolated", j],
        col = color2[j],
        lty = 1
      )
      lines(
        input.sub$m[input.sub$method[, j] == "extrapolated", j],
        input.sub$qD[input.sub$method[, j] == "extrapolated", j],
        col = color2[j],
        lty = 3
      )
    }
  }
  
  # add legend
  plot(
    0, 0,
    type = "n",
    axes = F,
    xlim = c(0, 5),
    ylim = c(1, ncol(input.sub$qD)),
    xlab = "",
    ylab = ""
  )
  points(
    x = rep(0.5, ncol(input.sub$qD)),
    y = 1:ncol(input.sub$qD),
    cex = 2,
    pch = 15,
    col = adjustcolor(color2, alpha.f = 0.3)
  )
  segments(
    x0 = rep(0, ncol(input.sub$qD)),
    y0 = 1:ncol(input.sub$qD),
    x1 = rep(1, ncol(input.sub$qD)),
    y1 = 1:ncol(input.sub$qD),
    lty = 1,
    col = color2
  )
  text(
    x = rep(1, ncol(input.sub$qD)),
    y = 1:ncol(input.sub$qD),
    pos = 4,
    labels = labels2,
    cex = 0.8
  )
}
