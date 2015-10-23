#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements:
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# plotting function for Alpha diversity indices (Hill numbers)
#
# input: 
# x - output from SumsampleNGS
# per.sample - should the lines be plotted per sample (T) or per condition (F)
# cond - grouping variable of class factor (if (per.sample == T) not required)
# color - vector of colors of class character for each sample (if per.sample == T) or for each of levels(cond) if per.sample == F
# method - should values and CI be based on mean and 95% confidence interval (mean) or median and interpquartile range (default: median)
# CI - plot confidence interval (default: T)
# raw - should lines for non-subsampled data also be plotted (default: F)
# select - which subset of samples to use, i.e. the column number of the sample by OTU table (default: all samples)
# label - sample names if (per.sample == T), defaults to the column names of the sample by OTU table 
# ylim - see plotting parameters
#
# output:
# plot with diversity index on y-axis and Hill number on x-axis
# 
# dependencies:
# SubsampleNGS.R
#=============================================================================
#documentation end

PlotHillNGS <- function(input, per.sample, cond, color,
                        method = "median", CI = T, raw = F,
                        select = NULL,
                        label = NULL,
                        ylim = NULL) {
  
  if (is.null(select)) {
    x <- input
  } else {
    x <- c(iterations = list(lapply(input$iterations, function(a) { a[,select] })),
           lapply(input[3:4], function(a) { a[,select] }))
  }
  
  if (is.null(label)) {
    label = colnames(x$summaryHill)
  }
  
  if (per.sample == T) {
    
    if (CI == T) {
      if (method == "mean") {
        N <- dim(x$iterations$nOTU)[1]
        S.0 <- apply(x$iterations$nOTU, 2, sd)
        S.1 <- apply(exp(x$iterations$shannon), 2, sd)
        S.2 <- apply(x$iterations$invS, 2, sd)
        CI.0 <- qt(0.975, df = (N - 1)) * (S.0 / sqrt(N))
        CI.1 <- qt(0.975, df = (N - 1)) * (S.1 / sqrt(N))
        CI.2 <- qt(0.975, df = (N - 1)) * (S.2 / sqrt(N))
        CI.Hill <- rbind(CI.0,CI.1,CI.2)
        
      }
      if (method == "median") {
        IQL.0 <- apply(x$iterations$nOTU, 2, function (a) { quantile(a)[2] })
        IQU.0 <- apply(x$iterations$nOTU, 2, function (a) { quantile(a)[4] })
        IQL.1 <- apply(exp(x$iterations$shannon), 2, function (a) { quantile(a)[2] })
        IQU.1 <- apply(exp(x$iterations$shannon), 2, function (a) { quantile(a)[4] })
        IQL.2 <- apply(x$iterations$invS, 2, function (a) { quantile(a)[2] })
        IQU.2 <- apply(x$iterations$invS, 2, function (a) { quantile(a)[4] })
        IQL.Hill <- rbind(IQL.0, IQL.1, IQL.2)
        IQU.Hill <- rbind(IQU.0, IQU.1, IQU.2)
      }
    }
    
    if (method == "mean") {
      input <- x$summaryHill
    }
    if (method == "median") {
      input <- rbind(apply(x$iterations$nOTU, 2, median),
                     exp(apply(x$iterations$shannon, 2, median)),
                     apply(x$iterations$invS, 2, median)
                     )
    }
    
    if (raw == T) {
      input2 <- x$summaryHillRaw
    } else {
      input2 <- input
    }
    
    plot(2, max(input2) + 0.1 * max(input2), 
         type = "n", 
         xlim = c(0, 2.5), 
         ylim = if (is.null(ylim)) { 
           c(0,max(input2) + 0.1 * max(input2)) 
           } else {
             ylim
           }, 
         axes = F,
         xlab = "Hill",
         ylab = "Diversity"
         )
    axis(1, at = c(0, 1, 2), las = 1, cex.axis = 0.7)
    axis(2, las = 1, cex.axis = 0.7)
    for (i in 1:ncol(input)) {
      if (raw == T) {
        lines(c(0, 1, 2), input2[, i],
              col = color[i],
              lty = 2
              )
        text(2.25, input2[3, i],
             labels = paste("all", label[i], collapse = "_"),
             col = "grey",
             pos = 4,
             cex = 0.7
             )
      }
      if (CI == T) {
        if (method == "mean") {
          polygon(x = c(0, 1, 2, 2, 1, 0),
                  y = c(input[, i] + CI.Hill[, i],rev(input[, i] - CI.Hill[, i])),
                  col = adjustcolor(color[i], alpha.f = 0.3),
                  border = NA
                  )
        }
        if (method == "median") {
          polygon(x = c(0, 1, 2, 2, 1, 0),
                  y = c(IQU.Hill[, i],rev(IQL.Hill[, i])),
                  col = adjustcolor(color[i], alpha.f = 0.3),
                  border = NA
                  )
        }
      }
      lines(c(0, 1, 2), input[, i],
            col = color[i],
            lty = 1,
            lwd = 1.5
            )
      text(2, input[3, i],
           labels = label[i],
           pos = 4,
           cex = 0.7
           )
    }
    
  } else {
    
    N <- c(table(cond))
    if (CI == T) {
      if (method == "mean") {
        S <- apply(x$summaryHill, 1, function(a) { c(by(a, cond, sd)) })
        CI.Hill <- apply(S, 2, function(a) { qt(0.975, df = (N - 1)) * (a / sqrt(N)) })
      }
      if (method == "median") {
        IQL.Hill <- apply(x$summaryHill, 1, function(a) { c(by(a, cond, function(b) { quantile(b)[2] })) })
        IQU.Hill <- apply(x$summaryHill, 1, function(a) { c(by(a, cond, function(b) { quantile(b)[4] })) })
      }
    }
    
    if (method == "mean") {
      input <- apply(x$summaryHill, 1, function(a) { c(by(a, cond, mean)) })
      if (raw == T) {
        S2 <- apply(x$summaryHillRaw, 1, function(a) { c(by(a, cond, sd)) })
        CI.Hill2 <- apply(S2, 2, function(a) { qt(0.975, df = (N - 1)) * (a / sqrt(N)) })
        input2 <- apply(x$summaryHillRaw, 1, function(a) { c(by(a, cond, mean)) })
      } else {
        input2 <- input
      }
    }
    if (method == "median") {
      input <- apply(x$summaryHill, 1, function(a) { c(by(a, cond, median)) })
      if (raw == T) {
        IQL.Hill2 <- apply(x$summaryHillRaw, 1, function(a) { c(by(a, cond, function(b) { quantile(b)[2] })) })
        IQU.Hill2 <- apply(x$summaryHillRaw, 1, function(a) { c(by(a, cond, function(b) { quantile(b)[4] })) })
        input2 <- apply(x$summaryHillRaw, 1, function(a) { c(by(a, cond, median)) })
      } else {
        input2 <- input
      }
    }
    
    plot(2, max(input2) + 0.1 * max(input2), 
         type = "n", 
         xlim = c(0, 2.5), 
         ylim = if (is.null(ylim)) { 
           c(0,max(input2) + 0.1 * max(input2)) 
         } else {
           ylim
         }, 
         axes = F,
         xlab = "Hill",
         ylab = "Diversity"
         )
    axis(1, at = c(0, 1, 2), las = 1, cex.axis = 0.7)
    axis(2, las = 1, cex.axis = 0.7)
    for (i in 1:nrow(input)) {
      if (raw == T) {
        if(CI == T) {
          if (method == "mean") {
            lines(c(0, 1, 2), c(input2[i, ] + CI.Hill2[i, ]),
                  col = adjustcolor(color[i], alpha.f = 0.5),
                  lty = 3
                  )
            lines(c(0, 1, 2), c(input2[i, ] - CI.Hill2[i, ]),
                  col = adjustcolor(color[i], alpha.f = 0.5),
                  lty = 3
                  )
          }
          if (method == "median") {
            lines(c(0, 1, 2), IQU.Hill2[i, ],
                  col = adjustcolor(color[i], alpha.f = 0.5),
                  lty = 3
                  )
            lines(c(0, 1, 2), IQL.Hill2[i, ],
                  col = adjustcolor(color[i], alpha.f = 0.5),
                  lty = 3
                  )
          }
        }      
        lines(c(0, 1, 2), input2[i, ],
              col = color[i],
              lty = 2
              )
        text(2.25, input2[i, 3],
             labels = paste("all", levels(cond)[i], collapse = "_"),
             col = "grey",
             pos = 4,
             cex = 0.7
             )
      }
      if (CI == T) {
        if (method == "mean") {
          polygon(x = c(0, 1, 2, 2, 1, 0),
                  y = c(input[i, ] + CI.Hill[i, ],rev(input[i, ] - CI.Hill[i, ])),
                  col = adjustcolor(color[i], alpha.f = 0.3),
                  border = NA
                  )
        }
        if (method == "median") {
          polygon(x = c(0, 1, 2, 2, 1, 0),
                  y = c(IQU.Hill[i, ],rev(IQL.Hill[i, ])),
                  col = adjustcolor(color[i], alpha.f = 0.3),
                  border = NA
                  )
        }
      }
      lines(c(0, 1, 2), input[i, ],
            col = color[i],
            lty = 1,
            lwd = 1.5
            )
      text(2, input[i, 3],
           labels = levels(cond)[i],
           pos = 4,
           cex = 0.7
           )
    }
  }
  
}