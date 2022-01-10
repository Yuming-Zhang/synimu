gg_color_hue <- function(n, alpha) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

#' @export

plot_results = function(singles, virtuals, yl, ci=T, legend, unit = "deg^s/s^2") {

  scales = singles[[1]]$scales

  plot(NA,
       ylim = yl,
       xlim = c(min(scales), max(scales)),
       ylab = paste("WV [", unit, "]", sep = ""),
       xlab = "scales",
       log = "xy",
       xaxt = "n",
       yaxt = "n")

  axis(1, at=scales, labels=log2(scales))
  ylab = floor(log10(yl[1])):ceiling(log10(yl[2]))
  axis(2, at=10^ylab, labels=parse(text=sprintf("10^%.0f",ylab)))
  # add grid
  abline(v=scales, col="grey")
  abline(h=10^ylab, col="grey")

  # # add polygon ci
  # if(ci ==T){
  #   polygon(x = c(scales, rev(scales)),
  #           y = c(wv_emph$ci_low, rev(wv_emph$ci_high)),
  #           col = "#E5EBF2",
  #           border = NA)
  # }

  n_cols = length(singles) + length(virtuals)

  cols = gg_color_hue(n_cols, 1)

  for (i in seq_along(singles)) {
    lines(singles[[i]]$scales, singles[[i]]$variance, type="l", col = cols[i], lt=2)
  }

  for (i in seq_along(virtuals)) {
    lines(virtuals[[i]]$scales, virtuals[[i]]$variance, type="l", col = cols[i+length(singles)], lw=3)
  }

  legend(
    "topright",
    legend=legend,
    col = cols,
    lt = c(rep(2, length(singles)), rep(1, length(virtuals))),
    lw = c(rep(1, length(singles)), rep(3, length(virtuals))),
    bty="n", bg="transparent",
    ncol = 2
  )
}
