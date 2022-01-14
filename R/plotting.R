gg_color_hue <- function(n, alpha) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

#' @title Wavelet Variance Plot for Virtual Gyroscopes
#' @description This function plots the wavelet variance of the obtained virtual gyroscopes with the provided coefficients
#' @export
#' @param Xt A \code{matrix} of dimension T by p, where T is the length of the time series and p is the number of processes.
#' @param ... Multiple \code{vector} of coefficients.
#' @param names A \code{vector} of strings to be used for the legend.
#' @param ration_ylim Custom y limits for the ratio plot
#' @importFrom wv wvar
#' @author Davide Antonio Cucci
plot_virtual_gyro <- function(Xt, ..., names, ratio_ylim = NA) {
  coeffs = list(...)

  names = c(names, "Equal")

  coeffs = list(...)
  coeffs[[length(coeffs)+1]] = rep(1, ncol(Xt))/ncol(Xt)

  # compute the virtual gyros

  vgs = list()
  for (j in seq_along(coeffs)) {
    vgs[[j]] = get_virtual_gyro(Xt, coeffs[[j]])
  }

  # compute wvs of each single gyroscope

  wvs = list()
  for (i in 1:ncol(Xt)) {
    wvs[[i]] = wvar(Xt[,i])
  }

  i_single = 1:i

  # compute wv for each virtual gyroscope
  for (j in seq_along(vgs)) {
    wvs[[i+j]] = wvar(vgs[[j]])
  }

  i_virtual = (i+1):(i+j)

  # compute ratios wrt equal coeffs

  i_eq = length(wvs)

  ratios = list()
  h = 1
  for (j in i_virtual) {
    ratios[[h]] = wvs[[j]]$variance/wvs[[i_eq]]$variance
    h = h + 1
  }

  # prepare plot

  # compute limits

  n_scales = length(wvs[[1]]$variance)
  yl = sapply(wvs, function (x) {c(max(x$ci_high), min(x$ci_low))})
  yl = c(min(yl[2,]),max(yl[1,]))

  xl = c(min(wvs[[1]]$scales), max(wvs[[1]]$scales))

  if (is.na(ratio_ylim[1])) {
    ylr = sapply(ratios, function (x) {c(max(x), min(x))})
    ylr = c(min(ylr[2,]),max(ylr[1,]))*c(1,2)
    ylr = c(10^floor(log10(ylr[1])), 10^ceiling(log10(ylr[2])))
  } else {
    ylr = ratio_ylim
  }

  cols = gg_color_hue(length(coeffs))
  colsa = gg_color_hue(length(coeffs), 0x20/256)

  layout(matrix(c(1,2), nrow=2, ncol=1), heights=c(2.5,4))
  par(mar = c(0.5, 4.5, 1, 0.5))

  # plot ratios

  plot(NA, xlim = xl, ylim = ylr, log = "xy", xlab = NA, ylab = "Ratio wrt Equal", xaxt = "n", yaxt="n")

  # fix axis

  xlogs = 0:30
  # axis(1, at=10^xlogs, label= sapply(xlogs, function(i) as.expression(bquote(10^ .(i)))))

  ylogs = -20:20
  ylogs2 = seq(-20,20,.1)
  axis(2, at=10^ylogs, label = sapply(ylogs, function(i) as.expression(bquote(10^ .(i)))))
  axis(2, at=10^ylogs2, label = FALSE, lwd=0.5)

  abline(h=10^ylogs, v=10^xlogs, col = "grey", lt=2)

  for (j in seq_along(ratios)) {
    lines(wvs[[1]]$scales, ratios[[j]], col = cols[j], lw = 2)
  }

  legend("top", horiz = T, legend = names, lw=2, col = cols, lt = 1)


  # plot wvs

  par(mar = c(4.5, 4.5, 1, 0.5))
  plot(NA, xlim = xl, ylim = yl, log = "xy", xlab = "Scales [samples]", ylab = "Wavelet Variance", xaxt = "n", yaxt = "n")

  # fix axis

  xlogs = 0:30
  axis(1, at=10^xlogs, label= sapply(xlogs, function(i) as.expression(bquote(10^ .(i)))))

  ylogs = -20:20
  axis(2, at=10^ylogs, label = sapply(ylogs, function(i) as.expression(bquote(10^ .(i)))))

  abline(h=10^ylogs, v=10^xlogs, col = "grey", lt=2)

  # plot single wvs

  for (i in i_single) {
    lines(wvs[[i]]$scales, wvs[[i]]$variance, col = "grey")
    polygon(c(wvs[[i]]$scales, rev(wvs[[i]]$scales)), c(wvs[[i]]$ci_high, rev(wvs[[i]]$ci_low)), border = NA, col = "#A0A0A020")
  }

  # plot virtual

  for (j in i_virtual) {
    lines(wvs[[j]]$scales, wvs[[j]]$variance, col = cols[j-i_virtual[1]+1], lw = 2)
    polygon(c(wvs[[j]]$scales, rev(wvs[[j]]$scales)), c(wvs[[j]]$ci_high, rev(wvs[[j]]$ci_low)), border = NA, col = colsa[j-i_virtual[1]+1])
  }

  legend("top", horiz = T, legend = names, lw=2, col = cols, lt = 1)
}
