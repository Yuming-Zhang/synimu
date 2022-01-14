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
#' @param freq The frequency of the data provided in Xt
#' @param ration_ylim Custom y limits for the ratio plot
#' @importFrom wv wvar
#' @author Davide Antonio Cucci
plot_virtual_gyro <- function(Xt, ..., names, freq = 1, ratio_ylim = NA) {
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

  xl = c(min(wvs[[1]]$scales), max(wvs[[1]]$scales))/freq

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

  xlogs = -10:30
  # axis(1, at=10^xlogs, label= sapply(xlogs, function(i) as.expression(bquote(10^ .(i)))))

  ylogs = -20:20
  # ylogs2 = seq(-20,20,.1)
  axis(2, at=10^ylogs, label = sapply(ylogs, function(i) as.expression(bquote(10^ .(i)))))
  # axis(2, at=10^ylogs2, label = FALSE, lwd=0.5)

  abline(h=10^ylogs, v=10^xlogs, col = "grey", lt=2)

  for (j in seq_along(ratios)) {
    lines(wvs[[1]]$scales/freq, ratios[[j]], col = cols[j], lw = 2)
  }

  legend("top", horiz = T, legend = names, lw=2, col = cols, lt = 1)


  # plot wvs

  par(mar = c(4.5, 4.5, 1, 0.5))
  plot(NA, xlim = xl, ylim = yl, log = "xy", xlab = "Averaging time [s]", ylab = expression(paste("Wavelet variance [", rad^2/s^2, "]")), xaxt = "n", yaxt = "n")

  # fix axis

  xlogs = -10:30
  axis(1, at=10^xlogs, label= sapply(xlogs, function(i) as.expression(bquote(10^ .(i)))))

  ylogs = -20:20
  axis(2, at=10^ylogs, label = sapply(ylogs, function(i) as.expression(bquote(10^ .(i)))))

  abline(h=10^ylogs, v=10^xlogs, col = "grey", lt=2)

  # plot single wvs

  for (i in i_single) {
    lines(wvs[[i]]$scales/freq, wvs[[i]]$variance, col = "grey")
    polygon(c(wvs[[i]]$scales, rev(wvs[[i]]$scales))/freq, c(wvs[[i]]$ci_high, rev(wvs[[i]]$ci_low)), border = NA, col = "#A0A0A020")
  }

  # plot virtual

  for (j in i_virtual) {
    lines(wvs[[j]]$scales/freq, wvs[[j]]$variance, col = cols[j-i_virtual[1]+1], lw = 2)
    polygon(c(wvs[[j]]$scales, rev(wvs[[j]]$scales))/freq, c(wvs[[j]]$ci_high, rev(wvs[[j]]$ci_low)), border = NA, col = colsa[j-i_virtual[1]+1])
  }

  legend("top", horiz = T, legend = names, lw=2, col = cols, lt = 1)
}

#' @title Plot Wavelet Variance
#' @description This function displays a plot of wavelet variance accounting for confidence interval values and supplied efficiency.
#' This function is a slight modification of the plot.wvar function provided in the wv package.
#' @export
#' @param x                A \code{wvar} object.
#' @param units            A \code{string} that specifies the units of time plotted on the x axis.
#' @param xlab             A \code{string} that gives a title for the x axis.
#' @param ylab             A \code{string} that gives a title for the y axis.
#' @param main             A \code{string} that gives an overall title for the plot.
#' @param col_wv           A \code{string} that specifies the color of the wavelet variance line.
#' @param col_ci           A \code{string} that specifies the color of the confidence interval polygon.
#' @param nb_ticks_x       An \code{integer} that specifies the maximum number of ticks for the x-axis.
#' @param nb_ticks_y       An \code{integer} that specifies the maximum number of ticks for the y-axis.
#' @param legend_position  A \code{string} that specifies the position of the legend (use \code{legend_position = NA} to remove legend).
#' @param ci_wv            A \code{boolean} that determines whether a confidence interval polygon will be drawn.
#' @param point_cex        A \code{double} that specifies the size of each symbol to be plotted.
#' @param point_pch        A \code{double} that specifies the symbol type to be plotted.
#' @param ...              Additional arguments affecting the plot.
#' @return Plot of wavelet variance and confidence interval for each scale.
#' @author Stephane Guerrier

plotwvar = function(x, units = NULL, xlab = NULL, ylab = NULL, main = NULL,
                    col_wv = NULL, col_ci = NULL, nb_ticks_x = NULL, nb_ticks_y = NULL,
                    legend_position = NULL, ci_wv = NULL, point_cex = NULL,
                    point_pch = NULL){

  # Labels
  if (is.null(xlab)){
    if (is.null(units)){
      xlab = expression(paste("Scale ", tau, sep =""))
    }else{
      xlab = bquote(paste("Scale ", tau, " [", .(units), "]", sep = " "))
    }
  }

  if (is.null(ylab)){
    ylab = expression(paste("Wavelet variance ", nu^2, sep = ""))
  }else{
    ylab = ylab
  }

  # Main Title
  if (is.null(main)){
    main = "Haar Wavelet Variance Representation"
  }

  # Line and CI colors
  if (is.null(col_wv)){
    col_wv = "darkblue"
  }

  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  # Range
  x_range = range(x$scales)
  x_low = floor(log10(x_range[1]))
  x_high = ceiling(log10(x_range[2]))

  y_range = range(c(x$ci_low, x$ci_high))
  y_low = floor(log10(y_range[1]))
  y_high = ceiling(log10(y_range[2]))

  # Axes
  if (is.null(nb_ticks_x)){
    nb_ticks_x = 6
  }

  if (is.null(nb_ticks_y)){
    nb_ticks_y = 5
  }

  x_ticks = seq(x_low, x_high, by = 1)
  if (length(x_ticks) > nb_ticks_x){
    x_ticks = x_low + ceiling((x_high - x_low)/(nb_ticks_x + 1))*(0:nb_ticks_x)
  }
  x_labels = sapply(x_ticks, function(i) as.expression(bquote(10^ .(i))))

  y_ticks <- seq(y_low, y_high, by = 1)
  if (length(y_ticks) > nb_ticks_y){
    y_ticks = y_low + ceiling((y_high - y_low)/(nb_ticks_y + 1))*(0:nb_ticks_y)
  }
  y_labels <- sapply(y_ticks, function(i) as.expression(bquote(10^ .(i))))

  # Legend Position
  if (is.null(legend_position)){
    if (which.min(abs(c(y_low, y_high) - log2(x$variance[1]))) == 1){
      legend_position = "topleft"
    }else{
      legend_position = "bottomleft"
    }
  }

  # Main Plot
  plot(NA, xlim = x_range, ylim = y_range, xlab = xlab, ylab = ylab,
       log = "xy", xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE)
  win_dim = par("usr")

  par(new = TRUE)
  plot(NA, xlim = x_range, ylim = 10^c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       xlab = xlab, ylab = ylab, log = "xy", xaxt = 'n', yaxt = 'n', bty = "n")
  win_dim = par("usr")

  # Add Grid
  abline(v = 10^x_ticks, lty = 1, col = "grey95")
  abline(h = 10^y_ticks, lty = 1, col = "grey95")

  # Add Title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = 10^c(win_dim[4], win_dim[4],
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
               win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = 10^(win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

  # Add Axes and Box
  lines(x_vec[1:2], rep(10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
  #y_ticks = y_ticks[(2^y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
  y_labels = y_labels[1:length(y_ticks)]
  box()
  axis(1, at = 10^x_ticks, labels = x_labels, padj = 0.3)
  axis(2, at = 10^y_ticks, labels = y_labels, padj = -0.2)

  # CI for WV
  if (ci_wv == TRUE || is.null(ci_wv)){
    polygon(c(x$scales, rev(x$scales)), c(x$ci_low, rev(x$ci_high)),
            border = NA, col = col_ci)
  }

  # Add legend
  CI_conf = 1 - x$alpha

  if (x$robust == TRUE){
    wv_title_part1 = "Empirical Robust WV "
  }else{
    wv_title_part1 = "Empirical WV "
  }

  if (!is.na(legend_position)){
    if (legend_position == "topleft"){
      legend_position = 10^c(1.1*win_dim[1], 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
      legend(x = legend_position[1], y = legend_position[2],
             legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))),
                        as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")")))),
             pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
    }else{
      if (legend_position == "topright"){
        legend_position = 10^c(0.7*win_dim[2]/freq, 0.98*(win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))
        legend(x = legend_position[1], y = legend_position[2],
               legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))),
                          as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }else{
        legend(legend_position,
               legend = c(as.expression(bquote(paste(.(wv_title_part1), hat(nu)^2))),
                          as.expression(bquote(paste("CI(",hat(nu)^2,", ",.(CI_conf),")")))),
               pch = c(16, 15), lty = c(1, NA), col = c(col_wv, col_ci), cex = 1, pt.cex = c(1.25, 3), bty = "n")
      }
    }
  }

  # Add WV
  lines(x$scales, x$variance, type = "l", col = col_wv, pch = 16)

  if (is.null(point_pch)){
    point_pch = 16
  }

  if (is.null(point_cex)){
    point_cex = 1.25
  }
  lines(x$scales, x$variance, type = "p", col = col_wv, pch = point_pch, cex = point_cex)
}
