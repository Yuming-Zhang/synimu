#' @title Simulate White Noise
#' @description This function simulates multiple white noise processes.
#' @export
#' @param N An \code{integer} indicating the length of the time series.
#' @param sigma_wn A \code{vector} of variances of the white noise processes.
#' @return A \code{matrix} of the simulated white noise processes.
#' @author Davide Antonio Cucci

simulate_wn = function(N, sigma_wn) {
  n_ts = length(sigma_wn)

  WN = matrix(0, nrow = n_ts, ncol = N)
  for (i in 1:n_ts) {
    WN[i,] = rnorm(N, mean = 0, sd = sqrt(sigma_wn[i]))
  }

  t(WN)
}

#' @title Simulate Correlated Random Walk
#' @description This function simulates multiple correlated random walk processes.
#' @export
#' @param N An \code{integer} indicating the length of the time series.
#' @param sigma_rw A \code{matrix} denoting the covariance of the random walk innovations.
#' @return A \code{matrix} of the simulated correlated random walk processes.
#' @author Davide Antonio Cucci

simulate_corr_rw = function(N, sigma_rw) {
  if (is.null(dim(sigma_rw))) {
    n_ts = 1
  } else {
    n_ts = dim(sigma_rw)[1]
  }

  INN = mvrnorm(N, rep(0,n_ts), sigma_rw)

  RW = apply(INN, 2, cumsum)

  RW
}


#' @title Simulate Correlated AR(1)
#' @description This function simulates multiple correlated AR(1) processes.
#' @export
#' @param N An \code{integer} indicating the length of the time series.
#' @param phi_ar1 A \code{vector} indicating the phi terms of the AR(1) processes.
#' @param sigma_ar1 A \code{matrix} denoting the covariance of the AR(1) innovations.
#' @return A \code{matrix} of the simulated correlated AR(1) processes.
#' @author Davide Antonio Cucci

simulate_corr_ar1 = function(N, phi_ar1, sigma_ar1) {
  n_ts = length(phi_ar1)

  INN_AR1 = t(mvrnorm(N, rep(0,n_ts), sigma_ar1))

  AR1 = matrix(0, nrow = n_ts, ncol = N)
  AR1[,1] = INN_AR1[,1]
  for (i in 2:N) {
    AR1[,i] = AR1[,i-1] * phi_ar1 + INN_AR1[,i]
  }

  t(AR1)
}

