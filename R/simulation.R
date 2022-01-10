#' @export

simulate_wn = function(N, sigma_wn) {
  n_ts = length(sigma_wn)

  WN = matrix(0, nrow = n_ts, ncol = N)
  for (i in 1:n_ts) {
    WN[i,] = rnorm(N, mean = 0, sd = sqrt(sigma_wn[i]))
  }

  t(WN)
}

#' @export

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

#' @export

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

#' Simulates a WN + coorelated RW
#'
#' @param N the number of samples
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#'
#' @export

simulate_wn_corr_rw = function(N, sigma_wn, sigma_rw) {

  n_ts = length(sigma_wn)

  WN = matrix(0, nrow = n_ts, ncol = N)
  for (i in 1:n_ts) {
    WN[i,] = rnorm(N, mean = 0, sd = sqrt(sigma_wn[i]))
  }

  INN = mvrnorm(N, rep(0,n_ts), sigma_rw)

  RW = apply(INN, 2, cumsum)

  X = WN + t(RW)

  return(t(X))
}

#' Simulates a WN + coorelated RW from a T-Student distribution
#'
#' @param N the number of samples
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#'
#' @export

simulate_wn_corr_rw_ts = function(N, sigma_wn, sigma_rw, df = 5) {

  n_ts = length(sigma_wn)

  WN = rmvt(N, sigma = diag(sigma_wn), df = df)

  INN = rmvt(N, sigma_rw, df = df)

  RW = apply(INN, 2, cumsum)

  X = WN + RW

  return(X)
}

#' Simulates a WN + AR1 + coorelated RW
#'
#' @param N the number of samples
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param phi_ar1 a vector of size \code{n}, the phi of the AR1s
#' @param sigma_ar1 a vector of size \code{n}, the variance of AR1s innovations
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#'
#' @export

simulate_wn_ar1_corr_rw = function(N, sigma_wn, phi_ar1, sigma_ar1, sigma_rw) {

  n_ts = length(sigma_wn)

  WN = matrix(0, nrow = n_ts, ncol = N)
  for (i in 1:n_ts) {
    WN[i,] = rnorm(N, mean = 0, sd = sqrt(sigma_wn[i]))
  }

  INN_AR1 = matrix(0, nrow = n_ts, ncol = N)
  for (i in 1:n_ts) {
    INN_AR1[i,] = rnorm(N, mean = 0, sd = sqrt(sigma_ar1[i]))
  }

  AR1 = matrix(0, nrow = n_ts, ncol = N)
  AR1[,1] = INN_AR1[,1]
  for (i in 2:N) {
    AR1[,i] = AR1[,i-1] * phi_ar1 + INN_AR1[,i]
  }

  INN = mvrnorm(N, rep(0,n_ts), sigma_rw)

  RW = apply(INN, 2, cumsum)

  X = WN + t(RW) + AR1

  return(t(X))
}

#' Simulates a WN + correlated AR1 + coorelated RW
#'
#' @param N the number of samples
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param phi_ar1 a vector of size \code{n}, the phi of the AR1s
#' @param sigma_ar1 a matrix of size \code{n}x\code{n}, the variance of AR1s innovations
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#'
#' @export

simulate_wn_corr_ar1_corr_rw = function(N, sigma_wn, phi_ar1, sigma_ar1, sigma_rw) {

  n_ts = length(sigma_wn)

  WN = matrix(0, nrow = n_ts, ncol = N)
  for (i in 1:n_ts) {
    WN[i,] = rnorm(N, mean = 0, sd = sqrt(sigma_wn[i]))
  }

  INN_AR1 = t(mvrnorm(N, rep(0,n_ts), sigma_ar1))

  AR1 = matrix(0, nrow = n_ts, ncol = N)
  AR1[,1] = INN_AR1[,1]
  for (i in 2:N) {
    AR1[,i] = AR1[,i-1] * phi_ar1 + INN_AR1[,i]
  }

  INN = mvrnorm(N, rep(0,n_ts), sigma_rw)

  RW = apply(INN, 2, cumsum)

  X = WN + t(RW) + AR1

  return(t(X))
}

#' Theoretical wavelet cross-covariance of a WN + coorelated RW
#'
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#' @param J maximum scale, so that tau = 2^(1:J)
#'
#' @export

theo_wccv_wn_corr_rw = function(sigma_wn, sigma_rw, J){
  num.ts = length(sigma_wn)
  wccv.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  counter = 0
  for(i in 1:(num.ts-1)){
    for(j in (i+1):num.ts){
      counter = counter+1
      for(k in 1:J){
        wccv.mat[counter,] = expect.wccv_RW(sigma_rw[i,j], J)
      }
    }
  }
  list(wccv = t(wccv.mat), wccv.cov = NA, ci_low = NA, ci_high = NA)
}

#' Theoretical wavelet cross-covariance of a WN + correlated AR1 + coorelated RW
#'
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param phi_ar1 a vector of size \code{n}, the phi of the AR1s
#' @param sigma_ar1 a matrix of size \code{n}x\code{n}, the variance of AR1s innovations
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#' @param J maximum scale, so that tau = 2^(1:J)
#'
#' @export

theo_wccv_wn_corr_ar1_corr_rw = function(sigma_wn, phi_ar1, sigma_ar1, sigma_rw, J){
  num.ts = length(sigma_wn)
  wccv.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  counter = 0
  for(i in 1:(num.ts-1)){
    for(j in (i+1):num.ts){
      counter = counter+1
      for(k in 1:J){
        wccv.mat[counter,] =
          expect.wccv_AR(phi_ar1[i], phi_ar1[j], sigma_ar1[i,j], J) +
          expect.wccv_RW(sigma_rw[i,j], J)
      }
    }
  }
  list(wccv = t(wccv.mat), wccv.cov = NA, ci_low = NA, ci_high = NA)
}

#' Theoretical wavelet cross-covariance of a WN + N correlated AR1 + coorelated RW
#'
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param phi_ar1 a list of vectors of size \code{n}, the phi of the AR1s
#' @param sigma_ar1 a list of matrices of size \code{n}x\code{n}, the variance of AR1s innovations
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#' @param J maximum scale, so that tau = 2^(1:J)
#'
#' @export

theo_wccv_wn_corr_Nar1_corr_rw = function(sigma_wn, phi_ar1, sigma_ar1, sigma_rw, J){
  num.ts = length(sigma_wn)
  wccv.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  counter = 0
  for(i in 1:(num.ts-1)){
    for(j in (i+1):num.ts){
      counter = counter+1
      for(k in 1:J){
        wccv.mat[counter,] = expect.wccv_RW(sigma_rw[i,j], J)

        for (m in seq_along(phi_ar1)) {
          wccv.mat[counter,] = wccv.mat[counter,] + expect.wccv_AR(phi_ar1[[m]][i], phi_ar1[[m]][j], sigma_ar1[[m]][i,j], J)
        }
      }
    }
  }
  list(wccv = t(wccv.mat), wccv.cov = NA, ci_low = NA, ci_high = NA)
}

#' Theoretical wavelet variance of a WN  + coorelated RW
#'
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#' @param J maximum scale, so that tau = 2^(1:J)
#'
#' @export

theo_wv_wn_corr_rw = function(sigma_wn, sigma_rw, J){
  num.ts = length(sigma_wn)
  wv.mat = matrix(NA, num.ts, J)
  for(i in 1:num.ts){
    wv.mat[i,] = expect.wv_WN(sigma_wn[i], J) + expect.wv_RW(sigma_rw[i,i], J)
  }
  t(wv.mat)
}

#' Theoretical wavelet variance of a WN + correlated AR1 + coorelated RW
#'
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param phi_ar1 a vector of size \code{n}, the phi of the AR1s
#' @param sigma_ar1 a matrix of size \code{n}x\code{n}, the variance of AR1s innovations
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#' @param J maximum scale, so that tau = 2^(1:J)
#'
#' @export

theo_wv_wn_corr_ar1_corr_rw = function(sigma_wn, phi_ar1, sigma_ar1, sigma_rw, J){
  num.ts = length(sigma_wn)
  wv.mat = matrix(NA, num.ts, J)
  for(i in 1:num.ts){
    wv.mat[i,] = expect.wv_WN(sigma_wn[i], J) +
                 expect.wv_AR(phi_ar1[i], sigma_ar1[i,i], J) +
                 expect.wv_RW(sigma_rw[i,i], J)
  }
  t(wv.mat)
}

#' Theoretical wavelet variance of a WN + N correlated AR1 + coorelated RW
#'
#' @param sigma_wn a vector of size \code{n} of WN variances for each time series
#' @param phi_ar1 a list of vectors of size \code{n}, the phi of the AR1s
#' @param sigma_ar1 a list of matrices of size \code{n}x\code{n}, the variance of AR1s innovations
#' @param sigma_rw a matrix of size \code{n}x\code{n}, the covariance of RW innovations
#' @param J maximum scale, so that tau = 2^(1:J)
#'
#' @export

theo_wv_wn_corr_Nar1_corr_rw = function(sigma_wn, phi_ar1, sigma_ar1, sigma_rw, J){
  num.ts = length(sigma_wn)
  wv.mat = matrix(NA, num.ts, J)
  for(i in 1:num.ts){
    wv.mat[i,] = expect.wv_WN(sigma_wn[i], J) + expect.wv_RW(sigma_rw[i,i], J)

    for (j in seq_along(phi_ar1)) {
      wv.mat[i,] = wv.mat[i,]+expect.wv_AR(phi_ar1[[j]][i], sigma_ar1[[j]][i,i], J)
    }
  }
  t(wv.mat)
}
