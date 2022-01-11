
# This function computes the theoretical wavelet variance of a white noise.
# R: variance of the white noise
# J: maximum number of scales
expect.wv_WN = function(R, J){
  wv = 1/2^(1:J) * R
  wv
}

# This function computes the theoretical wavelet variance of a random walk.
# Q: variance of the innovation
# J: maximum number of scales
expect.wv_RW = function(Q, J){
  wv = (2^(2*(1:J))+2)/(12*2^(1:J)) * Q
  wv
}

# This function computes the theoretical wavelet variance of an AR(1) process.
# phi: the phi term of the AR(1) process
# Z: variance of the AR(1) process
# J: maximum number of scales
expect.wv_AR = function(phi, Z, J){
  wv = ((phi^2-1)*2^(1:J) + 2*phi*(phi^(2^(1:J)) - 4*phi^(2^(0:(J-1))) + 3))/((phi-1)^3 * (phi+1) * 4^(1:J)) * Z
  wv
}

# This function computes the theoretical wavelet cross-covariance of two white noise processes.
# R.cross: covariance of the white noise processes.
# J: maximum number of scales
expect.wccv_WN = function(R.cross, J){
  wccv = 1/2^(1:J) * R.cross
  wccv
}

# This function computes the theoretical wavelet variance of two random walk processes.
# Q.cross: covariance of the innovations
# J: maximum number of scales
expect.wccv_RW = function(Q.cross, J){
  wv = (2^(2*(1:J))+2)/(12*2^(1:J)) * Q.cross
  wv
}

# This function computes the theoretical wavelet variance of two AR(1) processes.
# phi1, phi2: the phi terms of the AR(1) processes
# Z.cross: covariance of the AR(1) processes
# J: maximum number of scales
expect.wccv_AR = function(phi1, phi2, Z.cross, J){
  wccv = (2^(1:J) + 2*(phi2/(1-phi2)*(2^(1:J)/2 - 1) - (phi2/(1-phi2))^2*(1-phi2^(2^(1:J)/2 - 1))) + 2*(phi1/(1-phi1)*(2^(1:J)/2 - 1) - (phi1/(1-phi1))^2*(1-phi1^(2^(1:J)/2 - 1))) - ((1-phi1^(2^(1:J)/2))/(1-phi1))^2 * phi1 - ((1-phi2^(2^(1:J)/2))/(1-phi2))^2 * phi2) / (1-phi1*phi2) * Z.cross / 4^(1:J)
  wccv
}


#' @title TO DO
#' @description WHAT IT DOES
#' @export
#' @param R        what it does
#' @param J        what it does
#' @return something
#' @details
#' If needed
#' @author WHO WROTE IT
#' @examples
#' # TO BE ADDED
#'
d.wccv_RW = function(J) {
  dwv = (2^(2*(1:J))+2)/(12*2^(1:J))
  dwv
}

obj.fun.wv = function(theta_tr, wv){

  theta = exp(theta_tr)

  R = theta[1]
  Q = theta[2]

  J = floor(log2(tail(wv$scales,1)))

  theo.wv = expect.wv_WN(R, J) + expect.wv_RW(Q, J)

  omega = diag(1 / (wv$ci_low - wv$ci_high)^2)

  dif = wv$variance - theo.wv
  t(dif) %*% omega %*% dif
}

obj.fun.wccv = function(theta, wccv, wccv.cov) {
  Q = theta[1]

  J = length(wccv)

  theo.wccv = expect.wccv_RW(Q, J)

  omega = diag(1 / wccv.cov)

  dif = wccv - theo.wccv
  t(dif) %*% omega %*% dif
}

estimate.gmwm.wv = function(theta0, wv){
  theta0_tr = log(theta0)
  theta_hat_tr = optim(theta0_tr, obj.fun.wv, wv = wv)$par

  exp(theta_hat_tr)
}

estimate.gmwm.wccv = function(limits, wccv, wccv.cov){
  optimize(obj.fun.wccv, interval = limits, wccv = wccv, wccv.cov = wccv.cov)$minimum
}

estimate.smwm.wccv.cf = function(wccv, wccv.cov) {
  J = length(wccv)
  Jc = matrix(d.wccv_RW(J), ncol=1)
  omega = diag(1 / wccv.cov)

  solve(t(Jc) %*% omega %*% Jc) %*% t(Jc) %*% omega %*% wccv
}

#' @export

estimate_wn_corr_rw = function(X) {

  wvs = apply(X, 2, wvar)
  wccvs = wccv_local(X)

  n_ts = dim(X)[2]

  R = rep(NA, n_ts)
  Q = matrix(NA, nrow = n_ts, ncol = n_ts)

  for (i in 1:n_ts) {
    init = c(1,1)
    init = c(sigma_wn[i], sigma_rw[i,i])

    theta_hat = estimate.gmwm.wv(init, wvs[[i]])
    R[i] = theta_hat[1]
    Q[i,i] = theta_hat[2]

    # J = floor(log2(tail(wvs[[i]]$scales,1)))
    # theo.wv = expect.wv_WN(R[i], J) + expect.wv_RW(Q[i,i], J)
    # plot(wvs[[i]]$variance, type="l", log="y")
    # lines(1:J, theo.wv, col="red")
  }

  cnt = 1
  for (i in 1:(n_ts-1)) {
    for (j in (i+1):n_ts) {

      # iterative
      #Q[i,j] = estimate.gmwm.wccv(c(-1,1)*sqrt(Q[i,i])*sqrt(Q[j,j]), wccv = wccvs$wccv[,cnt], wccv.cov = wccvs$wccv.cov[,cnt])

      # closed form
      Q[i,j] = estimate.smwm.wccv.cf(wccv = wccvs$wccv[,cnt], wccv.cov = wccvs$wccv.cov[,cnt])

      Q[j,i] = Q[i,j]

      # J = floor(log2(tail(wvs[[i]]$scales,1)))
      # theo.wv = expect.wccv_RW(Q[i,j], J)
      # plot(wccvs$wccv[,cnt], type="l")
      # lines(1:J, theo.wv, col="red")

      cnt = cnt+1
    }
  }

  return(
    list(
      R = R,
      Q = Q
    )
  )
}

#' @export

estimate_wn_corr_rw_theo = function(wvs, wccvs) {

  n_ts = length(wvs)

  R = rep(NA, n_ts)
  Q = matrix(NA, nrow = n_ts, ncol = n_ts)

  for (i in 1:n_ts) {
    init = c(1,1)
    init = c(sigma_wn[i], sigma_rw[i,i])

    theta_hat = estimate.gmwm.wv(init, wvs[[i]])
    R[i] = theta_hat[1]
    Q[i,i] = theta_hat[2]

    # J = floor(log2(tail(wvs[[i]]$scales,1)))
    # theo.wv = expect.wv_WN(R[i], J) + expect.wv_RW(Q[i,i], J)
    # plot(wvs[[i]]$variance, type="l", log="y")
    # lines(1:J, theo.wv, col="red")
  }

  cnt = 1
  for (i in 1:(n_ts-1)) {
    for (j in (i+1):n_ts) {

      # iterative
      #Q[i,j] = estimate.gmwm.wccv(c(-1,1)*sqrt(Q[i,i])*sqrt(Q[j,j]), wccv = wccvs$wccv[,cnt], wccv.cov = wccvs$wccv.cov[,cnt])

      # closed form
      Q[i,j] = estimate.smwm.wccv.cf(wccv = wccvs$wccv[,cnt], wccv.cov = wccvs$wccv.cov[,cnt])

      Q[j,i] = Q[i,j]

      # J = floor(log2(tail(wvs[[i]]$scales,1)))
      # theo.wv = expect.wccv_RW(Q[i,j], J)
      # plot(wccvs$wccv[,cnt], type="l")
      # lines(1:J, theo.wv, col="red")

      cnt = cnt+1
    }
  }

  return(
    list(
      R = R,
      Q = Q
    )
  )
}

