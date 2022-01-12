#' @title Empirical Wavelet Cross-Covariance
#' @description This function provides an empirical estimate of the wavelet cross-covariance given multiple processes.
#' @export
#' @param Xt A \code{matrix} of dimension T by p, where T is the length of the time series and p is the number of processes.
#' @return A \code{list} with the following structure:
#' \itemize{
#'  \item wccv: A \code{matrix} of the estimated wavelet cross-covariance.
#'  \item wccv.cov: A \code{matrix} of the estimated covariance matrix of the estimated wavelet cross-covariance.
#'  \item ci_low: A \code{matrix} of the lower bound of the confidence interval of the estimated wavelet cross-covariance.
#'  \item ci_high: A \code{matrix} of the upper bound of the confidence interval of the estimated wavelet cross-covariance.
#' }
#' @importFrom wv wvar modwt
#' @author Haotian Xu

wccv_local = function(Xt){
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1
  wccv.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  cov.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  ci.low.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  ci.high.mat = matrix(NA, num.ts*(num.ts-1)/2, J)
  counter = 0
  for(i in 1:(num.ts-1)){
    coe1 = modwt(Xt[,i])
    for(j in (i+1):num.ts){
      counter = counter+1
      coe2 = modwt(Xt[,j])
      wv1 = wvar(Xt[,i])
      wv1_bar = (wv1$ci_high - wv1$ci_low)^2/4
      wv2 = wvar(Xt[,j])
      wv2_bar = (wv2$ci_high - wv2$ci_low)^2/4
      wv_min = apply(cbind(wv1_bar,wv2_bar), 1, min)
      for(k in 1:J){
        y = unlist(coe1[k])
        z = unlist(coe2[k])
        wccv.mat[counter,k] = mean((y) * (z))
        cov.mat[counter,k] = wv_min[k]
        ci.low.mat[counter,k] =  wccv.mat[counter,k] + qnorm(0.025) * sqrt(cov.mat[counter,k])
        ci.high.mat[counter,k] = wccv.mat[counter,k] + qnorm(1-0.025) * sqrt(cov.mat[counter,k])
      }
    }
  }
  list(wccv = t(wccv.mat), wccv.cov = t(cov.mat), ci_low = t(ci.low.mat), ci_high = t(ci.high.mat))
}



# Computes the matrix \hat{A} in Zhang et al. (2021)
get_A = function(Xt, scale_weights){
  # setting
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1

  # wccv
  wcov = wccv_local(Xt)
  wcov.mat = wcov$wccv
  wcov.cov.mat = wcov$wccv.cov

  # set up W array
  W = array(NA, c(J, num.ts, num.ts))

  counter = 1
  for (i in 1:num.ts) {
    for (j in i:num.ts) {
      if(i == j){
        W[,i,i] = wvar(Xt[,i])$variance
      }else{
        W[,i,j] = wcov.mat[,counter]
        W[,j,i] = wcov.mat[,counter]
        counter = counter + 1
      }
    }
  }

  # build A
  A = matrix(NA, nrow = num.ts, ncol = num.ts)
  for (i in 1:num.ts) {
    for (k in 1:num.ts) {
      A[i,k] = sum(scale_weights *W[,i,k])
    }
  }

  return(A)
}

#' @title Estimate Optimal Coefficients based on Scale-wise Variance Optimization
#' @description This function computes the estimated optimal coefficients based on the Scale-wise Variance Optimization approach.
#' The detailed definition can be found in Equation (7) in Zhang et al. (2021) (https://arxiv.org/abs/2106.15997).
#' @export
#' @param Xt A \code{matrix} of dimension T by p, where T is the length of the time series and p is the number of processes.
#' @param scale_weights A \code{vector} that denotes the weights on scales. All elements should be non-negative and sum to one.
#' @return A \code{vector} of the estimated optimal coefficients on the p individual processes.
#' @author Yuming Zhang

find_optimal_coefs = function(Xt, scale_weights){
  num.ts = ncol(Xt)
  ones = rep(1, num.ts)

  A = get_A(Xt = Xt, scale_weights = scale_weights)
  A_inv = solve(A)

  res = A_inv %*% ones / as.numeric(t(ones)%*%A_inv%*%ones)
  res = as.numeric(res)
  return(res)
}


#' @title Construct Virtual Gyroscope Signal
#' @description This function computes the virtual gyroscope signal by taking a linear combination of the individual gyroscope signals.
#' @export
#' @param Xt A \code{matrix} of dimension T by p, where T is the length of the time series and p is the number of processes.
#' @param weights A \code{vector} of coefficients on individual signals. All elements should be non-negative and sum to one.
#' @return A \code{vector} of the virtual gyroscope signal.
#' @author Yuming Zhang

get_virtual_gyro = function(Xt, weights){
  Xt %*% weights
}




# This function computes empirical WCCV and formulate it as a vector.
get_wccv = function(Xt){
  # setting
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1

  # wccv
  wcov = wccv_local(Xt)
  wcov.mat = wcov$wccv
  wcov.cov.mat = wcov$wccv.cov

  # set up W array
  W = array(NA, c(J, num.ts, num.ts))

  counter = 1
  for (i in 1:num.ts) {
    for (j in i:num.ts) {
      if(i == j){
        W[,i,i] = wvar(Xt[,i])$variance
      }else{
        W[,i,j] = wcov.mat[,counter]
        W[,j,i] = wcov.mat[,counter]
        counter = counter + 1
      }
    }
  }

  # wccv vector
  wccv_vec = c()
  for (i in 1:num.ts) {
    for (j in 1:num.ts) {
      wccv_vec = c(wccv_vec, W[,i,j])
    }
  }

  out = list(W = W,
             wccv_vec = wccv_vec)
  return(out)
}

# This function block bootstraps the wavelet coefficients.
# sB: constant related to the block size
block_boot_wave_coef = function(Xt, sB = 10){
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1

  # ----- block bootstrap on wavelet coefficients

  m = N - 2^1 + 1         # length of first scale
  h = sB*floor(m^(1/3))   # size of block
  nb_block = floor(m/h)
  missing = m - nb_block*h
  index_block = 1:h

  # fill all scales so that all scales have m wavelet coefficients

  wave_coef_mat = matrix(NA, nrow = m, ncol = num.ts*J)
  index_last_vec = rep(NA, J)

  for (j in 2:J) {
    index_last_vec[j] = sample(m-2*(2^j-2)+1, 1) # should use the same for all time series at each scale
  }

  for (i in 1:num.ts) {
    Xtmodwt = modwt(Xt[,i])
    wave_coef_mat[,1+(i-1)*J] = Xtmodwt[[1]]

    for (j in 2:J) {
      wave_coef_mat[1:(m-2^j+2), j+(i-1)*J] = Xtmodwt[[j]]
      index_last = index_last_vec[j]
      wave_coef_mat[(m-2^j+3):m, j+(i-1)*J] = wave_coef_mat[index_last+(0:(2^j-3)), j+(i-1)*J]
    }
  }

  # resample

  start_index = sample(0:missing, 1)
  index = sample(nb_block, nb_block, replace = TRUE)
  wave_coef_mat_star = matrix(NA, nrow = m, ncol = num.ts*J)

  for (l in 1:length(index)){
    wave_coef_mat_star[(l-1)*h + index_block, ] = wave_coef_mat[start_index + (index[l]-1)*h + index_block, ]
  }
  if (missing > 0){
    last_index = sample(m - missing, 1)
    wave_coef_mat_star[(nb_block*h + 1):m, ] = wave_coef_mat[(last_index + 1):(last_index + missing), ]
  }


  # ----- set up W array
  W = array(NA, c(J, num.ts, num.ts))

  for (i in 1:num.ts) {
    for (k in 1:num.ts) {
      for (j in 1:J) {
        coe1 = wave_coef_mat_star[1:(m - 2^j + 2), j+(i-1)*J]
        coe2 = wave_coef_mat_star[1:(m - 2^j + 2), j+(k-1)*J]
        W[j,i,k] = cov(coe1, coe2)
      }
    }
  }

  # ----- wccv vector
  wccv_vec = c()
  for (i in 1:num.ts) {
    for (j in 1:num.ts) {
      wccv_vec = c(wccv_vec, W[,i,j])
    }
  }

  out = list(W = W,
             wccv_vec = wccv_vec)
  return(out)

}


#' @title Empirical Covariance of Coefficients on Individual Signals
#' @description This function computes the estimated covariance matrix of the coefficients on individual signals using the Moving Block Bootstrap approach considered in Zhang et al. (2021).
#' The detailed definition can be found in Equation (8) in Zhang et al. (2021) (https://arxiv.org/abs/2106.15997).
#' @export
#' @param Xt A \code{matrix} of dimension T by p, where T is the length of the time series and p is the number of processes.
#' @param scale_weights A \code{vector} that denotes the weights on scales. All elements should be non-negative and sum to one.
#' @param sB A \code{double} that denotes the positive constant C associated with the block size, which is defined as floor(C*T^{1/3}). Default value is 10.
#' @param B An \code{integer} indicating the number of Monte-Carlo replications used in the Moving Block Bootstrap.
#' @return A \code{matrix} of the estimated covariance matrix of the coefficients on individual signals.
#' @author Yuming Zhang

est_cov = function(Xt, scale_weights, sB = 10, B){
  # setting
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1

  # ----- Sigma is the asymptotic covariance matrix of the wccv vector

  # --- block bootstrap on the time series

  # num.blocks = N/size.blocks
  # index = (0:(num.blocks-1))*size.blocks+1  # index of first element of each block
  #
  # wccv_mat = matrix(NA, nrow = B, ncol = J*num.ts^2)
  #
  # for (i in 1:B) {
  #   vec = sample(index, length(index), replace = T)
  #   vec2 = c()
  #
  #   for (j in 1:length(vec)) {
  #     vec2 = c(vec2, vec[j]:(vec[j]+size.blocks-1))
  #   }
  #
  #   Xt2 = Xt[vec2,]
  #   wccv_mat[i,] = get_wccv(Xt2)$wccv_vec
  #
  #   print(i)
  # }
  #
  # Sigma = var(wccv_mat)

  # --- block bootstrap on the wavelet coefficient

  wccv_mat = matrix(NA, nrow = B, ncol = J*num.ts^2)
  for (i in 1:B) {
    tmp = block_boot_wave_coef(Xt, sB)
    wccv_mat[i,] = tmp$wccv_vec
    # print(i)
  }

  Sigma = cov(wccv_mat)

  # --- derivative in front of Sigma (Delta method)
  tmp =  get_wccv(Xt)
  W = tmp$W
  gamma_hat = tmp$wccv_vec

  A = matrix(NA, nrow = num.ts, ncol = num.ts)
  for (i in 1:num.ts) {
    for (k in 1:num.ts) {
      A[i,k] = sum(scale_weights *W[,i,k])
    }
  }

  A_inv = solve(A)
  one = rep(1, num.ts)
  x = A_inv %*% one
  part1 = 1/(sum(x))*diag(1, num.ts) - 1/(sum(x)^2) * x %*% t(one)

  for (i in 1:num.ts) {
    for (k in 1:num.ts) {
      for (j in 1:J) {
        A_by_gamma_i = matrix(0, num.ts, num.ts)
        A_by_gamma_i[i,k] = scale_weights[j]
        part2 = -A_inv %*% A_by_gamma_i %*% A_inv %*% one

        if(i==1 & k==1 & j==1){
          my_deriv = part1 %*% part2
        }else{
          my_deriv = cbind(my_deriv, part1 %*% part2)
        }

      }
    }
  }

  my_deriv %*% Sigma %*% t(my_deriv)

}

