#' @export

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

### Transform Weights

inv_logit = function(x){
  exp(x)/(1 + exp(x))
}

transform_weights = function(kappa_R){
  n = length(kappa_R)
  kappa = rep(NA, n)
  sum_available = 1
  for (i in 1:n){
    kappa[i] = sum_available * inv_logit(kappa_R[i])
    sum_available = 1 - sum(kappa[1:i])
  }
  kappa
}

obj.fun.theo.wts = function(kappa_R, W, scale_weights){
  size = dim(W)
  kappa = transform_weights(kappa_R = kappa_R)
  weights = c(kappa, 1 - sum(kappa))
  inter = 0
  for (j in 1:size[1]){
    for (i in 1:size[2]){
      for (k in 1:size[3]){
        inter = inter + scale_weights[j]*weights[i]*weights[k]*W[j,i,k]
      }
    }
  }
  inter
}

#' iterative method for optimal coefs
#'
#' @export

find_optimal_coefs_iterative = function(Xt, scale_weights){
  # setting
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1

  # wccv
  wcov = wccv(Xt)
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

  # W[,1,1] = wvar(Xt[,1])$variance
  # W[,1,2] = wcov.mat[,1]
  # W[,1,3] = wcov.mat[,2]
  # W[,1,4] = wcov.mat[,3]
  #
  # W[,2,1] = wcov.mat[,1]
  # W[,2,2] = wvar(Xt[,2])$variance
  # W[,2,3] = wcov.mat[,4]
  # W[,2,4] = wcov.mat[,5]
  #
  # W[,3,1] = wcov.mat[,2]
  # W[,3,2] = wcov.mat[,4]
  # W[,3,3] = wvar(Xt[,3])$variance
  # W[,3,4] = wcov.mat[,6]
  #
  # W[,4,1] = wcov.mat[,3]
  # W[,4,2] = wcov.mat[,5]
  # W[,4,3] = wcov.mat[,6]
  # W[,4,4] = wvar(Xt[,4])$variance

  # compute coefficient
  # NOTE: c(-1.1, -0.7, 0) below corresponds to initial values of the iterative algo to find coefficients,
  # depending on your data, you may want to adjust it to avoid computation errors.
  res_kappa = transform_weights(optim(c(-1.1, -0.7, 0), obj.fun.theo.wts, W = W, scale_weights = scale_weights)$par)

  coefs = c(res_kappa, 1 - sum(res_kappa))
  return(coefs)
}

#' Construct the matrix A hat

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

#' closed-form solution to optimal coefs
#'
#' @export

find_optimal_coefs = function(Xt, scale_weights){
  num.ts = ncol(Xt)
  ones = rep(1, num.ts)

  A = get_A(Xt = Xt, scale_weights = scale_weights)
  A_inv = solve(A)

  res = A_inv %*% ones / as.numeric(t(ones)%*%A_inv%*%ones)
  res = as.numeric(res)
  return(res)
}


#' @export
get_virtual_gyro = function(Xt, weights){
  Xt %*% weights
}


get_A_theo = function(wccv.mat, wv.mat, scale_weights){
  # setting
  num.ts = dim(wv.mat)[2]
  J = dim(wv.mat)[1]

  # set up W array
  W = array(NA, c(J, num.ts, num.ts))

  counter = 1
  for (i in 1:num.ts) {
    for (j in i:num.ts) {
      if(i == j){
        W[,i,i] = wv.mat[,i]
      }else{
        W[,i,j] = wccv.mat[,counter]
        W[,j,i] = wccv.mat[,counter]
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

#' closed-form solution to optimal coefs from theoretical wvs and wccvs
#'
#' @export

find_optimal_coefs_theo = function(wccv.mat, wv.mat, scale_weights){
  num.ts = dim(wv.mat)[2]
  ones = rep(1, num.ts)

  A = get_A_theo(wccv.mat, wv.mat, scale_weights)
  A_inv = solve(A)

  res = A_inv %*% ones / as.numeric(t(ones)%*%A_inv%*%ones)
  res = as.numeric(res)
  return(res)
}


#' get estimated covariance matrix for the coefficients with block bootstrap
#'
#' @export

get_c_var = function(Xt, scale_weights, c_hat, B = 10^3, sB = 10){
  num.ts = ncol(Xt)
  N = nrow(Xt)
  c_hat_mat = matrix(NA, nrow = B, ncol = num.ts)

  for (i in 1:B) {
    c_hat_star = get_c_hat_star(Xt, scale_weights, sB)
    c_hat_mat[i,] = sqrt(N)*(c_hat_star - c_hat)
  }

  cov(c_hat_mat)
}

get_c_hat_star = function(Xt, scale_weights, sB){
  num.ts = ncol(Xt)
  N = nrow(Xt)
  J = floor(log2(N)) - 1

  # ----- block bootstrap on wavelet coefficients

  m = N - 2^1 + 1         # length of first scale
  h = sB*floor(m^(1/3))   # size of block  # range 0.1 0.5 1 5 10 50 100
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


  # ----- compute A
  A = matrix(NA, nrow = num.ts, ncol = num.ts)
  for (i in 1:num.ts) {
    for (k in 1:num.ts) {
      A[i,k] = sum(scale_weights *W[,i,k])
    }
  }

  A_inv = solve(A)

  # ----- compute c_hat_star
  one = rep(1, num.ts)
  c_hat = as.numeric(A_inv %*% one) / as.numeric(t(one)%*%A_inv%*%one)

  return(c_hat)
}

# ----------- second method for coefficient covariance

#' @export

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

#' @export

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


#' @export

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

