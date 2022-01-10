partial_inv <- function(A, n) {
  ret = matrix(0, nrow = dim(A)[1], ncol = dim(A)[2])

  ev = eigen(A, only.values = T)$values

  dec = svd(A)

  # for (i in which(ev > 0)) {
  for (i in 1:(length(dec$d)-n)) {
    ret = ret + dec$d[i]^-1 * dec$u[,i] %*% t(dec$v[,i])
  }

  t(ret)
}

#' @export
find_optimal_coefs_vaccaro = function(Q) {
  o = matrix(1, nrow=dim(Q)[1], ncol=1)

  c = ( solve(Q) %*% o ) / (t(o) %*% solve(Q) %*% o)[1]

  c
}

#' @export
find_optimal_coefs_vaccaro_pinv = function(Q, n) {
  o = matrix(1, nrow=dim(Q)[1], ncol=1)

  pinv = partial_inv(Q, n)

  c = ( pinv %*% o ) / (t(o) %*% pinv %*% o)[1]

  c
}


#' @export

av = function(X) {
  J = floor(log2(length(X))) - 3
  scales = 2^(1:J)

  variance = rep(0, length(scales))

  for (is in 1:length(scales)) {
    segs = seq(from=1, to=length(X), by=scales[is])

    segs_means = sapply( 1:(length(segs)-1),
                         function (iseg) {
                           mean(X[segs[iseg]:(segs[iseg+1]-1)])
                         })

    ds = diff(segs_means)^2

    variance[is] = mean(ds)/2
  }

  return(
    list(
      J = J,
      scales = scales,
      variance = variance
    )
  )
}

#' @export

acov = function(X) {
  n_ts = dim(X)[2]
  N = dim(X)[1]

  J = floor(log2(N)) - 3
  scales = 2^(1:J)

  covariance = list()

  for (is in 1:length(scales)) {
    covariance[[is]] = matrix(NA, nrow=n_ts, ncol=n_ts)

    segs_means = list()
    for (ix in 1:n_ts) {
      segs = seq(from=1, to=N, by=scales[is])

      segs_means[[ix]] = sapply( 1:(length(segs)-1),
                           function (iseg) {
                             mean(X[segs[iseg]:(segs[iseg+1]-1), ix])
                           })
    }

    diffs = lapply(segs_means, diff)

    for (ixa in 1:n_ts) {
      for (ixb in ixa:n_ts) {
        covariance[[is]][ixa,ixb] = mean(diffs[[ixa]]*diffs[[ixb]])/2
        covariance[[is]][ixb,ixa] = covariance[[is]][ixa,ixb]
      }
    }
  }

  return(
    list(
      J = J,
      scales = scales,
      covariance = covariance
    )
  )
}

CR = function(R, N, scales) {
  cr = matrix(NA, nrow = length(scales), ncol = length(scales))
  for (i_m1 in 1:length(scales)) {
    for (i_m2 in i_m1:length(scales)) {
      M1 = N/scales[i_m1]
      M2 = N/scales[i_m2]

      p = scales[i_m2]/scales[i_m1]

      cr[i_m1, i_m2] = (3*M2-4)/((M1-1)*(M2-1)*p^2)*R^2/scales[i_m1]^2
      cr[i_m2, i_m1] = cr[i_m1, i_m2]
    }
  }

  cr
}

CQ = function(Q, N, scales) {
  cq = matrix(NA, nrow = length(scales), ncol = length(scales))
  for (i_m1 in 1:length(scales)) {
    for (i_m2 in i_m1:length(scales)) {
      M1 = N/scales[i_m1]
      M2 = N/scales[i_m2]

      p = scales[i_m2]/scales[i_m1]

      cq[i_m1, i_m2] = (((12*p^3-6*p+3)*M2-2*(6*p^3-3*p+2))*Q^2*(scales[i_m1])^2)/(36*(M1-1)*(M2-1)*p^2)
      cq[i_m2, i_m1] = cq[i_m1, i_m2]
    }
  }

  cq
}

CCR = function(R1, R2, N, scales) {
  ccr = matrix(NA, nrow = length(scales), ncol = length(scales))
  for (i_m1 in 1:length(scales)) {
    for (i_m2 in i_m1:length(scales)) {
      M1 = N/scales[i_m1]
      M2 = N/scales[i_m2]

      p = scales[i_m2]/scales[i_m1]

      ccr[i_m1, i_m2] = (3*M2-4)/(2*(M1-1)*(M2-1)*p^2)*R1*R2/scales[i_m1]^2
      ccr[i_m2, i_m1] = ccr[i_m1, i_m2]
    }
  }

  ccr
}

CCQ = function(Q1, Q2, Q12, N, scales) {
  ccq = matrix(NA, nrow = length(scales), ncol = length(scales))
  for (i_m1 in 1:length(scales)) {
    for (i_m2 in i_m1:length(scales)) {
      M1 = N/scales[i_m1]
      M2 = N/scales[i_m2]

      p = scales[i_m2]/scales[i_m1]

      ccq[i_m1, i_m2] = (((12*p^3-6*p+3)*M2-2*(6*p^3-3*p+2))*(Q1*Q2+Q12^2)*(scales[i_m1])^2)/(72*(M1-1)*(M2-1)*p^2)
      ccq[i_m2, i_m1] = ccq[i_m1, i_m2]
    }
  }

  ccq
}

#' we assume that T = 1 everywhere
#'
#' @export

algorithm_1 = function(X) {
  allan = av(X)

  m0 = allan$scales[which.min(allan$variance)]

  i1 = which(allan$scales < m0/8)
  m1 = allan$scales[i1]
  N1 = tail(i1,1)

  Omega_Cm1 = solve(CR(1, length(X), m1))
  HR = matrix(1/2^(1:N1), nrow=N1, ncol=1)

  R0 = solve(t(HR) %*% Omega_Cm1 %*% HR) %*% t(HR) %*% Omega_Cm1 %*% allan$variance[i1]
  Q0 = 3*R0/m0^2

  Omega_C = solve(CR(R0, length(X), allan$scales) + CQ(Q0, length(X), allan$scales))

  HR = 2^(1:(allan$J))/3  # in the paper they say up to J-2 (non-comformable ...)
  HQ = 1/(2^(1:(allan$J)))
  H = cbind(HR, HQ)

  sol = solve(t(H) %*% Omega_C %*% H) %*% t(H) %*% Omega_C %*% allan$variance

  return(list
    (
      Q = sol[1],
      R = sol[2]
    ))
}

#' @export

algorithm_1_ms = function(X) {
  n_ts = dim(X)[2]

  R = rep(NA, n_ts)
  Q = rep(NA, n_ts)

  for (i in 1:n_ts) {
    ret = algorithm_1(X[,i])
    R[i] = ret$R
    Q[i] = ret$Q
  }

  return(list(
    R = R,
    Q = Q
  ))

}

#' @export

algorithm_2 = function(X) {
  if (log2(dim(X)[1]) %% 1 != 0) {
    X = X[1:(2^floor(log2(dim(X)[1]))),]
    warning(paste("truncating X to the first 2^", log2(dim(X)[1]),"samples"))
  }

  ret = algorithm_1_ms(X)

  R = ret$R
  Q = diag(ret$Q)

  allancov = acov(X)

  H = matrix(allancov$scales/3, nrow=length(allancov$scales), ncol=1)

  n_ts = dim(X)[2]
  N = dim(X)[1]

  for (ia in 1:n_ts) {
    for (ib in ia:n_ts) {
      if (ia == ib) {
        next
      }

      a = sapply(allancov$covariance, function(x) {x[ia,ib]})

      Omega_C = solve(CCR(R[ia], R[ib], N, allancov$scales) + CCQ(Q[ia,ia], Q[ib, ib], 0, N, allancov$scales))

      Q[ia, ib] = solve(t(H) %*% Omega_C %*% H) %*% t(H) %*% Omega_C %*% a
      Q[ib, ia] = Q[ia, ib]
    }
  }

  return(list(
    R = R,
    Q = Q
  ))
}

#' @export

make_healty = function(a) {
  e = eigen(a)

  lambda = e$values
  lambda[lambda < 0] = min(lambda[lambda > 0])

  e$vectors %*% diag(lambda) %*% t(e$vectors)
}
