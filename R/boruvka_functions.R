# boruvka simulation
library(MASS)

expit <- function(x) 1/(1+exp(-x))

create_corr_matrix <- function(n) {
  corr_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      corr_matrix[i, j] <- 0.5^(abs(i - j) / 2)
    }
  }
  return(corr_matrix)
}

ee <- function(param, Y, S, A, rho.hat, w) {
  alpha10 <- param[1]
  alpha11 <- param[2]
  beta1 <- param[3]
  U = 0
  for (i in 1:n) {
    for (t in 1:Time) {
      U <- U + ( Y[i,t+1] - (alpha10 + alpha11 * S[i,t]) - (A[i,t] - rho.hat) * beta1 ) * w[i,t] * c(1, S[i,t], (A[i,t] - rho.hat))
    }
  }
  return(U)
}

gee_ind <- function(param, Y, S, A, rho.hat, w) {
  alpha10 <- param[1]
  alpha11 <- param[2]
  beta1 <- param[3]
  U = 0
  for (i in 1:n) {
    for (t in 1:Time) {
      U <- U + ( Y[i,t+1] - (alpha10 + alpha11 * S[i,t]) - (A[i,t]) * beta1 ) * c(1, S[i,t], (A[i,t]))
    }
  }
  return(U)
}

simulate.boruvka <- function(At_1, At, St_1, St, xi, eta1, eta2, theta1, theta2, beta10, beta11) {
  E.St <- expit(xi * At_1) - (1 - expit(xi * At_1))
  #pt_1 <- expit(eta1 * At_2 + eta2 * St_1)
  pt <- expit(eta1 * At_1 + eta2 * St)
  return(theta1 * ( St - E.St ) +
         #theta2 * ( At_1 - pt_1 ) +
         ( At - pt ) * (beta10 + beta11 * St))
}

