# boruvka simulation
library(MASS)

#' expit function
#'
#' @param x numeric
#'
#' @return expit(x)
expit <- function(x) {1/(1+exp(-x))}

#' Get correlation matrix
#'
#' @param n integer
#'
#' @return (n x n) corr_matrix
create_corr_matrix <- function(n) {
  corr_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      corr_matrix[i, j] <- 0.5^(abs(i - j) / 2)
    }
  }
  return(corr_matrix)
}

#' Boruvka estimating equations
#'
#' @param param (numeric vector)
#' @param Y (numeric matrix)
#' @param S (numeric matrix)
#' @param A (numeric matrix)
#' @param rho.hat (numeric) Correlation
#' @param w (numeric matrix) Weights
#'
#' @return U (numeric)
#' @export
ee.boruvka <- function(param, Y, S, A, rho.hat, w) {
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

#' GEE independent estimating equation
#'
#' @param param (numeric vector)
#' @param Y (numeric matrix)
#' @param S (numeric matrix)
#' @param A (numeric matrix)
#' @param rho.hat (numeric) Correlation
#' @param w (numeric matrix) Weights
#'
#' @return U (numeric)
#' @export
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

#' Title
#'
#' @param At_1 (numeric) previous treatment
#' @param At (numeric) current treatment
#' @param St_1 (vector) previous covariate
#' @param St (vector) current covariate
#' @param xi (numeric)
#' @param eta1 (numeric) effect of baseline
#' @param eta2 (numeric) effect of time-varying covariate
#' @param theta1 (numeric) effect of something
#' @param theta2 (numeric) effect of previous treatment
#' @param beta10 (numeric) effect of treatment
#' @param beta11 (numeric) effect of interaction
#'
#' @return (numeric) simulated value
#' @export
simulate.boruvka <- function(At_1, At, St_1, St, xi, eta1, eta2, theta1, theta2, beta10, beta11) {
  E.St <- expit(xi * At_1) - (1 - expit(xi * At_1))
  #pt_1 <- expit(eta1 * At_2 + eta2 * St_1)
  pt <- expit(eta1 * At_1 + eta2 * St)
  return(theta1 * ( St - E.St ) +
         #theta2 * ( At_1 - pt_1 ) +
         ( At - pt ) * (beta10 + beta11 * St))
}

