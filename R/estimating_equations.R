# Functions ----

#' Title generate regimes
#'
#' @param time_points number of total randomization points within a day
#' @param total_doses number of doses given in a day
#'
#' @return a matrix of all possible dose combinations
generate_regimes <- function(time_points, total_doses) {
  if (total_doses == 0) {
    return(matrix(0, nrow = 1, ncol = time_points))
  }
  # Generate all combinations of time points taken `total_doses` at a time
  combinations <- combn(time_points, total_doses)
  n_combinations <- ncol(combinations)
  # Initialize the dose matrix
  dose_matrix <- matrix(0, nrow = n_combinations, ncol = time_points)
  # Fill the dose matrix based on combinations
  for (i in 1:n_combinations) {
    dose_matrix[i, combinations[, i]] <- 1
  }
  return(dose_matrix)
}

#' Title get regime probability
#'
#' @param a matrix of treatments (n x time points)
#' @param p probability of single treatment
#'
#' @return a list of cumulative dose, probability of a_t, and dose weight
get_p_a <- function(a, p) {
  # a: matrix of treatments
  # p: probability of single treatment
  # return:
  #   cumulative dose
  #   probability of a_t
  cum_d <- matrix(0, nrow = nrow(a), ncol = ncol(a))
  p_a <- matrix(nrow=nrow(a), ncol = ncol(a))
  d_w <- matrix(0, nrow=nrow(a), ncol = ncol(a))
  for (decision in 1:ncol(a)) {
    cum_d[,decision] <- if(decision==1) {as.matrix(a[,1:decision, drop=FALSE])} else {rowSums(a[,1:decision])}
    p_a[, decision] <- (p^cum_d[,decision]) * ((1-p)^(decision - cum_d[,decision]))
    d_w[, decision] <- (cum_d[,decision]==1) * 1/5 + (cum_d[,decision]==0) * (5-decision)/5
  }
  return(list(cum_d = cum_d, p_a = p_a, d_w = d_w))
}


ee6.4.improved <- function( beta, y, a, h, s, p_a, cum_d, dose ) {
  T.dp <- ncol(y)
  a_5 <- generate_regimes(ncol(a), dose)
  #U <- matrix(0, nrow = ncol(s), ncol = 1)
  U <- 0

  zeros <- lapply(1:5, rep, x = 0)

  for (decision in 1:T.dp) {
    df <- as.data.frame(cbind(y=y[,decision], h, a=a[,decision], cum_d=cum_d[,decision]))

    # Fit is too slow
    fit <- lm(y ~ -1 + h1 + h2 + a, data = df)
    summary(fit)

    df.m0 <- df
    df.m1 <- df

    for (regime in 1:nrow(a_5)) {
      df.m1$a <- a_5[regime, decision]
      df.m0$a <- 0

      m1 <- predict(fit, newdata = df.m1)
      m0 <- predict(fit, newdata = df.m0)

      # I mean, even faster here!
      I_at <- rowSums(a[, 1:decision, drop=FALSE] == matrix(a_5[regime, 1:decision], nrow=nrow(a), ncol=decision, byrow=TRUE)) == decision
      I_0t <- rowSums(a[, 1:decision, drop=FALSE] == matrix(zeros[[decision]], nrow=nrow(a), ncol=decision, byrow=TRUE)) == decision

      U <- U + t(s) %*% (I_at/p_a[,decision] * (y[,decision] - m1) -
                           I_0t/p_a[,decision] * (y[,decision] - m0) +
                           m1 - m0 - ( as.matrix(s) %*% beta )/T.dp)
    }
  }
  U <- U / nrow(y)
  return(U)
}


# Improved I_0t calculation by using zeros list
#' Title
#'
#' @param beta initial beta
#' @param y (n x time points) outcomes
#' @param a (n x time points) treatments
#' @param h (n x p) baseline covariates
#' @param s (n x q) effect modifying covariates
#' @param p_a (n x time points) probability of treatment
#' @param cum_d (n x time points) cumulative dose
#' @param dose (scalar) dose of interest
#'
#' @return U (scalar) estimating equation sum
#' @export
ee.cor.2 <- function( beta, y, a, h, s, p_a, cum_d, dose ) {
  T.dp <- ncol(y)
  a_5 <- generate_regimes(ncol(a), dose)
  #U <- matrix(0, nrow = ncol(s), ncol = 1)
  U <- 0
  zeros <- lapply(1:5, rep, x = 0)

  for (decision in 1:T.dp) {
    s.decision <- cbind(1, s[,decision])
    df <- data.frame("y"=y[,decision, drop=TRUE], "h"=h, "s"=s[,decision, drop=TRUE], "a"=a[,decision, drop=TRUE], "cum_d"=cum_d[,decision, drop=TRUE])
    fit <- lm(y ~  h + s + a, data = df)
    #X<-model.matrix(~h + a + s*a,data = df)
    #fit <- lm.fit(y =df$y,x=X) # included intercept

    df.m0 <- df
    df.m1 <- df

    df.m0$a <- 0

    a.interest <- a[, 1:decision, drop=FALSE]

    for (regime in 1:nrow(a_5)) {

      df.m1$a <- a_5[regime, decision]

      m1 <- predict(fit, newdata = df.m1)
      m0 <- predict(fit, newdata = df.m0)

      I_at <- rowSums(a.interest == matrix(a_5[regime, 1:decision], nrow=nrow(a), ncol=decision, byrow=TRUE)) == decision
      I_0t <- rowSums(a.interest == matrix(zeros[[decision]], nrow=nrow(a), ncol=decision, byrow=TRUE)) == decision

      U <- U + t(s.decision) %*% (I_at/p_a[,decision] * (y[,decision] - m1) -
                                    I_0t/p_a[,decision] * (y[,decision] - m0) +
                                    m1 - m0 - ( as.matrix(s.decision) %*% beta )/T.dp)
    }
  }
  U <- U / nrow(y)
  return(U)
}


#' Get Simulated Results
#'
#' @param m number of simulated datasets
#' @param dfs list of simulated datasets
#' @param print_progress print progress if true
#'
#' @return matrix with 2 rows: estimated mean, estimated standard deviation
#' @export
#'
#' @examples
#' df.corr <- simMhealth(m=100, n=37, time=5, days=42,
#' eta=1, rho=0.5,
#' theta1=1, theta2=0.8,
#' beta10=0.5, beta11=0, beta12=0.2,
#' p=0.5)
#' rslt <- get.simulated.rslts(m=100, dfs=df.corr, print_progress=TRUE)
get.sim.rslts <-function(m, dfs, dose, print_progress=FALSE) {
  ee.corr <- matrix(nrow = m, ncol = 2)
  for (rep in 1:m) {
    if ((rep %% 10 == 0) & (print_progress)) {
      cat(rep, "\n")
    }

    df.wide <- dfs[[rep]] |>
      dplyr::select(ID,DAY,SLOT,Y,A,H,S) |>
      dplyr::group_by(ID,DAY) |>
      tidyr::pivot_wider(names_from=SLOT, values_from=c(Y,A,S)) |>
      as.matrix()

    y <- df.wide[,c("Y_1", "Y_2", "Y_3", "Y_4", "Y_5")]
    a <- df.wide[,c("A_1", "A_2", "A_3", "A_4", "A_5")]
    h <- df.wide[,"H", drop=FALSE]
    s <- df.wide[,c("S_1", "S_2", "S_3", "S_4", "S_5")]
    p <- 0.5

    matrix <- get_p_a(a, p)
    cum_d <- matrix$cum_d
    p_a <- matrix$p_a

    init_beta <- c(0,0)
    ee.cor.2( init_beta, y, a, h, s, p_a, cum_d, dose )
    rslt <- nleqslv::nleqslv(init_beta, function(beta) ee.cor.2( beta, y, a, h, s, p_a, cum_d, dose ))
    ee.corr[rep,] <- rslt$x
  }
  return(rbind(colMeans(ee.corr), apply(ee.corr, 2, sd)))
}



# CeeDose <- function(id, y, trt, baseline, timevar, baselineInter, timevarInter, )












# # Boruvka ----
# boruvka <- function(beta, y, a, h, s, dose) {
#   beta1 <- beta[2]
#   alpha <- beta[1]
#   rho.hat = mean(a==dose)
#   a <- a==dose
#   p_a <- dbinom(dose, 5, 0.5)
#   p_0 <- dbinom(0, 5, 0.5)
#   p <- p_a / (p_a+p_0)
#   w = w = (rho.hat^a * (1-rho.hat)^(1-a)) / (p^a * (1-p)^(1-a))
#   U = 0
#   resid <- y - (a - rho.hat) * beta
#   resid.neg <- ifelse(resid <0, 1,0)
#   U <- U + sum((y - alpha - (a - rho.hat) * beta1) * w * cbind(rep(1,m), (a - rho.hat)))
#
#   #print(U)
#   U <- as.vector(U / length(y))
#   return(U)
# }
#
# boruvka.covar <- function(beta, y, a, h, s, p_a, dose) {
#   rho.hat = mean(a==dose)
#   a <- a==dose
#   p_a <- dbinom(dose, 5, 0.5)
#   w = w = (rho.hat^a * (1-rho.hat)^(1-a)) / (p_a^a * (1-p_a)^(1-a))
#   #U = c(0,0,0)
#   U = 0
#   for (i in 1:length(y)) {
#     #resid <- y[i] - (a[i] - rho.hat) * as.matrix(s[1,]) %*% beta
#     resid <- y[i] - (a[i] - rho.hat) * beta
#     weight <- w[i]
#     s_beta <- (a[i] - rho.hat) # * s[i,]
#     U <- U + as.matrix(resid * weight * s_beta)
#     #print(as.matrix(resid * weight * s_beta))
#     #U <- U + as.matrix(( y[i] - (a[i] - rho.hat) * as.matrix(s[1,]) %*% beta ) * w[i] * s[i,] * (a[i] - rho.hat))
#   }
#   #print(U)
#   U <- as.vector(U / length(y))
#   return(U)
# }
#







