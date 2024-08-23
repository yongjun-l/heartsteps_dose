# Functions ----
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
    cum_d[,decision] <- if(decision==1) {a[,1:decision]} else {rowSums(a[,1:decision])}
    p_a[, decision] <- (p^cum_d[,decision]) * ((1-p)^(decision - cum_d[,decision]))
    d_w[, decision] <- (cum_d[,decision]==1) * 1/5 + (cum_d[,decision]==0) * (5-decision)/5
  }
  return(list(cum_d = cum_d, p_a = p_a, d_w = d_w))
}

# Estimating Equations ----
ee0 <- function( beta, y, a, h, s, p_a, d_w, dose, p ) { # we don't consider any covariates
  #dose = 1
  #p = 0.5
  a_5 <- generate_regimes(ncol(a), dose)
  U <- matrix(0, nrow = nrow(a_5), ncol = ncol(y))
  for (decision in 1:ncol(y)) {
    #decision = 2
    for (regime in 1:nrow(a_5)) {
      #regime = 1
      I_at <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == a_5[regime,1:decision]))
      I_0t <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == rep(0, decision)))
      #w <- d_w[,decision] * I_at/p_a[,decision] - I_0t/p_a[,decision]
      w <- (1/5) * I_at/p_a[,decision] - I_0t/p_a[,decision]
      U[regime, decision] <- sum( w * y[,decision] - beta/5 )
      #U <- U + sum( (1/5) * w * ( y[,decision] ) )
    }
  }
  U <- sum(U) / nrow(y)
  return(U)
}

ee1 <- function( beta, y, a, h, s, p_a, cum_d, dose ) { # we don't consider any covariates

  T.dp <- ncol(y)
  a_5 <- generate_regimes(ncol(a), dose)
  U <- matrix(0, nrow = nrow(a_5), ncol = T.dp)
  for (decision in 1:T.dp) {
    for (regime in 1:nrow(a_5)) {
      #print(paste("Decision", decision, "Regime", regime))
      I_at <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == a_5[regime,1:decision]))
      I_0t <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == rep(0, decision)))
      p_a_noteq_0 <- 1 - mean(I_0t)
      w <- I_at/p_a[,decision] - I_0t/p_a[,decision]
      U[regime, decision] <- sum(w*y[,decision] - beta/T.dp)

      # print(cbind(I_at/p_a[,decision] * y[,decision],
      #             I_0t/p_a[,decision] * y[,decision],
      #             w*y[,decision]))
      # print(paste("sum(w*y[,decision]", sum(w*y[,decision])/10))
      # print(paste("beta/T.dp",beta/T.dp))
      # print(U/10)
      # print(paste("sum(U)",sum(U)))
      # print("\n\n")
    }
  }
  U <- sum(U) / nrow(y)
  return(U)
}

# we consider S
ee2 <- function( beta, y, a, h, s, p_a, cum_d, dose ) {
  #dose = 1
  #p = 0.5
  T.dp <- ncol(y)
  a_5 <- generate_regimes(ncol(a), dose)
  U <- matrix(0, nrow = ncol(s), ncol = 1)
  for (decision in 1:T.dp) {
    #decision = 1
    #cat("Decision: ", decision, "\n")
    for (regime in 1:nrow(a_5)) {
      #cat("Regime: ", regime, "\n")
      #regime = 1
      I_at <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == a_5[regime,1:decision]))
      I_0t <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == rep(0, decision)))
      p_a_noteq_0 <- 1 - mean(I_0t)
      w <- I_at/p_a[,decision] - I_0t/p_a[,decision]
      U <- U + t(s) %*% (w*y[,decision] - ( as.matrix(s) %*% beta )/T.dp)
    }
  }
  U <- U / nrow(y)
  return(U)
}

# We consider H
ee4 <- function( beta, y, a, h, s, p_a, cum_d, dose ) {
  #dose = 1
  #p = 0.5
  T.dp <- ncol(y)
  a_5 <- generate_regimes(ncol(a), dose)
  #U <- matrix(0, nrow = ncol(s), ncol = 1)
  U <- 0
  for (decision in 1:T.dp) {
    #decision = 1
    #cat("Decision: ", decision, "\n")

    df <- as.data.frame(cbind(y=y[,decision], h, a=a[,decision], cum_d=cum_d[,decision]))
    #df$cum_d <- factor(df$cum_d, levels = c(5,4,3,2,1,0))


    fit <- lm(y ~ -1 + h1 + h2 + h3 + cum_d, data = df)
    summary(fit)

    df.m0 <- df
    df.m1 <- df

    temp <- c(0,0,0,0)

    for (regime in 1:nrow(a_5)) {

      df.m1$cum_d <- sum(a_5[regime,1:decision])
      df.m0$cum_d <- 0

      m1 <- predict(fit, newdata = df.m1)
      m0 <- predict(fit, newdata = df.m0)

      I_at <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == a_5[regime,1:decision]))
      I_0t <- apply(a[, 1:decision, drop=FALSE], 1, function(x) all(x == rep(0, decision)))
      #print(paste("Regime", regime))
      #print(rbind(I_at/p_a[,decision] * y[,decision],
      #           I_0t/p_a[,decision] * y[,decision],
      #            I_at/p_a[,decision] * m1,
      #            m1,
      #            I_0t/p_a[,decision] * m0,₩
      #            m0))

      U <- U + sum(
        I_at/p_a[,decision] * (y[,decision] - m1) -
        I_0t/p_a[,decision] * (y[,decision] - m0) +
        m1 - m0 - beta / T.dp )

    }
  }
  U <- U / nrow(y)
  return(U)
}



# We consider H and S

boruvka <- function(beta, y, a, h, s, dose) {
  beta1 <- beta[2]
  alpha <- beta[1]
  rho.hat = mean(a==dose)
  a <- a==dose
  p_a <- dbinom(dose, 5, 0.5)
  p_0 <- dbinom(0, 5, 0.5)
  p <- p_a / (p_a+p_0)
  w = w = (rho.hat^a * (1-rho.hat)^(1-a)) / (p^a * (1-p)^(1-a))
  U = 0
  resid <- y - (a - rho.hat) * beta
  resid.neg <- ifelse(resid <0, 1,0)
  U <- U + sum((y - alpha - (a - rho.hat) * beta1) * w * cbind(rep(1,m), (a - rho.hat)))

  #print(U)
  U <- as.vector(U / length(y))
  return(U)
}

boruvka.covar <- function(beta, y, a, h, s, p_a, dose) {
  rho.hat = mean(a==dose)
  a <- a==dose
  p_a <- dbinom(dose, 5, 0.5)
  w = w = (rho.hat^a * (1-rho.hat)^(1-a)) / (p_a^a * (1-p_a)^(1-a))
  #U = c(0,0,0)
  U = 0
  for (i in 1:length(y)) {
    #resid <- y[i] - (a[i] - rho.hat) * as.matrix(s[1,]) %*% beta
    resid <- y[i] - (a[i] - rho.hat) * beta
    weight <- w[i]
    s_beta <- (a[i] - rho.hat) # * s[i,]
    U <- U + as.matrix(resid * weight * s_beta)
    #print(as.matrix(resid * weight * s_beta))
    #U <- U + as.matrix(( y[i] - (a[i] - rho.hat) * as.matrix(s[1,]) %*% beta ) * w[i] * s[i,] * (a[i] - rho.hat))
  }
  #print(U)
  U <- as.vector(U / length(y))
  return(U)
}








