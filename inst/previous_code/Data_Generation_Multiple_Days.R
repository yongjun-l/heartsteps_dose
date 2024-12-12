# parameters ----
p = 0.5
n = 100
days = 20
m = 200
t = 5
dose = 1



a <- matrix(nrow=n*days, ncol = t)
cum_d <- matrix(0, nrow = n*days, ncol = t)
p_a <- matrix(nrow=n*days, ncol = t)


# data generation ----
## Treatment randomization
dfs <- list()
for (rep in 1:m)  {
  #rep = 1
  for (decision in 1:t) {
    #decision = 1
    a[,decision] <- rbinom(n*days, 1, p)
    cum_d[,decision] <- if(decision==1) {a[,1:decision]} else {rowSums(a[,1:decision])}
    p_a[, decision] <- (p^cum_d[,decision]) * ((1-p)^(decision - cum_d[,decision]))
  }
  # First, Second setting ----
  ## Baseline Covariates
  'h0 <- rep(1, n) # intercept
  h1 <- rbinom(n, 1, 0.2) # binary covariate
  h2 <- rnorm(n) # continuous covariate # or ignore this variable
  h <- cbind(h0, h1, h2)

  ## Effect modifying Covariates
  s0 <- rep(1, n) # intercept
  s1 <- rbinom(n, 1, 0.5) # binary covariate
  s2 <- rnorm(n) # continuous covariate
  s <- cbind(s0, s1, s2)

  ## Outcome Generation
  beta <- c(100, 20, 0)
  eta <- c(50, 0, 0)
  beta_t <- outer(beta / 5, rep(1, t))
  eta_t <- outer(eta / 5, rep(1, t))

  eps <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  Y <- h %*% eta_t + s %*% beta_t * a + eps
  Y <- h %*% eta_t + beta_t * a + eps'


  # Fourth Setting: Incorporate H ----
  ## Baseline Covariates
  'setting="fourth"
  h0 <- rep(1, n) # intercept
  h1 <- rbinom(n, 1, 0.5) # binary covariate
  h2 <- rnorm(n) # continuous covariate # or ignore this variable
  h <- cbind(h0, h1, h2)

  ## Effect modifying Covariates
  s0 <- rep(1, n) # intercept
  s1 <- rbinom(n, 1, 0.5) # binary covariate
  s2 <- rnorm(n) # continuous covariate
  s <- cbind(s0, s1, s2)

  ## Outcome Generation
  beta <- c(20, 0, 0)
  eta <- c(30, 20, 10)
  beta_t <- outer(beta , rep(1, t))
  eta_t <- outer(eta , rep(1, t))

  eps <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  Y <- h %*% eta_t + (s %*% beta_t) * a + eps'

  # Fifth Setting: Misspecification----
  ## Baseline Covariates
  setting="seventh"
  h0 <- rep(1, n) # intercept
  #h1 <- rnorm(n, 1, 0.2) # binary covariate
  h1 <- rnorm(n) # continuous covariate
  h2 <- h1 ^ 2
  h <- cbind(h0, h1, h2)

  ## Effect modifying Covariates
  s0 <- rep(1, n) # intercept
  s1 <- rbinom(n, 1, 0.5) # binary covariate
  s2 <- rnorm(n) # continuous covariate
  s <- cbind(s0, s1, s2)

  #repeat each row 5 times
  h <- h[rep(1:n, each = days),]
  s <- s[rep(1:n, each = days),]

  ## Outcome Generation
  beta <- c(20, 10, 5)
  eta <- c(30, 20, 10)
  beta_t <- outer(beta, rep(1, t))
  eta_t <- outer(eta, rep(1, t))

  eps <- matrix(rnorm(n * days * 5), nrow = n * days, ncol = 5)
  Y <- h %*% eta_t + (s %*% beta_t) * a + eps


  # Put all data together into a data frame. ----
  colnames(Y) <- paste0("Y", 1:t)
  colnames(h) <- paste0("h", 1:ncol(h))
  colnames(s) <- paste0("s", 1:ncol(s))
  colnames(a) <- paste0("a", 1:ncol(a))

  data <- list(
    df = data.frame(Y, h, s, a),
    params = list(
      beta = beta,
      eta = eta,
      p = p
    )
  )
  dfs[[rep]] <- data
}


# get current script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
filename <- paste(setting,"_setting","_n",n,"_m",m,"_days",days,"_t",t,"_p",p,"_beta",beta[1],"_dose",dose,".rds", sep = "")
saveRDS(dfs, paste(script_dir, "simulated_data",filename, sep = "/"))

























