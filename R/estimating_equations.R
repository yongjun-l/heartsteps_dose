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


ee.cor.3 <- function( beta, df.wide, dose, y, trt, p_a, time, n.days,
                      baseline=NULL, timevar=NULL, b.prime=NULL, t.prime=NULL) {

  p <- 0.5; U <- 0; zeros <- lapply(1:5, rep, x = 0)
  a_5 <- generate_regimes(time, dose)

  if (!is.null(b.prime)) {
    ate <- cbind(1, df.wide[,b.prime]) |> as.matrix()
  } else {
    ate <- rep(1, n.days) |> as.matrix()
  }

  for (decision in 1:time) {
    y.name <- paste0(y,"_", decision)
    a.name <- paste0(trt,"_", decision)
    if (!is.null(timevar)) {s.name <- paste0(timevar,"_", decision)} else {s.name <- NULL}
    formula <- as.formula(paste(y.name, paste(c(baseline, s.name, a.name), collapse = " + "), sep = " ~ "))

    fit <- lm(formula, data = df.wide)

    if (!is.null(t.prime)) {
      ate.name <- c(b.prime,paste0(t.prime,"_", decision))
      ate <- cbind(1, df.wide[,ate.name]) |> as.matrix()
    }

    df.m0 <- df.wide[,c(y.name, baseline, s.name, a.name)]
    df.m1 <- df.wide[,c(y.name, baseline, s.name, a.name)]

    df.m0[,a.name] <- FALSE

    a.interest <- df.wide[, paste(trt, "_", 1:decision, sep="")] |> as.matrix()
    y.decision <- df.wide[,y.name] |> as.matrix()
    p.decision <- p_a[,decision]
    zeros.decision <- zeros[[decision]]
    for (regime in 1:nrow(a_5)) {

      df.m1[,a.name] <- as.logical(a_5[regime, decision])
      m1 <- predict(fit, newdata = df.m1)
      m0 <- predict(fit, newdata = df.m0)

      I_at <- rowSums(a.interest == matrix(a_5[regime, 1:decision], nrow=n.days, ncol=decision, byrow=TRUE)) == decision
      I_0t <- rowSums(a.interest == matrix(zeros.decision, nrow=n.days, ncol=decision, byrow=TRUE)) == decision

      U <- U + t(ate) %*% (I_at/p.decision * (y.decision - m1) -
                                    I_0t/p.decision * (y.decision - m0) +
                                    m1 - m0 - ( as.matrix(ate) %*% beta )/time)
    }
  }
  U <- ( U / n.days )
  return(U)
}

CeeDose <- function(df, id, day, slot, y, trt, dose, p,
                    baseline=NULL, timevar=NULL, b.prime=NULL, t.prime=NULL) {
  id <- rlang::sym(id)
  day <- rlang::sym(day)
  time <- length(unique(df[,slot]))
  df.wide <- df |>
    dplyr::select(!!id,!!day,slot,y,trt,baseline,timevar) |>
    dplyr::group_by(!!id,!!day) |>
    tidyr::pivot_wider(names_from={{slot}}, values_from=c(y,trt,timevar,t.prime))

  n.days <- nrow(df.wide);n.beta <- length(b.prime)+length(t.prime)+1
  a <- df.wide[, paste(trt, "_", 1:time, sep="")]
  matrix <- get_p_a(a, p)
  cum_d <- matrix$cum_d
  p_a <- matrix$p_a
  init_beta <- rep(0,n.beta)
  ee.cor.3( init_beta, df.wide, dose, y, trt, p_a, time, n.days,
            baseline, timevar, b.prime, t.prime)
  rslt <- nleqslv::nleqslv(init_beta, function(beta) ee.cor.3( beta, df.wide, dose, y, trt, p_a, time, n.days,
                                                               baseline, timevar, b.prime, t.prime))
  return(rslt$x)
}

#' Get Simulation Results
#'
#' @param m number of simulations
#' @param dfs list of dataframes
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param y (char) outcome variable
#' @param trt (char) treatment variable
#' @param dose (integer) dose variable
#' @param p (numeric) probability of single treatment
#' @param baseline (char vector) baseline variable
#' @param timevar (char vector) time-varying covariate
#' @param b.prime (char vector) baseline interaction
#' @param t.prime (char vector) time-varying interaction
#' @param print_progress (logical) print progress
#'
#' @return matrix of mean and standard deviation of estimated coefficients
#' @export
get.sim.rslts2 <- function(m, dfs, id, day, slot, y, trt, dose, p,
                           baseline=NULL, timevar=NULL, b.prime=NULL, t.prime=NULL,
                           print_progress=FALSE) {

  n.beta <- length(b.prime)+length(t.prime)+1
  ee.corr <- matrix(nrow = m, ncol = n.beta)
  for (rep in 1:m) {
    if ((rep %% 10 == 0) & (print_progress)) {
      cat(rep, "\n")
    }
    rslt <- CeeDose(dfs[[rep]], id, day, slot, y, trt, dose, p,
                    baseline, timevar, b.prime, t.prime)
    ee.corr[rep,] <- rslt
  }
  return(rbind(colMeans(ee.corr), apply(ee.corr, 2, sd)))
}



#' Get Empirical Results
#'
#' @param m number of simulations
#' @param dfs list of dataframes
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param y (char) outcome variable
#' @param trt (char) treatment variable
#' @param dose (integer) dose variable
#' @param p (numeric) probability of single treatment
#' @param baseline (char vector) baseline variable
#' @param timevar (char vector) time-varying covariate
#' @param b.prime (char vector) baseline interaction
#' @param t.prime (char vector) time-varying interaction
#' @param print_progress (logical) print progress
#'
#' @return matrix of mean and standard deviation of estimated coefficients
#' @export
mHealthDose <- function(m, df, id, day, slot, y, trt, dose, p,
                           baseline=NULL, timevar=NULL, b.prime=NULL, t.prime=NULL,
                           print_progress=FALSE) {
  
  n.beta <- length(b.prime)+length(t.prime)+1
  ee.corr <- matrix(nrow = m, ncol = n.beta)
  for (rep in 1:m) {
    if ((rep %% 10 == 0) & (print_progress)) {
      cat(rep, "\n")
    }
    rslt <- CeeDose(dfs[[rep]], id, day, slot, y, trt, dose, p,
                    baseline, timevar, b.prime, t.prime)
    ee.corr[rep,] <- rslt
  }
  return(rbind(colMeans(ee.corr), apply(ee.corr, 2, sd)))
}














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







