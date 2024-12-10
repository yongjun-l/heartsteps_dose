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

#' Estimating Equation
#'
#' @param beta (numeric) vector of coefficients
#' @param df.wide (data.frame) wide format of the dataset
#' @param p_a (matrix) probability of treatment regime
#' @param time (integer) number of decision points within a day
#' @param n.days (integer) number of total decision days
#' @param y (char) outcome variable
#' @param trt (char) treatment variable
#' @param dose (integer) dose variable
#' @param baseline (char vector) baseline variable
#' @param timevar (char vector) time-varying covariate
#' @param b.prime (char vector) baseline interaction
#' @param t.prime (char vector) time-varying interaction
#' @importFrom stats as.formula lm predict rbinom rnorm runif sd
#' @importFrom utils combn
#' @return a vector of estimating equation
ee.cor.3 <- function( beta, df.wide, p_a, time, n.days, dose, y, trt, 
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

CeeDose <- function(df, id, day, slot, p, ...) {
  args <- list(...)
  y <- args$y;trt <- args$trt;dose <- args$dose
  baseline <- if ("baseline" %in% names(args)) args$baseline else NULL
  timevar <- if ("timevar" %in% names(args)) args$timevar else NULL
  b.prime <- if ("b.prime" %in% names(args)) args$b.prime else NULL
  t.prime <- if ("t.prime" %in% names(args)) args$t.prime else NULL
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
  ee.cor.3( init_beta, df.wide,  p_a, time, n.days, ...)
  rslt <- nleqslv::nleqslv(init_beta, function(beta) ee.cor.3( beta, df.wide, p_a, time, n.days, ...))
  names(rslt$x) <- c("intercept", b.prime, t.prime)
  return(rslt$x)
}

#' Get Simulation Results
#'
#' @param m number of simulations
#' @param dfs list of dataframes
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param p (numeric) probability of single treatment
#' @param print_progress (logical) print progress
#' @param ... (char vector) model specification
#' variable names of \code{dose}, \code{y}, \code{trt}, \code{baseline}, \code{timevar}, \code{b.prime}, \code{t.prime}
#' 
#' @return matrix of mean and standard deviation of estimated coefficients
#' @export
get.sim.rslts2 <- function(m, dfs, id, day, slot, p, ...,
                           print_progress=FALSE) {
  
  ee.corr <- NULL
  for (rep in 1:m) {
    if ((rep %% 10 == 0) & (print_progress)) {
      cat(rep, "\n")
    }
    rslt <- CeeDose(dfs[[rep]], id, day, slot, p, ...)
    ee.corr <- rbind(ee.corr, rslt)
  }
  return(rbind(colMeans(ee.corr), apply(ee.corr, 2, sd)))
}

#' mHealthDose fit
#'
#' @param df dataset
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param p (numeric) probability of single treatment
#' @param se (logical) get standard error
#' @param boot (integer) number of bootstrap samples
#' @param cores (integer) number of cores to use for bootstrap
#' @param print_progress (logical) print progress
#' @param ... model specification
#' variable names of \code{dose}, \code{y}, \code{trt}, \code{baseline}, \code{timevar}, \code{b.prime}, \code{t.prime}
#' @param ncpus (integer) number of cores to use for bootstrap
#' @param cl (cluster) cluster object
#' @return list of coefficients and standard errors
#' @importFrom parallel mclapply
#' @import data.table
#' @export
mHealthDose <- function(df, id, day, slot, p, ...,
                        se = TRUE, boot = 100, cores = 1, ncpus = 1, cl = NULL, 
                        print_progress=FALSE) {
  point <- CeeDose(df, id, day, slot, p, ...)
  
  if (se==TRUE) {
    #v=matrix(rexp(n),nrow=n)
    # boot_coef <- parallel::mclapply(1:boot, function(rep) {
    #   boot_df <- df[sample(1:nrow(df), nrow(df), replace = TRUE),]
    #   return(CeeDose(boot_df, id, day, slot, p, ...))
    # }, mc.cores = cores)
    # std.error <- apply(do.call(rbind,boot_coef), 2, sd)
    #dta <- as.data.table(df)
    boot_dfs <- boot_samples(df, id, B=boot)
    
    rslts <- NULL
    for (rep in 1:boot) {
      boot_df <- cbind(df[boot_dfs[[rep]]$index,], boot_id = boot_dfs[[rep]]$no.repeat.id)
      id="boot_id"
      rslt <- CeeDose(boot_df, id, day, slot, p, ...)
      rslts <- rbind(rslts, rslt)
    }
    std.error <- apply(rslts, 2, sd)
    
  } else {
    std.error <- NULL
  }
  
  return(
    structure(
      list(
        coefficients = point,
        std.error = std.error, 
        ...
      ),
      class = "mHealthDose"
    )
  )
}


#' Print mHealthDose object
#'
#' @param x mHealthDose object
#' @param ... optional
#'
#' @return invisible
#' @export
print.mHealthDose <- function(x, ...) {
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\n")
  cat("Standard Errors:\n")
  print(x$std.error)
  invisible(x)
}

#' Get coefficients
#' @param object mHealthDose object
#' @return coefficients
#' @export
coef.mHealthDose <- function(object) {
  object$coefficients
}

#' Print confidence interval
#'
#' @param object mHealthDose object
#' @param param parameters of interest
#' @param level confidence level
#'
#' @return custom print function for class mHealthDose
#' @export
confint.mHealthDose <- function(object, param, level = 0.95) {
  cf <- coef.mHealthDose(object)
  ses <- object$std.error
  pnames <- names(ses)
  if (is.matrix(cf)) 
    cf <- setNames(as.vector(cf), pnames)
  if (missing(param)) 
    param <- pnames
  else if (is.numeric(param)) 
    param <- pnames[param]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  pct <- paste0(round(a * 100, 3), "%")
  ci <- array(NA_real_, dim = c(length(param), 2L), dimnames = list(param, 
                                                                    pct))
  ci[] <- cf[param] + ses[param] %o% qnorm(a)
  return(ci)
}

#' Get p-values
#'
#' @param object mHealthDose object
#' @param param parameters of interest
#'
#' @return p-values
#' @importFrom stats pnorm setNames qnorm
#' @export
pval <- function(object, param) {
  cf <- confint.mHealthDose(object)
  ses <- object$std.error
  pnames <- names(ses)
  
  if (is.matrix(cf)) 
    cf <- setNames(as.vector(cf), pnames)
  if (missing(param)) 
    param <- pnames
  else if (is.numeric(param)) 
    param <- pnames[param]
  pval <- 2 * (1 - qnorm(abs(cf[param]/ses[param])))
  pval <- setNames(pval, param)
  return(pval)
}


#' mHealthDose summary
#'
#' @param object mHealthDose object
#' @param ... summary
#'
#' @return the summary object
#' @importFrom stats as.formula
#' @export 
summary.mHealthDose <- function(object, ...) {
  return(
    structure(
      list(
        working = as.formula(paste0(object$y, " ~ ", object$trt, " + ", 
                                    paste0(object$baseline, collapse = " + "), " + ", 
                                    paste0(object$timevar, collapse = " + "))),
        treatment = as.formula(
          paste0(
            object$y, " ~ ", object$trt, " + ",
            paste(
              paste(
                object$trt, c(object$b.prime, object$t.prime), sep=":"), 
              collapse = " + "
            )
          )
        ),
        coefficients = cbind("Estimate"=coef.mHealthDose(object), 
                             "Std. Error"=object$std.error, 
                             "z value"=abs(coef.mHealthDose(object)/object$std.error), 
                             confint.mHealthDose(object), 
                             "p-val"=pval(object))
      ), class = "summary.mHealthDose"))
}

#' Print Summary
#'
#' @param x summary.mHealthDose object
#' @param ... optional
#'
#' @return invisible 
#' @export
print.summary.mHealthDose <- function(x, ...) {
  cat("Working Model:\n")
  cat(noquote(deparse(x$working)))
  cat("\n\n")
  cat("Treatment Model:\n")
  cat(noquote(deparse(x$treatment)))
  cat("\n\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  invisible(x)
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







