# Estimating Equation ----

#' Title generate regimes
#'
#' @param time_points number of total randomization points within a day
#' @param total_doses number of doses given in a day
#'
#' @return a matrix of all possible dose combinations
#' @examples 
#' # This is an internal function
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
#' @examples 
#' # This is an internal function
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
#' @examples
#' # This is an internal function
ee <- function( beta, df.wide, p_a, time, n.days, dose, y, trt, 
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

#' Get coefficients
#'
#' @param df dataset
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param p (numeric) probability of single treatment
#' @param dose (integer) dose
#' @param ... model specification
#'
#' @return a vector of coefficients
#' @examples
#' # This is an internal function
CeeDose <- function(df, id, day, slot, p, dose, ...) {
  args <- list(...)
  y <- args$y;trt <- args$trt
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
  ee( init_beta, df.wide,  p_a, time, n.days, dose, ...)
  rslt <- nleqslv::nleqslv(init_beta, function(beta) ee( beta, df.wide, p_a, time, n.days, dose, ...))
  return(rslt$x)
}

# Fit ----

#' mHealthDose fit
#'
#' @param df dataset
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param p (numeric) probability of single treatment
#' @param dose (integer) dose
#' @param boot (integer) number of bootstrap samples
#' @param cores (integer) number of cores to use for bootstrap
#' @param print_progress (logical) print progress
#' @param ... model specification
#' variable names of \code{dose}, \code{y}, \code{trt}, \code{baseline}, \code{timevar}, \code{b.prime}, \code{t.prime}
#' @param sim (logical) true if simulation setting, false if real data
#' @param cl (cluster) cluster object
#'
#' @return list of coefficients and standard errors
#' @importFrom parallel makeCluster stopCluster parSapply mclapply
#' @import data.table
#' @export
#' @examples
#' m=5; n=10; time=5; days=2; eta=c(3,2,1,4); rho=0.5; theta1=c(3,2,1); theta2=0; 
#' beta10=1; beta11=0.5; beta12=0.3; p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' mHealthDose(dfs, id="ID", day="DAY", slot="SLOT", p,
#'            dose=dose, y="Y", trt="A",
#'            baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'            b.prime=c("H1"), t.prime=c("S1"),
#'            boot = m, cores = 1, sim = TRUE, print_progress=FALSE)
mHealthDose <- function(df, id, day, slot, p, dose, ...,
                        boot = 100, cores = 1, cl = NULL, sim = FALSE,
                        print_progress=FALSE) {
  
  if (sim) { boot_dfs <- df } else { boot_dfs <- boot_samples(df, id, B=boot) }
  
  # Calling parallel
  if ((cores > 1) | length(cl)) {
    
    # Creating the cluster
    if (!length(cl)) {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl))
      
      # Loading R packages
      parallel::clusterEvalQ(cl, library(mrt.dose))
    }
    
    # Calling the function
    rslts <- parallel::parSapply(cl, X=seq_len(boot), function(i, CeeDose, df, id, day, slot, p, dose, sim, boot_dfs, ...) {
      if (sim) {
        boot_df <- df[[i]]
        boot_id <- id
      } else {
        boot_df <- cbind(df[boot_dfs[[i]]$index,], boot_id = boot_dfs[[i]]$no.repeat.id)
        boot_id <- "boot_id"
      }
      # Call the CeeDose function with all required arguments
      return(CeeDose(boot_df, boot_id, day, slot, p, dose, ...))
    }, CeeDose = mrt.dose:::CeeDose, df = df, id = id, day = day, slot = slot, p = p, dose = dose, sim = sim, boot_dfs = boot_dfs, ...)
    
    
  } else {
    
    # If no parallel apply
    rslts <- sapply(X = seq_len(boot), function(i, CeeDose, df, id, day, slot, p, dose, sim, boot_dfs, ...) {
      if (sim) {
        boot_df <- df[[i]]
        boot_id <- id
      } else {
        boot_df <- cbind(df[boot_dfs[[i]]$index,], boot_id = boot_dfs[[i]]$no.repeat.id)
        boot_id <- "boot_id"
      }
      # Call the CeeDose function with all required arguments
      return(CeeDose(boot_df, boot_id, day, slot, p, dose, ...))
    }, CeeDose = mrt.dose:::CeeDose, df = df, id = id, day = day, slot = slot, p = p, dose = dose, sim = sim, boot_dfs = boot_dfs, ...)
    
  }
  args <- list(...)
  trt <- args$trt; b.prime <- args$b.prime; t.prime <- args$t.prime
  rownames(rslts) <- c(paste0(trt, dose, sep=""), paste(paste0(trt, dose, sep=""), c(b.prime, t.prime), sep=":")) 
  
  std.error <- apply(rslts, 1, sd)
  if (sim) 
    point <- rowMeans(rslts)
  else 
    point <- CeeDose(df, id, day, slot, p, dose, ...)
  names(point) <- c(paste0(trt, dose, sep=""), paste(paste0(trt, dose, sep=""), c(b.prime, t.prime), sep=":"))
  
  return(structure(list(
    coefficients = point,
    std.error = std.error, 
    ...), class = "mHealthDose")
  )
}

#' wrapper for multiple doses
#'
#' @param df dataset
#' @param id (char) id variable
#' @param day (char) day variable
#' @param slot (char) slot variable
#' @param p (numeric) probability of single treatment
#' @param doses (integer vector) doses
#' @param ... model specification
#'
#' @return list of mHealthDose objects
#' @export
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fits <- mHealthDoses(dfs, id="ID", day="DAY", slot="SLOT", p,
#'                      doses=c(1,2,3,4,5), y="Y", trt="A",
#'                      baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                      b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                      boot = m, cores = 1, sim = TRUE)
mHealthDoses <- function(df, id, day, slot, p, doses, ...) {
  fits <- list()
  for (dose in doses) {
    fits[[dose]] <- mHealthDose(df, id, day, slot, p, dose=dose, ...)
  }
  return(structure(fits, class = "mHealthDoses"))
}

# Reporting Results ----

#' Print mHealthDose object
#'
#' @param x mHealthDose object
#' @param ... optional
#'
#' @return invisible
#' @export
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fit <- mHealthDose(dfs[[1]], id="ID", day="DAY", slot="SLOT", p,
#'                    dose=dose, y="Y", trt="A",
#'                    baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                    b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                    boot = m, cores = 1, print_progress=FALSE)
#' print(fit)
print.mHealthDose <- function(x, ...) {
  cat("Coefficients:\n")
  print(x$coefficients)
  cat("\n")
  cat("Standard Errors:\n")
  print(x$std.error)
  invisible(x)
}

#' Get coefficients
#'
#' @param param parameters of interest
#' @param object mHealthDose object
#'
#' @return coefficients
#' @export
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fit <- mHealthDose(dfs[[1]], id="ID", day="DAY", slot="SLOT", p,
#'                    dose=dose, y="Y", trt="A",
#'                    baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                    b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                    boot = m, cores = 1, print_progress=FALSE)
#' coef(fit)
coef.mHealthDose <- function(object, param=NULL) {
  if (is.null(param)) {
    object$coefficients
  } else {
    object$coefficients[param]
  }
}

#' Print confidence interval
#'
#' @param object mHealthDose object
#' @param param parameters of interest
#' @param level confidence level
#'
#' @return custom print function for class mHealthDose
#' @export
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fit <- mHealthDose(dfs[[1]], id="ID", day="DAY", slot="SLOT", p,
#'                    dose=dose, y="Y", trt="A",
#'                    baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                    b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                    boot = m, cores = 1, print_progress=FALSE)
#' confint(fit)
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
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fit <- mHealthDose(dfs[[1]], id="ID", day="DAY", slot="SLOT", p,
#'                    dose=dose, y="Y", trt="A",
#'                    baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                    b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                    boot = m, cores = 1, print_progress=FALSE)
#' pval(fit)
pval <- function(object, param) {
  cf <- coef.mHealthDose(object)
  ses <- object$std.error
  pnames <- names(ses)
  
  if (is.matrix(cf)) 
    cf <- setNames(as.vector(cf), pnames)
  if (missing(param)) 
    param <- pnames
  else if (is.numeric(param)) 
    param <- pnames[param]
  pval <- 2 * (1 - pnorm(abs(cf[param]/ses[param])))
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
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fit <- mHealthDose(dfs[[1]], id="ID", day="DAY", slot="SLOT", p,
#'                    dose=dose, y="Y", trt="A",
#'                    baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                    b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                    boot = m, cores = 1, print_progress=FALSE)
#' summary(fit)
summary.mHealthDose <- function(object, ...) {
  dose = object$dose
  working = as.formula(paste0(object$y, " ~ ", paste0(c(object$trt,object$baseline, object$timevar), collapse = " + ")))
  treatment = as.formula(
    paste0(
      object$y, " ~ ",
      paste(c(object$trt, paste(object$trt, c(object$b.prime, object$t.prime), sep=":")), collapse = " + ")
    )
  )
  coefficients = cbind("Estimate"=coef.mHealthDose(object), 
                       "Std. Error"=object$std.error, 
                       "z value"=abs(coef.mHealthDose(object)/object$std.error), 
                       confint.mHealthDose(object), 
                       "p-val"=pval(object))
  return(
    structure(
      list(
        dose = dose,
        working = working,
        treatment = treatment,
        coefficients = coefficients
      ), class = "summary.mHealthDose"))
}

#' Print Summary of mHealthDose
#'
#' @param x summary.mHealthDose object
#' @param ... optional
#'
#' @return invisible 
#' @export
#' @examples
#' m=5; n=10; time=5; days=5; eta=c(0,0,0,0); rho=0; theta1=c(0,0,0); theta2=0; 
#' beta10=5; beta11=c(0,0); beta12=c(0,0); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fit <- mHealthDose(dfs[[1]], id="ID", day="DAY", slot="SLOT", p,
#'                    dose=dose, y="Y", trt="A",
#'                    baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                    b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                    boot = m, cores = 1, print_progress=FALSE)
#' summary(fit)
print.summary.mHealthDose <- function(x, ...) {
  cat("Dose: ", x$dose, "\n\n")
  cat("Working Model:\n", noquote(deparse(x$working)), "\n\n")
  cat("Treatment Model:\n", noquote(deparse(x$treatment)),"\n\n")
  cat("Coefficients:\n")
  print(x$coefficients)
  invisible(x)
}

#' mHealthDoses summary
#'
#' @param object mHealthDoses object
#' @param ... not used
#'
#' @return plot dose response relationship
#' @export
#' @examples
#' m=5; n=10; time=5; days=10; eta=c(3,2,1,4); rho=0.5; theta1=c(3,2,1); theta2=0; 
#' beta10=5; beta11=c(3,2); beta12=c(2,1); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fits <- mHealthDoses(dfs, id="ID", day="DAY", slot="SLOT", p,
#'                      doses=c(1,2,3,4,5), y="Y", trt="A",
#'                      baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                      b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                      boot = m, cores = 1, sim = TRUE)
#' summary(fits)
summary.mHealthDoses <- function(object, ...) {
  allSummaries <- NULL
  for (dose in 1:length(object)) {
    allSummaries <- rbind(allSummaries, summary.mHealthDose(object[[dose]])$coefficients)
  }
  return(structure(allSummaries, class = c("mHealthDoses.summary", "matrix")))
}

#' Plot treatment effect
#'
#' @param x mHealthDoses object
#' @param y varlevels
#' @param ... not used
#'
#' @return ggplot of treatment effects
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr mutate filter slice n bind_rows transmute
#' @import ggplot2
#' @export
#'
#' @examples 
#' m=5; n=10; time=5; days=10; eta=c(3,2,1,4); rho=0.5; theta1=c(3,2,1); theta2=0; 
#' beta10=5; beta11=c(3,2); beta12=c(2,1); p=0.5; dose=1
#' dfs <- simMhealth(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p)
#' fits <- mHealthDoses(dfs, id="ID", day="DAY", slot="SLOT", p,
#'                      doses=c(1,2,3,4,5), y="Y", trt="A",
#'                      baseline=c("H1", "H2"), timevar=c("S1", "S2"),
#'                      b.prime=c("H1", "H2"), t.prime=c("S1", "S2"),
#'                      boot = m, cores = 1, sim = TRUE)
#' sum <- summary(fits)
#' varlevels <- rbind(H1 = c("Male", "Female"), 
#'                    S1 = c("Student", "Non-student"), 
#'                    H2 = c("Mean Age", "Mean Age+1"), 
#'                    S2 = c("Physical Low", "Physical High"))
#' colnames(varlevels) <- c("Main Effect", "Interaction Effect")
#' p <- plot(sum, varlevels)
#' p[[1]]
plot.mHealthDoses.summary <- function(x, y, ...) {
  x.df <- as.data.frame(x, check.names = FALSE)
  terms <- rownames(x.df)
  x.df2 <- x.df |>
    transmute(
      Treatment = sub(":.*", "", terms),  # Extract treatment
      Interaction = ifelse(grepl(":", terms), sub(".*:", "", terms), "None"), # Extract interaction or "None"
      EffectType = ifelse(Interaction == "None", "Main Effect", "Interaction Effect"),
      Estimate = Estimate,
      se = `Std. Error`, 
      lower = `2.5%`,
      upper = `97.5%`
    ) 
  
  # Duplicate main effects for each interaction type
  x.m <- x.df2 |>
    filter(Interaction == "None") |>
    slice(rep(1:n(), each = 4))  # Repeat rows for each interaction
  x.m$Interaction<-rep(c("H1", "H2", "S1", "S2"), times = nrow(x.m) / 4)
  
  # Combine main effects with interactions
  x.int <- x.df2 |>
    filter(Interaction != "None")
  x.int[,4:7] <- x.m[,4:7] + x.int[,4:7]
  
  x.f <- bind_rows(x.m, x.int)
  x.f$labels <- diag(y[x.f$Interaction, x.f$EffectType])
  
  plots <- list()
  for( var in 1:nrow(y)) {
    plot_data_var <- x.f |> filter(Interaction == rownames(y)[var])
    plots[[rownames(y)[var]]] <- ggplot(plot_data_var, aes(x=Treatment, y=Estimate, fill=labels)) + 
      geom_bar(
        aes(y=Estimate, group=labels),
        stat = "identity",
        position = position_dodge(width = 0.8),  # Offset the bars
        width = 0.3
      ) +
      geom_errorbar(
        aes(ymin = lower, ymax = upper, group = labels),
        position = position_dodge(width = 0.8),  # Match the dodge width
        width = 0.2
      ) +  # Narrow error bars
      geom_hline(yintercept = 0) +
      labs(
        title = "",
        x = "Dose",
        y = "Estimate",
        fill = varlevels[var]
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.title = element_text(size = 10),       # Adjust title size
        legend.text = element_text(size = 8),         # Adjust text size
        legend.key.size = unit(0.5, "cm"),            # Adjust key size
        legend.spacing.y = unit(0.1, "cm")  
      ) 
  }
  plots[["all"]] <- grid.arrange(grobs = plots, ncol = 2)
  plots[["all"]]
  invisible(plots)
}