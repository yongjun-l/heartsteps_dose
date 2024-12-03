

#' Title
#'
#' @param n number of observations (in HeartSteps, number of total randomization points)
#' @param rho correlation between adjacent observations
#' @return a correlation matrix with 0.5 correlation between adjacent observations
create_corr_matrix <- function(n, rho=0.5) {
  corr_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      corr_matrix[i, j] <- rho^(abs(i - j) / 2)
    }
  }
  return(corr_matrix)
}

#' Mobile Health Data Simulation
#'
#' \deqn{Y_{t+1} = \eta H + \theta_1 S_t + \theta_2 A_{t-1} + (\beta_{10} + \beta_{11}H' + \beta_{12}S_t')A_t}
#'
#' @param m number of replications
#' @param n number of observations (in HeartSteps, number of participants)
#' @param time number of time points (in HeartSteps, 5)
#' @param days number of days (in HeartSteps, 42)
#' @param eta effect of baseline covariate
#' @param rho correlation between adjacent observations
#' @param theta1 effect of time-varying covariate
#' @param theta2 effect of lagged treatment (lag 1)
#' @param beta10 effect of treatment at time t
#' @param beta11 interaction effect with some baseline covariate
#' @param beta12 interaction effect with some time-varying covariate
#' @param p probability of treatment
#'
#' @return A list of data frames that contains
#'   - ID: participant ID
#'   - OBS: observation number
#'   - SLOT: time slot
#'   - H: baseline covariate
#'   - S: time-varying covariate
#'   - A: treatment indicator
#'   - A.prev: treatment indicator at lag 1
#'   - eps: error term
#'   - Y: outcome
#' @export
simMhealth <- function(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p) {
  dfs <- list()
  n.h <- length(eta)
  n.s <- length(theta1)

  for (rep in 1:m) {
    if (rep %% 10 == 0) {
      cat(rep, "\n")
    }

    ID <- rep(1:n, each=time*days)
    OBS <- rep(1:(time*days), n)
    SLOT <- rep(1:time, n*days)
    H <- matrix(rep(rnorm(n*n.h, 0, 1), each=time*days), ncol=n.h) # baseline covariate. same for each person
    S <- matrix(rep((runif(n*days*n.s)>0.5), each=time), ncol=n.s) # time-varying covariate. differs for each day
    eps <- as.vector(t(MASS::mvrnorm(n=n, mu=rep(0, time*days), Sigma=create_corr_matrix(time*days, rho))))
    A <- matrix(runif(n*days*time)>p)
    A.prev <- c(0, (rep(c(rep(1, (time*days)-1),0),n)*A)[-time*days*n] )

    Y <- H %*% eta +
      S %*% theta1 +
      A.prev * theta2 +
      A * (beta10 + H[,1] %*% beta11 + S[,1] %*% beta12) +
      eps

    dfs[[rep]] <- data.frame(ID, OBS, SLOT, H=H, S=S, A, A.prev, eps, Y)

  }
  return(dfs)
}






