#' Create correlation matrix to simulate epsilon
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
#' \deqn{Y_{t+1} = \eta H + \theta_1 S_t + \theta_2 A_{t-1} + (\beta_{10} + \beta_{11}H' + \beta_{12}S_t')A_t + \epsilon}
#' \deqn{\begin{aligned}
#' H &= (\bar{h}_1, \bar{h}_2, \bar{h}_3),\\
#' h_1 &\sim Binomial(0.5), \\
#' h_2&\sim Normal(0,1), \\
#' h_3&= h_2^2, \\
#' h_{int} &= \bar{1}_n,\\\\
#' S &= (\bar{s}_1, \bar{s}_2, \bar{s}_3),\\
#' s_1 &\sim Binomial(0.5), \\
#' s_2&\sim Normal(0,1), \\
#' s_3&= s_2^2, 
#' \end{aligned}}
#' add text version of the document
#' h1=rep(rbinom(n, 1, 0.5),each=days*time)
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
#' @param print_progress (boolean) print progress
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
simMhealth <- function(m, n, time, days, eta, rho, theta1, theta2, beta10, beta11, beta12, p=0.5, print_progress=FALSE) {
  dfs <- list()
  n.h <- length(eta)
  n.s <- length(theta1)
  n.h.prime <- length(beta11)
  n.s.prime <- length(beta12)

  if ((n.h < n.h.prime)|(n.s < n.s.prime)) {
    stop("Interaction variables must be a subset of baseline covariates")
  }

  for (rep in 1:m) {
    if ((rep %% 10 == 0) & (print_progress)) {
      cat(rep, "\n")
    }

    ID <- rep(1:n, each=time*days)
    OBS <- rep(1:(time*days), n)
    DAY <- rep(rep(1:days, each=time), n)
    SLOT <- rep(1:time, n*days)


    h1=rep(rbinom(n, 1, 0.5),each=days*time)
    h2=rep(rnorm(n), each=days*time)
    h3=rep(h2^2)
    int=rep(1, n*days*time)

    s1=rep(rbinom(n*days, 1, 0.5),each=time)
    s2=rep(rnorm(n*days), each=time)
    s3=rep(s2^2)

    H <- cbind(H1=h1,H2=h2,H3=h3,HINT=int)
    S <- cbind(S1=s1,S2=s2,S3=s3)
    #H <- matrix(rep(rnorm(n*n.h, 0, 1), each=time*days), ncol=n.h) # baseline covariate. same for each person
    #S <- matrix(rep((runif(n*days*n.s)>0.5), each=time), ncol=n.s) # time-varying covariate. differs for each day

    eps <- as.vector(t(MASS::mvrnorm(n=n, mu=rep(0, time*days), Sigma=create_corr_matrix(time*days, rho))))
    A <- matrix(runif(n*days*time)>p)
    A.prev <- c(0, (rep(c(rep(1, (time*days)-1),0),n)*A)[-time*days*n] )

    Y <- H %*% eta +
      S %*% theta1 +
      A.prev * theta2 +
      A * (beta10 + H[,c(1:n.h.prime),drop=FALSE] %*% beta11 + S[,c(1:n.s.prime),drop=FALSE] %*% beta12) +
      eps

    dfs[[rep]] <- data.frame(ID, OBS, DAY, SLOT, H, S, A, A.prev, eps, Y)

  }
  return(dfs)
}






