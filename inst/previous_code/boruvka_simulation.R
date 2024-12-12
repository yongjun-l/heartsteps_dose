library(geepack)
library(MASS)
library(nleqslv)
library(gee)


# functions ----
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

# yt.generate <- function(At_1, At, St_1, St, xi, eta1, eta2, theta1, theta2, beta10, beta11) {
#   E.St <- expit(xi * At_1) - (1 - expit(xi * At_1))
#   #pt_1 <- expit(eta1 * At_2 + eta2 * St_1)
#   pt <- expit(eta1 * At_1 + eta2 * St)
#   return(theta1 * ( St - E.St ) +
#          #theta2 * ( At_1 - pt_1 ) +
#          ( At - pt ) * (beta10 + beta11 * St))
# }



# data generation ----
Time = 30
n = 30
N = 1000
theta1 = 0.8
theta2 = 0
eta1 = -0.8
eta2 = 0.8
xi = 0
beta10 = -0.2
beta11 = 0.2
A0 = 0

for (beta11 in c(0, 0.2, 0.5, 0.8)) {
  beta = c()
  beta2 = c()
  beta3 = c()
  for (rep in 1:N) {

    if (rep %% 100 == 0) {
      cat(rep, "\n")
    }
    A = matrix(0, nrow=n, ncol=Time)
    S = matrix(0, nrow=n, ncol=Time)
    for (i in 1:n) {
      S[i,] = rep(c(-1,1), length.out=Time)
    }

    p = matrix(0, nrow=n, ncol=Time)
    Pr.St = matrix(0, nrow=n, ncol=Time)
    Y = matrix(0, nrow=n, ncol=Time+1)
    for (i in 1:n) {
      #i=1
      for (t in 1:Time) {
        #t=1
        if (t == 1) {
          Pr.St[i,t] <- expit(0)
          #S[i,t] <- ifelse(Pr.St[i,t] > runif(1), 1, -1)
          S[i,t] <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
          p[i,t] <- expit(eta1 * A0 + eta2 * S[i,t])
        } else {
          Pr.St[i,t] <- expit(xi * A[i, t-1])
          #S[i,t] <- ifelse(Pr.St[i, t] > runif(1), 1, -1)
          S[i,t] <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
          p[i,t] <- expit(eta1 * A[i, t-1] + eta2 * S[i, t])
        }

        A[i, t] <- rbinom(1, 1, p[i, t])

        'if (t==1) {
      Y[i, t+1] <- yt.generate(A0, A[i, t], 0, S[i, t], xi, eta1, eta2, theta1, theta2, beta10, beta11)
    } else {
      Y[i, t+1] <- yt.generate(A[i, t-1], A[i, t], S[i, t-1], S[i, t], xi, eta1, eta2, theta1, theta2, beta10, beta11)
    }'

        if (t==1) {
          Y[i, t+1] <- theta1*S[i,t]  +  theta2*(0-1)  +  (A[i,t]-p[i,t]) * (beta10 + beta11 * S[i, t])
        } else {
          Y[i, t+1] <- theta1*S[i,t]  +  theta2*(A[i, t-1]-p[i, t-1])  +  (A[i,t]-p[i,t]) * (beta10 + beta11 * S[i, t])
        }


      }
      eps <- mvrnorm(n=1, mu=rep(0, Time), Sigma=create_corr_matrix(Time))
      Y[i,] <- Y[i,] + c(0,eps)
    }
    mean(rowMeans(S))
    mean(S)


    # estimating equation

    rho.hat = mean(rowMeans(A))
    w = (rho.hat^A * (1-rho.hat)^(1-A)) / (p^A * (1-p)^(1-A))
    #w = 1 / p
    param = c(1,1,1)
    out = nleqslv(param, ee, Y=Y, S=S, A=A, rho.hat=rho.hat, w=w)
    out2 = nleqslv(param, gee_ind, Y=Y, S=S, A=A, rho.hat=rho.hat, w=w)
    beta <- c(beta, out$x[3])
    beta2 <- c(beta2, out2$x[3])


    # put y in long format with id as row number
    Y.long <- data.frame(ID = rep(1:(Time), n), t = rep(1:n, each=(Time)), a = as.vector(A),
                         p.tilde = as.vector((rho.hat^A * (1-rho.hat)^(1-A))), p = as.vector(p),
                         w = as.vector(w),
                         s = as.vector(S), y = as.vector(Y[,2:31]) )
    # sort Y.long by ID
    Y.long <- Y.long[order(Y.long$ID),]

    fit.wcls <- geeglm(y ~ s + a, id=ID, data=Y.long, weights = w, corstr="independence")
    #fit.indep <- geeglm(y ~ s + a, id=ID, data=Y.long, corstr="independence")
    fit.ar1 <- geeglm(y ~ s + a, id=ID, data=Y.long, corstr="ar1")

    fit.wcls$coefficients
    #cat("wcls == estimating eqaution", round(fit.wcls$coefficients[3],5) == round(out$x[3],5),"\n")
    #cat("geeglm == estimating eqaution", round(fit.indep$coefficients[3],5) == round(out2$x[3],5),"\n")
    #fit.ar1$coefficients
    beta3 <- c(beta3, fit.ar1$coefficients[3])
    #s <- summary(fit)
    #s$coefficients[3,4]
  }



  cat("######################\n")
  cat("beta11", beta11, "\n")
  cat("WCLS\n")
  cat("mean",mean(beta),"\n")
  cat("sd",sd(beta),"\n")
  cat("RMSE", (mean((beta-beta10)^2))^0.5,"\n")
  cat("coverage",mean(ifelse(beta +2*sd(beta) > beta10 & beta -2*sd(beta) < beta10, 1, 0)),"\n")

  cat("\nGEE-IND\n")
  cat("mean",mean(beta2),"\n")
  cat("sd",sd(beta2),"\n")
  cat("RMSE", (mean((beta2-beta10)^2))^0.5,"\n")
  cat("coverage",mean(ifelse(beta2 +2*sd(beta2) > beta10 & beta2 -2*sd(beta2) < beta10, 1, 0)),"\n")

  cat("\nGEE-AR1\n")
  cat("mean",mean(beta3),"\n")
  cat("sd",sd(beta3),"\n")
  cat("RMSE", (mean((beta3-beta10)^2))^0.5,"\n")
  cat("coverage",mean(ifelse(beta3 +2*sd(beta3) > beta10 & beta3 -2*sd(beta3) < beta10, 1, 0)),"\n")
  cat("######################\n")
}
