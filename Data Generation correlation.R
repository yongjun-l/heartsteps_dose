# Data generation with correlation from boruvka
library(MASS)
create_corr_matrix <- function(n) {
  corr_matrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      corr_matrix[i, j] <- 0.5^(abs(i - j) / 2)
    }
  }
  return(corr_matrix)
}


# data generation ----
Time = 30
n = 100
m = 500
theta1 = 0.8
theta2 = 0
eta = 0.5
xi = 0
beta10 = -0.2
beta11 = 0.2
A0 = 0
p=0.5

dfs <- list()
for (beta11 in c(0, 0.2, 0.5, 0.8)) {
  #beta11=0.2
  for (rep in 1:m) { # Replication
    #rep=1
    if (rep %% 100 == 0) {
      cat(rep, "\n")
    }
    A = matrix(0, nrow=n, ncol=Time)
    S = matrix(0, nrow=n, ncol=Time)
    Y = matrix(0, nrow=n, ncol=Time)

    H = rnorm(n, 0, 1)

    for (i in 1:n) {
      #i=1
      for (t in 1:Time) {
        #t=1
        S[i, t] <- sample(c(-1, 1), size = 1, prob = c(0.5, 0.5))
        A[i, t] <- rbinom(1, 1, p)
        if (t == 1) {
          A[i, t] <- A0
          Y[i, t] <- eta * H[i] + theta1*S[i,t]  +  (A[i,t]) * (beta10 + beta11 * S[i, t])
        } else {
          Y[i, t] <- eta * H[i] + theta1*S[i,t]  +  theta2*(A[i, t-1])  +  (A[i,t]) * (beta10 + beta11 * S[i, t])
        }
      }
      eps <- mvrnorm(n=1, mu=rep(0, Time), Sigma=create_corr_matrix(Time))
      Y[i,] <- Y[i,] + c(eps)
    }
    mean(rowMeans(S))
    mean(S)


    Y <- t(matrix(t(Y), nrow=5, ncol=n*6))
    S <- t(matrix(t(S), nrow=5, ncol=n*6))
    A <- t(matrix(t(A), nrow=5, ncol=n*6))
    H <- rep(H, each=6)


    # Put all data together into a data frame. ----
    # colnames(Y) <- paste0("Y", 1:t)
    # colnames(H) <- paste0("h", 1:ncol(H))
    # colnames(S) <- paste0("s", 1:ncol(S))
    # colnames(A) <- paste0("a", 1:ncol(A))

    data <- list(
      df = list(Y=Y, H=H, S=S, A=A),
      params = list(
        beta10 = beta10,
        beta11 = beta11,
        theta1 = theta1,
        theta2 = theta2,
        p = p
      )
    )
    dfs[[rep]] <- data
  }

  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  filename <- paste("corr_setting_1","_n",n,"_m",m,"_t",t,"_beta", beta11,".rds", sep = "")
  saveRDS(dfs, paste(script_dir, "simulated_data",filename, sep = "/"))
}
