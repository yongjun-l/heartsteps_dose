source("~/Library/CloudStorage/GoogleDrive-yongjun.lee5@gmail.com/My Drive/1. UCI/2024-1 Winter/Tianchen_excursion/simulations/working scripts/estimating_equations.R")
library(rootSolve)
library(nleqslv)
script_dir <- getwd()


main.linear <-function(dfs, m=100) {
  ee.6.4 <- matrix(nrow = m, ncol = 3)
  for (rep in 1:m) {
    data <- dfs[[rep]]$df

    y <- data[,c("Y1", "Y2", "Y3", "Y4", "Y5")]
    a <- data[,c("a1", "a2", "a3", "a4", "a5")]
    h <- data[,c("h1", "h2", "h3")]
    s <- data[,c("s1", "s2", "s3")]
    p <- 0.5

    matrix <- get_p_a(a, p)
    cum_d <- matrix$cum_d
    p_a <- matrix$p_a


    init_beta <- c(0,0,0)
    #rslt2 <- nleqslv(init_beta, function(beta) ee2( beta, y, a, h, s, p_a, cum_d, dose ))
    #rslt6.4 <- nleqslv(init_beta, function(beta) ee6.4( beta, y, a, h, s, p_a, cum_d, dose ))
    rslt6.4 <- nleqslv(init_beta, function(beta) ee6.4.improved( beta, y, a, h, s, p_a, cum_d, dose ))
    #ee.2[rep,] <- rslt2$x
    ee.6.4[rep,] <- rslt6.4$x
  }
  return(ee.6.4)
}

library(parallel)
main.parallel <-function(dfs, m=100) {
  rslt.parallel <- mclapply(1:m, function(rep) {
    data <- dfs[[rep]]$df

    y <- data[,c("Y1", "Y2", "Y3", "Y4", "Y5")]
    a <- data[,c("a1", "a2", "a3", "a4", "a5")]
    h <- data[,c("h1", "h2", "h3")]
    s <- data[,c("s1", "s2", "s3")]
    p <- 0.5

    matrix <- get_p_a(a, p)
    cum_d <- matrix$cum_d
    p_a <- matrix$p_a

    init_beta <- c(0,0,0)
    rslt6.4 <- nleqslv(init_beta, function(beta) ee6.4.improved( beta, y, a, h, s, p_a, cum_d, dose ))
    return(rslt6.4$x)
  }, mc.cores = 4)

  ee.6.4 <- do.call(rbind, rslt.parallel)
  return(ee.6.4)
}
