rm(list = ls())
# Do you do Python?
# Also, programs this long should be split into stages.
# Functions with more than ~50 lines of code I think are too long.
main <-function(m=100, setting=8, dfs=dfs) {
  ee.corr <- matrix(nrow = m, ncol = 2)
  for (rep in 1:m) {
    #rep = 1
    if (rep %% 10 == 0) {
      cat(rep, "\n")
    }
    data <- dfs[[rep]]$df
    #dfs[[1]]$params

    # Non correlated
    # y <- data[,c("Y1", "Y2", "Y3", "Y4", "Y5")]
    # a <- data[,c("a1", "a2", "a3", "a4", "a5")]
    # h <- data[,c("h1", "h2", "h3")]
    # s <- data[,c("s1", "s2", "s3")]
    # p <- 0.5

    # Correlated
    y <- data$Y
    a <- data$A
    h <- data$H
    s <- data$S
    p <- 0.5

    matrix <- get_p_a(a, p)
    cum_d <- matrix$cum_d
    p_a <- matrix$p_a
    #d_w <- matrix$d_w

    # No corr ----

    ## Setting 1: no covariates ----
    'init_beta <- c(100)
  rslt1 <- nleqslv(init_beta, function(beta) ee1(beta, y, a, h, s, p_a, cum_d, dose=1))
  rslt2 <- nleqslv(init_beta, function(beta) ee1(beta, y, a, h, s, p_a, cum_d, dose=2))
  rslt3 <- nleqslv(init_beta, function(beta) ee1(beta, y, a, h, s, p_a, cum_d, dose=3))
  rslt4 <- nleqslv(init_beta, function(beta) ee1(beta, y, a, h, s, p_a, cum_d, dose=4))
  rslt5 <- nleqslv(init_beta, function(beta) ee1(beta, y, a, h, s, p_a, cum_d, dose=5))
  ee.first[rep] <- rslt1$x
  cat("mine ",rslt1$x)'
    #ee1 <- cbind(ee, c(rslt1$x, rslt2$x, rslt3$x, rslt4$x, rslt5$x))
    #print(rep)

    ## Setting 2: Covariates ----
    'init_beta <- matrix(0, nrow=ncol(s), ncol=1)
  rslt1 <- nleqslv(init_beta, function(beta) ee2(beta, y, a, h, s, p_a, cum_d, dose=1))
  ee.second[rep, 1:3] <- rslt1$x'

    ## Setting 3: Boruvka Comparison ----
    'y.agg <- rowSums(y)
  a.agg <- rowSums(a)

  dose = 1
  index <- which(ifelse(a.agg == dose | a.agg == 0, 1, 0)==1)
  init_beta = c(10,10)
  rslt1 <- nleqslv(init_beta, function(beta) boruvka(beta, y.agg[index], a.agg[index], h[index], s[index], dose=1))
  ee.boruvka[rep] <- rslt1$x
  cat("boruv",rslt1$x,"\n")'

    ## MSM comparison ----
    'y.msm <- rowSums(y)
  dose <- rowSums(a)
  dose_steps <- data.frame(y.msm, dose)
  # dummy variable coding for dose
  dose_steps$dose0 <- ifelse(dose_steps$dose==0, 1, 0)
  dose_steps$dose1 <- ifelse(dose_steps$dose==1, 1, 0)
  dose_steps$dose2 <- ifelse(dose_steps$dose==2, 1, 0)
  dose_steps$dose3 <- ifelse(dose_steps$dose==3, 1, 0)
  dose_steps$dose4 <- ifelse(dose_steps$dose==4, 1, 0)
  dose_steps$dose5 <- ifelse(dose_steps$dose==5, 1, 0)

  # MSM model..?
  dose_steps$trtprob <- dbinom(dose_steps$dose, 5, p)
  dose_steps$weight <- 1 / dose_steps$trtprob
  dose_steps$wy <- dose_steps$weight * dose_steps$y.msm
  fit1 <- lm(y.msm ~  dose1 + dose2 + dose3 + dose4 + dose5, weights = weight, data=dose_steps)
  #fit1 <- lm(y.msm ~  0 + dose0 + dose1 + dose2 + dose3 + dose4 + dose5, data=dose_steps)
  #fit1 <- lm(wy ~  0 + dose0 + dose1 + dose2 + dose3 + dose4 + dose5, data=dose_steps)

  coef(fit1)
  ee_msm <- cbind(ee_msm, coef(fit1)[2:6])'

    ## Setting 4: Robust EE ----
    if (setting==4) {
      init_beta <- 0
      rslt1 <- nleqslv(init_beta, function(beta) ee1( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt4.1 <- nleqslv(init_beta, function(beta) ee4.1( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt4.2 <- nleqslv(init_beta, function(beta) ee4.2( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt4.3 <- nleqslv(init_beta, function(beta) ee4.3( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt4.4 <- nleqslv(init_beta, function(beta) ee4.4( beta, y, a, h, s, p_a, cum_d, dose ))

      ee.1[rep] <- rslt1$x
      ee.4.1[rep] <- rslt4.1$x
      ee.4.2[rep] <- rslt4.2$x
      ee.4.3[rep] <- rslt4.3$x
      ee.4.4[rep] <- rslt4.4$x
    }

    ## Setting 5: Robust EE misspecified----
    if (setting==5) {
      init_beta <- 0
      rslt1 <- nleqslv(init_beta, function(beta) ee1( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt5.1 <- nleqslv(init_beta, function(beta) ee5.1( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt5.2 <- nleqslv(init_beta, function(beta) ee5.2( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt5.3 <- nleqslv(init_beta, function(beta) ee5.3( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt5.4 <- nleqslv(init_beta, function(beta) ee5.4( beta, y, a, h, s, p_a, cum_d, dose ))

      ee.1[rep] <- rslt1$x
      ee.5.1[rep] <- rslt5.1$x
      ee.5.2[rep] <- rslt5.2$x
      ee.5.3[rep] <- rslt5.3$x
      ee.5.4[rep] <- rslt5.4$x
    }

    if (setting == 6) {
      init_beta <- c(0,0,0)
      rslt2 <- nleqslv(init_beta, function(beta) ee2( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt6.4 <- nleqslv(init_beta, function(beta) ee6.4( beta, y, a, h, s, p_a, cum_d, dose ))
      ee.2[rep,] <- rslt2$x
      ee.6.4[rep,] <- rslt6.4$x
    }

    if (setting == 7) {
      init_beta <- c(0,0,0)
      #rslt2 <- nleqslv(init_beta, function(beta) ee2( beta, y, a, h, s, p_a, cum_d, dose ))
      #rslt6.4 <- nleqslv(init_beta, function(beta) ee6.4( beta, y, a, h, s, p_a, cum_d, dose ))
      rslt6.4 <- nleqslv(init_beta, function(beta) ee6.4.improved( beta, y, a, h, s, p_a, cum_d, dose ))
      #ee.2[rep,] <- rslt2$x
      ee.6.4[rep,] <- rslt6.4$x
    }

    # Correlation ----
    ##  Correlation Setting 1: ----
    ##    treatment affects immediate outcome
    ##    Time varying S,
    ##    baseline H,
    ##    correlated error per participant

    if (setting == 8) {
      init_beta <- c(0,0)
      rslt <- nleqslv(init_beta, function(beta) ee.cor.2( beta, y, a, h, s, p_a, cum_d, dose ))
      ee.corr[rep,] <- rslt$x

      # debugonce(ee.cor.2)
      # ee.cor.2( init_beta, y, a, h, s, p_a, cum_d, dose )
      # init_beta <- c(0,0)
      # profvis::profvis({
      #   ee.cor.1( init_beta, y, a, h, s, p_a, cum_d, dose )
      # })
      #
      # profvis::profvis({
      #   nleqslv(init_beta, function(beta) ee.cor.2( beta, y, a, h, s, p_a, cum_d, dose ))
      #})

    }
  }
  return(ee.corr)
}


source("estimating_equations.R")
library(nleqslv)


# Simulated data ----
# get current script directory

#script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_dir <- getwd()


# # first setting
# filename <- "first_setting.rds"
#
# # second setting beta = (100, 20, 0) eta = (50, 0, 0)
# filename <- "second_setting_n150_m1000_t5_p0.5_dose1.rds"
# filename <- "second_setting_n500_m1000_t5_p0.5_dose1.rds"
#
# # fourth setting
# filename <- "fourth_setting_n100_m1000_t5_p0.5_beta20_dose1.rds"
# filename <- "fourth_setting_n300_m1000_t5_p0.5_beta20_dose1.rds"
#
# # fifth setting - check robustness with misspecification
# filename <- "fifth_setting_n100_m1000_t5_p0.5_beta20_dose1.rds"
# filename <- "fifth_setting_n500_m1000_t5_p0.5_beta20_dose1.rds"
# filename <- "fifth_setting_n300_m1000_t5_p0.5_beta20_dose1.rds"
#
# # sitxh setting - check robustness with misspecification and non null S
# filename <- "sixth_setting_n100_m1000_t5_p0.5_beta20_dose1.rds"
# filename <- "sixth_setting_n300_m1000_t5_p0.5_beta20_dose1.rds"
#
# # seventh setting - extend sixth setting to multiple days
# filename <- "seventh_setting_n100_m1000_days10_t5_p0.5_beta20_dose1.rds"
# filename <- "seventh_setting_n300_m1000_days10_t5_p0.5_beta20_dose1.rds"
# filename <- "seventh_setting_n50_m1000_days10_t5_p0.5_beta20_dose1.rds"

# Lines like this could probably be optimized for clarity,
# i.e, a for-loop iterating through the various combinations.
# eighth setting - correlation setting 1
filename <- "corr_setting_1_n100_m100_t30_beta0.2_theta20.5.rds"
filename <- "corr_setting_1_n200_m100_t30_beta0.2_theta20.5.rds"
filename <- "corr_setting_1_n300_m100_t30_beta0.2_theta20.5.rds"
filename <- "corr_setting_1_n300_m100_t30_beta0.2_theta20.rds"
filename <- "corr_setting_2_n300_m400_t30_beta0.2_theta20.rds"
filename <- "corr_setting_2_n300_m400_t30_beta0.2_theta20.5.rds"
filename <- "corr_setting_2_n300_m400_t30_beta0.2_theta20.8.rds"
filename <- "corr_setting_3_n300_m400_t30_beta0.2_theta20.8.rds"
filename <- "corr_setting_4_n300_m400_t30_beta0.2_theta20.8.rds"
filename <- "corr_setting_5_n300_m400_t30_beta0.2_theta20.8.rds"

#dfs <- readRDS(paste(script_dir, "simulated_data",filename, sep = "/"))
dfs <- readRDS(paste(script_dir, "simulated_data_original/correlation",filename, sep = "/"))

m <- length(dfs)
#dfs[[1]]$params
dose=1

# result place holders
# ee.1 <- rep(0, m)
#
# ee.2 <- matrix(nrow = m, ncol = 3)
#
# ee.4.1 <- rep(0, m)
# ee.4.2 <- rep(0, m)
# ee.4.3 <- rep(0, m)
# ee.4.4 <- rep(0, m)
#
# ee.5.1 <- rep(0, m)
# ee.5.2 <- rep(0, m)
# ee.5.3 <- rep(0, m)
# ee.5.4 <- rep(0, m)
#
# ee.6.4 <- matrix(nrow = m, ncol = 3)
# ee.6.4.improved <- matrix(nrow = m, ncol = 3)

ee.corr <- matrix(nrow = m, ncol = 2)

# ee.boruvka <- rep(0, m)
setting = 8
#debugonce(main)
ee.corr <- main(m, setting, dfs)

cat("\n\n")
colMeans(ee.corr, na.rm = TRUE)
apply(ee.corr, 2, median, na.rm = TRUE)
apply(ee.corr, 2, sd, na.rm = TRUE)

hist(ee.corr[,2])

#
# profvis::profvis({
#   main(m, setting, dfs)
# })
# # benchmark ee.cor.1, and ee.cor.2
# data <- dfs[[rep]]$df
# y <- data[,c("Y1", "Y2", "Y3", "Y4", "Y5")]
# a <- data[,c("a1", "a2", "a3", "a4", "a5")]
# h <- data[,c("h1", "h2", "h3")]
# s <- data[,c("s1", "s2", "s3")]
# p <- 0.5
# matrix <- get_p_a(a, p)
# cum_d <- matrix$cum_d
# p_a <- matrix$p_a
#
# bench <- bench::mark(
#   ee6.4( c(0,0,0), y, a, h, s, p_a, cum_d, dose ),
#   ee6.4.improved( c(0,0,0), y, a, h, s, p_a, cum_d, dose ), relative = TRUE
# )

# So much commented code!!!!
# truncate extreme values of ee
# par(mfrow=c(1,2))
# ee <- ee[abs(ee) < 1000]
# median(ee[1,])
# abs(rowMeans(ee) - c(20,40,60,80,100))
# apply(ee, 1, sd)
# cov <- rowSums(abs(ee - c(20,40,60,80,100)) < 2*apply(ee, 1, sd))/1000
#
# median(ee_msm)
# abs(rowMeans(ee_msm, na.rm = TRUE) - c(20,40,60,80,100))
# apply(ee_msm, 1, sd, na.rm = TRUE)
# cov <- (ee_msm - c(20,40,60,80,100))
# cov_msm <- rowSums(abs(ee_msm - c(20,40,60,80,100)) < 2*apply(ee, 1, sd, na.rm = TRUE), na.rm = TRUE)/1000
# table(dose_steps$dose)
#
# par(mfrow=c(1,3))
# hist(ee[,1], breaks = 20, main = "beta1")
# hist(ee[,2], breaks = 20, main = "beta2")
# hist(ee[,3], breaks = 20, main = "beta3")

# # setting 4 rslts
# mean(ee.1[1:m])
# median(ee.1[1:m])
# var(ee.1[1:m])
# cat("\n\n")
# mean(ee.4.1[1:m])
# median(ee.4.1[1:m])
# var(ee.4.1[1:m])
# cat("\n\n")
# mean(ee.4.2[1:m])
# median(ee.4.2[1:m])
# var(ee.4.2[1:m])
# cat("\n\n")
# mean(ee.4.3[1:m])
# median(ee.4.3[1:m])
# var(ee.4.3[1:m])
# cat("\n\n")
# mean(ee.4.4[1:m])
# median(ee.4.4[1:m])
# var(ee.4.4[1:m])
# cat("\n\n")
# cat("\n\n")
#
# # setting 5 rslts
# mean(ee.1[1:m])
# median(ee.1[1:m])
# sd(ee.1[1:m])
# cat("\n\n")
# mean(ee.5.1[1:m])
# median(ee.5.1[1:m])
# sd(ee.5.1[1:m])
# cat("\n\n")
# mean(ee.5.2[1:m])
# median(ee.5.2[1:m])
# sd(ee.5.2[1:m])
# cat("\n\n")
# mean(ee.5.3[1:m])
# median(ee.5.3[1:m])
# sd(ee.5.3[1:m])
# cat("\n\n")
# mean(ee.5.4[1:m])
# median(ee.5.4[1:m])
# sd(ee.5.4[1:m])

# setting 6 rslts
# cat("\n\n")
# colMeans(ee.2, na.rm = TRUE)
# apply(ee.2, 2, median, na.rm = TRUE)
# apply(ee.2, 2, sd, na.rm = TRUE)
#

#
#
#

