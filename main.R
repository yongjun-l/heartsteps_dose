source("~/Library/CloudStorage/GoogleDrive-yongjun.lee5@gmail.com/My Drive/1. UCI/2024-1 Winter/Tianchen DIS/simulations/working scripts/estimating_equations.R")
library(rootSolve)
library(nleqslv)


# Simulated data ----
# get current script directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
script_dir <- getwd()


dfs <- readRDS(paste(script_dir, "simulated_data","first_setting.rds", sep = "/"))

# second setting beta = (100, 20, 0) eta = (50, 0, 0)
dfs <- readRDS(paste(script_dir, "simulated_data","second_setting.rds", sep = "/"))
dfs <- readRDS(paste(script_dir, "simulated_data","third_setting.rds", sep = "/"))

filename <- "second_setting_n150_m1000_t5_p0.5_dose1.rds"
filename <- "second_setting_n500_m1000_t5_p0.5_dose1.rds"
# fourth setting
filename <- "fourth_setting_n100_m1000_t5_p0.5_dose1.rds"
# fourth setting debug
filename <- "fourth_setting_n100_m1000_t5_p0.5_beta20_dose1.rds"
filename <- "fourth_setting_n300_m1000_t5_p0.5_beta20_dose1.rds"

# fifth setting - check robustness with misspecification
filename <- "fifth_setting_n100_m1000_t5_p0.5_beta20_dose1.rds"
filename <- "fifth_setting_n300_m1000_t5_p0.5_beta20_dose1.rds"

dfs <- readRDS(paste(script_dir, "simulated_data",filename, sep = "/"))

m <- length(dfs)

dfs[[1]]$params

dose=1
ee <- c(0,0,0,0,0)
ee_msm <- c(0,0,0,0,0)
ee <- matrix(0, nrow = m, ncol = 5*3 )
ee.1 <- rep(0, m)

ee.2 <- matrix(0, nrow = m, ncol = 3)

ee.4.2 <- rep(0, m)
ee.4.3 <- rep(0, m)
ee.4.4 <- rep(0, m)

ee.5.2 <- rep(0, m)
ee.5.3 <- rep(0, m)
ee.5.4 <- rep(0, m)

ee.boruvka <- rep(0, m)

#m=100

for (rep in 1:m) {
  #rep = 1
  if (rep %% 1 == 0) {
  cat(rep, "\n")
  }
  data <- dfs[[rep]]$df
  #dfs[[1]]$params
  y <- data[,c("Y1", "Y2", "Y3", "Y4", "Y5")]
  a <- data[,c("a1", "a2", "a3", "a4", "a5")]
  h <- data[,c("h1", "h2", "h3")]
  s <- data[,c("s1", "s2", "s3")]
  p <- 0.5

  matrix <- get_p_a(a, p)
  cum_d <- matrix$cum_d
  p_a <- matrix$p_a
  #d_w <- matrix$d_w

  # Proposed model ----

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



  # Setting 2: Covariates ----
  'init_beta <- matrix(0, nrow=ncol(s), ncol=1)
  rslt1 <- nleqslv(init_beta, function(beta) ee2(beta, y, a, h, s, p_a, cum_d, dose=1))
  ee.second[rep, 1:3] <- rslt1$x'



  # Setting 3: Boruvka Comparison ----
  'y.agg <- rowSums(y)
  a.agg <- rowSums(a)

  dose = 1
  index <- which(ifelse(a.agg == dose | a.agg == 0, 1, 0)==1)
  init_beta = c(10,10)
  rslt1 <- nleqslv(init_beta, function(beta) boruvka(beta, y.agg[index], a.agg[index], h[index], s[index], dose=1))
  ee.boruvka[rep] <- rslt1$x
  cat("boruv",rslt1$x,"\n")'

  # MSM comparison ----
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

  # Setting 4: Robust EE ----
  init_beta <- 0
  rslt1 <- nleqslv(init_beta, function(beta) ee1( beta, y, a, h, s, p_a, cum_d, dose ))
  rslt4.2 <- nleqslv(init_beta, function(beta) ee4.2( beta, y, a, h, s, p_a, cum_d, dose ))
  rslt4.3 <- nleqslv(init_beta, function(beta) ee4.3( beta, y, a, h, s, p_a, cum_d, dose ))
  rslt4.4 <- nleqslv(init_beta, function(beta) ee4.4( beta, y, a, h, s, p_a, cum_d, dose ))

  ee.1[rep] <- rslt1$x
  ee.4.2[rep] <- rslt4.2$x
  ee.4.3[rep] <- rslt4.3$x
  ee.4.4[rep] <- rslt4.4$x


  # Setting 5: Robust EE misspecified----
  init_beta <- 0
  rslt1 <- nleqslv(init_beta, function(beta) ee1( beta, y, a, h, s, p_a, cum_d, dose ))
  rslt5.2 <- nleqslv(init_beta, function(beta) ee5.2( beta, y, a, h, s, p_a, cum_d, dose ))
  rslt5.3 <- nleqslv(init_beta, function(beta) ee5.3( beta, y, a, h, s, p_a, cum_d, dose ))
  rslt5.4 <- nleqslv(init_beta, function(beta) ee5.4( beta, y, a, h, s, p_a, cum_d, dose ))

  ee.1[rep] <- rslt1$x
  ee.5.2[rep] <- rslt5.2$x
  ee.5.3[rep] <- rslt5.3$x
  ee.5.4[rep] <- rslt5.4$x

}


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

mean(ee.1[1:m])
median(ee.1[1:m])
var(ee.1[1:m])
mean(ee.4.2[1:m])
median(ee.4.2[1:m])
var(ee.4.2[1:m])
mean(ee.4.3[1:m])
median(ee.4.3[1:m])
var(ee.4.3[1:m])
mean(ee.4.4[1:m])
median(ee.4.4[1:m])
var(ee.4.4[1:m])


mean(ee.1[1:m])
median(ee.1[1:m])
var(ee.1[1:m])
mean(ee.5.2[1:m])
median(ee.5.2[1:m])
var(ee.5.2[1:m])
mean(ee.5.3[1:m])
median(ee.5.3[1:m])
var(ee.5.3[1:m])
mean(ee.5.4[1:m])
median(ee.5.4[1:m])
var(ee.5.4[1:m])




hist(ee.first[1:m], breaks=20)
hist(ee.fourth[1:m], breaks=30)
#
#
# mean(ee.second[,1])
# mean(ee.second[,2])
# mean(ee.second[,3])
# #
# sd(ee.second[,1])
# sd(ee.second[,2])
# sd(ee.second[,3])
#
# hist(ee[2,], breaks = 20)
# hist(ee[3,], breaks = 20)
# hist(ee[4,], breaks = 20)
# hist(ee[5,], breaks = 20)

#max(ee.boruvka)
#hist(ee.boruvka, breaks = 20)
#hist(ee.first, breaks = 20)

# cat("\n\n")
# mean(ee.first)
# median(ee.first)
# sd(ee.first)
# cat("\n\n")
# mean(ee.boruvka)
# median(ee.boruvka)
# sd(ee.boruvka)











