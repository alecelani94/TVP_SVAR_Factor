cat("\014")
rm(list = ls())

source("Data/prepare_data.R")
source("Models/MCMC_BVAR_Factor.R")
source("simulate_DGP.R")

set.seed(123456)

#################
#### setting ####
#################

## real data (for DGP calibration)
d      <- prepare_data("Data/data_for_MC.xlsx")
Data   <- d$Data
vnames <- d$names

## dimensions
Tt <- 1000
N  <- ncol(Data)     # 14
P  <- 6              # estimation lags
R  <- 5              # number of factors

## MCMC
MCMC <- 10000
Burn <- 2500
Thin <- 1
disp <- 2500          # print every disp iterations

################
#### priors ####
################

## Extended Minnesota (fixed, non-hierarchical)
##   delta[1] : intercept variance scale
##   delta[2] : own-lag variance
##   delta[3] : cross-lag variance
##   delta[4] : factor loading variance (factor var = 1, scales by sigma_i)
##   delta[5] : lag-decay exponent
delta   <- c(10^2, 0.2^2, 0.1^2, 0.4^2, 2)

## Idiosyncratic variances: d_i ~ Ga(alpha_d, alpha_d / v_d)
##   E(d_i) = v_d = sigma_hat_i^2 (AR scale, set inside sampler)
##   alpha_d = 0.5 corresponds to sqrt(d_i) ~ N(0, v_d)
alpha_d <- 0.5

###############################
#### true parameters (OLS) ####
###############################

dgp <- get_param_realGDP(Data, P_dgp = P, R_gdp = R)

## true covariance
Sigma_true <- dgp$Lambda %*% t(dgp$Lambda) + dgp$D

########################
#### simulate & run ####
########################

## simulate from VAR(P) DGP
sim <- simulate_VAR(dgp$B, dgp$Lambda, dgp$D, P_dgp = P, T_sim = Tt)

## estimate on simulated data
out <- MCMC_BVAR_Factor(sim$Y_sim, P = P, R = R,
                        delta = delta, alpha_d = alpha_d,
                        MCMC = MCMC, Burn = Burn, Thin = Thin, disp = disp)

cat(sprintf("done in %.1fs\n", out$time))

######################################
#### compare Sigma = LL' + D chain ####
######################################

## compute Sigma_s = L_s L_s' + D_s for each stored draw
nkeep <- dim(out$L)[1]
Sigma_chain <- array(0, c(nkeep, N, N))
for (s in 1:nkeep) {
  Ls <- out$L[s, , ]
  Sigma_chain[s, , ] <- Ls %*% t(Ls) + diag(out$d[s, ])
}

## plot lower-triangular + diagonal elements as trace plots with true value
## N diagonal + N*(N-1)/2 lower = N*(N+1)/2 = 105 elements
## arrange in N x N grid, only lower triangle + diagonal
par(mfrow = c(N, N), mar = c(1, 1, 1, 0.5), oma = c(2, 2, 2, 0))
for (i in 1:N) {
  for (j in 1:N) {
    if (j > i) {
      plot.new()
    } else {
      plot(1:nkeep, Sigma_chain[, i, j], type = "l", col = "grey50",
           xlab = "", ylab = "", xaxt = "n", yaxt = "n")
      abline(h = Sigma_true[i, j], col = "red", lwd = 2)
      if (i == j) {
        mtext(vnames[i], side = 3, line = 0, cex = 0.4)
      }
    }
  }
}
mtext("Sigma chain (grey) vs true (red)", outer = TRUE, side = 3, line = 0.5, cex = 0.8)

##############################
#### MCMC efficiency (ESS) ####
##############################

library(coda)
K     <- N * P + 1
nkeep <- dim(out$B)[1]

## ESS % for each parameter class
ess_B <- numeric(K * N)
idx <- 0
for (i in 1:N) {
  for (k in 1:K) {
    idx <- idx + 1
    ess_B[idx] <- 100 * effectiveSize(as.mcmc(out$B[, k, i])) / nkeep
  }
}

ess_L <- numeric(N * R)
idx <- 0
for (i in 1:N) {
  for (j in 1:R) {
    idx <- idx + 1
    ess_L[idx] <- 100 * effectiveSize(as.mcmc(out$L[, i, j])) / nkeep
  }
}

ess_d <- numeric(N)
for (i in 1:N) {
  ess_d[i] <- 100 * effectiveSize(as.mcmc(out$d[, i])) / nkeep
}

Tt_est <- dim(out$Ft)[2]
ess_Ft <- numeric(Tt_est * R)
idx <- 0
for (j in 1:R) {
  for (tt in seq(1, Tt_est, length.out = min(Tt_est, 50))) {  # subsample for speed
    idx <- idx + 1
    ess_Ft[idx] <- 100 * effectiveSize(as.mcmc(out$Ft[, round(tt), j])) / nkeep
  }
}
ess_Ft <- ess_Ft[1:idx]

## ESS % for lower-triangular + diagonal of Sigma = LL' + D
ess_Sig <- numeric(N * (N + 1) / 2)
idx <- 0
for (i in 1:N) {
  for (j in 1:i) {
    idx <- idx + 1
    ess_Sig[idx] <- 100 * effectiveSize(as.mcmc(Sigma_chain[, i, j])) / nkeep
  }
}

## boxplot
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
boxplot(list(B = ess_B, Lambda = ess_L, Factors = ess_Ft, D = ess_d, Sigma = ess_Sig),
        main = "ESS % by parameter class",
        ylab = "ESS %", ylim = c(0, 100),
        col = c("steelblue", "tomato", "goldenrod", "seagreen", "orchid"))
abline(h = c(10, 50), col = "grey", lty = 2)
