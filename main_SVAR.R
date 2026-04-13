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
R  <- 3              # number of factors / structural shocks

## MCMC
MCMC <- 10000
Burn <- 2500
Thin <- 1
disp <- 2500

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
alpha_d <- 0.5

## Sign restrictions on Lambda — set after extracting true parameters below

###############################
#### true parameters (OLS) ####
###############################

dgp <- get_param_realGDP(Data, P_dgp = P, R_gdp = R)

## true covariance
Sigma_true <- dgp$Lambda %*% t(dgp$Lambda) + dgp$D

## sign restrictions from true Lambda: sg[i,j] = sign(Lambda[i,j])
sg <- sign(dgp$Lambda)


## bounds for hdtg
lb <- matrix(-Inf, N, R)
ub <- matrix( Inf, N, R)
for (i in 1:N) {
  for (j in 1:R) {
    if (sg[i, j] ==  1) { lb[i, j] <- 0 }
    if (sg[i, j] == -1) { ub[i, j] <- 0 }
  }
}

## display
rownames(dgp$Lambda) <- vnames
colnames(dgp$Lambda) <- paste0("S", 1:R)
cat("True Lambda:\n")
print(round(dgp$Lambda, 3))

cat("\nTrue diag(D):\n")
names(diag(dgp$D)) <- vnames
print(round(diag(dgp$D), 4))

########################
#### simulate & run ####
########################

## simulate from VAR(P) DGP
sim <- simulate_VAR(dgp$B, dgp$Lambda, dgp$D, P_dgp = P, T_sim = Tt)

## estimate on simulated data (TODO: replace with MCMC_BSVAR_Factor using hdtg)
# out <- MCMC_BSVAR_Factor(sim$Y_sim, P = P, R = R,
#                          delta = delta, alpha_d = alpha_d,
#                          sg = sg, lb = lb, ub = ub,
#                          MCMC = MCMC, Burn = Burn, Thin = Thin, disp = disp)
