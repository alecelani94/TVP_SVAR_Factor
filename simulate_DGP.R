## get_param_realGDP: extract true VAR parameters from real US data via OLS
##
## Steps:
##   1. OLS VAR(P_dgp) -> B, Sigma
##   2. PCA on residuals (raw, not standardized) -> first R eigenvectors
##   3. Normalize factors to unit variance: var(F_j) = 1
##   4. OLS loadings: Lambda = (F'F)^{-1} F' resid, transposed to N x R
##   5. Sign flip: ensure Lambda[j,j] > 0 for j = 1,...,R
##   6. Idiosyncratic: D = diag(var(resid - F Lambda'))
##
## Returns list with:
##   B, Sigma, Lambda, D, Ft

get_param_realGDP <- function(Data, P_dgp = 2, R_gdp = 3) {

  T0 <- nrow(Data)
  N  <- ncol(Data)
  Tt <- T0 - P_dgp

  ## OLS VAR(P_dgp)
  X <- matrix(1, Tt, 1)
  for (l in 1:P_dgp) {
    X <- cbind(X, Data[(P_dgp + 1 - l):(T0 - l), ])
  }
  Y     <- Data[(P_dgp + 1):T0, ]
  B     <- solve(crossprod(X), crossprod(X, Y))
  resid <- Y - X %*% B
  Sigma <- crossprod(resid) / (Tt - P_dgp - 1)

  ## PCA on residuals
  ev       <- eigen(crossprod(resid), symmetric = TRUE)
  pct_var  <- ev$values / sum(ev$values)
  Ft       <- resid %*% ev$vectors[, 1:R_gdp]

  ## unit-variance factors
  for (j in 1:R_gdp) {
    Ft[, j] <- Ft[, j] / sd(Ft[, j])
  }

  ## OLS loadings (N x R)
  Lambda <- t(solve(crossprod(Ft), crossprod(Ft, resid)))

  ## sign: Lambda[j,j] > 0
  for (j in 1:R_gdp) {
    if (Lambda[j, j] < 0) {
      Lambda[, j] <- -Lambda[, j]
      Ft[, j]     <- -Ft[, j]
    }
  }

  ## idiosyncratic variance
  D <- diag(colMeans((resid - Ft %*% t(Lambda))^2))

  list(B = B, Sigma = Sigma, Lambda = Lambda, D = D, Ft = Ft, pct_var = pct_var)
}

## simulate_VAR: generate data from VAR(P) with factor structure on errors
##   u_t = Lambda * f_t + e_t,  f_t ~ N(0,I_R),  e_t ~ N(0,D)
##   y_t = B x_t + u_t
##
## Returns list with:
##   Y_sim  : T_sim x N
##   Ft_sim : T_sim x R  (true factors)

simulate_VAR <- function(B, Lambda, D, P_dgp, T_sim = 1000, T_burn = 1000) {

  N       <- ncol(B)
  R_gdp   <- ncol(Lambda)
  T_total <- T_sim + T_burn
  cholD   <- sqrt(diag(D))        # D is diagonal, so chol = sqrt of diagonal

  y  <- matrix(0, P_dgp + T_total, N)
  ft <- matrix(0, P_dgp + T_total, R_gdp)

  for (tt in (P_dgp + 1):(P_dgp + T_total)) {
    xlag <- 1
    for (l in 1:P_dgp) {
      xlag <- c(xlag, y[tt - l, ])
    }
    ft[tt, ] <- rnorm(R_gdp)
    e_t      <- cholD * rnorm(N)
    y[tt, ]  <- xlag %*% B + as.numeric(Lambda %*% ft[tt, ]) + e_t
  }

  idx <- (P_dgp + T_burn + 1):(P_dgp + T_total)
  list(Y_sim  = y[idx, ],
       Ft_sim = ft[idx, ])
}
