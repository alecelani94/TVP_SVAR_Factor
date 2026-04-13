## AR_scales: OLS residual variance from AR(P) for each variable
## Used to construct the Minnesota prior scaling factors sigma_hat_i^2
AR_scales <- function(Data, P) {

  T0 <- nrow(Data)
  N  <- ncol(Data)

  X <- matrix(1, T0 - P, 1)
  for (l in 1:P) {
    X <- cbind(X, Data[(P + 1 - l):(T0 - l), ])
  }
  Y <- Data[(P + 1):T0, ]

  s <- numeric(N)
  for (i in 1:N) {
    cols <- c(1, 1 + i + N * (0:(P - 1)))
    b    <- solve(crossprod(X[, cols]), crossprod(X[, cols], Y[, i]))
    s[i] <- var(Y[, i] - X[, cols] %*% b)
  }
  s
}

## Minnesota: extended prior variance for [b_i, lambda_i]
## Returns:
##   C0 : (K+R) x N  constant scaling (sigma ratios + lag decay, no delta_2/3/4)
##   V0 : (K+R) x N  full prior variance (C0 * delta_2/3/4 applied)
## C0 is kept separate so delta_2, delta_3, delta_4 can be re-estimated
## without rebuilding the scaling structure
Minnesota <- function(scales, P, R, delta) {

  N <- length(scales)
  K <- N * P + 1

  # C0: constant part (intercept includes delta_1, lags include decay delta_5)
  C0 <- matrix(0, K + R, N)
  for (i in 1:N) {
    C0[1, i] <- delta[1] * scales[i]
    for (l in 1:P) {
      for (j in 1:N) {
        if (i == j) {
          C0[1 + (l - 1) * N + j, i] <- 1 / l^delta[5]
        } else {
          C0[1 + (l - 1) * N + j, i] <- scales[i] / scales[j] / l^delta[5]
        }
      }
    }
    for (j in 1:R) {
      C0[K + j, i] <- scales[i]
    }
  }

  # V0: multiply C0 by delta_2 (own), delta_3 (cross), delta_4 (loadings)
  V0 <- C0
  for (i in 1:N) {
    for (l in 1:P) {
      for (j in 1:N) {
        row <- 1 + (l - 1) * N + j
        if (i == j) {
          V0[row, i] <- delta[2] * C0[row, i]
        } else {
          V0[row, i] <- delta[3] * C0[row, i]
        }
      }
    }
    for (j in 1:R) {
      V0[K + j, i] <- delta[4] * C0[K + j, i]
    }
  }
  list(V0 = V0, C0 = C0)
}

MCMC_BVAR_Factor <- function(Data, P = 12, R = 3,
                             delta = c(100, 0.2^2, 0.1^2, 0.2^2, 2),
                             alpha_d = 0.5,
                             MCMC = 5000, Burn = 1000, Thin = 1,
                             disp = 1000, verbose = TRUE) {

  ## dimensions
  T0 <- nrow(Data)
  N  <- ncol(Data)
  Tt <- T0 - P
  K  <- N * P + 1

  ## VAR matrices: X = [1, Y_{-1}, ..., Y_{-P}], Y = Y_{P+1:T0}
  X <- matrix(1, Tt, 1)
  for (l in 1:P) {
    X <- cbind(X, Data[(P + 1 - l):(T0 - l), ])
  }
  Y <- Data[(P + 1):T0, ]

  ## prior construction
  scales <- AR_scales(Data, P)
  mn     <- Minnesota(scales, P, R, delta)
  V0     <- mn$V0
  C0     <- mn$C0
  v_d    <- scales

  ## storage
  nkeep     <- MCMC %/% Thin
  B_save <- array(0, c(nkeep, K, N))
  L_save    <- array(0, c(nkeep, N, R))
  F_save    <- array(0, c(nkeep, Tt, R))
  d_save    <- matrix(0, nkeep, N)

  ## initialization: OLS for b_i, PCA for L and Ft
  ## (normalization: L'L/N = I, factors not unit-variance; doesn't matter,
  ##  overwritten after first Gibbs sweep)
  B <- matrix(0, K, N)
  d     <- numeric(N)
  for (i in 1:N) {
    B[, i] <- solve(crossprod(X), crossprod(X, Y[, i]))
    d[i]       <- var(Y[, i] - X %*% B[, i])
  }
  v_hat <- Y - X %*% B
  ev    <- eigen(crossprod(scale(v_hat)), symmetric = TRUE)
  L     <- sqrt(N) * ev$vectors[, 1:R]
  Ft    <- scale(v_hat) %*% L / N

  ## precompute: X'X and X'Y (don't change across sweeps)
  XtX <- crossprod(X)                    # K x K
  XtY <- crossprod(X, Y)                 # K x N

  ## precompute: prior precisions (diagonal, per equation)
  iV0 <- 1 / V0[1:(K + R), ]            # (K+R) x N

  ## precompute: GIG parameter that doesn't change
  gig_lam <- alpha_d - Tt / 2
  gig_psi <- 2 * alpha_d / v_d          # length N

  ## MCMC loop
  Nsim <- Burn + Thin * MCMC
  k  <- 0
  t0 <- Sys.time()
  M  <- K + R                            # dimension of theta_i

  for (it in 1:Nsim) {
    if (verbose && it %% disp == 0) {
      message(sprintf("iter %d / %d  (%.1fs)", it, Nsim,
                      as.numeric(Sys.time() - t0, units = "secs")))
    }

    ## Step 1: factors f_t | rest  (joint Normal, t = 1,...,Tt)
    v_hat  <- Y - X %*% B
    iD     <- 1 / d
    LtiD   <- t(L) * rep(iD, each = R)   # R x N: Lambda' D^{-1} row-scaled
    Vf     <- solve(diag(R) + LtiD %*% L)
    Cf     <- t(chol(Vf))
    Ft_mean <- v_hat %*% t(LtiD) %*% t(Vf)  # Tt x R: each row = Vf Lambda' D^{-1} u_hat_t
    Ft     <- Ft_mean + matrix(rnorm(Tt * R), Tt, R) %*% t(Cf)
    Ft     <- Ft - rep(1, Tt) %*% t(colMeans(Ft))

    ## Step 2: theta_i = [b_i, lambda_i] | rest  (equation by equation)
    ## precompute blocks involving Ft (change every sweep, but same for all i)
    FtF  <- crossprod(Ft)                # R x R
    XtF  <- crossprod(X, Ft)             # K x R
    FtY  <- crossprod(Ft, Y)             # R x N
    ZtZ  <- rbind(cbind(XtX, XtF),
                  cbind(t(XtF), FtF))    # M x M
    ZtYi <- rbind(XtY, FtY)             # M x N

    for (i in 1:N) {
      # posterior precision = Z'Z / d_i + diag(1/v0_i)
      iV <- ZtZ / d[i]
      diag(iV) <- diag(iV) + iV0[, i]

      # Cholesky solve for posterior mean, then draw
      Ch <- chol(iV)
      m  <- backsolve(Ch, forwardsolve(t(Ch), ZtYi[, i] / d[i]))
      th <- m + backsolve(Ch, rnorm(M))

      B[, i] <- th[1:K]
      L[i, ] <- th[(K + 1):M]
    }

    ## Step 3: d_i | rest  (GIG draw, Gamma prior)
    e   <- Y - X %*% B - Ft %*% t(L)
    sse <- colSums(e^2)
    for (i in 1:N) {
      d[i] <- GIGrvg::rgig(1, lambda = gig_lam, chi = sse[i], psi = gig_psi[i])
    }

    ## store
    if (it > Burn && (it - Burn) %% Thin == 0) {
      k <- k + 1
      B_save[k, , ] <- B
      L_save[k, , ] <- L
      F_save[k, , ] <- Ft
      d_save[k, ]   <- d
    }
  }

  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  if (verbose) {
    message(sprintf("done: %d draws in %.1fs (%.0f draws/sec)", nkeep, elapsed, nkeep / elapsed))
  }

  list(B = B_save, L = L_save, Ft = F_save, d = d_save, time = elapsed)
}
