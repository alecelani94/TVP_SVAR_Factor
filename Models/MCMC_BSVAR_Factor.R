# MCMC_BSVAR_Factor.R
#
# R port of Korobilis (2022)'s Gibbs_FSR.m: a Gibbs sampler for a Bayesian
# VAR with a static factor structure on the reduced-form residuals,
#
#       y_t = c + B1 y_{t-1} + ... + Bp y_{t-p} + v_t
#       v_t = L f_t + e_t,    f_t ~ N(0, I_N),   e_t ~ N(0, R)
#
# with sign restrictions on the loadings L (matrix `sg`, entries in
# {-1, 0, 1, NA}) and a horseshoe prior on each VAR equation's coefficients.
#
# This file is a *faithful* port of the original MATLAB sampler. In
# particular, the loadings L are still drawn one element at a time using a
# univariate truncated normal — exactly the inefficient step Korobilis flags
# in the paper. A second version of this sampler will replace that block
# with a joint draw via the `hdtg` package (sampling efficient
# high-dimensional truncated Gaussians); that is the whole point of the
# project.
#
# Inputs
#   Y      : T x M data matrix (already stationary / standardized)
#   p      : VAR lag length
#   sg     : M x N matrix of sign restrictions on the loadings
#               1   -> positive   (L_ij > 0)
#              -1   -> negative   (L_ij < 0)
#               0   -> exact zero
#              NA   -> unrestricted
#   nhor   : IRF horizon
#   ngibbs : number of post-burn-in draws to keep (before thinning)
#   nburn  : burn-in iterations
#   nthin  : thinning factor
#
# Output: list with
#   beta   : (ngibbs/nthin) x KK x M     VAR coefficients (KK = M*p+1)
#   SIGMA  : (ngibbs/nthin) x M  x M     reduced-form covariance L L' + R
#   L      : (ngibbs/nthin) x M  x N     factor loadings
#   F      : (ngibbs/nthin) x T  x N     factors
#   R      : (ngibbs/nthin) x M          idiosyncratic variances
#   irf    : (ngibbs/nthin) x M  x N x nhor   structural IRFs
#   DIC, DIC2 : scalar deviance information criteria

suppressPackageStartupMessages({
  library(MASS)
})

# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------

# build lagged regressor matrix Ylag of dimension T x (M*p) (rows 1..p have NA)
.mlag2 <- function(Y, p) {
  T <- nrow(Y); M <- ncol(Y)
  out <- matrix(NA_real_, T, M * p)
  for (i in seq_len(p)) out[(i + 1):T, ((i - 1) * M + 1):(i * M)] <- Y[1:(T - i), ]
  out
}

# Y, X, dimensions for a VAR(p) with intercept
.prepare_BVAR_matrices <- function(Y, p) {
  T <- nrow(Y); M <- ncol(Y)
  Ylag <- .mlag2(Y, p)[(p + 1):T, , drop = FALSE]
  Yp   <- Y[(p + 1):T, , drop = FALSE]
  X    <- cbind(1, Ylag)
  list(y = Yp, x = X, M = M, T = T - p, KK = M * p + 1)
}

# PCA factor extraction matching extract.m: lam'lam/n = I, fac = data*lam/n
.extract_factors <- function(data, k) {
  n <- ncol(data)
  xx <- crossprod(data)
  ev <- eigen(xx, symmetric = TRUE)
  lam <- sqrt(n) * ev$vectors[, seq_len(k), drop = FALSE]
  fac <- data %*% lam / n
  list(F = fac, L = lam)
}

# set_lambda.m: write the free elements `vecL` of an M x N loading matrix
# whose unit-diagonal elements are fixed to 1
.set_lambda <- function(vecL, n, r) {
  vecL1 <- rep(1, n * r)
  if (r > 1) {
    for (jj in seq_len(r - 1)) {
      idx <- (jj * (n + 1) - (n - 1)):(jj * (n + 1))
      vecL1[idx] <- vecL[(1 + n * (jj - 1)):(jj * n)]
    }
  }
  idx <- (r * (n + 1) - (n - 1)):length(vecL1)
  vecL1[idx] <- vecL[(1 + n * (r - 1)):length(vecL)]
  matrix(vecL1, n, r)
}

# horseshoe.m: one Gibbs sweep updating (beta, lambda, tau) for a single
# regression equation. type = 1 (Bhattacharya et al., p > n) or 2 (Rue, n > p)
.horseshoe <- function(y, X, lambda, tau, sigma_sq, type, ieq) {
  n <- nrow(X); p <- ncol(X)

  # ---- step 1: sample beta ----
  lambda_star <- tau * lambda
  if (type == 1) {
    U <- (lambda_star^2) * t(X)            # p x n
    u <- rnorm(p, 0, lambda_star)
    v <- as.numeric(X %*% u + rnorm(n))
    M <- X %*% U + diag(n)
    v_star <- solve(M, y / sqrt(sigma_sq) - v)
    Beta <- sqrt(sigma_sq) * (u + U %*% v_star)
    Beta <- as.numeric(Beta)
  } else {
    Q_star <- crossprod(X)
    A <- (1 / sigma_sq) * (Q_star + diag(1 / lambda_star^2))
    L_chol <- t(chol(A))                   # lower-triangular
    rhs <- as.numeric(crossprod(X, y)) / sigma_sq
    v_tmp <- forwardsolve(L_chol, rhs)
    mu <- backsolve(t(L_chol), v_tmp)
    u  <- backsolve(t(L_chol), rnorm(p))
    Beta <- mu + u
  }

  # ---- update lambda_j (slice sampler) ----
  eta    <- 1 / lambda^2
  upsi   <- runif(p, 0, 1 / (1 + eta))
  tempps <- Beta^2 / (2 * sigma_sq * tau^2)
  ub     <- (1 - upsi) / upsi
  Fub    <- 1 - exp(-tempps * ub)
  Fub[Fub < 1e-4] <- 1e-4
  up     <- runif(p, 0, Fub)
  eta    <- -log(1 - up) / tempps
  lambda <- 1 / sqrt(eta)
  lambda[1]   <- 1e10  # don't shrink intercept
  lambda[ieq] <- 1e2   # don't shrink own AR(1)

  # ---- update tau (slice sampler) ----
  tempt <- sum((Beta / lambda)^2) / (2 * sigma_sq)
  et    <- 1 / tau^2
  utau  <- runif(1, 0, 1 / (1 + et))
  ubt   <- (1 - utau) / utau
  Fubt  <- pgamma(ubt, shape = (p + 1) / 2, scale = 1 / tempt)
  Fubt  <- max(Fubt, 1e-8)
  ut    <- runif(1, 0, Fubt)
  et    <- qgamma(ut, shape = (p + 1) / 2, scale = 1 / tempt)
  tau   <- 1 / sqrt(et)

  list(Beta = Beta, lambda = lambda, tau = tau)
}

# univariate truncated standard normal on [a, b], inverse-CDF method.
# This is the slot we will replace with `hdtg` in version 2.
.trandn <- function(a, b) {
  pa <- pnorm(a); pb <- pnorm(b)
  u  <- runif(1, pa, pb)
  qnorm(min(max(u, 1e-15), 1 - 1e-15))
}

# ---------------------------------------------------------------------------
# main sampler
# ---------------------------------------------------------------------------

MCMC_BSVAR_Factor <- function(Y, p, sg, nhor = 60,
                              ngibbs = 5000, nburn = 1000, nthin = 1,
                              verbose = TRUE) {

  pm <- .prepare_BVAR_matrices(Y, p)
  y <- pm$y; x <- pm$x; M <- pm$M; T <- pm$T; KK <- pm$KK
  N <- ncol(sg)
  stopifnot(nrow(sg) == M)

  # OLS init
  beta_OLS <- solve(crossprod(x), crossprod(x, y))

  # algorithm choice for horseshoe (Korobilis: type=2 if T>KK, else 1)
  est_alg <- if (nrow(x) > ncol(x)) 2 else 1

  # horseshoe state per equation
  psi      <- replicate(M, rep(1, KK), simplify = FALSE)
  tau      <- rep(1, M)
  sigma_sq <- rep(1, M)

  # prior mean / variance on L (loose, matches Korobilis defaults)
  cL <- c(1, 0.5, 0, -0.5, -0.7, 0, -0.7, -0.6, 1, 0.8, 1.2, 0.7, 0.6, -1.5,
          0, 0, 0, 0, -0.5, -0.5, 0.5, -1, -1, 1, 0.5, 0.5, -0.5, -1.5, -1,
          -3, -2, -0.44, -2, 0.2, -0.5, -0.5, -2, -2.5, -1.5)
  if (M * N < length(cL) + N) {
    Ltemp <- .set_lambda(cL, 14, 3)[1:M, 1:N, drop = FALSE]
    Ltemp[diag(nrow(Ltemp)) == 1] <- NA
    cL <- as.numeric(t(Ltemp)); cL <- cL[!is.na(cL)]
  } else if (M * N > length(cL) + N) {
    Ltemp <- numeric(M * N - N); Ltemp[seq_len(min(39, length(Ltemp)))] <- cL[seq_len(min(39, length(Ltemp)))]
    cL <- Ltemp
  }
  lam_prmean <- .set_lambda(cL, M, N)
  lam_iprvar <- matrix(1 / 4, M, N)
  L_bar      <- matrix(0, M, N)
  index_vec  <- seq_len(N)

  # init factors via PCA on standardized OLS residuals
  yhat   <- y - x %*% beta_OLS
  yhat_z <- scale(yhat)
  fl     <- .extract_factors(yhat_z, N)
  F      <- fl$F
  L      <- fl$L

  # idiosyncratic state
  R   <- diag(M)
  iR  <- rep(1, M)
  beta_mat <- matrix(0, KK, M)

  # storage
  nkeep      <- ngibbs %/% nthin
  beta_save  <- array(0, c(nkeep, KK, M))
  SIGMA_save <- array(0, c(nkeep, M,  M))
  L_save     <- array(0, c(nkeep, M,  N))
  F_save     <- array(0, c(nkeep, T,  N))
  R_save     <- array(0, c(nkeep, M))
  irf_save   <- array(0, c(nkeep, M, N, nhor))
  D_theta    <- numeric(nkeep)
  D_theta2   <- numeric(nkeep)

  counter <- 0L
  total   <- ngibbs + nburn
  t0 <- Sys.time()
  for (iter in seq_len(total)) {
    if (verbose && iter %% 500 == 0)
      message(sprintf("iter %d / %d  (%.1f%%)", iter, total, 100 * iter / total))

    y_til <- y - F %*% t(L)

    # ----- STEP 1: VAR coefficients (one equation at a time) -----
    sample_betas <- function(y_til) {
      for (ieq in seq_len(M)) {
        hs <- .horseshoe(y_til[, ieq], x, psi[[ieq]], tau[ieq],
                         R[ieq, ieq], est_alg, ieq)
        psi[[ieq]] <<- hs$lambda
        tau[ieq]   <<- hs$tau
        beta_mat[, ieq] <<- hs$Beta
      }
    }
    sample_betas(y_til)
    Bcomp <- rbind(t(beta_mat[-1, , drop = FALSE]),
                   cbind(diag(M * (p - 1)), matrix(0, M * (p - 1), M)))

    # stationarity check (cap at 100 retries; fall back to y as in MATLAB)
    rej <- 0L
    while (max(abs(eigen(Bcomp, only.values = TRUE)$values)) > 0.999) {
      rej <- rej + 1L
      if (rej > 100) y_til <- y
      sample_betas(y_til)
      Bcomp <- rbind(t(beta_mat[-1, , drop = FALSE]),
                     cbind(diag(M * (p - 1)), matrix(0, M * (p - 1), M)))
    }

    # ----- STEP 2: factors and loadings -----
    yhat <- y - x %*% beta_mat

    # sample F_t | rest, t = 1..T  (closed-form Normal)
    F_var  <- solve(diag(N) + t(L) %*% (iR * L))
    F_chol <- t(chol(F_var))
    LtiR   <- t(L) * matrix(iR, N, M, byrow = TRUE)   # N x M
    F_mean <- yhat %*% t(LtiR) %*% F_var              # T x N
    F      <- F_mean + matrix(rnorm(T * N), T, N) %*% t(F_chol)
    F      <- sweep(F, 2, colMeans(F), "-")           # demean

    # sample L element-by-element via univariate truncated normal
    # (THE block we will replace with hdtg)
    FtF <- crossprod(F)
    for (i in seq_len(M)) {
      Lvar <- diag(lam_iprvar[i, ]) + FtF * iR[i]     # N x N (precision)
      rhs  <- diag(lam_iprvar[i, ]) %*% lam_prmean[i, ] +
              as.numeric(crossprod(F, yhat[, i])) * iR[i]
      L_bar[i, ] <- as.numeric(solve(Lvar, rhs))
      for (j in seq_len(N)) {
        whole <- Lvar[j, ] * (L[i, ] - L_bar[i, ])
        Lpostvar  <- 1 / Lvar[j, j]
        ss        <- sqrt(Lpostvar)
        Lpostmean <- L_bar[i, j] - sum(whole[index_vec != j]) * Lpostvar
        s_ij <- sg[i, j]
        if (!is.na(s_ij) && s_ij == 1) {
          z <- .trandn((0 - Lpostmean) / ss, (10000 - Lpostmean) / ss)
          L[i, j] <- Lpostmean + ss * z
        } else if (!is.na(s_ij) && s_ij == -1) {
          z <- .trandn((-10000 - Lpostmean) / ss, (0 - Lpostmean) / ss)
          L[i, j] <- Lpostmean + ss * z
        } else if (!is.na(s_ij) && s_ij == 0) {
          L[i, j] <- 0
        } else {
          L[i, j] <- Lpostmean + ss * rnorm(1)
        }
        if (i == j) L[i, j] <- 1
      }
    }

    # ----- STEP 3: idiosyncratic variances -----
    sse2 <- colSums((yhat - F %*% t(L))^2)
    iR   <- rgamma(M, shape = (1 + T) / 2, rate = (0.1 + sse2) / 2)
    R    <- diag(1 / iR)

    SIGMA <- L %*% t(L) + R

    # ----- store -----
    if (iter > nburn && (iter %% nthin) == 0) {
      counter <- counter + 1L
      beta_save[counter, , ]  <- beta_mat
      SIGMA_save[counter, , ] <- SIGMA
      L_save[counter, , ]     <- L
      F_save[counter, , ]     <- F
      R_save[counter, ]       <- 1 / iR

      # IRFs
      Bcomp <- rbind(t(beta_mat[-1, , drop = FALSE]),
                     cbind(diag(M * (p - 1)), matrix(0, M * (p - 1), M)))
      Bpow <- diag(nrow(Bcomp))
      IRF  <- array(0, c(M, N, nhor))
      for (h in seq_len(nhor)) {
        IRF[, , h] <- Bpow[1:M, 1:M] %*% L
        Bpow <- Bpow %*% Bcomp
      }
      irf_save[counter, , , ] <- IRF

      # DIC contributions
      cnst <- -T * M / 2 * log(2 * pi)
      u    <- as.numeric(t(y - x %*% beta_mat - F %*% t(L)))
      D_theta[counter] <-
        cnst - T / 2 * log(det(R)) -
        0.5 * sum(u^2 * rep(iR, T))
      m8 <- min(8, M)
      u8 <- as.numeric(t((y - x %*% beta_mat - F %*% t(L))[, seq_len(m8)]))
      D_theta2[counter] <-
        -T * m8 / 2 * log(2 * pi) -
        T / 2 * log(det(R[seq_len(m8), seq_len(m8)])) -
        0.5 * sum(u8^2 * rep(iR[seq_len(m8)], T))
    }
  }
  if (verbose) message(sprintf("done in %.1fs", as.numeric(Sys.time() - t0, units = "secs")))

  # posterior medians for DIC point estimate
  med <- function(a) apply(a, seq_along(dim(a))[-1], median)
  BB  <- med(beta_save)
  LL  <- med(L_save)
  RRv <- apply(R_save, 2, median); iRR <- 1 / RRv
  FF  <- med(F_save)
  cnst <- -T * M / 2 * log(2 * pi)
  u  <- as.numeric(t(y - x %*% BB - FF %*% t(LL)))
  D_theta_bar  <- cnst - T / 2 * log(prod(RRv)) - 0.5 * sum(u^2 * rep(iRR, T))
  m8 <- min(8, M)
  u8 <- as.numeric(t((y - x %*% BB - FF %*% t(LL))[, seq_len(m8)]))
  D_theta_bar2 <- -T * m8 / 2 * log(2 * pi) -
                  T / 2 * log(prod(RRv[seq_len(m8)])) -
                  0.5 * sum(u8^2 * rep(iRR[seq_len(m8)], T))

  list(
    beta  = beta_save,
    SIGMA = SIGMA_save,
    L     = L_save,
    F     = F_save,
    R     = R_save,
    irf   = irf_save,
    DIC   = -4 * mean(D_theta)  + 2 * D_theta_bar,
    DIC2  = -4 * mean(D_theta2) + 2 * D_theta_bar2
  )
}
