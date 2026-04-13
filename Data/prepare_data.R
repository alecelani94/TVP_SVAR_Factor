# prepare_data.R
#
# R port of prepare_data.m: load and transform the 14-variable monthly
# dataset of Korobilis (2022). Bit-for-bit faithful to realVARDGP.m's
# in-line preprocessing pipeline:
#
#   transx (per tcode)  ->  adjout(., 4.5, 4)  ->  drop row 1  ->  zscore
#
# Returns a list with:
#   Data  : T x 14 numeric matrix, ready to feed into a VAR routine
#   dates : T-length Date vector aligned with Data
#   names : length-14 character vector of variable names
#   tcode : length-14 integer vector of transformation codes
#
# Usage:
#   source("Data/prepare_data.R")
#   d <- prepare_data()        # default path: Data/data_for_MC.xlsx
#   Y <- d$Data

suppressPackageStartupMessages(library(readxl))

# ---- transx: Stock & Watson stationarity transformations -------------------
# Mirrors transx.m. Only the codes actually present in data_for_MC.xlsx
# (1 = level, 5 = log first difference) are exercised, but we implement the
# full set used by Korobilis for completeness.
transx <- function(x, tcode) {
  small <- 1e-40
  n <- length(x)
  y <- rep(NA_real_, n)
  needs_log <- tcode %in% c(4, 5, 6, 16, 17)
  if (needs_log && min(x, na.rm = TRUE) < small) return(rep(NA_real_, n))

  switch(as.character(tcode),
    "1"  = { y <- x },
    "2"  = { y[2:n] <- x[2:n] - x[1:(n - 1)] },
    "3"  = { y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)] },
    "4"  = { y <- log(x) },
    "5"  = { lx <- log(x); y[2:n] <- lx[2:n] - lx[1:(n - 1)] },
    "6"  = { lx <- log(x); y[3:n] <- lx[3:n] - 2 * lx[2:(n - 1)] + lx[1:(n - 2)] },
    "16" = { lx <- log(x); y[3:n] <- lx[3:n] - 2 * lx[2:(n - 1)] + lx[1:(n - 2)] },
    "17" = { lx <- log(x)
             y[14:n] <- lx[14:n] - lx[13:(n - 1)] - lx[2:(n - 12)] + lx[1:(n - 13)] },
    stop(sprintf("transx: tcode %s not supported in this port", tcode))
  )
  y
}

# ---- adjout: outlier adjustment, tflag = 4 (one-sided rolling median) ------
# Mirrors adjout.m for the (thr = 4.5, tflag = 4) call used by Korobilis.
adjout <- function(y, thr = 4.5, tflag = 4) {
  stopifnot(tflag == 4)               # only branch we need
  small <- 1e-6
  n <- length(y)
  z <- y[!is.na(y)]
  q <- quantile(z, c(0.25, 0.50, 0.75), names = FALSE, type = 7)
  zm  <- q[2]
  iqr <- q[3] - q[1]
  if (iqr < small) return(rep(NA_real_, n))

  ya  <- abs(y - zm)
  iya <- !is.na(ya) & ya >  thr * iqr  # outlier flag (NA -> FALSE, like MATLAB)
  iyb <- !is.na(ya) & ya <= thr * iqr

  # one-sided rolling median over the 5 preceding obs (window of 6 incl. self)
  iwin <- 5
  ymvec <- vapply(seq_len(n), function(i) {
    j1 <- max(1L, i - iwin)
    median(y[j1:i], na.rm = TRUE)
  }, numeric(1))

  out <- y
  out[iya] <- ymvec[iya]
  out[iyb] <- y[iyb]
  out
}

# ---- main entry point ------------------------------------------------------
prepare_data <- function(xlsx_file = NULL) {
  if (is.null(xlsx_file)) {
    xlsx_file <- file.path(
      dirname(sys.frame(1)$ofile %||% "Data"),  # script-relative if sourced
      "data_for_MC.xlsx"
    )
    if (!file.exists(xlsx_file)) xlsx_file <- "Data/data_for_MC.xlsx"
  }

  raw <- suppressMessages(read_excel(
    xlsx_file, sheet = "data", col_names = FALSE, range = "A1:O518"
  ))
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  vnames <- as.character(unlist(raw[1, 2:15]))
  tcode  <- as.integer(unlist(raw[2, 2:15]))
  data_mat <- as.matrix(sapply(raw[3:nrow(raw), 2:15], as.numeric))
  colnames(data_mat) <- vnames

  date_serials <- as.numeric(unlist(raw[3:nrow(raw), 1]))
  dates_full   <- as.Date(date_serials, origin = "1899-12-30")

  n <- ncol(data_mat)
  Dt <- matrix(NA_real_, nrow = nrow(data_mat), ncol = n)
  for (i in seq_len(n)) {
    Dt[, i] <- transx(data_mat[, i], tcode[i])
    Dt[, i] <- adjout(Dt[, i], thr = 4.5, tflag = 4)
  }

  Data  <- scale(Dt[-1, , drop = FALSE])           # zscore: mean 0, sd 1
  attr(Data, "scaled:center") <- NULL
  attr(Data, "scaled:scale")  <- NULL
  colnames(Data) <- vnames
  dates <- dates_full[-1]

  list(Data = Data, dates = dates, names = vnames, tcode = tcode)
}

# small null-coalescing helper used above
`%||%` <- function(a, b) if (is.null(a)) b else a
