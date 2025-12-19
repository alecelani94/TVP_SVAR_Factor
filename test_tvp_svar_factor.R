cat("\014")  
rm(list=ls())

library("hdtg")
library("Matrix")
library("coda")

set.seed(123456)

#################
#### setting ####
#################

zigzag <- FALSE # sampler to use

# MCMC setting
MCMC <- 5000
burn <- 1000

# dimensions
N  <- 2
Tt <- 300

K    <- N*Tt # tot. parameters to sample
time <- 1:Tt

sigma_u   <- .05       # measurement error st. dev.
sigma_lam <- c(2, .05) # states st. devs

#########################
#### true parameters ####
#########################

x <- seq(0, pi, length.out = Tt)
lambda_11 <- .25+sin(x)
lambda_22 <- lambda_11

x <- seq(-1*pi, 1.5*pi, length.out = Tt)
lambda_12 <- -(.9+sin(x))/2

x <- seq(-2*pi, .5*pi, length.out = Tt)
lambda_21 <- (.9+sin(x))/2

Lambda <- array(NA, dim = c(2, 2, Tt))

Lambda[1, 1, ] <- lambda_11
Lambda[1, 2, ] <- lambda_12
Lambda[2, 1, ] <- lambda_21
Lambda[2, 2, ] <- lambda_22

###################
#### simulate #####
###################

E  <- matrix(NA, N, Tt)
Ff <- matrix(rnorm(Tt * N), N, Tt)
U  <- matrix(rnorm(Tt * N, mean = 0, sd = sigma_u), N, Tt)
for (t in 1:Tt){
  
     E[,t] <- Lambda[,,t] %*% Ff[,t] + U[,t]
}
E  <- t(E)
Ff <- t(Ff)
U  <- t(U)

###########################
#### posterior objects ####
###########################

# factors and measurement errors given

D <- diag(1, K) - cbind(rbind(matrix(0, N, N*(Tt-1)),diag(1, N*(Tt - 1))),matrix(0, K, N))

D_sparse <- Matrix(D, sparse = TRUE)

browser()

iV <- diag(c(sigma_ss[1]^(-2), rep(sigma_ss[2]^(-2), Tt - 1))) 

# prior precision
Omega <- t(D) %*% iV %*% D

post_prec <- Omega + diag(rep(sigma_e^(-2),Tt))
post_mean <- solve(post_prec) %*% (y_t / sigma_e^(2))

chol_post_prec <- chol(post_prec)

##################
#### estimate ####
##################

for (i in 1:2){
  
  t0 <- Sys.time()
  
  if (zigzag){
    
    # zig-zag HMC
    draws_mu <- zigzagHMC(nSample = MCMC, burnin = burn, 
                          mean = post_mean, prec = post_prec, 
                          lowerBounds = lb[,i], upperBounds = ub) 
  } else {
    
    # harmonic HMC
    draws_mu <- harmonicHMC(nSample = MCMC, burnin = burn, 
                            mean = post_mean, choleskyFactor = chol_post_prec, precFlg = TRUE,
                            constrainDirec = Fh, constrainBound = gh[,i], 
                            init = rep(mean(y_t), Tt))  
  }
  
  
  t1 <- Sys.time()
  
  print(t1 - t0)
  
  q05 <- apply(draws_mu, 2, quantile, probs = 0.025)
  q50 <- apply(draws_mu, 2, quantile, probs = 0.50)
  q95 <- apply(draws_mu, 2, quantile, probs = 0.975)
  
  ##############
  #### plot ####
  ##############
  
  # data and true signal
  plot(time, y_t, type = "l", lwd = 1, xlab = "", ylab = "", main = title_plot[i])
  lines(time, x_t, col = "red", lwd = 2, lty = 1)
  
  # 95% band shaded
  polygon(c(time, rev(time)), c(q05,  rev(q95)),
          border = NA, col = rgb(0.6, 0.8, 1.0, alpha = 0.4))
  
  # posterior median 
  lines(time, q50, col = "blue", lwd = 3, lty = 2)
  abline(h = 0, col = "black", lwd = 1, lty = 2)
  
  legend("topright",legend = c("data","true signal","95% c.i.","median"),
         col    = c("black","red",rgb(0.6, 0.8, 1.0, alpha = 0.4),"blue"),
         lwd    = c(1, 2, NA, 3),lty    = c(1, 1, NA, 2),pch    = c(NA, NA, 15, NA),
         pt.cex = c(NA, NA, 2, NA), bty    = "o")
}


par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))  # margins

plot(lambda_11, type = "l", col = "red", lty = 2, lwd = 3,
     main = expression(lambda[11]), xlab = "", ylab = "")

plot(lambda_12, type = "l", col = "red", lty = 2, lwd = 3,
     main = expression(lambda[12]), xlab = "", ylab = "")
abline(h = 0, col = "black", lwd = 1, lty = 2)

plot(lambda_21, type = "l", col = "red", lty = 2, lwd = 3,
     main = expression(lambda[21]), xlab = "", ylab = "")
abline(h = 0, col = "black", lwd = 1, lty = 2)

plot(lambda_22, type = "l", col = "red", lty = 2, lwd = 3,
     main = expression(lambda[22]), xlab = "", ylab = "")

