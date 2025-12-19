cat("\014")  
rm(list=ls())

library("hdtg")
library("Matrix")

set.seed(123456)

#################
#### setting ####
#################

# MCMC setting #
MCMC <- 5000
burn <- 1000

# time dimension
Tt   <- 250
time <- 1:Tt

# process noise st. dev.
sigma <- 1 

# states noise st. devs.
ss <- c(1, .2)

# harmonic HMC constraints
Fh <- diag(1, Tt) - cbind(rbind(rep(0, Tt - 1),diag(1, Tt - 1)),rep(0, Tt))
Fh <- Fh[-1, ]

gh <- cbind(rep(Inf,Tt-1),rep(0,Tt-1))

title_plot = c("unconstrained","truncated")

##################
#### simulate ####
##################

eps <- rnorm(Tt, 0, sigma)

x <- log(time)

y <- x + eps

y <- numeric(Tt)
y[1] <- 0

for (i in 2:Tt) {
     
     y[i] <- mu + y[i-1] + mu + eps[i]
}

###########################
#### posterior objects ####
###########################

D <- diag(1, Tt) - cbind(rbind(rep(0, Tt - 1),diag(1, Tt - 1)),rep(0, Tt))

iV <- diag(c(ss[1]^(-2), rep(ss[2]^(-2), Tt - 1))) 

# prior precision
Omega <- t(D) %*% iV %*% D

post_prec <- Omega + diag(rep(sigma^(-2),Tt))
post_mean <- solve(post_prec) %*% (y / sigma^(2))

chol_post_prec <- chol(post_prec)

t <- 1:length(y)        # time trend
mod <- lm(y ~ t)

y_hat <- fitted(mod)

##################
#### estimate ####
##################

for (i in 1:2){
  
     t0 <- Sys.time()
     
     # harmonic HMC
     draws_mu <- harmonicHMC(nSample = MCMC, burnin = burn, 
                             mean = post_mean, choleskyFactor = chol_post_prec, precFlg = TRUE,
                             constrainDirec = Fh, constrainBound = gh[,i], 
                             init = y_hat)  
     
     
     t1 <- Sys.time()
     
     print(t1 - t0)
     
     q05 <- apply(draws_mu, 2, quantile, probs = 0.025)
     q50 <- apply(draws_mu, 2, quantile, probs = 0.50)
     q95 <- apply(draws_mu, 2, quantile, probs = 0.975)
  
     ##############
     #### plot ####
     ##############
     
     # data and true signal
     plot(time, y, type = "l", lwd = 1, xlab = "", ylab = "", main = title_plot[i])
     # lines(time, x_t, col = "red", lwd = 2, lty = 1)
     
     # 95% band shaded
     polygon(c(time, rev(time)), c(q05,  rev(q95)),
             border = NA, col = rgb(0.6, 0.8, 1.0, alpha = 0.4))
     
     # posterior median 
     lines(time, q50, col = "blue", lwd = 3, lty = 2)
     #abline(h = 0, col = "black", lwd = 1, lty = 2)
     
     #legend("topright",legend = c("data","true signal","95% c.i.","median"),
     #       col    = c("black","red",rgb(0.6, 0.8, 1.0, alpha = 0.4),"blue"),
     #       lwd    = c(1, 2, NA, 3),lty    = c(1, 1, NA, 2),pch    = c(NA, NA, 15, NA),
     #       pt.cex = c(NA, NA, 2, NA), bty    = "o")
}
