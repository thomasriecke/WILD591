# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(latex2exp)
library(jagsUI)
library(MuMIn)
library(loo)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample size (n)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 100
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# covariate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x <- rnorm(n, 0, 1)
hist(x, main = NULL, xlab = "x")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate some data (y)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y <- rnorm(n, x, 1)
plot(y ~ x, pch = 21, bg = 'dodgerblue', las = 1, cex.lab = 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here's a standard lm() model with the intercept fixed to 0
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m0F <- lm(y ~ 1))
summary(m1F <- lm(y ~ 0 + x))
AICc(m0F)
AICc(m1F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("lm.jags")
cat("
      model {

      beta ~ dnorm(0, 1)
      sigma ~ dgamma(1,1)
      tau = 1/(sigma * sigma)
      
      for (i in 1:n){
      
        # model
        y[i] ~ dnorm(beta * x[i], tau)
        
        # calculate log-likelihood
        logL[i] <- log(dnorm(y[i], beta * x[i], tau))
        
        # generate new data for PPC
        new[i] ~ dnorm(beta * x[i], tau)
        
      }

      
      }
      ",fill = TRUE)
sink()

jags.data <- list(y = y, x = x, n = n)
inits <- function(){list()}  
parameters <- c('sigma','beta','logL','new')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "lm.jags", 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
          parallel = T)
Sys.time()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# ok, let's derive a few things here...
# we've thoughtfully saved our parameter (k = 2) estimates (beta and sigma)
# as well as the log-likelihood for each data point
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First let's see if we can derive the deviance (which JAGS saves for us)
# using the data and our parameter estimates, and calculate DIC
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of iterations, likelihood, and expected value
n.iter <- nc * (ni-nb)/nt
lik <- NULL
E <- NULL
for (i in 1:n.iter){
  
  # here E is just the expected value of each data point at each time step
  # we don't use it here, but could plug it in at the end of the next line of code
  # e.g., sum((y-E)^2) rather than sum((y-m$sims.list$beta[i]*x)^2)
  E <- m$sims.list$beta[i] * x
  
  # this is the log-likelihood function for normal linear regression
  # it's in many places on the internet (sometimes with typos, watch out!)
  lik[i] <- -((n/2) * log(2*pi)) - ((n/2) * log(m$sims.list$sigma[i]^2)) -
            (1/(2 * m$sims.list$sigma[i]^2)) * sum((y-m$sims.list$beta[i]*x)^2)
}
hist(lik)

# What JAGS calls 'Deviance' equals -2 * the log likelihood (lik)
D <- -2 * lik
plot(m$sims.list$deviance ~ D) # nice!

# the number of effective parameters (as measured for DIC) = the variance of the Deviance/2
pD <- var(D)/2
m$pD

# JAGS calculates DIC as the mean of the Deviance + the number of effective parameters (pD)
DIC <- mean(D) + pD
DIC
m$DIC
# well that's handy dandy! i'm assuming the 0.07 difference is a result of rounding error somewhere
# and can live with it, but it may be due to some tiny mistake in the code...




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WAIC
# Gelman, Hwang, and Vehtari (2013) Statistics and Computing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here rather than a log-likelihood for each iteration (see above), we'll save
# the log-likelihood for each data point at each iteration!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.iter <- nc * (ni-nb)/nt

# note the change in the dimensions of 'lik' from above
lik <- matrix(NA, n.iter, n)
E <- NULL

for (i in 1:n.iter){
  E <- m$sims.list$beta[i] * x
  for (j in 1:n){
    # note that we substitute 1 for n from the equation above, as 
    # this is for each single data point...
    lik[i,j] <- -((1/2) * log(2*pi)) - ((1/2) * log(m$sims.list$sigma[i]^2)) -
      (1/(2 * m$sims.list$sigma[i]^2)) * (y[j]-E[j])^2
  }
}

hist(rowSums(lik))
# cool, now we have a log-likelihood for each data point!
plot(lik[,1] ~ m$sims.list$logL[,1])

# now we can calculate WAIC!
waic(lik)
waic(m$sims.list$logL)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Below is a method to calculate WAIC by hand, this is from Gelman et al.
#
# See here for more details:
#
# Gelman, Hwang, and Vehtari (2013)
# http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf
# 
# Vehtari, Gelman and Gabry (2016)
# http://www.stat.columbia.edu/~gelman/research/unpublished/loo_stan.pdf
#
# Note that if you can define your model in Stan (and brms can often help
# you do that quite easily!), Stan will auto-calculate WAIC for you
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lppd <- sum(log(colMeans(exp(lik))))
p.waic <- sum(apply(lik, 2, var))
waic <- -2*lppd + 2*p.waic
waic; p.waic
waic(m$sims.list$logL)
# ok... that was nice! let's do it a bunch of times, eh?

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ok, finally let's explore posterior predictive checks via
# Bayesian p-values (an assessment of model fit)
# Here we'll compare the residual sum of squares of the real data
# to the RSS of the generated data
# Here we'll pull directly from chapter 8 in Kery (2010) Academic Press
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.iter <- nc * (ni-nb)/nt
fit.data <- NULL
fit.newd <- NULL
for (i in 1:n.iter){
  fit.data[i] <- sum((y - (m$sims.list$beta[i] * x))^2)
  fit.newd[i] <- sum((m$sims.list$new[i,] - (m$sims.list$beta[i] * x))^2)  
}

Bp <- mean(fit.newd > fit.data)

smoothScatter(fit.newd ~ fit.data, ylim = c(0,250), xlim = c(0,250),
     ylab = 'Residual sum of squares (simulated data)', xlab = 'Residual sum of squares (actual data)')
points(fit.newd ~ fit.data, cex = 0.25, col = 'navy')
abline(0,1,col = 'grey20', lty = 2)
Bp














# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# [THIS WILL TAKE 30+ minutes even on a fast computer, don't run in class]
# a run with 100 simulations will take > 3minutes, you might try that by
# adjusting the simulations run below...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 
# Note that, we've already established that a lot of these things are accurately
# calculated in JAGS (and via loo, etc., not that we needed too, JAGS and Stan
# are extremely reliable in my experience), so we're going to cut a lot of the comments
# and run two models:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sims <- 1000

# sample size and vectors to save results of interest
n <- 100

# mean and f-value (proportion of posterior that is on same side of zero as mean)
# from Bayes model as well as mean and p-value for covariate effect
muBeta <- NULL
fBeta <- NULL
lmBeta <- NULL
pBeta <- NULL

# model ICs (different information criteria scores)
aicc <- NULL
waic <- NULL
dic <- NULL

# effective number of parameters
pd <- NULL
pwaic <- NULL

# Bayesian p-value
Bayesp <- NULL

# our JAGS model
sink("lm.jags")
cat("
      model {

      beta ~ dnorm(0, 1)
      sigma ~ dgamma(1,1)
      tau = 1/(sigma * sigma)
      
      for (i in 1:n){
        y[i] ~ dnorm(beta * x[i], tau)
        logL[i] <- log(dnorm(y[i], beta * x[i], tau))
        E[i] ~ dnorm(beta * x[i], tau)
      }
      
      
      
      
      
      }
      ",fill = TRUE)
sink()

# simulation parameters
nc <- 4
nt <- 10
ni <- 25000
nb <- 10000
n.iter <- nc * (ni-nb)/nt

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# run simulation
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (ii in 1:sims){
  
  # help forecast total runtime
  print(Sys.time())
  print(ii)
  
  # simulate a covariate (x) and response (y)
  x <- rnorm(n, 0, 1)
  y <- rnorm(n, x, 1)
  
  summary(mF <- lm(y ~ 0 + x))
  aicc[ii] <- AICc(mF)
  lmBeta[ii] <- mF$coefficients[1]
  pBeta[ii] <- summary(mF)$coefficients[1,4]
  
  jags.data <- list(y = y, x = x, n = n)
  inits <- function(){list()}  
  parameters <- c('sigma','beta','logL','E')
  m <- jags(jags.data, inits, parameters, "lm.jags", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
            parallel = T)
  
  dic[ii] <- m$DIC
  muBeta[ii] <- m$q50$beta
  fBeta[ii] <- m$f$beta
  pd[ii] <- m$pD
  
  
  lppd <- sum(log(colMeans(exp(m$sims.list$logL))))
  pwaic[ii] <- sum(apply(m$sims.list$logL, 2, var))
  waic[ii] <- -2*lppd + 2*pwaic[ii]
  
  
  fit.data <- NULL
  fit.newd <- NULL
  for (i in 1:n.iter){
    fit.data[i] <- sum((y - (m$sims.list$beta[i] * x))^2)
    fit.newd[i] <- sum((m$sims.list$E[i,] - (m$sims.list$beta[i] * x))^2)  
  }
  
  Bayesp[ii] <- mean(fit.newd > fit.data)
  
 
}



hist(Bayesp, breaks = 100, main = '', xlab = 'Bayes p-values')
smoothScatter(dic~waic, ylab = 'DIC', xlab = 'WAIC')
points(dic ~ waic, cex = 0.25, col = 'navy')

smoothScatter(pD~p.waic, ylab = TeX("pD"), xlab = TeX("p[WAIC]"))
points(pD~p.waic, cex = 0.25, col = 'navy')

smoothScatter(waic ~ aicc, xlab = 'AICc', ylab = 'WAIC')
points(waic~aicc, cex = 0.25, col = 'navy')

smoothScatter(dic ~ aicc, xlab = 'AICc', ylab = 'DIC')
points(dic ~ aicc, cex = 0.25, col = 'navy')

library(vioplot)
smoothScatter(lmBeta ~ muBeta, xlab = TeX("$\\beta_\\JAGS$"), ylab = TeX("$\\beta_\\glm$"))
points(lmBeta ~ muBeta, cex = 0.25, col = 'navy')

