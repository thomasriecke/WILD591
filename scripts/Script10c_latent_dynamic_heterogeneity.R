library(jagsUI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of alligators released per year 
# and number of years of study
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r <- 50         # number of releases per year
t <- 10         # number of years
n <- r * (t-1)  # total number of individuals released

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate variation in landscape wetness
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
w <- rnorm(t, 0, 1)
par(mfrow = c(1,1), oma = rep(0,4), mar = c(5,5,2,2))
plot(w, xlab = 'Year', ylab = 'Wetness (w)', type = 'b', cex = 2, cex.lab = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate individual variation...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eta <- rnorm(n, 0, 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll simulate three latent processes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z: a matrix indicating whether or not an individual
#    was alive
# b: a matrix indicating whether or not an individual was
#    breeding
# x: a matrix containing the 'size' of each individual during
#    each time step. This will be zero-centered (at release)
#    it might be easier to think of this as the amount of nutrient
#    reserves available for breeding
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- matrix(NA, n, t)
b <- matrix(NA, n, t)
x <- matrix(NA, n, t)

phi <- matrix(NA, n, t)     # survival probability matrix
psi <- matrix(NA, n, t)     # breeding probability matrix
delta <- matrix(NA, n, t)   # growth matrix

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(2, -0.7, 0.5)
beta <- c(1, 0.25, 0.25)
gamma <- c(0.5, -0.5, 0.15)
lambda <- c(0.25, 0.25, 0.05)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# entry into the marked population (first encounter; f)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f <- rep(1:(t-1), each = r)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's simulate some alligator life-histories!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in 1:n){
  
  # release new individuals and give them a random 'size'
  z[i,f[i]] <- 1
  x[i,f[i]] <- rnorm(1, 0, 1)
  
  for (j in (f[i]):(t-1)){
    # breeding (if alive)
    if (z[i, j] == 1 & !is.na(z[i, j] == 1)){
      psi[i, j] <- plogis(beta[1] + beta[2] * x[i, j] + beta[3] * w[j] + lambda[2] * eta[i])
      b[i, j] <- rbinom(1, z[i, j], psi[i, j])
    }

    # survival
    if (z[i, j] == 1 & !is.na(z[i, j] == 1)){
      phi[i, j] <- plogis(alpha[1] + alpha[2] * b[i, j] + alpha[3] * w[j] + lambda[1] * eta[i])
      z[i, j + 1] <- rbinom(1, z[i, j], phi[i, j])
    }

    # growth (if alive)
    if (z[i, j + 1] == 1 & !is.na(z[i, j + 1])){
      delta[i, j] <- gamma[1] + gamma[2] * b[i, j] + gamma[3] * w[j] + lambda[3] * eta[i]
      x[i, j + 1] <- rnorm(1, x[i, j] + delta[i, j], 0.1)
    }

  }
  
  # if alive in the last year, did it breed?
  if (z[i, t] == 1 & !is.na(z[i, t] == 1)){
    psi[i, t] <- plogis(beta[1] + beta[2] * x[i,j] + beta[3] * w[t])
    b[i, t] <- rbinom(1, 1, psi[i, t])
  }
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# figure for size
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplot(x[1:r, 1], boxwex = 0.25, col = 'lightgreen',
        las = 1, ylab = 'Condition (x)', cex.lab = 2)
points(x[1:r, 1] ~ jitter(rep(0.8, r), 1), col = 'darkgreen', cex = 2)
vioplot::vioplot(x[1:r, 1], at = 1.2, side = 'right', add = T, 
                 col = 'lightgreen', wex = 0.25)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot of breeding propensity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(psi ~ x, cex.lab = 2, las = 1, xlab = 'Condition (x)',
     col = 'darkgreen',
     ylab = latex2exp::TeX("Breeding propensity ($\\psi$)"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot of growth
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(delta ~ jitter(b,0.1), cex.lab = 1.5, las = 1,
     col = 'darkgreen',
     ylab = latex2exp::TeX("Growth from t to t+1 ($\\delta$)"),
     xlab = "Nested in t-1")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot of survival
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(phi ~ jitter(b,0.1), cex.lab = 1.5, las = 1,
     col = 'darkgreen',
     ylab = latex2exp::TeX("Survival from t to t+1 ($\\delta$)"),
     xlab = "Nested in t-1")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ok... let's see some alligators change condition
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1))
plot(x[1,], type = 'b', ylim = c(-5,5))
for (j in 2:r){
  points(x[j,], type = 'b') # fun!
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Make a quick plot of pop-level parameter values...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(2,2), mar = c(2,5,2,2), oma = c(3,0,0,0))
boxplot(phi, ylab = 'Survival')
boxplot(psi, ylab = 'Breeding propensity')
boxplot(delta, ylab = 'Change in condition')
plot(w, ylab = 'Available habitat', type = 'b', cex = 2, xlab = '')
mtext('Year', side = 1, line = 1, outer = T, cex = 1.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate capture-recapture data (y)
#
# This will be a little different than the other examples
# note that our three states:
# 1) alive, breeding
# 2) alive, non-breeding
# 3) dead
# could also be resulted as the outcome of two sequential Bernoulli trials
# 1) are you alive? yes (1) or no (0)
# 2) if you're alive, are you breeding? yes (1) or no (0)
#
# We'll assume different p's for different states
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# we'll have different detection probabilities for breeders (p[1]) and 
# non-breeders (p[2])
# Y will be our detection/non-detection data
# knownB will be our observations of breeding effort
p <- c(0.9, 0.5)
Y <- array(0, dim = c(n, t, 2))
for (i in 1:n){
  
  # we'll observe everyone during their release
  Y[i,f[i],1] <- Y[i,f[i],2] <- 1
  
  for (j in (f[i]+1):t){
    
    # if it's alive and breeding
    if (z[i,j] == 1 & b[i,j] == 1 & !is.na(z[i,j]) & !is.na(b[i,j])){
      tmp <- rbinom(1, 1, p[1])
      if (tmp == 1){
        Y[i,j,1] <- 1
      }
      tmp <- rbinom(1, 1, p[1])
      if (tmp == 1){
        Y[i,j,2] <- 1
      }
    }
    
    # if it's alive and non-breeding
    if (z[i,j] == 1 & b[i,j] == 0 & !is.na(z[i,j]) & !is.na(b[i,j])){
      tmp <- rbinom(1, 1, p[2])
      if (tmp == 1){
        Y[i,j,1] <- 1
      }
      tmp <- rbinom(1, 1, p[2])
      if (tmp == 1){
        Y[i,j,2] <- 1
      }
    }
    
  }
}

knownB <- b
# sum number of encounters per occasion across secondary occasions
e <- apply(Y, c(1,2), sum)
# if not seen, then breeding status is unknown
knownB[e == 0] <- NA

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# y will be 'multi-state' capture-recapture data
# 1) breeder
# 2) non-breeder
# 3) not seen
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# another thing to think about. we can only measure condition when 
# we recapture these things, and size changes over time, i.e.,
# using the same size covariate for 15 years wouldn't be 
# very intelligent
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knownX <- x
knownX[e == 0] <- NA



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 5: almost the 'truth'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(Y[,,1], Y[,,2])
table(Y[,,1], knownB)
table(Y[,,2], knownB)
table(e, knownB)


sink("ms_model5.jags")
cat("
model{


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # priors
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  alpha[1] ~ dlogis(0,1)
  alpha[2] ~ dnorm(0, 0.1)  
  alpha[3] ~ dnorm(0, 0.1)  
  beta[1] ~ dlogis(0,1)
  beta[2] ~ dnorm(0, 0.1) 
  beta[3] ~ dnorm(0, 0.1) 
  gamma[1] ~ dnorm(0, 0.1) 
  gamma[2] ~ dnorm(0, 0.1) 
  gamma[3] ~ dnorm(0, 0.1)  
  lambda[1] = 1   
  lambda[2] ~ dnorm(0, 0.1) 
  lambda[3] ~ dnorm(0, 0.1)    
  p[1] ~ dbeta(1, 1)
  p[2] ~ dbeta(1,1)
    
  for (j in 1:3){
    sigma[j] ~ dgamma(1, 1)
    tau[j] = pow(sigma[j], -2)
  }

  for (i in 1:n){
    eta.star[i] ~ dnorm(0, 1)
    eta[i] = eta.star[i] * sigma[1]
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # glms for survival and breeding propensity
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  for (i in 1:n){
    for (j in (f[i]):(t-1)){
    
      # survivals differ by previous breeding status
      logit(phi[i,j]) = beta[1] + 
                        beta[2] * b[i,j] + 
                        beta[3] * w[j] + 
                        lambda[1] * eta[i]

      # note that wetness corresponds to next year's wetness for breeding propensity
      logit(psi[i,j]) = alpha[1] + 
                        alpha[2] * x[i, j + 1] + 
                        alpha[3] * w[j + 1] +
                        lambda[2] * eta[i]
      
      # here growth will be a function of wetness and individual het.                  
      delta[i,j] = gamma[1] + 
                   gamma[2] * b[i, j]  +
                   gamma[3] * w[j + 1] +
                   lambda[3] * eta[i]
    }
  }


  for (i in 1:n){
    
    # size at first occasion
    x[i,f[i]] ~ dnorm(0, tau[3])

    # animals transition, grow and are detected following release
    for (j in (f[i] + 1):t){
    
      # state-transition (survival and breeding propensity)    
      z[i,j] ~ dbern(z[i,j-1] * phi[i,j-1])
      b[i,j] ~ dbern(z[i,j] * psi[i,j-1])
      
      # robust observation process (detection)
      # note that detection is conditional on breeding, so that
      # if b[i,j] = 1 then we use p[2], otherwise p[1]
      Y[i,j,1] ~ dbern(z[i,j] * p[(b[i,j] + 1)])
      Y[i,j,2] ~ dbern(z[i,j] * p[(b[i,j] + 1)])
      
      # growth process
      x[i,j] ~ dnorm(x[i,j-1] + delta[i,j-1], tau[2])
    }
  }

}
",fill = TRUE)
sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knownZ <- matrix(NA, n, t)
getL <- function(x){max(which(x > 0))}
l <- apply(e, 1, getL)

for (i in 1:n){
  knownZ[i,f[i]:l[i]] <- 1
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(Y = Y, n = n, w = w, b = knownB,
           t = t, f = f, x = knownX, z = knownZ)
pars = c('alpha','beta','gamma','p','sigma','lambda')
inits <- function(){list(p = c(0.5,0.9))}  
Sys.time()
# this will take at least a few minutes depending on your machine
m5 <- jags(data = dat, inits=inits, model.file = "ms_model5.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()

print(summary(m5), digits = 3)
m5$mcmc.info$elapsed.mins



plot(m5$sims.list$alpha[,1])
plot(m5$sims.list$alpha[,2])
plot(m5$sims.list$alpha[,3])

plot(m5$sims.list$beta[,1])
plot(m5$sims.list$beta[,2])
plot(m5$sims.list$beta[,3])

plot(m5$sims.list$gamma[,1])
plot(m5$sims.list$gamma[,2])
plot(m5$sims.list$gamma[,3])

plot(m5$sims.list$sigma[,1])
plot(m5$sims.list$sigma[,2])
plot(m5$sims.list$sigma[,3])