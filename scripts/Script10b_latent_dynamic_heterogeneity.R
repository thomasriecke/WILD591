library(jagsUI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# number of alligators released per year 
# and number of years of study
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
r <- 25         # number of releases per year
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
gamma <- c(0.3, -0.3, 0.15)
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
boxplot(delta, ylab = 'Growth rate')
plot(w, ylab = 'Available habitat', type = 'b', cex = 2, xlab = '')
mtext('Year', side = 1, line = 1, outer = T, cex = 1.5)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate capture-recapture data (y)
#
# There will be three states
# 1) alive, breeding
# 2) alive, non-breeding
# 3) dead
#
# We'll assume different p's for different states
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# we'll have different detection probabilities for breeders (p[1]) and 
# non-breeders (p[2])
p <- c(0.9, 0.5)
y <- matrix(3, n, t)
for (i in 1:n){
  
  # we'll observe everyone during their release
  if (b[i,f[i]] == 1){
    y[i,f[i]] <- 1
  }
  if (b[i,f[i]] == 0){
    y[i,f[i]] <- 2
  }
  
  for (j in (f[i]+1):t){
    
    # if it's alive and breeding
    if (z[i,j] == 1 & b[i,j] == 1 & !is.na(z[i,j]) & !is.na(b[i,j])){
      tmp <- rbinom(1, 1, p[1])
      if (tmp == 1){
        y[i,j] <- 1
      }
    }
    
    # if it's alive and non-breeding
    if (z[i,j] == 1 & b[i,j] == 0 & !is.na(z[i,j]) & !is.na(b[i,j])){
      tmp <- rbinom(1, 1, p[2])
      if (tmp == 1){
        y[i,j] <- 2
      }
    }
    
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# y will be 'multi-state' capture-recapture data
# 1) breeder
# 2) non-breeder
# 3) not seen
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# another thing to think about. we can only measure size when 
# we recapture these things, and size changes over time, i.e.,
# using the same size covariate for 15 years wouldn't be 
# very intelligent
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knownX <- x
knownX[y == 3] <- NA




sink("ms_model3.jags")
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
  delta ~ dnorm(0, 1)  
  p ~ dbeta(1, 1)
    
  for (j in 1:2){
    sigma[j] ~ dgamma(1, 1)
    tau[j] = pow(sigma[j], -2)
  }



  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # glms for survival and breeding propensity
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (j in 1:(t-1)){
  
    # survivals differ by previous breeding status
    logit(phi[j,1]) = beta[1] + beta[2] + beta[3] * w[j]
    logit(phi[j,2]) = beta[1] + beta[3] * w[j]
    
  }
  
  for (i in 1:n){
    for (j in (f[i]):(t-1)){
      # note that wetness corresponds to next year's wetness for breeding propensity
      logit(psi[i,j]) = alpha[1] + alpha[2] * w[j + 1] + alpha[3] * x[i, j + 1]
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # state-transition and observation matrices
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i in 1:n){
    for (j in f[i]:(t-1)){
      zeta[1,i,j,1] <- phi[j,1] * psi[i,j]
      zeta[1,i,j,2] <- phi[j,1] * (1 - psi[i,j])
      zeta[1,i,j,3] <- 1 - phi[j,1]
    
      zeta[2,i,j,1] <- phi[j,2] * psi[i,j]
      zeta[2,i,j,2] <- phi[j,2] * (1 - psi[i,j])
      zeta[2,i,j,3] <- 1 - phi[j,2]
    
      zeta[3,i,j,1] <- 0
      zeta[3,i,j,2] <- 0
      zeta[3,i,j,3] <- 1
    }
  }


    omega[1,1] <- p
    omega[1,2] <- 0
    omega[1,3] <- 1 - p

    omega[2,1] <- 0
    omega[2,2] <- p
    omega[2,3] <- 1 - p
    
    omega[3,1] <- 0
    omega[3,2] <- 0
    omega[3,3] <- 1    


  for (i in 1:n){
    
    # size at first occasion
    x[i,f[i]] ~ dnorm(0, tau[1])

    # animals transition, grow and are detected following release
    for (j in (f[i] + 1):t){
      z[i,j] ~ dcat(zeta[z[i, j - 1], i, j - 1, ])
      y[i,j] ~ dcat(omega[z[i, j],  ])
      x[i,j] ~ dnorm(x[i,j-1] + delta, tau[2])
    }
  }

}
",fill = TRUE)
sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knownZ <- y
getL <- function(x){max(which(x < 3))}
l <- apply(y, 1, getL)
for (i in 1:n){
  for (j in (f[i] + 1) : l[i]){
    if (knownZ[i,j] == 3){
      knownZ[i,j] <- NA
    }
  }
  if (f[i] > 1){
    knownZ[i, 1:(f[i] - 1)] <- NA
  }
  if (l[i] < t){
    knownZ[i, (l[i]+1):t] <- NA
  }
}

guess <- matrix(NA, n, t)
for (i in 1:n){
  for (j in (f[i] + 1) : l[i]){
    if (is.na(knownZ[i,j])){
      guess[i,j] <- which(rmultinom(1, 1, c(0.8,0.2)) == 1)
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(y = y, n = n, w = w,
           t = t, f = f, x = knownX, z = knownZ)
pars = c('alpha','beta','p','sigma','delta')
inits <- function(){list(z = guess)}  
Sys.time()
# this will take at least a few minutes depending on your machine
m3 <- jags(data = dat, inits=inits, model.file = "ms_model3.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()

print(summary(m3), digits = 3)
m3$mcmc.info$elapsed.mins













# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 4: almost the 'truth'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we'll simulate 'robust' capture-recapture data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y <- array(3, dim = c(n, t, 2))
# we'll have different detection probabilities for breeders (p[1]) and 
# non-breeders (p[2])
p <- c(0.9, 0.5)

for (i in 1:n){
  
  # we'll observe everyone during their release
  if (b[i,f[i]] == 1){
    Y[i,f[i],1:2] <- 1
  }
  if (b[i,f[i]] == 0){
    Y[i,f[i],1:2] <- 2
  }
  
  for (j in (f[i]+1):t){
    
    # if it's alive and breeding
    if (z[i,j] == 1 & b[i,j] == 1 & !is.na(z[i,j]) & !is.na(b[i,j])){
      # first obseration occasion
      tmp <- rbinom(1, 1, p[1])
      if (tmp == 1){
        Y[i,j,1] <- 1
      }
      # second observation occasion
      tmp <- rbinom(1, 1, p[1])
      if (tmp == 1){
        Y[i,j,2] <- 1
      }
    }
    
    # if it's alive and non-breeding
    if (z[i,j] == 1 & b[i,j] == 0 & !is.na(z[i,j]) & !is.na(b[i,j])){
      # first observation occasion
      tmp <- rbinom(1, 1, p[2])
      if (tmp == 1){
        Y[i,j,1] <- 2
      }
      # second observation occasion      
      tmp <- rbinom(1, 1, p[2])
      if (tmp == 1){
        Y[i,j,2] <- 2
      }
    }
    
  }
}

table(Y[,,1], Y[,,2])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# new measurement data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knownX <- x
knownX[Y[,,1] == 3 & Y[,,2] == 3] <- NA



sink("ms_model4.jags")
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
  # lambda[1] = 1 
  # lambda[2] ~ dnorm(0, 0.1) 
  # lambda[3] ~ dnorm(0, 0.1)     
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
      logit(phi[i,j,1]) = beta[1] + beta[2] + beta[3] * w[j]
      logit(phi[i,j,2]) = beta[1] + beta[3] * w[j] + lambda[1] * eta[i]
    
      # note that wetness corresponds to next year's wetness for breeding propensity
      logit(psi[i,j]) = alpha[1] + alpha[2] * w[j + 1] + alpha[3] * x[i, j + 1]
      delta[i,j] = gamma[1] + gamma[2] * w[j + 1]
    }
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # state-transition and observation matrices
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (i in 1:n){
    for (j in f[i]:(t-1)){
      zeta[1,i,j,1] <- phi[i,j,1] * psi[i,j]
      zeta[1,i,j,2] <- phi[i,j,1] * (1 - psi[i,j])
      zeta[1,i,j,3] <- 1 - phi[i,j,1]
    
      zeta[2,i,j,1] <- phi[i,j,2] * psi[i,j]
      zeta[2,i,j,2] <- phi[i,j,2] * (1 - psi[i,j])
      zeta[2,i,j,3] <- 1 - phi[i,j,2]
    
      zeta[3,i,j,1] <- 0
      zeta[3,i,j,2] <- 0
      zeta[3,i,j,3] <- 1
    }
  }


    omega[1,1] <- p[1]
    omega[1,2] <- 0
    omega[1,3] <- 1 - p[1]

    omega[2,1] <- 0
    omega[2,2] <- p[2]
    omega[2,3] <- 1 - p[2]
    
    omega[3,1] <- 0
    omega[3,2] <- 0
    omega[3,3] <- 1    


  for (i in 1:n){
    
    # size at first occasion
    x[i,f[i]] ~ dnorm(0, tau[3])

    # animals transition, grow and are detected following release
    for (j in (f[i] + 1):t){
    
      # state-transition (survival and breeding propensity)    
      z[i,j] ~ dcat(zeta[z[i, j - 1], i, j - 1, ])
      
      # robust observation process (detection)
      Y[i,j,1] ~ dcat(omega[z[i, j],  ])
      Y[i,j,2] ~ dcat(omega[z[i, j],  ])
      
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
getL <- function(x){max(which(x < 3))}
l <- apply(Y, c(1,3), getL)
l <- apply(l, 1, max)
for (i in 1:n){
  knownZ[i,f[i]] <- Y[i,f[i],1]
  for (j in (f[i] + 1) : l[i]){
    
    # known...
    if (Y[i,j,1] == 3 & Y[i,j,2] == 3){
      knownZ[i,j] <- NA
    }
    if (Y[i,j,1] == 1 | Y[i,j,2] == 1){
      knownZ[i,j] <- 1
    }
    if (Y[i,j,1] == 2 | Y[i,j,2] == 2){
      knownZ[i,j] <- 2
    }
    
  }
  if (f[i] > 1){
    knownZ[i, 1:(f[i] - 1)] <- NA
  }
  if (l[i] < t){
    knownZ[i, (l[i]+1):t] <- NA
  }
}

guess <- matrix(NA, n, t)
for (i in 1:n){
  for (j in (f[i] + 1) : l[i]){
    if (is.na(knownZ[i,j])){
      guess[i,j] <- which(rmultinom(1, 1, c(0.8,0.2)) == 1)
    }
  }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(Y = Y, n = n, w = w,
           t = t, f = f, x = knownX, z = knownZ)
pars = c('alpha','beta','gamma','p','sigma','lambda')
inits <- function(){list(z = guess, p = c(0.9,0.5))}  
Sys.time()
# this will take at least a few minutes depending on your machine
m4 <- jags(data = dat, inits=inits, model.file = "ms_model4.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()

print(summary(m4), digits = 3)
m4$mcmc.info$elapsed.mins

plot(m4$sims.list$lambda[,2])
