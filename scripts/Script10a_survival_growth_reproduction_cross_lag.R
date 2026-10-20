
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
# now we'll simulate three latent processes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# z: a matrix indicating whether or not an individual
#    was alive
# b: a matrix indicating whether or not an individual was
#    breeding
# x: a matrix containing the 'body condition' of each individual during
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
      psi[i, j] <- plogis(beta[1] + beta[2] * x[i, j] + beta[3] * w[j])
      b[i, j] <- rbinom(1, z[i, j], psi[i, j])
    }

    # survival
    if (z[i, j] == 1 & !is.na(z[i, j] == 1)){
      phi[i, j] <- plogis(alpha[1] + alpha[2] * b[i, j] + alpha[3] * w[j])
      z[i, j + 1] <- rbinom(1, z[i, j], phi[i, j])
    }

    # growth (if alive)
    if (z[i, j + 1] == 1 & !is.na(z[i, j + 1])){
      delta[i, j] <- gamma[1] + gamma[2] * b[i, j] + gamma[3] * w[j]
      x[i, j + 1] <- rnorm(1, x[i, j] + delta[i, j], 0.1)
    }

  }
  
  # if alive in the last year of study, did it breed?
  if (z[i, t] == 1 & !is.na(z[i, t] == 1)){
    psi[i, t] <- plogis(beta[1] + beta[2] * x[i,j] + beta[3] * w[t])
    b[i, t] <- rbinom(1, 1, psi[i, t])
  }
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# figure for size
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), oma = rep(0,4), mar = c(5,5,2,2))
boxplot(x[1:r, 1], boxwex = 0.25, col = 'lightgreen',
        las = 1, ylab = 'Body condition (x)', cex.lab = 2)
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
# ok... let's see some alligators grow
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
# m will be simpler capture-recapture data
# 1) seen
# 2) not seen
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m <- y
m[m == 3] <- 0
m[m != 0] <- 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# another thing to think about. we can only measure size when 
# we recapture these things, and size changes over time, i.e.,
# using the same size covariate for 15 years wouldn't be 
# very intelligent
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
knownX <- x
knownX[y == 3] <- NA


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# specify that individuals were alive b/w first and last sighting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getL <- function(x){max(which(x == 1))}
l <- apply(m, 1, getL)
knownZ <- matrix(NA, n, t)
for (i in 1:n){
  knownZ[i,f[i]:l[i]] <- 1
}


sink("cjs_model1.jags")
cat("
model{

  p ~ dbeta(1,1)
  alpha[1] ~ dlogis(0,1)
  alpha[2] ~ dnorm(0, 0.1)

  for (j in 1:(t-1)){
    logit(phi[j]) = alpha[1] + alpha[2] * w[j]
  }
  
  for (i in 1:n){
    # animals survive and are detected following release
    for (j in (f[i] + 1):t){
      z[i,j] ~ dbern(z[i, j-1] * phi[j-1])
      m[i,j] ~ dbern(z[i,j] * p)
    }
  }

}
",fill = TRUE)
sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(m = m, n = n, w = w,
           t = t, f = f, x = knownX,
           z = knownZ)
pars = c('alpha','p')
inits <- function(){list()}  
Sys.time()
m1 <- jags(data = dat, inits=inits, model.file = "cjs_model1.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()

print(summary(m1), digits = 3)
m1$mcmc.info$elapsed.mins


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract and plot median, 80%, and 95% credible intervals for expected
# values of phi given wetness
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iter <- m1$mcmc.info$n.samples
res <- 50
predw <- seq(min(w), max(w), length.out = res)
predS <- matrix(NA, iter, res)
qS <- matrix(NA, 5, res)
for (j in 1:res){
  predS[,j] <- plogis(m1$sims.list$alpha[,1] + m1$sims.list$alpha[,2] * predw[j])
  qS[,j] <- quantile(predS[,j], c(0.025,0.1,0.5,0.9,0.975))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2), oma = rep(0,4))
plot(qS[3,] ~ predw, ylim = c(0.5,1), type = 'l', lwd = 3,
     xlab = 'Wetness (w)',
     ylab = latex2exp::TeX("Survival ($\\phi$)"))
lines(qS[1,] ~ predw, lty = 2, lwd = 2)
lines(qS[5,] ~ predw, lty = 2, lwd = 2)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# add the actual values, note that they're replicated for many individuals
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (i in 1:n){
  points(phi[i,] ~ w, col = 'red')
}




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 2:
# a multi-state model
# now we'll use the data that are formatted as 
# 1) breeder
# 2) non-breeder
# 3) not seen
# We'll also explicitly test a 'cost-of-breeding' hypothesis on survival,
# albeit at the population, not indivdual level.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("ms_model2.jags")
cat("
model{


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # priors
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  p ~ dbeta(1,1)
  alpha[1] ~ dlogis(0,1)
  alpha[2] ~ dnorm(0, 0.1)  
  beta[1] ~ dlogis(0,1)
  beta[2] ~ dnorm(0, 0.1) 
  beta[3] ~ dnorm(0, 0.1) 

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # glms for survival and breeding propensity
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (j in 1:(t-1)){
  
    # survivals differ by previous breeding status
    logit(phi[j,1]) = beta[1] + beta[2] + beta[3] * w[j]
    logit(phi[j,2]) = beta[1] + beta[3] * w[j]
    
    # note that wetness corresponds to next year's wetness for breeding propensity
    logit(psi[j]) = alpha[1] + alpha[2] * w[j + 1]
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # state-transition and observation matrices
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for (j in 1:(t-1)){
    zeta[1,j,1] <- phi[j,1] * psi[j]
    zeta[1,j,2] <- phi[j,1] * (1 - psi[j])
    zeta[1,j,3] <- 1 - phi[j,1]
    
    zeta[2,j,1] <- phi[j,2] * psi[j]
    zeta[2,j,2] <- phi[j,2] * (1 - psi[j])
    zeta[2,j,3] <- 1 - phi[j,2]
    
    zeta[3,j,1] <- 0
    zeta[3,j,2] <- 0
    zeta[3,j,3] <- 1
    
    omega[1,j,1] <- p
    omega[1,j,2] <- 0
    omega[1,j,3] <- 1 - p

    omega[2,j,1] <- 0
    omega[2,j,2] <- p
    omega[2,j,3] <- 1 - p
    
    omega[3,j,1] <- 0
    omega[3,j,2] <- 0
    omega[3,j,3] <- 1    
  }

  for (i in 1:n){
    # animals transition and are detected following release
    for (j in (f[i] + 1):t){
      z[i,j] ~ dcat(zeta[z[i, j - 1], j - 1, ])
      y[i,j] ~ dcat(omega[z[i, j], j - 1, ])
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
pars = c('alpha','beta','p')
inits <- function(){list(z = guess)}  
Sys.time()
# this will take at least a few minutes depending on your machine
m2 <- jags(data = dat, inits=inits, model.file = "ms_model2.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()

print(summary(m2), digits = 3)
m2$mcmc.info$elapsed.mins







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 3:
# a multi-state model
# now we'll use the data that are formatted as 
# 1) breeder
# 2) non-breeder
# 3) not seen
# We'll also explicitly test a 'cost-of-breeding' hypothesis on survival,
# albeit at the population, not indivdual level.
# and we'll model the population level average change in size over time...
# as well as the effect of size on breeding propensity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
# HOMEWORK:
# 1) Here we simulated a somewhat complex known truth, and we
#    used models to test increasingly complex hypotheses that more
#    accurately represented 'truth'
#
#    a. How did our hypotheses change from models 1 to 3?
#    b. How did our inference change from models 1 to 3?
#    c. What are we still missing when we compare model 3 to 'truth'
#
# 2) In script 10b, we'll add an additional piece to 'truth.' Put
#    simply, alligators will be different (beyond differences
#    in size), where some alligators will simply be 'better' than others.
#    This may represent underlying genetic or epigenetic differences,
#    differences in territory quality, etc. 
#    We'll also add an additional piece to the observation process
#    where we will look for alligators using multiple observation methods
#    This is called a 'robust design', and is similar to the approach
#    we use in occupancy models (i.e., look more than once).
# 
#    These are fairly complex capture-recapture models. m4 is an 
#    an extension of Riecke et al. [2025] in Ecology,
#    Riecke et al. [2018] in Methods in Ecol & Evol, and 
#    van de Pol & Brouwer [2021] in Journal of Animal Ecology,
#    cool, we all learned something...
#    at the same time, it still barely explains the 'toy' reality 
#    that we simulated. That reality is missing simple components
#    such as alligator age, real measurements of alligator size (length,
#    mass). Take a minute, think about your real system, and think about
#    what your model is missing. It shouldn't take too long. Write down
#    a brief list, and consider ways that you might address some of these
#    shortcomings, &| discuss how you've already tried to address them.
#
#    In other words, reconsider Box's 'all models are wrong, but some are useful'
#    aphorism, and think about how it applies to your system.
#    
# 3a) How does inference change in script 10b b/w model 3 (which does not 
#     incorporate individual heterogeneity) and m4 (which does)? Do
#     you think individual heterogeneity exists in real populations?
#     What implications does this have for inference?
#
# 3b) How does inference change for model 3 b/w scripts 10a and 10b? i.e.,
#     How does the same model perform when unmodeled individual heterogeneity is 
#     present? beta[2] is the key parameter here...
#
# 4) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



