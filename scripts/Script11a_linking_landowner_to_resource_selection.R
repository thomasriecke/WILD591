# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) simulate a landscape
# 2) simulate variation in habitat quality and landowner attitude
# 3) simulate wolf use
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
library(MASS)
library(geoR)
library(sf)
library(jagsUI)
library(latex2exp)
library(reshape2)
set.seed(1234) 

par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(rnorm(100000, 0, 1), xlab = 'Landowner attitude',
     main = NULL, breaks = 1000, ylab = '', yaxt = 'n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up an empty raster and pull coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.res <- 25
n <- grid.res^2
r=rast(ncol=grid.res,nrow=grid.res,extent=c(-grid.res,grid.res,-grid.res,grid.res))                      # make an empty raster
r[]=0                                                                # cero; 0
s=crds(r)                                                          # pull coordinates

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first, we'll simulate random variation in landowner attitudes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z=grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10))         # gaussian (auto-correlated) random field
d = r                                                               # make dist a spatial raster
attitude <- z$data
d$a = attitude                                                      # assign auto-correlated attitude values
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(d['a'], main='Landowner attitude')
rm(z)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# second, we'll simulate random variation in habitat quality
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
w=grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10))         # gaussian (auto-correlated) random field
habitat <- w$data
d$h = habitat                                                      # assign auto-correlated habitat values
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(d['h'], main='Habitat quality')
rm(w,r)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some quick plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,3))
plot(d['a'], main='Landowner attitude')
plot(d['h'], main='Habitat quality')
plot(habitat ~ attitude)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate resource selection
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(0, 1, 1)
mu = exp(alpha[1] + alpha[2] * attitude + alpha[3] * habitat)
y <- rpois(n, mu)
d$r <- y
par(mfrow = c(1,1), mar = c(5,5,2,2))

str(d)
plot(d$a)

plot(y ~ habitat)
plot(y ~ attitude)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate survey (s) responses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta <- c(-2, -3, -3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create matrices for relative probabilities of each outcome (p)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5 outcomes, so five probabilities
outcomes <- 5
questions <- 3
q <- p <- array(NA, dim = c(n, outcomes, 3))
# we'll also need to specify 'cutpoints'
cut <- matrix(NA, outcomes - 1, questions)
cut[,1] <- seq(-5, 5, length.out = 4)
cut[,2] <- seq(-5, 5, length.out = 4)
cut[,3] <- seq(-5, 5, length.out = 4)

# continuous scores
epsilon <- matrix(NA, n, 3)

for (k in 1:3){
  epsilon[,k] <- beta[k] * attitude
  
  q[,1,k] <- plogis(cut[1,1] - epsilon[,k])
  p[,1,k] <- q[,1,k]
  
  for (j in 2:(outcomes-1)){
    q[,j,k] <- plogis(cut[j,k] - epsilon[,k])
    p[,j,k] <- q[,j,k] - q[,j-1,k]
  }
  p[,outcomes,k] <- 1 - q[,outcomes-1,k]
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create survey responses (S)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S <- matrix(NA, n, 3)

for (i in 1:n){
  for (k in 1:questions){
    S[i,k] <- which(rmultinom(1, 1, p[i,,k]) == 1)
  }
}

table(S[,1])
table(S[,2])
table(S[,3])

# correlations
cor(S[,1], S[,2])
cor(S[,1], S[,3])
cor(S[,2], S[,3])

par(mfrow = c(1,3))
plot(jitter(S[,1]) ~ jitter(S[,2]),
     ylab = 'I feel threatend by bears',
     xlab = 'Bears neg. affect my economic well-being')

plot(jitter(S[,1]) ~ jitter(S[,3]),
     ylab = 'I feel threatend by bears',
     xlab = 'Bears do not belong near people')

plot(jitter(S[,2]) ~ jitter(S[,3]),
     xlab = 'Bears do not belong near people',
     ylab = 'Bears neg. affect my economic well-being')





sink("ordinal_regression.jags")
cat("
  model {

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # cut-points
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   for (j in 1:K){
     for (k in 1:C){
       theta.star[j,k] ~ dnorm(0, 0.01)
     }
     theta[j,1:C] <- sort(theta.star[j,1:C])
   }
   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # ordered logit model for survey responses
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   beta[1] ~ dnorm(0, 0.1)
   beta[2] = -1
   beta[3] ~ dnorm(0, 0.1)
   sigma ~ dgamma(3,1)
   tau = pow(sigma, -2)
   
   for (i in 1:n){
     eta[i] ~ dnorm(0, tau)
     
     for (j in 1:K){
       mu[i,j] <- beta[j] * eta[i]
       S[i,j] ~ dordered.logit(mu[i,j], theta[j,1:C])
     }

   } 


   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # inhomogeneous point process
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   alpha[1] ~ dnorm(0, 1)
   alpha[2] ~ dnorm(0, 1) T(0,)
   alpha[3] ~ dnorm(0, 1)
   for (i in 1:n){
     log(lambda[i]) = alpha[1] + alpha[2] * eta[i] + alpha[3] * x[i]
     y[i] ~ dpois(lambda[i])
   }


      
  }
  ",fill = TRUE)
sink()

jags.data <- list(y = y, S = S, x = habitat, n = n, K = 3, C = 4)

init.theta <- matrix(NA, 3, 4)
init.theta[1,] <- seq(-3,3,length.out = 4)
init.theta[2,] <- seq(-3,3,length.out = 4)
init.theta[3,] <- seq(-3,3,length.out = 4)

inits <- function(){list(theta.star = init.theta,
                         beta = c(-1, NA, -1))}  
parameters <- c('alpha','beta','theta','sigma')

nc <- 4
nt <- 10
ni <- 10000
nb <- 5000


Sys.time()
m1 <- jags(jags.data, inits, parameters, "ordinal_regression.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
print(m1, digits = 3)

plot(m1$sims.list$sigma, type = 'l')
plot(m1$sims.list$beta[,1], type = 'l')
plot(m1$sims.list$beta[,2], type = 'l')
plot(m1$sims.list$beta[,3], type = 'l')


plot(m1$sims.list$theta[,1,1], type = 'l')
plot(m1$sims.list$theta[,1,2], type = 'l')
plot(m1$sims.list$theta[,1,3], type = 'l')
plot(m1$sims.list$theta[,1,4], type = 'l')
boxplot(m1$sims.list$theta[,1,])

plot(m1$sims.list$theta[,2,1], type = 'l')
plot(m1$sims.list$theta[,2,2], type = 'l')
plot(m1$sims.list$theta[,2,3], type = 'l')
plot(m1$sims.list$theta[,2,4], type = 'l')
boxplot(m1$sims.list$theta[,2,])

plot(m1$sims.list$alpha[,1], type = 'l')
plot(m1$sims.list$alpha[,2], type = 'l')
plot(m1$sims.list$alpha[,3], type = 'l')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualizing how this works...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5,5,2,2))
hist(rnorm(100000, 0, 2), xlab = 'Landowner attitude',
     main = NULL, breaks = 1000, ylab = '', yaxt = 'n')
abline(v = c(-4,-1,1,4), col = 'red', lwd = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate probabilities of outcomes for different values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- 100
test <- seq(-4, 4, length.out = tmp) 
q <- p <- matrix(NA, tmp, outcomes)
breaks <- c(-4,-1,1,4)

q[,1] <- plogis(breaks[1] - test)
p[,1] <- q[,1]

for (j in 2:(5-1)){
  q[,j] <- plogis(breaks[j] - mu)
  p[,j] <- q[,j] - q[,j-1]
}
p[,5] <- 1 - q[,5-1]
rowSums(p)


plogis(-4)
hist(tst <- rlogis(100000, 0, 1), xlim = c(-10,10),
     breaks = 500, xlab = 'Logistic(4, 1)', main = NULL)
abline(v = breaks, col = 'red')
length(which(tst < -4))
plogis(4) - plogis(1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's imagine a second model system where landowner attitudes also affect
# habitat quality?!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up an empty raster and pull coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.res <- 25
n <- grid.res^2
r=rast(ncol=grid.res,nrow=grid.res,extent=c(-grid.res,grid.res,-grid.res,grid.res))                      # make an empty raster
r[]=0                                                                # cero; 0
s=crds(r)                                                          # pull coordinates

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first, we'll simulate random variation in landowner attitudes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z=grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10))         # gaussian (auto-correlated) random field
d = r                                                               # make dist a spatial raster
attitude <- z$data
d$a = attitude                                                      # assign auto-correlated attitude values
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(d['a'], main='Landowner attitude')
rm(z)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# second, we'll simulate random variation in habitat quality given attitudes
# which have influenced quality via a history of land management
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
w=grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10))         # gaussian (auto-correlated) random field
habitat <- w$data
d$h = habitat + attitude                                            # assign auto-correlated habitat values
par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(d['h'], main='Habitat quality')
rm(w,r)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate resource selection
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- c(-1, 1, 1)
mu = exp(alpha[1] + alpha[2] * attitude + alpha[3] * habitat)
y <- rpois(n, mu)
d$r <- y

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate survey (s) responses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta <- c(-2, -3, -3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create matrices for relative probabilities of each outcome (p)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5 outcomes, so five probabilities
outcomes <- 5
questions <- 3
q <- p <- array(NA, dim = c(n, outcomes, 3))
# we'll also need to specify 'cutpoints'
cut <- matrix(NA, outcomes - 1, questions)
cut[,1] <- seq(-5, 5, length.out = 4)
cut[,2] <- seq(-5, 5, length.out = 4)
cut[,3] <- seq(-5, 5, length.out = 4)

# continuous scores
epsilon <- matrix(NA, n, 3)

for (k in 1:3){
  epsilon[,k] <- beta[k] * attitude
  
  q[,1,k] <- plogis(cut[1,1] - epsilon[,k])
  p[,1,k] <- q[,1,k]
  
  for (j in 2:(outcomes-1)){
    q[,j,k] <- plogis(cut[j,k] - epsilon[,k])
    p[,j,k] <- q[,j,k] - q[,j-1,k]
  }
  p[,outcomes,k] <- 1 - q[,outcomes-1,k]
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create survey responses (S)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S <- matrix(NA, n, 3)

for (i in 1:n){
  for (k in 1:questions){
    S[i,k] <- which(rmultinom(1, 1, p[i,,k]) == 1)
  }
}

table(S[,1])
table(S[,2])
table(S[,3])

par(mfrow = c(1,3))
plot(jitter(S[,1]) ~ jitter(S[,2]),
     ylab = 'I feel threatend by bears',
     xlab = 'Bears neg. affect my economic well-being')

plot(jitter(S[,1]) ~ jitter(S[,3]),
     ylab = 'I feel threatend by bears',
     xlab = 'Bears do not belong near people')

plot(jitter(S[,2]) ~ jitter(S[,3]),
     xlab = 'Bears do not belong near people',
     ylab = 'Bears neg. affect my economic well-being')



sink("ordinal_regression2.jags")
cat("
  model {

   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # cut-points
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   for (j in 1:K){
     for (k in 1:C){
       theta.star[j,k] ~ dnorm(0, 0.01)
     }
     theta[j,1:C] <- sort(theta.star[j,])
   }
   
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # ordered logit model for survey responses
   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   beta[1] = -1
   beta[2] ~ dnorm(0, 0.01)
   beta[3] ~ dnorm(0, 0.01)
   sigma[1] ~ dgamma(1,1)
   tau[1] = pow(sigma[1], -2)
   for (i in 1:n){
     eta[i] ~ dnorm(0, tau[1])
     
     for (j in 1:K){
       mu[i,j] <- beta[j] * eta[i]
       S[i,j] ~ dordered.logit(mu[i,j], theta[j,])
     }

   } 

   alpha[1] ~ dnorm(0, 1)
   alpha[2] ~ dnorm(0, 1) T(0,)
   alpha[3] ~ dnorm(0, 1)
   sigma[2] ~ dgamma(1,1)
   tau[2] = pow(sigma[2], -2)
   psi ~ dnorm(1, 1)
   for (i in 1:n){
     x[i] ~ dnorm(psi * eta[i], tau[2])
     log(lambda[i]) = alpha[1] + alpha[2] * eta[i] + alpha[3] * x[i]
     y[i] ~ dpois(lambda[i])
   }


      
  }
  ",fill = TRUE)
sink()

jags.data <- list(y = y, x = habitat, n = n, K = 3, C = 4, S = S)

init.theta <- matrix(NA, 3, 4)
init.theta[1,] <- seq(-3,3,length.out = 4)
init.theta[2,] <- seq(-3,3,length.out = 4)
init.theta[3,] <- seq(-3,3,length.out = 4)

inits <- function(){list(theta.star = init.theta)}  
parameters <- c('alpha','beta','theta','psi')

nc <- 4
nt <- 10
ni <- 2500
nb <- 1000


Sys.time()
m2 <- jags(jags.data, inits, parameters, "ordinal_regression2.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)
Sys.time()
print(m2, digits = 3)















