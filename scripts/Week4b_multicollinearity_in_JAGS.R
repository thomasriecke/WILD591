# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Simulation 1:
#
# Let's imagine a pretty simple (and absurd) dataset. 
# Our response variable will be counts (y) of a bird species (yellow-footed weeble-wobbles).
#
# our two covariates will be canopy height (ch, in meters) 
# and sub-canopy height (sch, in meters). These will be driven by a latent
# variable ('forest maturity')
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
library(vioplot)
set.seed(1234)
# first we'll define the number of sites
n <- 200

# then we'll simulate forest 'maturity' or age as a z-standardized covariate
# large positive values will represent 'old-growth' or 'over-mature' forests,
# very negative values will represent doghair pole-sapling thickets
maturity <- rnorm(n, 0, 1)

# then we'll simulate canopy height
canopy <- rlnorm(n, 3.65 + 0.25 * maturity, 0.05)
hist(canopy, breaks = 50, main = NULL, xlab = 'Canopy height (m)')

# similarly, we'll simulate sub-canopy height
subcan <- rlnorm(n, 2 + 0.5 * maturity, 0.1) 

# extreme multicollinearity
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(canopy ~ subcan, ylab = 'Canopy height (m)',
     xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)

# how many warblers?
warblers <- rpois(n, exp(0.5 + 0.75 * maturity))


# so what's the signal
plot(jitter(warblers) ~ canopy,
     ylab = 'Warbler abundance', xlab = 'Canopy height (m)', cex.lab = 2, las = 1)


plot(jitter(warblers) ~ subcan,
     ylab = 'Warbler abundance', xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)


# format data for JAGS
y <- warblers
c <- as.numeric(scale(canopy))
s <- as.numeric(scale(subcan))
log(mean(y))


sink("m_both.jags")
cat("
    model {

    alpha ~ dnorm(1, 0.1)
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)

    for (i in 1:n){
      log(psi[i]) = alpha + beta[1] * c[i] + beta[2] * s[i]    # log-link
      y[i] ~ dpois(psi[i])                                  # Poisson distribution
    }

    }
    ",fill = TRUE)
sink()

sink("m_sub.jags")
cat("
    model {

    alpha ~ dnorm(1, 0.1)
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)

    for (i in 1:n){
      log(psi[i]) = alpha + beta[2] * s[i]    # log-link
      y[i] ~ dpois(psi[i])                                  # Poisson distribution
    }

    }
    ",fill = TRUE)
sink()

sink("m_can.jags")
cat("
    model {

    alpha ~ dnorm(1, 0.1)
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)

    for (i in 1:n){
      log(psi[i]) = alpha + beta[1] * c[i]                  # log-link
      y[i] ~ dpois(psi[i])                                  # Poisson distribution
    }

    }
    ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here is where we provide JAGS with the data it will need to run the model,
# in this instance it needs y (our response variable, mass), x (our covariate, length
# in feet), and n (our sample size)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(y = y, c = c, s = s, n = n)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we provide initial values. This is where JAGS will begin sampling
# to build a posterior distribution, in this case, I won't provide any 
# initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list()}  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is where we tell JAGS which parameters we wish to monitor
# i.e., which posterior distributions we want to save, 
# generally we'll save all of them, but 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c('alpha','beta','WAIC')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000


# Call JAGS from R 
# 50k iterations takes 3 seconds on an i9
library(jagsUI)
# library(R2jags)
Sys.time()
m.b <- jags(jags.data, inits, parameters, "m_both.jags", parallel = T, 
            modules = c('glm','dic'),
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

m.c <- jags(jags.data, inits, parameters, "m_can.jags", parallel = T,
            modules = c('glm','dic'),
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

m.s <- jags(jags.data, inits, parameters, "m_sub.jags", parallel = T, 
            modules = c('glm','dic'),
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Sys.time()



print(m.b, digits = 3)
print(m.c, digits = 3) # we can ignore beta[2], that's just the prior
print(m.s, digits = 3) # we can ignore beta[1], that's just the prior





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 'HOMEWORK'
# Whew! Now things are really getting 'fun.'
# Let's take a minute to look at the parameter estimates of the effects of 
# canopy height (beta[1]) and sub-canopy height (beta[2]) from the models
# that only included single covariates (m.c and m.s)
#
par(mfrow = c(1,2))
vioplot(c(m.c$sims.list$beta[,1]), c(m.s$sims.list$beta[,2]),
        ylim = c(-0.5,1),
        las = 1, ylab = expression(beta))
#
# Ok. They look fairly similar! Now let's take a peak at the effects of these 
# covariates from a model when both covariates were included
#
vioplot(m.b$sims.list$beta, add = T, col = adjustcolor('red',alpha = 0.5),
        names = expression(beta[1],beta[2]), ylim = c(-0.5,1),
        las = 1, ylab = expression(beta), wex = 0.5)
# huh... what happened?!
# 
# let's take a look at the posterior distributions plotted against each other...
plot(m.b$sims.list$beta[,1] ~ m.b$sims.list$beta[,2],
     ylab = 'Effect of canopy', xlab = 'Effect of subcanopy')
# this is where the term 'variance inflation' comes from.
#
# Now this is starting to get interesting! What is happening here?
# here are some clues to get you started, and we'll revisit this as a group
#
par(mfrow = c(1,1))
plot(warblers ~ c, ylab = 'Counts', xlab = 'Z-standardized canopy')
points(0, exp(m.c$q50$alpha), cex = 2, pch = 19, col = 'red') # this is our intercept
lines(exp(c(m.c$q50$alpha + m.c$q50$beta[1]*seq(-3,3,length.out = 100))) ~ seq(-3,3,length.out = 100), col = 'red')
#
joint <- rowSums(m.b$sims.list$beta) # the effect of both covariates from the model with both
vioplot(m.c$sims.list$beta[,1], m.s$sims.list$beta[,2], joint,
        ylim = c(-0.5,1),
        las = 1, ylab = expression(beta))
vioplot(m.b$sims.list$beta, add = T, col = adjustcolor('red',alpha = 0.5),
        ylim = c(-0.5,1),
        las = 1, ylab = expression(beta), wex = 0.5)
# this plot shows the effect of canopy height from a stand-alone model (leftmost violin)
# the effect of subcanopy height from a stand-alone model (center violin)
# and the joint effect of both covariates (rightmost violin)
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation 2:
#
# Let's imagine a second dataset. Now our response variable will be the 
# amount of 'browse' available on the landscape (y)
#
# at each site, we'll also have estimates of wolf (w) and ungulate (u)
# abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12345)
# 100 sites
n <- 100

# wolf density
w <- rlnorm(n, 0.35, 0.75)
hist(w, xlab = 'Wolf density', breaks = 15)

# ungulate density will be a function of wolf density (not on the same scale!)
u <- rlnorm(n, 3.5 + -0.25 * w, 0.25)
plot(u ~ w, ylab = 'Ungulate density', xlab = 'Wolf density') # these are on different scales!

# veg density
v <- rlnorm(n, 6 + -0.15 * u, 0.25)
plot(v ~ u, ylab = 'Browse/vegetation density', xlab = 'Ungulate density')

plot(v ~ w, ylab = 'Browse/vegetation density', xlab = 'Wolf density')

zw <- as.numeric(scale(w))
zu <- as.numeric(scale(u))

sink("m_both.jags")
cat("
    model {

    alpha ~ dnorm(1, 0.1)
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)
    sigma ~ dgamma(1,1)
    tau = pow(sigma,-2)  # equivalent to tau = 1/(sigma^2) or tau = 1/(sigma*sigma)
    
    for (i in 1:n){
      psi[i] = alpha + beta[1] * w[i] + beta[2] * u[i]    # log-link
      v[i] ~ dlnorm(psi[i], tau)                                  # Poisson distribution
    }

    }
    ",fill = TRUE)
sink()

sink("m_wolf.jags")
cat("
    model {

    alpha ~ dnorm(1, 0.1)
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)
    sigma ~ dgamma(1,1)
    tau = pow(sigma,-2)  # equivalent to tau = 1/(sigma^2) or tau = 1/(sigma*sigma)
    
    for (i in 1:n){
      psi[i] = alpha + beta[1] * w[i]                             # log-link
      v[i] ~ dlnorm(psi[i], tau)                                  # Poisson distribution
    }

    }
    ",fill = TRUE)
sink()

sink("m_deer.jags")
cat("
    model {

    alpha ~ dnorm(1, 0.1)
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)
    sigma ~ dgamma(1,1)
    tau = pow(sigma,-2)  # equivalent to tau = 1/(sigma^2) or tau = 1/(sigma*sigma)    

    for (i in 1:n){
      psi[i] = alpha + beta[2] * u[i]    # log-link
      v[i] ~ dlnorm(psi[i], tau)                                  # Poisson distribution
    }

    }
    ",fill = TRUE)
sink()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here is where we provide JAGS with the data it will need to run the model,
# in this instance it needs y (our response variable, mass), x (our covariate, length
# in feet), and n (our sample size)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
jags.data <- list(v = v, w = zw, u = zu, n = n)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we provide initial values. This is where JAGS will begin sampling
# to build a posterior distribution, in this case, I won't provide any 
# initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list()}  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is where we tell JAGS which parameters we wish to monitor
# i.e., which posterior distributions we want to save, 
# generally we'll save all of them, but 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c('alpha','beta')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000


# Call JAGS from R 
# 50k iterations takes 3 seconds on an i9
library(jagsUI)
Sys.time()
m.b <- jags(jags.data, inits, parameters, "m_both.jags", parallel = T, 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

m.w <- jags(jags.data, inits, parameters, "m_wolf.jags", parallel = T, 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

m.u <- jags(jags.data, inits, parameters, "m_deer.jags", parallel = T, 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Sys.time()

print(m.b, digits = 3)
print(m.w, digits = 3) # note that beta[2] is simply the prior! (i.e., meaningless)
print(m.u, digits = 3) # note that beta[1] is simply the prior! (i.e., meaningless)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 'HOMEWORK'
# Make some similar plots as those above, exploring the effects of wolf
# and ungulate density on vegetation density, e.g.,
plot(m.b$sims.list$beta[,1] ~ m.b$sims.list$beta[,2], 
     ylab = 'Effect of wolves', xlab = 'Effect of ungulates')
# 
vioplot(m.w$sims.list$beta[,1], m.u$sims.list$beta[,2], 
        names = expression(beta[1],beta[2]), ylim = c(-4,4),
        las = 1, ylab = expression(beta))
vioplot(m.b$sims.list$beta, add = T, col = adjustcolor('red',alpha = 0.5),
        wex = 0.5)
# what's happening here with our parameter estimates?
# next week, we'll use the exact same system to explore the use of path
# analysis, or sequences of linear models that allow us to estimate 
# causal relationships among measured values.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~