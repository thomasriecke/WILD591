# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Welcome to latent variables. This type of modelling is one of the most
# intuitive things I've ever done conceptually. There are also statistical components of
# this process that are incredibly non-intuitive :)
#
# Our goal this week is to familiarize ourselves with the concept (i.e.,
# we're observing multiple measurements of an underlying 'latent' process), and 
# to begin to familiarize ourselves with implementing these types of models
# in lavaan/blavaan and JAGS. We will take our time. This is a critical concept
# for understanding SEMs. Ask questions as they arise!
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# our first simulated dataset will be counts of warblers (y; Setophoga townsendi)
# in different age-classes of forests. Counts will be greater in older forests
# 
# We will (initially) measure two components of forest age, canopy height (c) and 
# sub-canopy height (s).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load packages
library(jagsUI)
library(vioplot)
library(lavaan)
library(blavaan)
library(lavaan)
library(latex2exp)

set.seed(1234)
n <- 200
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Let's simulate forest age (i.e., maturity) as a z-standardized covariate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
maturity <- rnorm(n, 0, 1)
hist(maturity, breaks = 50, main = NULL, xlab = 'Forest maturity', cex.lab =2)
# large positive values indicate 'old growth', very negative numbers 
# are doghair thickets (dense, short, regen following logging or stand-replacing
# fires)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) let's simulate canopy height and sub-canopy height
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
canopy <- rlnorm(n, 3.65 + 0.25 * maturity, 0.05)
hist(canopy, breaks = 50, main = NULL, xlab = 'Canopy height (m)')
plot(canopy ~ maturity, ylab = 'Canopy height (m)',
     xlab = 'Forest maturity', cex.lab = 2)


# similarly, we'll simulate sub-canopy height
subcan <- rlnorm(n, 2 + 0.5 * maturity, 0.1) 
hist(subcan, breaks = 50, main = NULL, xlab = 'Sub-canopy height (s)')
plot(subcan ~ maturity, ylab = 'Sub-canopy height (m)',
     xlab = 'Forest maturity', cex.lab = 2)

# Note the extreme multicollinearity between the two covariates... this
# is because they're arising from the same underlying 'latent' process.
# As forest stand mature, trees get taller, similarly, as stands mature and 
# gaps develop, there is more room for regeneration and second-growth
plot(canopy ~ subcan, ylab = 'Canopy height (m)',
     xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) let's simulate warblers, more warblers in older forests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
warblers <- rpois(n, exp(0.5 + 0.75 * maturity))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Let's visualize our data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,3), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
plot(jitter(warblers) ~ canopy,
     ylab = 'Warbler abundance', xlab = 'Canopy height (m)', cex.lab = 2, las = 1)
plot(jitter(warblers) ~ subcan,
     ylab = 'Warbler abundance', xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)
plot(canopy ~ subcan, ylab = 'Canopy height (m)',
     xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Let's format our data for JAGS and write a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format data for JAGS
y <- warblers
c <- canopy
s <- subcan

log(mean(s))
log(mean(c))
log(mean(y))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# checking some vague priors
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hist(exp(rnorm(10000, log(mean(c)), 1)))  # canopy height is somewhere between 0 and 1500m tall...
hist(exp(rnorm(10000, log(mean(s)), 1)))  # subcanopy height is somewhere between 0 and 300m tall...
hist(exp(rnorm(10000, log(mean(y)), 1)))  # somewhere between 0 and 100 wrablers






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# using lavaan for the first time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create a data.frame (d) with relevant data
# warbler counts, canopy, and subcanopy heights
d <- data.frame(y,c,s)

sem1 <- sem('m =~ c + s
             c ~ 1
             s ~ 1
             y ~ m', data = d)
summary(sem1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's explore some things... first, where do our intercepts come from?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(c) # average canopy height is our intercept for c
mean(s) # average subcanopy height is our intercept for s
mean(y) # average counts are our intercept for y

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# derive estimates of latent maturity...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m.hat <- predict(sem1)


par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(as.numeric(m.hat) ~ maturity, cex.lab = 2, las = 1,
     ylab = TeX("$\\hat{m}"), xlab = 'Simulated latent maturity (m)')











# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Our JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("our_first_latent_variable.jags")
cat("
    model {

    # intercepts (alpha[1])
    alpha[1] ~ dnorm(2.09, 1)
    alpha[2] ~ dnorm(3.64, 1)
    alpha[3] ~ dnorm(0.75, 1)
    
    beta[1] ~ dnorm(0, 0.1)
    beta[2] ~ dnorm(0, 0.1)
    beta[3] ~ dnorm(0, 0.1)

    for (j in 1:3){
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }


    for (i in 1:n){
    
      # modeling latent maturity as an unobserved, centered random variable
      m[i] ~ dnorm(0, tau[3])
      
      # modeling measured canopy and subcanopy heights
      c[i] ~ dlnorm(alpha[1] + beta[1] * m[i], tau[1])
      s[i] ~ dlnorm(alpha[2] + beta[2] * m[i], tau[2])
      
      # modeling counts of birds
      log(psi[i]) = alpha[3] + beta[3] * m[i]
      y[i] ~ dpois(psi[i]) 
      
      
      # generating new data (to test goodness of fit)
      # not delving into that yet just showing you how easy it is
      # to generate data from a model in JAGS...
      # y.new[i] ~ dpois(psi[i])
      
      
    }

    }
    ",fill = TRUE)
sink()

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
parameters <- c('alpha','beta','sigma','m','y.new')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000


# Call JAGS from R 
# 50k iterations takes 3 seconds on an i9
library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "our_first_latent_variable.jags", parallel = T, 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m)


plot(m$sims.list$sigma[,3])
par(mfrow = c(1,2))

par(mfrow = c(1,3))
plot(m$sims.list$beta[,1])
plot(m$sims.list$m[,4])

huh <- m$sims.list$m[,4] * m$sims.list$beta[,1]
plot(huh)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 0) there is code below to run similar models in blavaan
#
# 1) Make a plot of the relationship between the rate of the poisson
#    distribution (psi) and forest maturity
#
# 2) Make a plot of the relationship between the rate of the poisson
#    distribution (psi) and canopy height (this might be tricky)
#
# 3) Make a plot of the expected counts of warblers drawn from a Poisson distribution
#    and forest maturity
# 
# 4) Change which loading is fixed, i.e., rather than fixing the relationship
#    between canopy and maturity to 1, fix the relationship between subcanopy and maturity
#    What happens? Pay particular attention to sigma[3] and beta. Does
#    that change our predicted number of warblers?
#
# 5) If you're feeling 'wild,' don't fix anything to 1. Try to estimate every 
#    parameter after assigning a reasonable prior. What happens? Why?!
#    We'll discuss this on Thursday morning. See the bottom of the script for a
#    brief explanation (lines 250ish)
#
# 6) If you literally cannot stop yourself from playing with this code, change
#    a few of the canopy OR sub-canopy heights to NAs to simulate somebody
#    not writing everything down on a datasheet (--it happens). Something like this will
#    do the trick: 
#    c[c(1,3,7)] <- NA
#    Run the model again while saving the estimates for 'c'
#    what happens? why? that's cool, huh?!
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# note we're log-transforming canopy and subcanopy height for lavaan and blavaan
# I'd do the same for the counts but we need to draw them from a Poisson (
# and could do that via modifying some JAGS code created by blavaan)
can <- log(canopy)
sub <- log(subcan)
warb <- y
dat <- data.frame(cbind(can,sub,warb))


library(lavaan)
ml <- sem(model = '
          maturity =~ can + sub
          warb ~ 1 + maturity',
           data = dat)
summary(ml)



library(blavaan)
mb <- bsem(model = '
          maturity =~ can + sub
          warb ~ 1 + maturity',
          data = dat,
          mcmcfile = T,
          target = 'jags')
summary(mb)
fitmeasures(mb)
# ppp: posterior predictive p-value is horrible (we're shooting for ~0.5 here)
# because our counts are from a poisson log-link and we're modeling them as linear here






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fixing things to 1...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# What does fixing something to 1 do?! That feels weird (it still does to me...
# reviewers also seem to really $*#&$&@ hate it).
#
# Remember way up in the code when we specified that the variance and sd of
# forest maturity equaled 1? Let's look at our estimate from the model
vioplot(m$sims.list$sigma[,3], names = 'sigma_maturity')
# it's about 0.25, or one-fourth of what we simulated.
#
# We also specified a linear relationship between maturity and canopy height,
# that beta was 0.25... but in the model we fixed it to 1 (or 4 times greater
# than 0.25)
#
# Similarly, we specified a relationship between subcanopy height and maturity.
# that was 0.5, but our model estimated it at 2 (4 times greater than simulated)
vioplot(m$sims.list$beta[,2], names = 'sigma_maturity')
#
# So what the heck is happening? Well, let's look at our estimates of maturity
plot(m$q50$m ~ maturity, ylab = 'Estimated maturity', xlab = 'Simulated maturity')
arrows(maturity, m$q2.5$m, maturity, m$q97.5$m, length = 0, lty = 2)
points(m$q50$m ~ maturity, pch = 21, bg = 'grey50', cex = 2)
# Now that is interesting! Simulated values of -3 correspond to estimated values
# of -0.75 (1/4). Similarly, a simulated value of 2 corresponds to an estimate 
# of about 0.5 on the y-axis...
#
# Now let's think about the relationship between warblers and maturity, we simulated
# this to be 0.75... and estimated it to be... 3!! That's 4x greater.
vioplot(m$sims.list$beta[,3], names = 'sigma_maturity')
#
# So we are recovering our parameter estimates (!!), they're just on a different scale.
# if a one-unit change in maturity leads to a 1 unit change in canopy height on
# the log-scale (that's what we predetermined when we fixed that loading to 1),
# then the standard deviation of canopy is just going to have to shrink!!
# 
# Go back up to your model and fix the loading between canopy height and maturity
# to 0.25... what happens to your estimate of sigma?
#
# We'll discuss this thoroughly on Thursday as well!
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

