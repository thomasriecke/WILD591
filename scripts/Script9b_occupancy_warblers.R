# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script includes an example of SEM/occupancy models:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(latex2exp)
library(reshape2)
library(jagsUI)

set.seed(1234)
n <- 300
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
# 3) let's simulate warbler presence, 
#    warblers are more likely to occur in older forests
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha3 <- 0; beta3 <- 3
psi <- plogis(alpha3 + beta3 * maturity)
xm <- seq(-3,3,length.out = 100)
plot(psi ~ maturity, ylab = TeX("Occupancy probability ($\\psi$)"))
warblers <- rbinom(n, 1, psi)


plot(jitter(warblers,0.25) ~ maturity, cex.lab = 1.5,
     las = 1, ylab = TeX("Presence (z)"), xlab = 'Forest maturity (m)')
lines(plogis(alpha3 + beta3 * xm) ~ xm, lwd = 2)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4) Let's visualize our data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
plot(jitter(warblers, 0.25) ~ canopy,
     ylab = 'Warbler presence', xlab = 'Canopy height (m)', cex.lab = 2, las = 1)
plot(jitter(warblers, 0.25) ~ subcan,
     ylab = 'Warbler presence', xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)
plot(canopy ~ subcan, ylab = 'Canopy height (m)',
     xlab = 'Subcanopy height (m)', cex.lab = 2, las = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create detection/non-detection data (y)
# we'll visit and survey each site 3 times
# we'll summarize z, or whether or not a site is occupied, as a function of
# the number of warblers
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p = 0.35
z <- warblers
visits <- 5
y <- rbinom(n, z*visits, p)
table(y,z) # if warblers aren't present (z = 0) then we'll never see them (y = 0)
           # if they are, we can miss them (y = 0), or see them once, twice, or three times


# we can also calculate the probability of seeing warblers at least once (p*)
# as a function of the probability of missing them every time (1 - p)
p.star <- 1 - (1 - p)^visits

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5) Let's format our data for JAGS and write a model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# format data for JAGS
c <- as.numeric(scale(log(canopy)))
s <- as.numeric(scale(log(subcan)))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Our JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("townsends_occupancy.jags")
cat("
    model {

    # intercepts (alphas)
    alpha ~ dlogis(0, 1)
    
    # slopes (betas)
    beta[1] = 1
    beta[2] ~ dnorm(0, 0.1)
    beta[3] ~ dnorm(0, 0.1)

    # standard deviations (sigma) and precisions (tau)
    for (j in 1:3){
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }

    # detection probability
    p ~ dbeta(1,1)

    for (i in 1:n){
    
      # modeling latent maturity as an unobserved, centered random variable
      m.star[i] ~ dnorm(0, 1)
      m[i] = m.star[i] * sigma[3]
      # m[i] ~ dnorm(0, tau[3])
      
      
      # modeling z-standardized, log-transformed canopy and subcanopy heights
      # note, we no longer need an intercept (the mean of each is 0,
      # and sites with average canopy and sub-canopy heights will thus have 
      # 'maturity' (m) values of 0)
      can[i] ~ dnorm(beta[1] * m[i], tau[1])
      sub[i] ~ dnorm(beta[2] * m[i], tau[2])
      
      # modeling occupancy
      logit(psi[i]) = alpha + beta[3] * m[i]
      z[i] ~ dbern(psi[i])
      y[i] ~ dbin(p, v * z[i])
      

    }

    }
    ",fill = TRUE)
sink()


z.dat <- rep(NA,n)
z.dat[y > 0] <- 1
jags.data <- list(y = y, can = c, sub = s, n = n, 
                  v = visits, z = z.dat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we provide initial values. This is where JAGS will begin sampling
# to build a posterior distribution, in this case, I won't provide any 
# initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list(beta = c(NA, 1, 3))}  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is where we tell JAGS which parameters we wish to monitor
# i.e., which posterior distributions we want to save, 
# generally we'll save all of them, but 
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parameters <- c('alpha','p','beta','sigma','m',
                'elpd_waic','p_waic')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 35000
nb <- 25000


# Call JAGS from R 
# 50k iterations takes 3 seconds on an i9
library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "townsends_occupancy.jags", parallel = T, 
          modules = c('glm','dic'),
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m)
m$sims.list$



plot(m$sims.list$alpha, type = 'l', ylab = expression(alpha), las = 1)
plot(m$sims.list$p, type = 'l', ylab = expression(p), las = 1)
plot(m$sims.list$sigma[,1], type = 'l', ylab = expression(sigma[1]), las = 1)

###
### plot standard deviations of parameters
###
plot(m$sims.list$sigma[,1], type = 'l', ylab = expression(sigma[1]), las = 1)
plot(m$sims.list$sigma[,2], type = 'l', ylab = expression(sigma[2]), las = 1)
plot(m$sims.list$sigma[,3], type = 'l', ylab = expression(sigma[3]), las = 1)

###
### plot regression relationships
###
plot(m$sims.list$beta[,1], type = 'l', ylab = expression(beta[2]), las = 1)
plot(m$sims.list$beta[,2], type = 'l', ylab = expression(beta[2]), las = 1)
plot(m$sims.list$beta[,3], type = 'l', ylab = expression(beta[2]), las = 1)


###
### examine traceplots of 'maturity' estimates for specific sites
###
plot(m$sims.list$m[,2], type = 'l', ylab = expression(m[2]), las = 1)
plot(m$sims.list$m[,1], type = 'l', ylab = expression(m[1]), las = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
#
# 1) Change the number of visits to 10. What happens to the precision of your
#    estimates? Now change it to 2. Same question (precision?). Why?
#    hint: examine the value of p* we calculate above on line 90
#
# 2) make a plot of expected warbler occupancy (psi) regressed against 
#    expected canopy height (start by simulating a range of values for 
#    forest maturity)
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
