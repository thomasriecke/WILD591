# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-lags!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'serif')
set.seed(12345)


# number of individuals
n <- 100
# number of breeding Seasons
S <- 7
# mean clutch size and body condition
mu <- c(30,0)


# create empty matrices for condition (X) and clutch size (Y)
y <- matrix(NA, n, S)
x <- matrix(NA, n, S)

# simulate correlated individual random effects using a latent
# individual or territory quality (q)
q <- rnorm(n, 0, 1)
hist(q, xlab = 'Latent quality', main = NULL, cex.lab = 2)

# individual condition and clutch size intercepts
alpha <- matrix(NA, n, 2)
alpha[,1] <- log(mu[1]) + 0.1 * q
alpha[,2] <- mu[2] + 1 * q

# plot individual intercepts...
plot(alpha[,2] ~ exp(alpha[,1]), xlab = 'Mean individual clutch size',
     ylab = 'Mean individual condition', cex.lab = 2)

# define regression parameters
beta <- c(0.1, -0.1)

hist(rnorm(10000,0,1), xlab = '')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make mean effects plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
pxt <- seq(-3,3, length.out = res)
pyt <- exp(log(mu[1]) + beta[1] * pxt)
plot(pyt ~ pxt, ylab = 'Expected clutch size in t', xlab = 'Body condition in t',
     las = 1, cex.axis = 1.5, cex.lab = 2, type = 'l', lwd = 2)
points(0, 30, cex = 2, pch = 19)

# predict condition in t+1 as a function of clutch size in t
pxtplus1 <- mu[2] + beta[2] * (pyt - mu[1])
plot(pxtplus1 ~ pyt, xlab = 'Clutch size in t', ylab = 'Body condition in t+1',
     las = 1, cex.axis = 1.5, cex.lab = 2, type = 'l', lwd = 2)
points(30, 0, cex = 2, pch = 19)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign each individual's initial quality (X[,1]) as the intercept
# model clutch size as a function of that, then
# simulate clutch size and body condition moving forward through
# time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x[,1] <- alpha[,2]
y[,1] <- rpois(n, exp(alpha[,1] + beta[1] * x[,1]))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# there are many ways to simulate this! I've chosen one that depends on:
# 1) each individual's average condition (alpha[,2])
# 2) whether they laid more or fewer eggs than average (alpha[,1])
# here are some additional potential options to consider:

# auto-regressive given difference between previous clutch and 'average' for that individual
# here average is the individual intercept (alpha[,1])
# x[,t] <- rnorm(n, x[,t-1] + beta[2] * (y[,t-1] - exp(alpha[,1])), 0.05)

# condition given number of eggs above or below what they 'should have optimally laid' adjusted from intercept
# here what the 'should have' laid is the expected clutch size in t-1
# x[,t] <- rnorm(n, alpha[,2] + beta[2] * (y[,t-1] - exp(alpha[,1] + beta[1] * x[,t-1])), 0.05)  

# condition given # eggs above or below what they 'should have optimally laid' adjusted from condition in t-1
# x[,t] <- rnorm(n, x[,t-1] + beta[2] * (y[,t-1] - exp(alpha[,1] + beta[1] * x[,t-1])), 0.05)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (t in 2:S){
  # condition given average condition and clutch size relative to individual mean
  x[,t] <- rnorm(n, alpha[,2] + beta[2] * (y[,t-1] - exp(alpha[,1])), 0.05)  
  
  # clutch size is a function of condition
  y[,t] <- rpois(n, exp(alpha[,1] + beta[1] * x[,t]))
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's make a plot of observed clutch size regressed against observed
# body condition! THIS SHOULD BE A STRONG POSITIVE EFFECT!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(jitter(y[]) ~ x[], ylab = 'Clutch size', xlab = 'Body condition prior to breeding',
     cex.lab = 2, las = 1, cex.axis = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now let's make a plot of observed condition as a function of the previous
# year's clutch size. Remember that laying more eggs means that alligator 
# body condition should decline
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(x[,2:S] ~ jitter(y[,1:(S-1)]), ylab = 'Body condition during t',
     xlab = 'Clutch size during t-1', cex.lab = 2, las = 1, cex.axis = 1.5)
cor.test(x[,2:S], y[,1:(S-1)])
text(x = 45, y = -2.9, 'r = 0.238; p < 0.0001', cex = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Goshdarnit! Jiminy crickets!! 
# What the heck is going on here?!
# Why is there a significant, weak positive correlation between clutch size
# and body condition? We just simulated that laying more eggs means losing condition
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's add one individual at a time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2))
plot(x[,2:S] ~ jitter(y[,1:(S-1)]), ylab = 'Body condition during t',
     xlim = c(10,60), ylim = c(-3,3),
     xlab = 'Clutch size during t-1', cex.lab = 2, las = 1, cex.axis = 1.5)
# cor.test(x[,2:S], y[,1:(S-1)])
# text(x = 45, y = -2.9, 'r = 0.238; p < 0.0001', cex = 1.5)

# grab the best and worst alligators
loq <- order(alpha[,1], decreasing = F)[1:3] # get the three 'worst' alligators
hiq <- order(alpha[,1], decreasing = T)[1:3]     # get the three 'best' alligators
med <- order(alpha[,1], decreasing = T)[n/2]

# this is the absolute 'worst' alligator,
# it's always in poor condition, it never lays that many eggs, and if it attempts
# to even lay the population average number of eggs, it's condition will go 
plot(x[loq[1],2:S] ~ jitter(y[loq[1],1:(S-1)]), ylab = 'Body condition during t',
       xlim = c(10,60), ylim = c(-3,3), col = 'red',
       xlab = 'Clutch size during t-1', cex.lab = 2, las = 1, cex.axis = 1.5)

# add their data to the plot...
points(x[loq[1:2],2:S] ~ jitter(y[loq[1:2],1:(S-1)]), col = 'red')
points(x[loq[1:3],2:S] ~ jitter(y[loq[1:3],1:(S-1)]), col = 'red')

# 3 high quality individuals
points(x[hiq[1:3],2:S] ~ jitter(y[hiq[1:3],1:(S-1)]), col = 'blue')
# median individual
points(x[med,2:S] ~ jitter(y[med,1:(S-1)]), col = 'black')





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cross-lag model
# note, this isn't exactly the data-generating model (we didn't include the latent
# 'quality' variable here), but it's close enough to obtain accurate inference
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sink("cross_lag.jags")
cat("
    model {

    mu[1] ~ dnorm(3.5, 1)
    mu[2] ~ dnorm(0,   1)
    
    for (j in 1:3){
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }

    for (j in 1:2){
      beta[j] ~ dnorm(0, 1)
    }

    for (i in 1:n){
    
      # random individual intercepts
      alpha[i,1] ~ dnorm(mu[1], tau[1])
      alpha[i,2] ~ dnorm(mu[2], tau[2])

      # model for clutch size in t as a function of condition in t
      for (t in 1:S){
        y[i,t] ~ dpois(exp(alpha[i,1] + beta[1] * x[i,t]))
      }
      
      # model for condition in t as a function of clutch size in t-1
      for (t in 2:S){
        x[i,t] ~ dnorm(alpha[i,2] + beta[2] * (y[i,t-1] - exp(alpha[i,1])), tau[3])
      }
    }

    }
    ",fill = TRUE)
sink()


jags.data <- list(y = y, x = x, S = S, n = n)

inits <- function(){list()}  

# Parameters monitored
parameters <- c('alpha','beta','sigma','mu')

nc <- 4
nt <- 25
ni <- 50000
nb <- 25000


# Call JAGS from R 
# 7s for 50k iterations with an Intel i9-10900 10-core processor with some other computing going on as well...
library(jagsUI)
Sys.time()
m1 <- jags(jags.data, inits, parameters, "cross_lag.jags", parallel = T, 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Sys.time()

vioplot(m1$sims.list$beta, drawRect = F)
arrows(1:2, m1$q2.5$beta, 1:2, m1$q97.5$beta, length = 0, col = 'white', lwd = 3)
points(m1$q50$beta ~ c(1:2), cex = 3, col = 'white', pch = 19)

# plot chains to visually assess convergence
plot(m1$sims.list$mu[,1])
plot(m1$sims.list$mu[,2])
plot(m1$sims.list$beta[,1])
plot(m1$sims.list$beta[,2])
plot(m1$sims.list$sigma[,1])
plot(m1$sims.list$sigma[,2])
plot(m1$sims.list$sigma[,3])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now let's try it without individual random effects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEFORE YOU RUN THE MODEL...
# what do you anticipate will happen? Will the effects change? how?
# This model is identical to the model above, but individual effects have 
# been replaced with population mens
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("cross_lag_no_individual.jags")
cat("
    model {

    mu[1] ~ dnorm(3.5, 1)
    mu[2] ~ dnorm(0,   1)
    
    for (j in 1:3){
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }

    for (j in 1:2){
      beta[j] ~ dnorm(0, 1)
    }

    for (i in 1:n){
    
      # model for clutch size in t as a function of condition in t
      for (t in 1:S){
        y[i,t] ~ dpois(exp(mu[1] + beta[1] * x[i,t]))
      }
      
      # model for condition in t as a function of clutch size in t-1
      for (t in 2:S){
        x[i,t] ~ dnorm(mu[2] + beta[2] * (y[i,t-1] - exp(mu[1])), tau[3])
      }
    }

    }
    ",fill = TRUE)
sink()


jags.data <- list(y = y, x = x, S = S, n = n)

inits <- function(){list()}  

# Parameters monitored
parameters <- c('alpha','beta','sigma','mu')

nc <- 4
nt <- 25
ni <- 50000
nb <- 25000


# Call JAGS from R 
# 7s for 50k iterations with an Intel i9-10900 10-core processor with some other computing going on as well...
library(jagsUI)
Sys.time()
m2 <- jags(jags.data, inits, parameters, "cross_lag_no_individual.jags", parallel = T, 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Sys.time()


vioplot(m2$sims.list$beta)
# malarkey


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) make a plot or two from the first model, e.g., the anticipated number
#    of eggs given body condition for an 'average' individual (estimates
#    of mu will be useful here)
#
# 2) Consider another hypothetical example (beyond our simulated alligators
#    or pintails) where this type of pattern might cause problems for inference.
#    Draw a diagram (similar to the diagrams provided on slides) demonstrating
#    hypothesized paths and how they might cause problems. 
#    Bonus points if it's related to your research, but
#    feel free to let your mind wander around, even beyond Ecology and Evolution!
#
# 3) Ensure that you clearly understand how this issue is arising. Try clearly 
#    explaining it to your 'class neighbor' (teaching is the best way to learn!)
#    or tell somebody that's not in this class!
#
# 4) If you're feeling 'wild,' consider an additional confounding factor! In this
#    example, we tracked all of these alligators across seven years with no
#    mortalities (i.e., these are magical immortal alligators!). What if alligators
#    died during the study? Specifically, what if the 'weakest' or 'lowest quality'
#    alligators  (those in poorer condition laying smaller clutches)
#    were more likely to die. This is often called 'survivorship bias.
#    How might that affect inference?! Would it?
#
# 5) If you can't stop yourself: Try to write the data-generating model, i.e.,
#    add in a latent quality variable (q in the simulation above) 
#    and link it to the individual level intercepts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~











sink("cross_lag_latent.jags")
cat("
    model {

    mu[1] ~ dnorm(3.5, 1)
    mu[2] ~ dnorm(0,   1)
    
    for (j in 1:4){
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }

    for (j in 1:2){
      beta[j] ~ dnorm(0, 1)
    }

    zeta[1] ~ dnorm(0,1)
    zeta[2] = 1

    for (i in 1:n){
     
      q.star[i] ~ dnorm(0, 1)
      q[i] = q.star[i] * sigma[4]
      
      # random individual intercepts
      alpha[i,1] ~ dnorm(mu[1] + zeta[1] * q[i], tau[1])
      alpha[i,2] ~ dnorm(mu[2] + zeta[2] * q[i], tau[2])

      # model for clutch size in t as a function of condition in t
      for (t in 1:S){
        y[i,t] ~ dpois(exp(alpha[i,1] + beta[1] * x[i,t]))
      }
      
      # model for condition in t as a function of clutch size in t-1
      for (t in 2:S){
        x[i,t] ~ dnorm(alpha[i,2] + beta[2] * (y[i,t-1] - exp(alpha[i,1])), tau[3])
      }
    }

    }
    ",fill = TRUE)
sink()


jags.data <- list(y = y, x = x, S = S, n = n)

inits <- function(){list()}  

# Parameters monitored
parameters <- c('beta','sigma','mu','zeta','alpha','q')

nc <- 4
nt <- 25
ni <- 50000
nb <- 25000


# Call JAGS from R 
# 7s for 50k iterations with an Intel i9-10900 10-core processor with some other computing going on as well...
library(jagsUI)
Sys.time()
m3 <- jags(jags.data, inits, parameters, "cross_lag_latent.jags", parallel = T, 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

Sys.time()



hist(m3$sims.list$zeta[,1])


