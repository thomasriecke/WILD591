# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation:
#
# Let's imagine a dataset. Our response variable will be the 
# amount of 'browse' available on the landscape (y)
#
# At each site, we'll also have estimates of wolf (w) and ungulate (u)
# abundance. Our hypothesis is that wolf abundance affects ungulate abundance,
# and ungulate abundance affects browse. Obviously any real system
# would be more complex (we're leaving out human drivers of wolf abundance,
# and I assume this would all be regulated via bottom-up processes as well),
# but let's keep it simple for the moment.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
library(lavaan)
library(vioplot)



set.seed(12345)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
# 100 sites
n <- 100

# wolf density
w <- rlnorm(n, 0.35, 0.75)
hist(w, xlab = 'Wolf density', 
     breaks = 25, main = NULL, cex.lab = 2)

# ungulate density will be a function of wolf density (not on the same scale!)
alpha0 = 3.75
alpha1 = -0.225
u <- rlnorm(n, alpha0 + alpha1 * w, 0.25)
plot(u ~ w, ylab = 'Ungulate density', xlab = 'Wolf density',
     cex.lab = 2) # these are on different scales!
lines(exp(alpha0 + alpha1*seq(0,7,length.out = 100)) ~ seq(0,7,length.out = 100))

# veg density
beta0 <- 0
beta1 <- -0.125
v <- rlnorm(n, beta0 + beta1 * u, 0.25)
plot(v ~ u, ylab = 'Browse/vegetation density', xlab = 'Ungulate density',
     cex.lab = 2)
lines(exp(beta0 + beta1*seq(0,70,length.out = 100)) ~ seq(0,70,length.out = 100))

# after: wolves -> ungulates AND ungulates -> vegetation
# we're left with this relationship!
plot(v ~ w, ylab = 'Browse/vegetation density', xlab = 'Wolf density', cex.lab = 2)
# how would we formally estimate this??!!

# let's look at some quick ML estimates from lm()
summary(m.uw <- lm(u ~ w))
summary(m.vu <- lm(v ~ u))






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's go back to the board and explore this!?
# what's the effect on ungulates of going from 2 wolves to NO wolves
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# first we'll pull our regression coefficients from our first model; lm(u ~ w)
alpha0 <- as.numeric(m.uw$coefficients[1])
alpha1 <- as.numeric(m.uw$coefficients[2])

# number of wolves at beginning and end
start.wolf <- 2         # initial number of wolves
end.wolf <- 0           # number of wolves after treatment

# make a plot
plot(u ~ w, ylab = 'Ungulate density', xlab = 'Wolf density',
     cex.lab = 2)

# add lines and points at specific values
lines(alpha0 + alpha1*seq(0,7,length.out = 100) ~ seq(0,7,length.out = 100))
points(start.wolf, alpha0 + alpha1*start.wolf, col = 'red', pch = 19, cex = 2); 
points(end.wolf, alpha0 + alpha1*end.wolf, col = 'royalblue4', pch = 19, cex = 2); 39.2352 + -4.2195*0

# calculate values on y-axis
start.deer <- alpha0 + alpha1*start.wolf
end.deer <- alpha0 + alpha1*end.wolf
end.deer - start.deer


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# how does that change in ungulates affect veg?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pull our regression coefficients from our second model; lm(v ~ u)
beta0 <- as.numeric(m.vu$coefficients[1])
beta1 <- as.numeric(m.vu$coefficients[2])


plot(v ~ u, ylab = 'Browse/vegetation density', xlab = 'Ungulate density',
     cex.lab = 2)
lines(beta0 + beta1*seq(0,60,length.out = 100) ~ seq(0,60,length.out = 100))
points(start.deer, beta0 + beta1*start.deer, 
       col = 'red', pch = 19, cex = 2)
points(end.deer, beta0 + beta1*end.deer, 
       col = 'royalblue4', pch = 19, cex = 2)

start.veg <- beta0 + beta1*start.deer
end.veg <- beta0 + beta1*end.deer

# the change in veg = alpha1 * beta1 * change in wolves (-2)
end.veg - start.veg
alpha1 * beta1 * -2

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's calculate the total change
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
change <- end.veg - start.veg      # total change in wolves
indirect <- alpha1 * beta1   # indirect effect = alpha_1 * beta_1
indirect * -2                      # indirect effect times loss of two wolves!?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# woah!?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's calculate and plot the predicted effect of wolves on vegtation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(v ~ w, ylab = 'Browse/vegetation density', xlab = 'Wolf density',
     cex.lab = 2)
intercept <- beta0 + beta1 * alpha0  # amount of veg as a function of deer density given no wolves
slope <- alpha1 * beta1
lines(intercept + slope * seq(0,10,length.out = 100) ~ seq(0,10,length.out = 100))

# if no wolves
points(0, beta0 + beta1*alpha0, cex = 2, pch = 19, col = 'red')
# if two wolves (where start.wolf <- 2)
points(2, beta0 + beta1*alpha0 + alpha1*beta1*start.wolf, cex = 2, pch = 19, col = 'royalblue4')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK (change the stop and start values for wolves on Lines 70 and 71)
# 1) What would the total change in vegetation be if we moved from 3 wolves to 0 wolves?
# 2) What about if we moved from 2 wolves to 4 wolves?
#
# You can calculate this simply using the 'indirect effect' above, or by working through the
# system of equations (I recommend both!)
# 
# there is a JAGS implementation of this model (normal distributions) at the bottom of the script.
# But for now, let's avoid negative values for vegetation!?
# To do that, we'll incorporate link functions into our JAGS model. This can be done using 
# 'PiecewiseSEM', but is far more flexible in the BUGS (or any other fully customizable) language
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m_path.jags")
cat("
    model {

    # priors for everything assigned in a loop to minimize unnecessary coding mistakes
    for (j in 1:2){
      alpha[j] ~ dnorm(0, 0.1)
      beta[j] ~ dnorm(0, 0.1)
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }

    for (i in 1:n){
      u[i] ~ dlnorm(alpha[1] + alpha[2] * w[i], tau[1])
      v[i] ~ dlnorm(beta[1] + beta[2] * u[i], tau[2])
    }

    }
    ",fill = TRUE)
sink()

jags.data <- list(n = n, w = w, u = u, v = v)
inits <- function(){list()}  
parameters <- c('alpha','beta','sigma')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000



library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "m_path.jags", parallel = T, 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()


vioplot(m$sims.list$alpha[,2], m$sims.list$beta[,2],
        # names = c(expression(alpha[1]),expression(beta[1]))
        names = c('alpha1','beta1'),
        las = 1)
# note that multiplying alpha_1 * beta_1 will be fairly meaningless
# on the log-link!



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's make a plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iters <- length(m$sims.list$alpha[,1]) # we saved 4000 iterations
res <- 100

xw <- seq(0,10, length.out = res)
xu <- seq(0,60, length.out = res)
Eu <- matrix(NA, iters, res)
Evw <- matrix(NA, iters, res) 
Evu <- matrix(NA, iters, res)


qEu <- matrix(NA, res, 5)
qEvw <- matrix(NA, res, 5)
qEvu <- matrix(NA, res, 5)


for (j in 1:res){
  Eu[,j] <- exp(m$sims.list$alpha[,1] + m$sims.list$alpha[,2] * xw[j])  # expected value of ungulates given values of wolves (xw)
  Evw[,j] <- exp(m$sims.list$beta[,1] + m$sims.list$beta[,2] * Eu[,j])  # expected value of veg given ungulates estimated from values of wolves
  Evu[,j] <- exp(m$sims.list$beta[,1] + m$sims.list$beta[,2] * xu[j])   # expected value of veg given range of ungulates (xu)
  
  qEu[j, ] <- quantile(Eu[,j], c(0.025,0.05,0.5,0.95,0.975))
  qEvw[j, ] <- quantile(Evw[,j], c(0.025,0.05,0.5,0.95,0.975))
  qEvu[j, ] <- quantile(Evu[,j], c(0.025,0.05,0.5,0.95,0.975))
}

mEu <- melt(Eu); names(mEu) <- c('iter','x','u')
mEvw <- melt(Evw); names(mEvw) <- c('iter','x','v')
mEvu <- melt(Evu); names(mEvu) <- c('iter','x','v')

# plot of expected ungulate density as a function of wolves
smoothScatter(mEu$u ~ xw[mEu$x], las = 1, nrpoints = 0,
              ylab = 'Expected density of ungulates',
              xlab = 'Density of wolves')
lines(qEu[,3] ~ xw, lty = 1, lwd = 3, col = 'white')
lines(qEu[,1] ~ xw, lty = 2, lwd = 3, col = 'white')
lines(qEu[,5] ~ xw, lty = 2, lwd = 3, col = 'white')


# plot of expected vegetation density as a function of wolves
smoothScatter(mEvw$v ~ xw[mEvw$x], las = 1, nrpoints = 0,
              ylab = 'Expected density of browse',
              xlab = 'Density of wolves')
lines(qEvw[,3] ~ xw, lty = 1, lwd = 3, col = 'white')
lines(qEvw[,1] ~ xw, lty = 2, lwd = 3, col = 'white')
lines(qEvw[,5] ~ xw, lty = 2, lwd = 3, col = 'white')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) above we made a plot demonstrating the predicted effects of changes in wolves
#    on changes in vegetation. Could you use components we made to generate that plot to plot the effect
#    of wolves on ungulates, and ungulates on vegetation?
#
# 2) What is the probability that an increase in density from two wolves to four would lead to 
#    an increase in browse density? 
# 
#    Note: to do that we might have to incorporate our estimates of
#    uncertainty as well! 
#
#    If you're uncertain about how to do this, just work on obtaining predictions of 
#    expected values of vegetation at 2 wolves and 4 wolves above...? Some suggested code
#    for incorporating uncertainty (i.e. prediction intervals) below
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xw <- seq(0,10, length.out = res)
xu <- seq(0,60, length.out = res)
Eu <- matrix(NA, iters, res)
Evw <- matrix(NA, iters, res) 
Evu <- matrix(NA, iters, res)


qEu <- matrix(NA, res, 5)
qEvw <- matrix(NA, res, 5)
qEvu <- matrix(NA, res, 5)


for (j in 1:res){
  Eu[,j] <- rlnorm(iters, m$sims.list$alpha[,1] + m$sims.list$alpha[,2] * xw[j], m$sims.list$sigma[,1])  # expected value of ungulates given values of wolves (xw)
  Evw[,j] <- rlnorm(iters, m$sims.list$beta[,1] + m$sims.list$beta[,2] * Eu[,j], m$sims.list$sigma[,2])  # expected value of veg given ungulates estimated from values of wolves
  Evu[,j] <- rlnorm(iters, m$sims.list$beta[,1] + m$sims.list$beta[,2] * xu[j], m$sims.list$sigma[,2])   # expected value of veg given range of ungulates (xu)
  
  qEu[j, ] <- quantile(Eu[,j], c(0.025,0.05,0.5,0.95,0.975))
  qEvw[j, ] <- quantile(Evw[,j], c(0.025,0.05,0.5,0.95,0.975))
  qEvu[j, ] <- quantile(Evu[,j], c(0.025,0.05,0.5,0.95,0.975))
}

mEu <- melt(Eu); names(mEu) <- c('iter','x','u')
mEvw <- melt(Evw); names(mEvw) <- c('iter','x','v')
mEvu <- melt(Evu); names(mEvu) <- c('iter','x','v')



# plot of expected ungulate density as a function of wolves
smoothScatter(mEu$u ~ xw[mEu$x], las = 1, nrpoints = 0,
              ylab = 'Expected density of ungulates',
              xlab = 'Density of wolves')
lines(qEu[,3] ~ xw, lty = 1, lwd = 3, col = 'white')
lines(qEu[,1] ~ xw, lty = 2, lwd = 3, col = 'white')
lines(qEu[,5] ~ xw, lty = 2, lwd = 3, col = 'white')


# plot of expected vegetation density as a function of wolves
smoothScatter(mEvw$v ~ xw[mEvw$x], las = 1, nrpoints = 0,
              ylab = 'Expected density of browse',
              xlab = 'Density of wolves')
lines(qEvw[,3] ~ xw, lty = 1, lwd = 3, col = 'white')
lines(qEvw[,1] ~ xw, lty = 2, lwd = 3, col = 'white')
lines(qEvw[,5] ~ xw, lty = 2, lwd = 3, col = 'white')

# interesting... what just happened to our uncertainty intervals?! :) Why?!

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we find values on our vector from 0-10 wolves 
# we use for prediction that roughly correspond to 2 and 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
xw[21]
xw[41]
# note that we could actually include the values 2 and 4 in this gradient

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# then we test which of the veg values with four wolves are greater than
# the simulated veg values with 2 wolves
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(which(Evw[,41] > Evw[,21]))/iters
# we conclude that 87.6% of our sites with 4 wolves will have more veg 
# than sites with 2 wolves...








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's run our first path analysis in a ML framework, we're just doing this to
# demonstrate that the lavaan package exists (i.e., you don't have to use Bayesian programs to 
# analyze data with SEMs). We'll also use lavaan's "cousin", 'blavaan' 
# to help us write models in the future (blavaan is short for Bayesian lavaan... go figure)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to do that we'll use 'lavaan', a package that allows for the construction of
# structural equation models. Lavaan can obtain ML parameter estimates for path coefficients
# produce cool figures, etc., but we're going to stick in BUGS-code land because lavaan
# CAN'T handle: 1) link functions, 2) joint analysis with other (e.g., capture-mark-recapture)
# data types, or 3) missing data.
data <- data.frame(w,u,v)
path1 <- sem(model = 
               'u ~ 1 + a1 * w   
              v ~ 1 + b1 * u
              vw := a1*b1',   # this allows us to estimate the indirect effect of wolves on veg
             # also note we can specify parameter names!
             data = data)
print(summary(path1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we can dig into the lavaan object (path1) and pull out some parameter estimates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path1@ParTable$label
# this is the product of our two path coefficients (alpha[1] & beta[1])
path1@ParTable$est[2] * path1@ParTable$est[4]
# this is our estimate of the indirect effect?! 
path1@ParTable$est[9]

# this is all the same as the code above, lavaan can just do it for us















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here is a JAGS model for the structural equation model parameterized with normal distributions:
# wolves -> ungulates -> vegetation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m_path_identity.jags")
cat("
    model {

    # priors for everything assigned in a loop to minimize unnecessary coding mistakes

    alpha[1] ~ dnorm(40, 0.1)
    alpha[2] ~ dnorm(0, 0.01)
    for (j in 1:2){
      beta[j] ~ dnorm(0, 0.01)
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])
    }

    for (i in 1:n){
      u[i] ~ dnorm(alpha[1] + alpha[2] * w[i], tau[1])
      v[i] ~ dnorm(beta[1] + beta[2] * u[i], tau[2])
    }

    }
    ",fill = TRUE)
sink()

jags.data <- list(n = n, w = w, u = u, v = v)
inits <- function(){list()}  
parameters <- c('alpha','beta','sigma')

# number of chains (nc), thinning rate (nt), number of iterations (ni), and number to burn-in
nc <- 4
nt <- 10
ni <- 20000
nb <- 10000



library(jagsUI)
Sys.time()
m <- jags(jags.data, inits, parameters, "m_path_identity.jags", parallel = T, 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
Sys.time()

print(m)

