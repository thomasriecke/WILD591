# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulation:
#
# Let's imagine a dataset. Our ultimate response variable will be ungulate density (u)
# in valley bottoms across the state of Montana. Please note these data
# are completely simulated!
#
# At each site, we'll also have estimates of wolf (w) and human (u)
# density. Our hypotheses are that:
# 1) wolf density will be negatively affected by humans (alpha_1)
# 3) ungulate density will be negatively affected by 
# both humans (beta_1) and wolves (beta_2)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
library(lavaan)
library(vioplot)
library(blavaan)
library(piecewiseSEM)


set.seed(12345)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
# 100 sites
n <- 100
# resolution for plots
res <- 100




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate human density (h)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h <- rlnorm(n, 1, 0.75)
hist(h, xlab = 'Human density (h)', breaks = 100, main = NULL, cex.lab = 2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate wolf density (w) as a function of humans (h)
#
# note that wolf density is high but variable in areas of low human density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha0 <- 2
alpha1 <- -0.75
w <- rlnorm(n, alpha0 + alpha1 * h, 0.25)

# histogram
hist(w, xlab = 'Wolf density', 
     breaks = 25, main = NULL, cex.lab = 2)

# scatterplot and P(Wolf | human)
plot(w ~ h, ylab = 'Wolf density (w)', xlab = 'Human density (h)', cex.lab = 2)
lines(exp(alpha0 + alpha1 * seq(min(h),max(h),length.out = res)) ~ seq(min(h),max(h),length.out = res), lwd = 2)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate ungulate density (u) as a function of human density (h), wolf density (w),
# and net primary productivity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta0 <- 3
beta1 <- -0.15
beta2 <- -0.05

u <- rlnorm(n, beta0 + beta1 * w + beta2 * h, 0.1)
hist(u, xlab = 'Ungulate density', 
     breaks = 25, main = NULL, cex.lab = 2)

# now let's plot ungulate density as a function of wolves
plot(u ~ w, ylab = 'Ungulate density', xlab = 'Wolf density', cex.lab = 2)

# now let's plot ungulate density as a function of humans
plot(u ~ h, ylab = 'Ungulate density', xlab = 'Human density', cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# predict wolf abundance as a function of human abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a sequence of values for human density
seq.h <- seq(0,25,length.out = res)
# predict wolves as a function of humans
pw <- exp(alpha0 + alpha1 * seq.h) 
# predict ungulates as a function of wolves given humans (beta1) and humans (beta2)
pu <- exp(beta0 + beta1*pw + beta2 * seq.h)

# now let's plot ungulate density as a function of humans (which will also be a function of wolves)
plot(u ~ h, ylab = 'Ungulate density', xlab = 'Human density', cex.lab = 2)
lines(pu ~ seq.h)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK PART I
#
# 1) The last few plots (regressing ungulate density against human and wolf density)
#    appear quadratic? i.e., ungulate density is highest at intermediate levels
#    of wolf or human density?
#
#    Why? We didn't simulate any quadratic relationships. What's happening here?!
#
#    Make sure you understand this clearly. If it's not immediately intuitive, no worries,
#    let's talk about it!
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First let's apply multivariate regression
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor(h,w) # below |r| > 0.7, proceed as you wish :)
summary(m.glm <- glm(log(u) ~ w + h))
# here we get the correct answers for the direct effects. but there's no indirect effect
# what effect does that have on our inference?

plot(u ~ h, ylab = 'Ungulate density', xlab = 'Human density', cex.lab = 2, ylim = c(0,25))
lines(exp(m.glm$coefficients[1] + m.glm$coefficients[3] * seq.h) ~ seq.h)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Second, let's explore lavaan syntax
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- data.frame(h,w,u)
m.lavaan <- sem(model = '
                w ~ 1 + h
                u ~ 1 + w + h', data = data)
summary(m.lavaan)


plot(w ~ h, ylab = 'Wolf density (w)', xlab = 'Human density (h)', cex.lab = 2)
lines(c(2.284 + -0.247 * seq.h) ~ seq.h)
# here we're able to account for direct and indirect effects, but lavaan still
# doesn't incorporate 'logit' link functions, etc.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Third, let's explore piecewiseSEM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m.psem <- psem(
  glm(log(w) ~ h, data = data),
  glm(log(u) ~ h + w, data = data),
  data = data
)
summary(m.psem)
# we recover the data-generating parameter values. This package is incredible,
# and we'll use it a bit more as we progress, but retaining a focus on Bayesian methods
# so that we can link SEMs to other types of analyses (e.g. capture-recapture, occupancy)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fourth, let's explore blavaan!
# note you can specify the number of chains, burn-in, even priors, but the real
# value of this (from my perspective) is that it will spit out a working JAGS
# or Stan model if you specify mcmcfile = T
# https://rdrr.io/cran/blavaan/man/bsem.html
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(runjags)
m.blavaan <- bsem(model = '
                w ~ 1 + h
                u ~ 1 + w + h', data = data, 
                target = 'jags',
                mcmcfile = T)
summary(m.blavaan)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if you don't have experience building models in the BUGS or Stan language, using these
# types of packages is an absolutely wonderful way to immediately build working code
# using 'lm()-like' model syntax
#
# you can then carefully alter variable names, adjust priors, and be on your merry way!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# search for lavExport on your computer (this will be a folder in your wd() or Documents)
getwd()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# JAGS model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(jagsUI)
sink("m_path.jags")
cat("
    model {

    # priors for everything assigned in a loop to minimize unnecessary coding mistakes
    for (j in 1:2){
      alpha[j] ~ dnorm(0, 0.1)
      sigma[j] ~ dgamma(1,1)
      tau[j] = 1/(sigma[j] * sigma[j])      
    }

    for (j in 1:3){
      beta[j] ~ dnorm(0, 0.1)
    }

    # the model is nearly identical to Tuesday's but now 
    for (i in 1:n){
      w[i] ~ dlnorm(alpha[1] + alpha[2] * h[i], tau[1])
      u[i] ~ dlnorm(beta[1] + beta[2] * w[i] + beta[3]*h[i], tau[2])
    }

    }
    ",fill = TRUE)
sink()

jags.data <- list(n = n, w = w, u = u, h = h)
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

print(m)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's make a plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iters <- length(m$sims.list$alpha[,1]) # we saved 4000 iterations
res <- 100

xh <- seq(0,20, length.out = res)
Ew <- matrix(NA, iters, res)
Eu <- matrix(NA, iters, res) 

qEw <- matrix(NA, res, 5)
qEu <- matrix(NA, res, 5)


for (j in 1:res){
  Ew[,j] <- exp(m$sims.list$alpha[,1] + m$sims.list$alpha[,2] * xh[j])  # expected value of wolves given humans
  Eu[,j] <- exp(m$sims.list$beta[,1] + m$sims.list$beta[,2] * Ew[,j] + m$sims.list$beta[,3] * xh[j])   # expected value of ungulates given wolves and humans
  
  qEw[j, ] <- quantile(Ew[,j], c(0.025,0.05,0.5,0.95,0.975))
  qEu[j, ] <- quantile(Eu[,j], c(0.025,0.05,0.5,0.95,0.975))
}


library(reshape2)
mEw <- melt(Ew); names(mEw) <- c('iter','x','w')
mEu <- melt(Eu); names(mEu) <- c('iter','x','u')


# plot of expected ungulate density as a function of human density
smoothScatter(mEu$u ~ xh[mEu$x], las = 1, nrpoints = 0,
              ylab = 'Expected density of ungulates',
              xlab = 'Density of humans', ylim = c(5,20))
lines(qEu[,3] ~ xh, lty = 1, lwd = 3, col = 'white')
lines(qEu[,1] ~ xh, lty = 2, lwd = 3, col = 'white')
lines(qEu[,5] ~ xh, lty = 2, lwd = 3, col = 'white')
points(u ~ h)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK PART II
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) above we made a plot demonstrating the predicted effects of changes in humans
#    on changes in ungulates. Could you use components we made to generate that plot to plot the effect
#    of humans on wolves?
#
# 2) What is the probability that an increase in density from 2 to 5 humans per [area] would lead to 
#    an increase in ungulates? 
# 
#    Note: to do that we might have to incorporate our estimates of
#    uncertainty as well! Reexamine your script from Tuesday to 
#    draw the expected values from a lognormal distribution
#
#  3) What would predicted densities of ungulates be if wolves were extirpated from Montana?
#     (at least in the short-term)
#
#  4) Below is some script (accompanied by slides in the .pdf of the lecture) that adds an
#    effect of gross primary productivity on human density, as well as ungulate density.
#    Run the script, and examine the relationships. Try changing a few of them, and perhaps
#    computing the total effect of GPP (see lines 79-91 above)
#
# 5) If you're really feeling 'wild:' delete all the script in the JAGS model.
#    Recreate the model. We'll start doing pieces of models in the next week or so, but 
#    there's no time like the present to break some code.
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some supplementary script (see additional slides) to envision how pathways 
# can work with more variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we'll add one additional variable, gross primary productivity, that will inform
# both human abundance and deer abundance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set the seed
# set some base plotting rules
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12345)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'sans')
# 100 sites
n <- 100
# resolution for plots
res <- 100

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gross primary productivity (x) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we might expect that this would affect deer density, and even human density
# (i.e., humans tend to live in nice beautiful valleys near agua). I assume
# it would also affect wolf density indirectly via it's effects on prey abundance
x <- rnorm(n, 0, 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate human density (h) as a function of NPP (x)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha0 <- 2
alpha1 <- 0.5
h <- rlnorm(n, alpha0 + alpha1 * x, 0.5)

# scatterplot and P(human | GPP)
plot(h ~ x, xlab = 'Gross primary productivity (x)', ylab = 'Human density (h)', cex.lab = 2)
lines(exp(alpha0 + alpha1 * seq(min(x),max(x),length.out = res)) ~ seq(min(x),max(x),length.out = res), lwd = 2)
hist(h, xlab = 'Human density', breaks = 100, main = NULL, cex.lab = 2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate wolf density (w) as a function of humans (h)
#
# note that wolf density is high but variable in areas of low human density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta0 <- 3
beta1 <- -0.75
w <- rlnorm(n, beta0 + beta1 * h, 0.25)

# histogram
hist(w, xlab = 'Wolf density', 
     breaks = 25, main = NULL, cex.lab = 2)

# scatterplot and P(Wolf | human)
plot(w ~ h, ylab = 'Wolf density', xlab = 'Human density', cex.lab = 2)
lines(exp(beta0 + beta1 * seq(min(h),max(h),length.out = res)) ~ seq(min(h),max(h),length.out = res), lwd = 2)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate ungulate density (u) as a function of human density (h), wolf density (w),
# and net primary productivity
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gamma0 <- 4
gamma1 <- 0.25
gamma2 <- -0.15
gamma3 <- -0.15
u <- rlnorm(n, gamma0 + gamma1 * x + gamma2 * h + gamma3 * w, 0.1)
hist(u, xlab = 'Ungulate density', 
     breaks = 25, main = NULL, cex.lab = 2)

# now let's plot ungulate density as a function of wolves
plot(u ~ w, ylab = 'Ungulate density', xlab = 'Wolf density', cex.lab = 2)

# now let's plot ungulate density as a function of humans
plot(u ~ h, ylab = 'Ungulate density', xlab = 'Human density', cex.lab = 2)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# 1) The last two plots (regressing ungulate density against human and wolf density)
#    appear somewhat quadratic? i.e., ungulate density is highest at intermediate levels
#    of wolf or human density?
#
#    Why? We didn't simulate any quadratic relationships. What's happening here...
#
# 2) Let's look at two more plots, where we'll regress ungulate abundance against net primary
#    productivity (x), as well as wolf abundance against net primary productivity.
#    note that you'll have to simulate some variation in x (NPP) to have these work effectively
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# now let's plot ungulate density as a function of net primary productivity
plot(u ~ x, ylab = 'Ungulate density', xlab = 'z(GPP)', cex.lab = 2)


# now let's plot ungulate density as a function of net primary productivity
plot(w ~ x, ylab = 'Wolf density', xlab = 'z(GPP)', cex.lab = 2)
# from this simple univariate relationship we might conclude that net primary productivity is bad for wolves!?
# why?! do you think that's ecologically accurate? what is driving this?















