# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# custom sampler (!)
# jump-in-and-let-your-feet-stick-out
# P(theta | y) = [P(y | theta) * P(theta)] / P(y)
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# (1) SIMULATE DATA
# (2) CUSTOM SAMPLER ONE (default-ish prior)
# (3) CUSTOM SAMPLER TWO (biologically informed prior)
# (4) CUSTOM SAMPLER THREE (playing with tuning parameter)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(1234) # so the numbers in the lecture are consistent
library(latex2exp)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (1) simulate data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# alpha is the mean pack size (mu = 6) on the log scale
# n is the number of independent counts
# y is the counts of wolf pack size
#
# please note that we didn't use truncated distributions here, (i.e.,
# there are some 0's and 1's, which are biologically impossible). it's to keep
# the math simpler.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
alpha <- log(6)
n <- 500
y <- rpois(n, exp(alpha))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ML estimate and different ways to approximate conf. intervals
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(ml <- glm(y ~ 1, family = 'poisson'))
exp(ml$coefficients)
exp(confint(ml))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
table(y)
par(mfrow = c(1,1), mar = c(5,5,2,2))
barplot(table(y), main = NULL, xlab = 'Counts (y)', las = 1, cex.lab = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summary statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mean(y)
var(y)
sd(y)
sd(y)/sqrt(n)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log link
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
a <- seq(-3, 3, length.out = 100)
plot(exp(a) ~ a, ylab = TeX("$\\e^{a}$"), cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (2) write a custom MCMC sampler! with a 'bad' prior
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize our prior
par(mar = c(5,5,2,2), mfrow = c(1,1))
hist(exp(rnorm(100000, 0, 1)), breaks = 100,
     main = TeX("$e^{N(0,1)}$"), xlab = '')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-likelihood
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LL <- function(theta, y) {
  tmp <- theta
  LLi <- dpois(y, exp(tmp), log=TRUE) # LL contribution of each datum i
  LL <- sum(LLi) # LL for all observations in the data vector y
  return(LL)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-posterior function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.posterior <- function(theta, y){
  loglike <- LL(theta, y)
  logprior <- dnorm(theta, mean = 0, sd = 1, log = TRUE)
  return(loglike + logprior)
}


# log-likelihoods of 1 and 2
LL(2, y)
LL(1, y)

# log-priors of 1 and 2
dnorm(2, mean = 0, sd = 1, log = TRUE)
dnorm(1, mean = 0, sd = 1, log = TRUE)

# log-posteriors of 1 and 2
LL(2, y) + dnorm(2, mean = 0, sd = 1, log = TRUE)
LL(1, y) + dnorm(1, mean = 0, sd = 1, log = TRUE)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampling objects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp <- NULL
potential.theta <- seq(0,3, length.out = n)
for (i in 1:n){
  tmp[i] <- LL(potential.theta[i], y)
}
plot(tmp ~ exp(potential.theta), ylab = TeX("l($\\theta$ | y)"), 
     xlab = TeX("$e^{\\theta}$"), cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampling objects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iterations <- 10000
tuning <- 0.05
# initial value for theta
theta <- 1
acc <- rep(0, iterations)
posterior <- rep(NA, iterations)
candidate <- rep(NA, iterations)
logpost.candidate <- rep(NA, iterations)
likelihood <- rep(NA, iterations)
log.prior <- rep(NA, iterations)
logpost.current <- log.posterior(theta, y)


LL(theta, y) + dnorm(theta, 0, 1, log = T)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampler
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (i in 1:iterations){
  # draw a candidate value from the starting value given tuning parameter
  candidate[i] <- rnorm(1, theta, tuning)
  # evaluate log(posterior) for candidate value
  logpost.candidate[i] <- log.posterior(candidate[i], y)
  # compute metropolis acceptance ratio r
  r <- exp(logpost.candidate[i] - logpost.current)
  # if we accept, adjust current value, log-posterior, and note acceptance
  if (runif(1) < r){
    theta <- candidate[i]
    logpost.current <- logpost.candidate[i]
    acc[i] <- 1
  }
  posterior[i] <- theta
  # likelihood 
  likelihood[i] <- LL(posterior[i], y)
  log.prior[i] <- dnorm(theta, mean = 0, sd = 1, log = TRUE)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# entire chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(posterior)
plot(posterior[2000:2050] ~ seq(2000,2050), ylab = TeX("$\\theta$"), xlab = 'Iteration',
     cex.lab = 2, las = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# burn post-sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
burn <- 1000
posterior <- posterior[(burn+1):iterations]
acc <- acc[(burn+1):iterations]
logpost.candidate <- logpost.candidate[(burn+1):iterations]
candidate <- candidate[(burn+1):iterations]
likelihood <- likelihood[(burn+1):iterations]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summary stats
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basic summary
summary(posterior)
# 2.5%, 25%, 50%, 75%, and 97.5% Bayesian credible intervals
quantile(posterior, c(0.025,0.25,0.5,0.75,0.975))
# what proportion of draws were accepted? (25-40% is a good rule; Gelman et al. 2014)
summary(acc)

plot(posterior, ylab = TeX("$\\theta$"), cex.lab = 1.5, type = 'l')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot all proposed values (black) vs. 
# accepted values (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(exp(candidate), ylab = TeX("$e^{\\theta}$"), cex.lab = 1.5)
points(exp(posterior), col = 'red')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot likelihood vs. parameter estimate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(likelihood ~ posterior)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# posterior after burnin
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(posterior, type = 'l', 
     ylab = TeX("$\\theta$"), xlab = 'iteration', cex.lab = 1.5, las = 1)

hist(posterior, breaks = 250, main = '',
     xlab = TeX("$\\theta$"), cex.lab = 1.5)










# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (3) But wait, there's more!!!!!!!!!!!!!!!!
# we could change our prior to be more biologically reasonable given previous
# research, e.g., Sells et al. (2022; JWM) found that mean
# wolf pack size in Montana was b/w 5 & 6 wolves
# https://wildlife.onlinelibrary.wiley.com/doi/pdf/10.1002/jwmg.22193
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# visualize our prior, that's better
par(mar = c(5,5,2,2), mfrow = c(1,1))
hist(exp(rnorm(100000, log(5.5), 0.5)), breaks = 100,
     main = TeX("$e^{N(ln(5.5), 0.5)}$"), xlab = '')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-likelihood
# this didn't change at all..
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LL <- function(theta, y) {
  tmp <- theta
  LLi <- dpois(y, exp(tmp), log=TRUE) # LL contribution of each datum i
  LL <- sum(LLi) # LL for all observations in the data vector y
  return(LL)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-posterior function, 
# note that we changed the prior given Sells et al. (2022)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.posterior <- function(theta, y){
  loglike <- LL(theta, y)
  logprior <- dnorm(theta, mean = log(5.5), sd = 0.5, log = TRUE)
  return(loglike + logprior)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's visualize the two priors we've used so far,
# one of these seems more reasonable than the other...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prior1 <- rnorm(10000, 0, 1)
prior2 <- rnorm(10000, log(5.5), 0.5)
par(mfrow = c(1,2), mar = c(5,5,2,2))
hist(exp(prior1), breaks = 250, xlab = '', main = TeX("$e^{N(0,1)}$"))
hist(exp(prior2), breaks = 250, xlab = '', main = TeX("$e^{N(ln(5.5),0.5)}$"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampling objects
# none of this changes...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iterations <- 10000
tuning <- 0.05
theta <- rnorm(1)
acc <- rep(0, iterations)
posterior <- rep(NA, iterations)
candidate <- rep(NA, iterations)
logpost.candidate <- rep(NA, iterations)
likelihood <- rep(NA, iterations)
log.prior <- rep(NA, iterations)
logpost.current <- log.posterior(theta, y)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampler
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (i in 1:iterations){
  # draw a candidate value from the starting value given tuning parameter
  candidate[i] <- rnorm(1, theta, tuning)
  # evaluate log(posterior) for candidate value
  logpost.candidate[i] <- log.posterior(candidate[i], y)
  # compute metropolis acceptance ratio r
  r <- exp(logpost.candidate[i] - logpost.current)
  # if we accept, adjust current value, log-posterior, and note acceptance
  if (runif(1) < r){
    theta <- candidate[i]
    logpost.current <- logpost.candidate[i]
    acc[i] <- 1
  }
  posterior[i] <- theta
  # likelihood 
  likelihood[i] <- LL(posterior[i], y)
  log.prior[i] <- dnorm(theta, mean = log(5.5), sd = 0.5, log = TRUE)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# entire chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(posterior)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# burn post-sampling
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
burn <- 1000
posterior <- posterior[(burn+1):iterations]
acc <- acc[(burn+1):iterations]
logpost.candidate <- logpost.candidate[(burn+1):iterations]
candidate <- candidate[(burn+1):iterations]
likelihood <- likelihood[(burn+1):iterations]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# summary stats
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# basic summary
summary(posterior)
# 2.5%, 25%, 50%, 75%, and 97.5% Bayesian credible intervals
quantile(posterior, c(0.025,0.25,0.5,0.75,0.975))
# what proportion of draws were accepted? (25-40% is a good rule; Gelman et al. 2014)
summary(acc)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot all proposed values (black) vs. 
# accepted values (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(candidate)
points(posterior, col = 'red')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot likelihood vs. parameter estimate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(likelihood ~ posterior)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(posterior, type = 'l', 
     ylab = TeX("$\\theta$"), xlab = 'iteration', cex.lab = 1.5, las = 1)

hist(posterior, breaks = 250, main = '',
     xlab = TeX("$\\theta$"), cex.lab = 1.5)









# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# (4) write a custom MCMC sampler! and change the tuning parameters?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-likelihood
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
LL <- function(theta, y) {
  tmp <- theta
  LLi <- dpois(y, exp(tmp), log=TRUE) # LL contribution of each datum i
  LL <- sum(LLi) # LL for all observations in the data vector y
  return(LL)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# log-posterior function
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
log.posterior <- function(theta, y){
  loglike <- LL(theta, y)
  logprior <- dnorm(theta, mean = 0, sd = 1, log = TRUE)
  return(loglike + logprior)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampling objects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iterations <- 10000
tuning <- 0.15
# initial value for theta
theta <- 1
acc <- rep(0, iterations)
posterior <- rep(NA, iterations)
candidate <- rep(NA, iterations)
logpost.candidate <- rep(NA, iterations)
likelihood <- rep(NA, iterations)
log.prior <- rep(NA, iterations)
logpost.current <- log.posterior(theta, y)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sampler
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
for (i in 1:iterations){
  # draw a candidate value from the starting value given tuning parameter
  candidate[i] <- rnorm(1, theta, tuning)
  # evaluate log(posterior) for candidate value
  logpost.candidate[i] <- log.posterior(candidate[i], y)
  # compute metropolis acceptance ratio r
  r <- exp(logpost.candidate[i] - logpost.current)
  # if we accept, adjust current value, log-posterior, and note acceptance
  if (runif(1) < r){
    theta <- candidate[i]
    logpost.current <- logpost.candidate[i]
    acc[i] <- 1
  }
  posterior[i] <- theta
  # likelihood 
  likelihood[i] <- LL(posterior[i], y)
  log.prior[i] <- dnorm(theta, mean = 0, sd = 1, log = TRUE)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# entire chain
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(posterior, ylab = TeX("$\\theta$"), cex.lab = 1.5, type = 'l')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot all proposed values (black) vs. 
# accepted values (red)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(exp(candidate), ylab = TeX("$e^{\\theta}$"), cex.lab = 1.5)
points(exp(posterior), col = 'red')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot rolling mean of acceptance rate
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

chunk <- 100
acc.rate <- rep(NA, iterations - chunk)

for (i in (chunk+1):(iterations)){
  acc.rate[i - chunk] <- mean(acc[(i-chunk):(i-1)])
}


par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(acc.rate, type = 'l', ylab = 'Mean acceptance rate for previous 100 samples')

plot(exp(candidate), ylab = TeX("$e^{\\theta}$"), cex.lab = 1.5, ylim = c(3,8))
points(exp(posterior), col = 'red', type = 'l')

