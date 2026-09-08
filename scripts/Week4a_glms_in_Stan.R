# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# First, we'll use exactly the same simulation code we built to analyze
# the alligator mass ~ legnth regression relationship, but now we'll use a 
# log-normal distribution in our model
#
# Second, we'll use a logit-link to include a covariate for the probability
# of band-recovery probability using our duck banding example given a declining 
# number of duck hunters...
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First, let's simulate some alligator mass and length data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(reshape2)
library(vioplot)
library(jagsUI)
library(rstan)

set.seed(1234)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# base plot settings
# (1x1 grid, 'times new roman-ish' font, large margins on x- and y-axes)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,1), family = 'serif', mar = c(5.1,5.1,2.1,2.1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate data for linear model demonstration
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# sample size (n)
n <- 200

# simulate length (x)
x <- rlnorm(n, 1.75, 0.15)
hist(x, breaks = 50, main = NULL, xlab = 'Alligator length (ft)',
     col = 'forestgreen', border = 'forestgreen', cex.lab = 2)

# true regression parameters on the log scale
beta0 <- 1.6
beta1 <- 0.375

# simulate mass (y) as a function of length (x)
y <- rlnorm(n, beta0 + beta1 * x, 0.1)
plot(y ~ x, ylab = 'Mass (lbs)', xlab = 'Length (ft)',
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)

# z-standardize mass as a 
z <- as.numeric(scale(x))



m1.stan <- stan_model(model_code = 
                        "                        
data {
  int n;            
  vector[n] y;
  vector[n] x;
}

parameters {
  real<lower=0> sigma;
  real beta0;
  real beta1;
}

transformed parameters {
  vector[n] E;
  E = beta0 + beta1 * x;
}

model {
  beta0 ~ normal(0, 100);
  beta1 ~ normal(0, 100);

  y ~ lognormal(E, sigma);


}
")




#############################################################################################
# Bundle data
#############################################################################################


stan_data <- list(y = y, n = length(y), x = x)

inits <- function(){list()}  

nc <- 4
nt <- 25
ni <- 5000
nb <- 2500


Sys.time()
m1 <- sampling(m1.stan, 
               data = stan_data,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)


Sys.time()

print(m1)
m1 <- extract(m1)

lci <- NULL
uci <- NULL
for (j in 1:n){
  lci[j] <- exp(quantile(m1$E[,j], c(0.025)))
  uci[j] <- exp(quantile(m1$E[,j], c(0.975)))
}
plot(colMeans(exp(m1$E)) ~ x, ylim = c(0,200), xlim = c(0,10),
     ylab = 'Expected mass (lbs)',
     xlab = 'Length', cex.lab = 2)
arrows(x, lci, x, uci, length = 0, lty = 2, lwd = 1.5)
points(exp(colMeans(m1$E)) ~ x, pch = 21, bg = 'forestgreen', cex = 1.25)









# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework PART I
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) What does our prior for beta0 represent here? How would this look
#    on the scale of our data? i.e.,
#    is this a good prior? why or why not? how much of our prior support is less than 1?, i.e.,
#    what might a better prior look like? What does beta0 represent in your own words?
#
# First, let's imagine we do have a pretty good idea what an alligator of average (5.7 ft)
# weighs, and that's about 45 lbs...
# we could back-transform this value to the real (-Inf, Inf) scale that our model 
# will work on quite easily
# so there's that number! 3.806662
# 
# now let's imagine we want some uncertainty around that, perhaps anywhere from 
# 25 to 100 pounds as 95% credible intervals on the prior?
hist(rlnorm(10000, log(45), 1)) # too broad
hist(rlnorm(10000, log(45), 0.01)) # way too informative :)
hist(rlnorm(10000, log(45), 0.1))  # maybe too informative?
hist(rlnorm(10000, log(45), 0.25)) # not so bad
# we could even use what we know about the median and variance of a log-normal
# to try to pick values that would correspond to our desired prior?
# https://en.wikipedia.org/wiki/Log-normal_distribution
#
#
# 2) What about our prior for beta1? What does that represent, are all of these values
#    reasonable? How might we modify this prior to be more reflective of our 
#    expectations and prior knowledge rather than an incredibly vague statement
#    such as normal(0, sigma^2 = 100000)
#
# 3) There is code below to make a plot of our expectations of mass given 
#    different values of length. It's very similar to the plot we used in
#    week 3 when reviewing linear models. How (and why) does it differ?
#    Review the use of the melt function...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Second, let's simulate some ducks being marked, released, and eventually
# shot, recovered, and reported by hunters...
# 
# We'll have a single covariate, the number of hunters. We might reasonably
# assume that there is a causal linkage between the number of hunters and 
# the proportion of ducks that are shot (note: if we don't at least test this 
# hypothesis we should probably just all throw in the towel on this whole 
# 'applied science' thing!) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# number of time periods, or years
nT <- 20
# number of releases each year [mark release 1000 ducks every year]
nR <- rep(1000, nT)
# number of hunters
nH <- rpois(nT, exp(7 + -0.125 * 1:nT))
plot(nH)

# z-standardize the number of hunters
z <- as.numeric(scale(nH))

# build logit-link model to estimate
beta0 <- -2
beta1 <- 0.25
f <- plogis(beta0 + beta1 * z)
plot(f ~ nH)

# simulate the number of actual recoveries and plot 
# number of recoveries divided by number of releases...
y <- rbinom(nT, nR, f)
plot(c(y/nR), ylab = 'ML band-recovery probability', xlab = 'Year', cex.lab = 2)



m2.stan <- stan_model(model_code = "                        
data {
  int n;            
  array[n] int y;
  array[n] int R;
  vector[n] z;
}

parameters {
  real beta0;
  real beta1;
}

transformed parameters {
  vector[n] f;
  f = inv_logit(beta0 + beta1 * z);
}

model {
  beta0 ~ logistic(0,1);
  beta1 ~ normal(0, 10);

  y ~ binomial(R, f);


}
")

stan_data <- list(y = y, n = length(y), z = z, R = nR)
inits <- function(){list()}  
nc <- 4
nt <- 25
ni <- 5000
nb <- 2500

Sys.time()
m2 <- sampling(m2.stan, 
               data = stan_data,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)


Sys.time()

print(m2)
summary(m2)
m2 <- extract(m2)


lci <- NULL
uci <- NULL
for (j in 1:n){
  lci[j] <- (quantile(m2$f[,j], c(0.025)))
  uci[j] <- (quantile(m2$f[,j], c(0.975)))
}

plot(colMeans(m2$f) ~ seq(1,nT), xlab = 'Year', ylab = 'Band-recovery probability',
     las = 1, cex.axis = 1.5, cex.lab = 2, ylim = c(0,0.2))
arrows(seq(1,nT), lci, seq(1,nT), uci, length = 0, lty = 2)
points(colMeans(m2$f) ~ seq(1,nT), pch = 21, bg = 'forestgreen', cex = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework PART II
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) What does our prior for beta0 represent here? How would this look
#    on the scale of our data? i.e.,
#    is this a good prior? why or why not? 
#    what might a better prior look like? What does beta0 represent in your own words?
#
# 2) What about our prior for beta1? What does that represent, are all of these values
#    reasonable? How might we modify this prior to be more reflective of our 
#    expectations and prior knowledge rather than an incredibly vague statement
#    such as normal(0, sigma^2 = 100000)
#
# 3) attempt to write code below to make a density plot 
#    along with credible intervals of expected values of band-recovery
#    probability (or the expected number of recoveries... or both?!) as a function
#    of the number of hunters
#
# 4) If you're really feeling wild, try to put a random effect on band-recovery
#    probability around the effect of the long-term decline in hunters?!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




