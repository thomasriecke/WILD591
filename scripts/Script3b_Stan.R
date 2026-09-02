# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# An introduction to categorical covariates, 
# in which we explore fixed and random effects, play with 'Bergmann's rule' a bit,
# and think about how random effects behave
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(1234)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1), family = 'serif')
library(jagsUI)
library(vioplot)
library(lme4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# randomly sample 'l' locations, 
# note to fix: pull the actual Red fox home range and generate realistic coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
l <- 30
# spatial coordinates (unused at at the moment, but may rework to demonstrate
# spatially correlated random effects)
s <- matrix(NA, l, 2) # number of rows = l, lat+long = 2 columns
s[,1] <- runif(l,-120,-80)
s[,2] <- runif(l, 35, 55)

# let's extract latitude from the spatial coordinates (x)
x <- as.numeric(scale(s[,2]))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# randomly sample 'n' foxes from each population and assign a group id
# and latitude to each fox...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lambda will represent the average number of samples
lambda <- 14
# we'll simulate the actual number of samples from a poisson distribution
n <- rpois(l, lambda)


# now let's generate the average size of an adult female red fox (mu.star)
# as well as the variance among populations (sigma[1])
mu.star <- 12
sigma <- NULL
sigma[1] <- 0.5

# now we'll simulate each populations mean as a function of mu.star, 
# sigma[1], and an effect of latitude on population size (beta)
beta <- 0.5
mu <- rnorm(l, mu.star + beta * x, sigma[1])

# let's plot these population level means as a function of z-standardized latitude
plot(mu ~ x, xlab = 'Z-standardized latitude', 
     ylab = 'Population mean mass (lbs)', cex.lab = 2)

# let's plot these population level means as a function of actual latitude
plot(mu ~ s[,2], xlab = "Latitude (N\u00B0)", 
     ylab = 'Population mean mass (lbs)', cex.lab = 2)


# simulate a population identifier for each fox
sum(n)
p <- rep(seq(1:l), times = n)
p

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's simulate each foxes mass
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sigma[2] <- 0.5
y <- rnorm(length(p), mu[p], sigma[2])

# make a boxplot of fox mass as a function of population
boxplot(y ~ p, ylab = 'Female fox mass (lbs)', xlab = 'Population id', cex.lab = 2)

# regress fox mass against latitude
plot(y ~ s[,2][p], xlab = "Latitude (N\u00B0)",
     ylab = 'Female fox mass (lbs)', pch = 21, bg = 'red2', cex.lab = 2, cex = 1.5)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's run two ML models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m1 <- lm(y ~ as.factor(p)))       # fixed effects
summary(m2 <- lmer(y ~ (1|as.factor(p)))) # random effects















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Specify random effects model in Stan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rstan)
m3.stan <- stan_model(model_code = 
"                        
data {
  int n;            
  int l;               
  int p[n];         
  vector[n] y;           
}

parameters {
  real<lower=-100, upper = 100> mu_star; 
  real<lower=0> sigma[2];
  real mu[l];
}

model {

  mu ~ normal(mu_star, sigma[1]);
  sigma ~ gamma(1,1);
  mu_star ~ normal(12, 3.16);
  for (i in 1:n){
    y[i] ~ normal(mu[p[i]], sigma[2]);
  }

}
")




#############################################################################################
# Bundle data
#############################################################################################


stan_data <- list(y = y, l = l, n = length(p), p = p)

inits <- function(){list()}  

nc <- 4
nt <- 25
ni <- 5000
nb <- 2500


Sys.time()
m3 <- sampling(m3.stan, 
               data = stan_data,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)


Sys.time()

print(m3)
m3 <- extract(m3)

# this is the mean of the posterior distribution of the grand mean of female fox mass
mean(m3$mu_star)
# this is the median of the posterior distribution of of the grand mean of female fox mass
quantile(m3$mu_star, c(0.5))
# this is the lower 95% credible interval of the posterior distribution of of the grand mean of female fox mass
quantile(m3$mu_star, c(0.025))
# this is the upper 95% credible interval of the posterior distribution of of the grand mean of female fox mass
quantile(m3$mu_star, c(0.975))


lc.mu <- NULL
uc.mu <- NULL
for (j in 1:l){
  lc.mu[j] <- quantile(m3$mu[,j], c(0.025))
  uc.mu[j] <- quantile(m3$mu[,j], c(0.975))
}
plot(colMeans(m3$mu), ylim = c(9,15), xlim = c(0,32),
     ylab = 'Mean mass (lbs)',
     xlab = 'Population', cex.lab = 2)
arrows(1:l, lc.mu, 1:l, uc.mu, length = 0, lty = 2, lwd = 1.5)
points(colMeans(m3$mu), pch = 21, bg = 'red2', cex = 2.5)

points(y ~ p, pch = 1, cex = 0.5)
points(mu, pch = 1, cex = 1.5) # the actual means

































# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Specify random effects model in Stan
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(rstan)
m4.stan <- stan_model(model_code = 
                        "                        
data {
  int n;            
  int l;               
  int p[n];         
  vector[n] y;           
}

parameters {
  // real<lower=-100, upper = 100> mu_star; 
  real<lower=0> sigma[2];
  real mu[l];
}

model {

  mu ~ normal(12, 3.16);
  sigma ~ gamma(1,1);
  for (i in 1:n){
    y[i] ~ normal(mu[p[i]], sigma[2]);
  }

}
")



stan_data <- list(y = y, l = l, n = length(p), p = p)
inits <- function(){list()}  
nc <- 4
nt <- 25
ni <- 5000
nb <- 2500

Sys.time()
m4 <- sampling(m4.stan, 
               data = stan_data,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)
Sys.time()

print(m4)
m4 <- extract(m4)


lc <- NULL
uc <- NULL
for (j in 1:l){
  lc[j] <- quantile(m4$mu[,j], c(0.025))
  uc[j] <- quantile(m4$mu[,j], c(0.975))
}

plot(colMeans(m4$mu), ylim = c(9,15),
     ylab = 'Mean mass (lbs)',
     xlab = 'Population', cex.lab = 2)
arrows(1:l, lc, 1:l, uc, length = 0, lty = 2)
points(colMeans(m4$mu), pch = 21, bg = 'red2', cex = 2)
points(y ~ p, pch = 1, cex = 0.5)


plot(colMeans(m3$mu) ~ colMeans(m4$mu), 
     ylab = 'Random effect estimates of pop. mass',
     xlab = 'Fixed effect estimates of pop. mass')
abline(0,1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework 1:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# There's a slight signal here [in the last plot above], do you see it?
# The estimates are the same near the mean, but the residuals are slightly
# above the line at low weights, and below the line at high weights.
# 
# This is called 'shrinkage.' Put simply, the random effects are shrinking ever so
# slighlty towards the mean. Go back up to the beginning of the script where
# we define lambda, or the average number of foxes sampled per population. Change
# lambda to ~4 and run both models again. What happens? Why?
#
# Now change lambda to 35 (we'll measure more foxes) and run the 
# entire script again. There should be less shrinkage. What do you predict would
# happen if we measured ~1000 foxes in each population (don't do this it will take
# hours!)
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~












# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model 5:
#
# Now let's imagine we want to estimate variation among populations 
# around a trend-line that we estimate as a function of latitude?
#
# How would we do that?
# 
# First, we can't with fixed effects because the estimates of independent population
# variation and the effect of latitude would be confounded (don't take my word
# for it, give it a shot and see what happens)
#
# We can do this with random effects. We will estimate the means of each population
# as a random draw around the trend-line estimated using an effect of latitude,
# i.e., random residuals.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m5.stan <- stan_model(model_code = 
                        "                        
data {
  int n;            
  int l;               
  int p[n];         
  vector[n] y;
  vector[l] x;
}

parameters {
  real<lower=-100, upper = 100> mu_star; 
  real<lower=0> sigma[2];
  real mu[l];
  real beta;
}

transformed parameters {
  vector[l] E;
  E = mu_star + beta * x;
}

model {
  mu ~ normal(E, sigma[1]);
  sigma ~ gamma(1,1);
  beta ~ normal(0, 3.16);
  mu_star ~ normal(12, 3.16);
  
  for (i in 1:n){
    y[i] ~ normal(mu[p[i]], sigma[2]);
  }

}
")




#############################################################################################
# Bundle data
#############################################################################################


stan_data <- list(y = y, l = l, n = length(p), p = p, x = x)

inits <- function(){list()}  

nc <- 4
nt <- 25
ni <- 5000
nb <- 2500


Sys.time()
m5 <- sampling(m5.stan, 
               data = stan_data,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)


Sys.time()

print(m5)
m5 <- extract(m5)


lci <- NULL
uci <- NULL
for (j in 1:l){
  lci[j] <- quantile(m5$mu[,j], c(0.025))
  uci[j] <- quantile(m5$mu[,j], c(0.975))
}
plot(colMeans(m5$mu) ~ x, ylim = c(9,15), xlim = c(-3,3),
     ylab = 'Mean mass (lbs)',
     xlab = 'Truth', cex.lab = 2)
arrows(x, lci, x, uci, length = 0, lty = 2, lwd = 1.5)
points(colMeans(m5$mu) ~ x, pch = 21, bg = 'red2', cex = 2.5)

vioplot(m5$beta)
plot(m5$beta, type = 'l')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Now let's compare estimates from this model to our random effects only model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(colMeans(m5$mu) ~ colMeans(m3$mu),
     ylab = 'Estimates from m(random + latitude)',
     xlab = 'Estimates from m(random)', cex.lab = 2)
abline(0,1)

boxplot(m3$sigma[,1], m5$sigma[,1])
# less variance because the covariate (latitude) explains much
# of the among population variance!



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework Part II:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Rerun the entire script, but change the number of sampled populations
#    to 4? What happens to our estimates of the grand mean and among-group
#    variance? Why?
plot(m3$sigma[,1] ~ m3$mu.star, cex.lab = 1,
     xlab = expression(mu*'*'), ylab = expression(sigma[1]))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



