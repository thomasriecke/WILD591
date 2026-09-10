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
maturity = rnorm(n, 0, 1)

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


# format data for Stan
y <- warblers
c <- as.numeric(scale(canopy))
s <- as.numeric(scale(subcan))
log(mean(y))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# write STan models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m.both.stan <- stan_model(model_code = 
"                        
data {
  int n;            
  array[n] int y;            
  vector[n] c;
  vector[n] s;
}

parameters {
  real alpha;
  array[2] real beta;
}

transformed parameters {
  vector[n] psi;
  for (i in 1:n){
    psi[i] = exp(alpha + beta[1] * c[i] + beta[2] * s[i]);
  }
}

model {
  alpha ~ normal(1, 1);
  beta ~ normal(0, 10);
  for (i in 1:n){
    y[i] ~ poisson(psi[i]);
  }

}
")


m.can.stan <- stan_model(model_code = 
"                        
data {
  int n;            
  array[n] int y;                
  vector[n] c;
  // vector[n] s;
}

parameters {
  real alpha;
  array[2] real beta;
}

transformed parameters {
  vector[n] psi;
  for (i in 1:n){
    psi[i] = exp(alpha + beta[1] * c[i]);
  }
}

model {
  alpha ~ normal(1, 1);
  beta ~ normal(0, 10);
  for (i in 1:n){
    y[i] ~ poisson(psi[i]);
  }
}
")


m.sub.stan <- stan_model(model_code = 
                           "                        
data {
  int n;            
  array[n] int y;           
  // vector[n] c;
  vector[n] s;
}

parameters {
  real alpha;
  array[2] real beta;
}

transformed parameters {
  vector[n] psi;
  for (i in 1:n){
    psi[i] = exp(alpha + beta[2] * s[i]);
  }
}

model {
  alpha ~ normal(1, 1);
  beta ~ normal(0, 10);
  for (i in 1:n){
    y[i] ~ poisson(psi[i]);
  }
}
")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# provide  data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stan_data_both <- list(y = y, n = n, c = c, s = s)
stan_data_can <- list(y = y, n = n, c = c)
stan_data_sub <- list(y = y, n = n, s = s)

inits <- function(){list()}  
nc <- 4
nt <- 25
ni <- 5000
nb <- 2500

Sys.time()
m.b <- sampling(m.both.stan, 
               data = stan_data_both,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)

m.c <- sampling(m.can.stan, 
                data = stan_data_can,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                open_progress = FALSE, refresh = -1, cores = 4)

m.s <- sampling(m.sub.stan, 
                data = stan_data_sub,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                open_progress = FALSE, refresh = -1, cores = 4)
Sys.time()

m.b <- extract(m.b)
m.c <- extract(m.c)
m.s <- extract(m.s)


# estimates from model with both covariates (grey) and 'stand-alone' models 'red'
vioplot(m.b$beta)
vioplot(m.c$beta[,1], at = 1, add = T, col = 'red')
vioplot(m.s$beta[,2], at = 2, add = T, col = 'red')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 'HOMEWORK'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Whew! Now things are really getting 'fun.'
# Let's take a minute to look at the parameter estimates of the effects of 
# canopy height (beta[1]) and sub-canopy height (beta[2]) from the models
# that only included single covariates (m.c and m.s)
#
# Ok. They look fairly similar! Now let's take a peak at the effects of these 
# covariates from a model when both covariates were included
#
# huh... what happened?!
# 
# let's take a look at the posterior distributions plotted against each other...
# this is where the term 'variance inflation' comes from.
#
# Now this is starting to get interesting! What is happening here?
# here are some clues to get you started, and we'll revisit this as a group
#
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


m.both.stan <- stan_model(model_code = 
"                        
data {
  int n;            
  vector[n] v;         
  vector[n] w;
  vector[n] u;
}

parameters {
  real alpha;
  array[2] real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n] psi;
  for (i in 1:n){
    psi[i] = alpha + beta[1] * w[i] + beta[2] * u[i];
  }
}

model {
  alpha ~ normal(1, 1);
  beta ~ normal(0, 10);
  sigma ~ gamma(1,1);
  v ~ lognormal(psi, sigma);
}
")

m.w.stan <- stan_model(model_code = 
                            "                        
data {
  int n;            
  vector[n] v;         
  vector[n] w;
  // vector[n] u;
}

parameters {
  real alpha;
  array[2] real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n] psi;
  for (i in 1:n){
    psi[i] = alpha + beta[1] * w[i];
  }
}

model {
  alpha ~ normal(1, 1);
  beta ~ normal(0, 10);
  sigma ~ gamma(1,1);
  v ~ lognormal(psi, sigma);
}
")

m.u.stan <- stan_model(model_code = 
"                        
data {
  int n;            
  vector[n] v;         
  // vector[n] w;
  vector[n] u;
}

parameters {
  real alpha;
  array[2] real beta;
  real<lower=0> sigma;
}

transformed parameters {
  vector[n] psi;
  for (i in 1:n){
    psi[i] = alpha + beta[2] * u[i];
  }
}

model {
  alpha ~ normal(1, 1);
  beta ~ normal(0, 10);
  sigma ~ gamma(1,1);
  v ~ lognormal(psi, sigma);
}
")


stan_data_both <- list(v = v, n = n, w = zw, u = zu)
stan_data_wolf <- list(v = v, n = n, w = zw)
stan_data_ungulate <- list(v = v, n = n, u = zu)

inits <- function(){list()}  
nc <- 4
nt <- 25
ni <- 5000
nb <- 2500

Sys.time()
m.b <- sampling(m.both.stan, 
                data = stan_data_both,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                open_progress = FALSE, refresh = -1, cores = 4)

m.w <- sampling(m.w.stan, 
                data = stan_data_wolf,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                open_progress = FALSE, refresh = -1, cores = 4)

m.u <- sampling(m.u.stan, 
                data = stan_data_ungulate,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                open_progress = FALSE, refresh = -1, cores = 4)
Sys.time()

m.b <- extract(m.b)
m.w <- extract(m.w)
m.u <- extract(m.u)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 'HOMEWORK'
# Make some similar plots as those above, exploring the effects of wolf
# and ungulate density on vegetation density, e.g.,

# what's happening here with our parameter estimates?
# next week, we'll use the exact same system to explore the use of path
# analysis, or sequences of linear models that allow us to estimate 
# causal relationships among measured values.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~