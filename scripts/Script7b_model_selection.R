# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(MuMIn)
library(latex2exp)
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BEFORE Example 0:
# What is likelihood...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 50
x <- rnorm(n, 0, 1)
y <- rnorm(n, x, 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# silly figure, which line is more likely? why?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(y ~ x, cex.lab = 2, las = 1); abline(0,1)
plot(y ~ x, cex.lab = 2, las = 1); abline(0,0)
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(y ~ x, cex.lab = 2, las = 1); abline(0,1)
plot(y ~ x, cex.lab = 2, las = 1); abline(0,0.5)
par(mfrow = c(1,1), mar = c(5,5,2,2))


LL <- function(y, x, intercept, slope){
  # sample size
  n = length(y)
  # expected mass (e) given tarsus, intercept, and slope
  e = intercept + slope * x
  # residual standard deviation (avg. distance of sample from expected)
  v = (1/n) * sum((y - intercept - slope*x)^2)
  s = sqrt(v)
  ll <- -(n/2) * log(2*pi) - (n/2) * log(v) -(1/(2*v)) * sum((y - e)^2)
  return(ll)
  
}

# extract logLik using R functions and our own function
summary(model <- lm(y ~ x))
logLik(model)
LL(y, x, model$coefficients[1], model$coefficients[2])

# think about it
plot(y ~ x, cex.lab = 2, las = 1, cex = 1.5)
abline(model$coefficients[1], model$coefficients[2], lwd = 2)


# extract AIC
aic.fxn <- function(y, x, intercept, slope){
  # sample size
  n = length(y)
  # expected mass (e) given tarsus, intercept, and slope
  e = intercept + slope * x
  # residual standard deviation (avg. distance of sample from expected)
  v = (1/n) * sum((y - intercept - slope*x)^2)
  s = sqrt(v)
  ll <- -(n/2) * log(2*pi) - (n/2) * log(v) -(1/(2*v)) * sum((y - e)^2)
  # return(ll)
  k <- 3
  aic = (2 * k) + -2*ll
  return(aic)
}

aic.fxn(y, x, model$coefficients[1], model$coefficients[2])
LL(y, x, model$coefficients[1], model$coefficients[2])
logLik(model)
AIC(model)

# an alternative way to compute AIC
# 2k + n * ln(RSS/n)
e <- model$coefficients[1] + model$coefficients[2]*x
RSS <- sum((y - e)^2)
n * log(RSS/n) + 6

# an alternative way to compute AIC
6 + n * (log(2 * pi * (RSS/n)) + 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMPLE 0:
# What are residuals, what is deviance, how is AIC and AICc calculated,
# and what does that all mean?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)
n <- 25
x <- rnorm(n, 0, 1)
y <- rnorm(n, x, 1)
m1 <- glm(y ~ x)
summary(m1)
m1$coefficients
# here eta is the expected value of y
eta <- m1$coefficients[1] + m1$coefficients[2] * x

# create a line of expected y values given x
res <- 100
px <- seq(-5,5, length.out = res)
peta <- m1$coefficients[1] + m1$coefficients[2] * px


plot(y ~ x, las = 1, cex.lab = 2, pch = 21, bg = 'grey50', cex = 2,
     xlim = c(-4,4), ylim = c(-4,4))
points(peta ~ px, col = 'red', type = 'l', lwd = 3)
points(eta ~ x, col = 'red', pch = 19, cex = 1.5)
arrows(x, eta, x, y, lty = 2, col = 'red', length = 0)
points(y ~ x, pch = 21, bg = 'grey50', cex = 2)

m0 <- glm(y ~ 1)
eta0 <- rep(m0$coefficients, n)


# sum of squared error
deviance.hand <- sum((y-eta)^2)
deviance.comp <- m1$deviance

# number of parameters = 3 (intercept, slope, variance)
k <- 3
aic <- 2*k + n * (log(2 * pi * (deviance.hand/n)) + 1)
aicc <- aic + (2 * k^2 + 2*k)/(n - k - 1)
AICc(m1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# exploring the extra penalty in AICc
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cpenalty <- function(n,k){
  cpen <- (2 * k^2 + 2*k)/(n - k - 1)
  print(cpen)
}
cpenalty(10,2)
cpenalty(1000,2)
cpenalty(10000,2)

rm(list=ls())


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMPLE 1a
# here we'll simulate a single covariate (x) that drives
# all of the variation in our response (y). 
#
# then we'll compare model 1 (y ~ x) to model 0 (y ~ 1). We can
# also think of y ~ 1 as an intercept only model!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)
n <- 25
x <- rnorm(n, 0, 1)
y <- rnorm(n, x, 1)

m0 <- glm(y ~ 1)
m1 <- glm(y ~ x)
summary(m0)
summary(m1)

# here eta is the expected value of y
eta0 <- rep(m0$coefficients, n)
eta1 <- m1$coefficients[1] + m1$coefficients[2] * x


# create a line of expected y values given x
res <- 100
px <- seq(-5,5, length.out = res)
peta0 <- rep(m0$coefficients, res)
peta1 <- m1$coefficients[1] + m1$coefficients[2] * px

# residual sum of squares
rss0 <- sum((y - eta0)^2)
rss1 <- sum((y - eta1)^2)

# R can do this for us :)
deviance(m0)
deviance(m1)

k0 <- 2
k1 <- 3
aic0 <- 2*k0 + n * (log(2 * pi * (rss0/n)) + 1)
aicc0 <- aic0 + (2 * k0^2 + 2*k0)/(n - k0 - 1)

aic1 <- 2*k1 + n * (log(2 * pi * (rss1/n)) + 1)
aicc1 <- aic1 + (2 * k1^2 + 2*k1)/(n - k1 - 1)

# R can do this for us too :)
AICc(m0)
AICc(m1)


plot(y ~ x, las = 1, cex.lab = 2, pch = 21, bg = 'grey50', cex = 2,
     xlim = c(-4,4), ylim = c(-4,4))
points(peta0 ~ px, col = 'red', type = 'l', lwd = 3)
arrows(x, eta0, x, y, lty = 2, col = 'red', length = 0)

plot(y ~ x, las = 1, cex.lab = 2, pch = 21, bg = 'grey50', cex = 2,
     xlim = c(-4,4), ylim = c(-4,4))
points(peta1 ~ px, col = 'red', type = 'l', lwd = 3)
arrows(x, eta1, x, y, lty = 2, col = 'red', length = 0)





# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fun. let's do that 10000 times just to be really confident!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### number of simulations (sims)
### deltaAIC (delta)
### correlation between x and m (r)
### relative model weights (w.x, w.m)
sims <- 10000
delta <- NULL
w0 <- NULL
w1 <- NULL
n <- 25
px <- NULL
sx <- NULL
mx <- NULL
lx <- NULL # lower 85% confidence interval

Sys.time()
for (i in 1:sims){
  
  x <- rnorm(n, 0, 1)
  y <- rnorm(n, x, 1)
  
  m0 <- glm(y ~ 1)  
  m1 <- glm(y ~ x)

  aic0 <- AICc(m0)
  aic1 <- AICc(m1)

  tmp <- model.sel(list(m0,m1))
  delta[i] <- aic1 - aic0

  mx[i] <- summary(m1)$coefficients[2,1]  
  px[i] <- summary(m1)$coefficients[2,4]
  sx[i] <- summary(m1)$coefficients[2,2]
  lx[i] <- mx[i] - sx[i]*1.44
  w0[i] <- tmp$weight[is.na(tmp$x)]     
  w1[i] <- tmp$weight[!is.na(tmp$x)]
  
}
Sys.time()


hist(delta, breaks = 500, xlab = TeX("$\\delta$AICc (m1)"), cex.lab = 2, main = NULL)

plot(delta ~ lx, las = 1,
     ylab = TeX("$\\delta$AICc"),
     xlab = TeX("Lower 85% confidence interval of $\\beta_1$"))
points(0,0, pch = 19, col = 'red')


mean(delta)
median(delta)
length(which(delta > 0))/sims
length(which(delta > -2))/sims
length(which(delta > -7))/sims




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMPLE 1b
# here we'll simulate a single covariate (x) that drives
# all of the variation in our response (y). 
#
# we'll also simualte a meaningless (m) covariate m ~ normal(0, 1)
# and compare results from two models
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sims <- 10000
delta <- NULL
w1 <- NULL
w2 <- NULL
n <- 25
px <- NULL
pm <- NULL
sx <- NULL
mx <- NULL
lx <- NULL # lower 85% confidence interval
r <- NULL

Sys.time()
for (i in 1:sims){
  
  m <- rnorm(n, 0, 1)
  x <- rnorm(n, 0, 1)
  y <- rnorm(n, x, 1)
  
  m1 <- glm(y ~ x)  
  m2 <- glm(y ~ x + m)
  
  aic1 <- AICc(m1)
  aic2 <- AICc(m2)
  
  tmp <- model.sel(list(m1,m2))
  r[i] <- cor(m,x)
  delta[i] <- aic1 - aic2
  
  pm[i] <- summary(m2)$coefficients[3,4] 
  
  mx[i] <- summary(m1)$coefficients[2,1]  
  px[i] <- summary(m1)$coefficients[2,4]
  sx[i] <- summary(m1)$coefficients[2,2]
  lx[i] <- mx[i] - sx[i]*1.44
  w1[i] <- tmp$weight[is.na(tmp$m)]     
  w2[i] <- tmp$weight[!is.na(tmp$m)]
  
}
Sys.time()




eta1 <- predict(m1)
eta2 <- predict(m2)

m1$coefficients
plot(y ~ x, cex.lab = 2, las = 1)
points(eta1 ~ x, col = 'red', pch = 19)

m2$coefficients
plot(y ~ x, cex.lab = 2, las = 1)
points(eta2 ~ x, col = 'purple', pch = 19)

deviance(m1)
deviance(m2)

hist(delta, breaks = 500, main = '', cex.lab = 2,
     xlab = TeX("$\\delta$AIC"))


length(which(delta > 0))/sims
length(which(delta > -2))/sims
length(which(delta > -4))/sims

mean(delta)
median(delta)


plot(delta ~ r, ylab = TeX("$\\delta$AIC"), xlab = TeX("$r_{xm}"), cex.lab = 2)

plot(delta ~ pm, ylab = TeX("$\\delta$AIC"), xlab = TeX("$p_{m}"), cex.lab = 2)
length(which(pm < 0.1))/sims

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXAMPLE 1c
# here we'll simulate a single covariate (x) that drives
# all of the variation in our response (y). We'll also simulate 
# a multi-collinear covariate (m) as a function of x
# that has no impact on y...
#
# then we'll compare model 1 (y ~ x) to model 2 (y ~ m)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


### number of simulations (sims)
### deltaAIC (delta)
### correlation between x and m (r)
### relative model weights (w.x, w.m)
sims <- 10000
delta <- NULL
r <- NULL
w1 <- NULL
w2 <- NULL
n <- 25
# run 10k simulations of the same process
# takes about 45 seconds on a:
# 12th Gen Intel(R) Core(TM) i9-12900K   3.20 GHz
# that's juggling a few other tasks at the moment
Sys.time()
for (i in 1:sims){
  
  x <- rnorm(n, 0, 1)
  m <- rnorm(n, x, 0.5)
  y <- rnorm(n, x, 1)
  
  r[i] <- cor(x,m)
  
  m1 <- glm(y ~ x)
  m2 <- glm(y ~ m)

  aic1 <- AICc(m1)
  aic2 <- AICc(m2)
 
  tmp <- model.sel(list(m1,m2))
  
  delta[i] <- aic1 - aic2
  
  w1[i] <- tmp$weight[!is.na(tmp$x)]
  w2[i] <- tmp$weight[!is.na(tmp$m)]    

}
Sys.time()


# plot distribution of deltaAIC
hist(delta, breaks = 500, main = '', 
     xlab = TeX("$\\delta$AIC"), las = 1, cex.lab = 1.5)

mean(delta)
median(delta)
# These tell us the proportion of simulations in which
# the 'wrong' model (y ~ m) was ranked better than
# the data generating (y ~ x) model [delta > 0]
length(which(delta > 0))/sims
length(which(delta > -2))/sims
length(which(delta > -7))/sims

# the data generating (y ~ x) model [delta > 0]
length(which(w2 > 0.25))/sims
length(which(w2 > 0.5))/sims


# scatter plot 
par(mfrow = c(1,1))
smoothScatter(delta ~ r, nrpoints = 0, xlim = c(0.5,1),
              ylab = TeX("$\\delta$AIC"), xlab = TeX("$r_{xm}$"), cex.lab = 2, las = 1) 
points(delta ~ r, cex = 0.25, col = 'navy')
cor.test(delta,r)














# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Final vignette
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 10000
x <- rnorm(n, 0, 1)
w <- rnorm(n, 0, 0.025)
e = x + w
y <- rnorm(n, e, 1)

par(mfrow = c(1,2))
smoothScatter(y ~ x, las = 1, cex.lab = 2)
points(y ~ x, col = 'navy', cex = 0.25)
points(c(x) ~ x, col = 'red', cex = 0.25)

smoothScatter(y ~ x, las = 1, cex.lab = 2)
points(y ~ x, col = 'navy', cex = 0.25)
points(c(x + w) ~ x, col = 'red', cex = 0.25)


m1 <- lm(y ~ x)
m2 <- lm(y ~ x + w)
AICc(m1)
AICc(m2)

