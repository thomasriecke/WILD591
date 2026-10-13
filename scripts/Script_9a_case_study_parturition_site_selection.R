# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# some spatial and statistical packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(terra)
library(MASS)
library(geoR)
library(sf)
library(jagsUI)
library(latex2exp)
library(reshape2)

set.seed(23)
inv.logit=plogis # redefine a function
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up an empty raster and pull coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.res <- 50
n.cells <- grid.res^2
r=rast(ncol=grid.res,nrow=grid.res,extent=c(-grid.res,grid.res,-grid.res,grid.res))                      # make an empty raster
r[]=0                                                                # cero; 0
s.r=crds(r)                                                          # pull coordinates

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first, we'll simulate random variation in disturbance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rf=grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10))         # gaussian (auto-correlated) random field
dist=r                                                               # make dist a spatial raster
dist[]=rf$data                                                       # assign auto-correlated disturbance values
plot(dist,main='simulated disturbance across a homogenous habitat')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note that there is spatial auto-correlation in disturbance
# i.e., some cells are almost completely undisturbed (-; dark blue) 
# while other cells have high disturbance levels (+; yellow)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2nd, we'll simulate predator activity (0 -> Inf) as a function of disturbance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pred = r                         # here we use r to assign class/structure to pred
alpha1 <- 1
beta1 <- -0.35
pred[] = exp(alpha1 + beta1 * dist[] + rnorm(n.cells, 0, 0.25))
plot(pred,main='simulated predator activity given disturbance')

hist(pred, breaks = 100, main = NULL, xlab = 'Predator activity')
plot(pred[] ~ dist[], ylab = 'Predator activity', xlab = "Disturbance",
     las = 1, cex.lab = 2, cex.axis = 1.5)
cor(pred[], dist[])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# third, we will simulate cover (centered) as a function of disturbance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cov = r                         # here we use r to assign class/structure to cov
alpha2 <- 0
beta2 <- 0.5
cov[] = alpha2 + beta2 * dist[] + rnorm(n.cells, 0, 0.5)
plot(cov,main='simulated cover as a function of disturbance')

hist(cov, breaks = 100, main = NULL, xlab = 'Cover')
plot(cov[] ~ dist[], 
     ylab = 'Cover', xlab = "Disturbance",
     las = 1, cex.lab = 2, cex.axis = 1.5)
cor(cov[], dist[])

plot(cov[] ~ pred[], ylab = 'Cover (c)', xlab = 'Predator activity (p)')
cor(cov[],pred[])
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Finally, we need to simulate parturition sites and random sites....
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mu <- 0
gamma1 <- -0.5
gamma2 <- 1.5

# make probability of selection a spatial raster
rsf <- r
rsf[] <- inv.logit(mu + gamma1*pred[] + gamma2*cov[])
plot(rsf)
plot(rsf[] ~ cov[])
plot(rsf[] ~ pred[])
plot(rsf, main = 'Relative probability of parturition')


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# there are certainly less hacky ways to do this using spatial packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.collar <- 100
n.random <- 250
p <- r
p[] <- rsf[]/sum(rsf[]) # make probability of selection relative (i.e., sum to 1)
plot(p, main = 'probability of parturition site selection')

sites <- rmultinom(1, n.collar, p[])
random <- rmultinom(1, n.random, rep(1/grid.res^2, grid.res^2))


plot(dist,main='simulated disturbance across a homogenous habitat')
points(s.r[which(sites > 0),1], s.r[which(sites > 0),2],
       col = 'red', pch = 19)

plot(dist,main='simulated disturbance across a homogenous habitat')
points(s.r[which(random > 0),1], s.r[which(random > 0),2],
       col = 'red', pch = 1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create data for JAGS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y <- c(rep(1, n.collar), rep(0, n.random))
c <- NULL
p <- NULL
d <- NULL

# pull covariate values for parturition sites... there are less hacky ways to do this :)
for (i in 1:length(sites)){
  if(sites[i] != 0){
    # note here we're using c as a variable name and as a function
    # that's poor practice on my part
    c <- c(c, rep(cov[][i], sites[i]))
    p <- c(p, rep(pred[][i], sites[i]))
    d <- c(d, rep(dist[][i], sites[i]))
  }
}

# pull covariate values for random points... there are less hacky ways to do this :)
for (i in 1:length(random)){
  if(random[i] != 0){
    # note here we're using c as a variable name and as a function
    # that's poor practice on my part
    c <- c(c, rep(cov[][i], random[i]))
    p <- c(p, rep(pred[][i], random[i]))
    d <- c(d, rep(dist[][i], random[i]))
  }
}

plot(jitter(y, 0.1) ~ c)
plot(jitter(y, 0.1) ~ p)
plot(jitter(y, 0.1) ~ d)

summary(glm(y ~ p + c + d))
summary(glm(y ~ d))
summary(glm(y ~ d + c))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data-generating model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m_path.jags")
cat("
model{

  alpha[1] ~ dnorm(1, 0.1)
  alpha[2] ~ dnorm(0, 0.1)
  alpha[3] ~ dlogis(0,1)

  for (k in 1:2){
    beta[k] ~ dnorm(0, 0.1)
    gamma[k] ~ dnorm(0, 0.1)
    sigma[k] ~ dgamma(1,1)
    tau[k] = pow(sigma[k], -2)   # easier to write than 1/(sigma[k] * sigma[k])
  }

  for (i in 1:n){
    p[i] ~ dlnorm(alpha[1] + beta[1] * d[i], tau[1])
    c[i] ~ dnorm(alpha[2] + beta[2] * d[i], tau[2])
  
    logit(psi[i]) = alpha[3] + gamma[1] * p[i] + gamma[2] * c[i]
    y[i] ~ dbern(psi[i])
  }

}
",fill = TRUE)
sink()



dat = list(y = y, d = d, c = c, p = p, n = length(y))
pars = c('alpha','beta','sigma','gamma')
inits <- function(){list()}  
m1 <- jags(data = dat, inits=inits, model.file = "m_path.jags",
           parameters.to.save = pars, n.chains = 5, n.iter = 15000, 
           n.thin = 10, n.burnin = 5000, parallel = T)







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data-generating model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m_latent.jags")
cat("
model{

  alpha[1] ~ dnorm(1, 0.1)
  alpha[2] ~ dnorm(0, 0.1)
  alpha[3] ~ dnorm(0, 0.1)
  alpha[4] ~ dlogis(0,1)

  beta[1] ~ dnorm(0,0.1)
  beta[2] ~ dnorm(0,0.1)
  beta[3] = 1
  gamma[1] ~ dnorm(0, 0.1)

  for (k in 1:4){
    sigma[k] ~ dgamma(1,1)
    tau[k] = pow(sigma[k], -2)   # easier to write than 1/(sigma[k] * sigma[k])
  }

  for (i in 1:n){
  
    ## this is called a 'non-centered' parameterization
    ## it's mathematically identical to: 
    ## eta[i] ~ dnorm(0, tau[4])
    ## but sometimes samples more efficiently
    eta.star[i] ~ dnorm(0, 1)
    eta[i] = eta.star[i] * sigma[4]
    
    p[i] ~ dlnorm(alpha[1] + beta[1] * eta[i], tau[1])
    c[i] ~ dnorm(alpha[2] + beta[2] * eta[i], tau[2])
    d[i] ~ dnorm(alpha[3] + beta[3] * eta[i], tau[3])
    
    logit(psi[i]) = alpha[4] + gamma[1] * eta[i]
    y[i] ~ dbern(psi[i])
  }

}
",fill = TRUE)
sink()



dat = list(y = y, d = d, c = c, p = p, n = length(y))
pars = c('alpha','beta','sigma','gamma')
inits <- function(){list()}  
m2 <- jags(data = dat, inits=inits, model.file = "m_latent.jags",
           parameters.to.save = pars, n.chains = 5, n.iter = 15000, 
           n.thin = 10, n.burnin = 5000, parallel = T)

print(m2)






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 0) add credible intervals to the plots (see plotting code below)
# 
# 1) How might the y-axis of 'selection' (psi) change if we added
#    more random points (e.g., 500 rather than 250)? Why?
#    if you're unfamiliar with RSFs, ask someone who is :)
#
# 2) Which of the two models explicitly allows you to modify values
#    of disturbance and observe direct and indirect effects on other
#    variables of interest? Explicitly predict probability of 
#    parturition site selection at -1 disturbance using that model
# 
# 3) compare plots from m1 and m2, specifically assessing the relationship
#    between selection and disturbance... how much does inference vary
#    as a function of the different parameterizations?
#
# 4) How would the interpretation of the relationship between
#    predator density and selection change across the two models?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's make a plot or two using the results from Model 1
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplot(m1$sims.list$beta, names = c(TeX("$\\beta_1$"), TeX("$\\beta_2$")))
boxplot(m1$sims.list$gamma, names = c(TeX("$\\gamma_1$"), TeX("$\\gamma_2$")))

n.iter <- length(m1$sims.list$alpha[,1])
res <- 100
xd <- seq(min(dist[]), max(dist[]), length.out = res)
xd
pc <- matrix(NA, n.iter, res)
pp <- matrix(NA, n.iter, res)
prsf <- matrix(NA, n.iter, res)

qc <- matrix(NA, res, 5)
qp <- matrix(NA, res, 5)
qrsf <- matrix(NA, res, 5)

for (j in 1:res){
  # expected parameter values
  pp[,j] <- exp(m1$sims.list$alpha[,1] + m1$sims.list$beta[,1] * xd[j])
  pc[,j] <- m1$sims.list$alpha[,2] + m1$sims.list$beta[,2] * xd[j]
  prsf[,j] <- plogis(m1$sims.list$alpha[,3] + m1$sims.list$gamma[,1] * pp[,j] + m1$sims.list$gamma[,2] * pc[,j])
  
  # quantiles of derived quantities  
  qp[j,] <- quantile(pp[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))
  qc[j,] <- quantile(pc[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))
  qrsf[j,] <- quantile(prsf[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))  
}

mpp <- melt(pp); names(mpp) = c('row','j','p')
mpc <- melt(pc); names(mpc) = c('row','j','c')
mprsf <- melt(prsf); names(mprsf) = c('row','j','rsf')

# Plot parturition site selection vs. disturbance
smoothScatter(mprsf$rsf ~ xd[mprsf$j], nrpoints = 0, las = 1,
              ylab = TeX("Selection ($\\psi$)"),
              xlab = 'Disturbance (d)')
# Plot parturition site selection vs. cover
smoothScatter(mprsf$rsf ~ mpc$c, nrpoints = 0, las = 1,
              ylab = TeX("Selection ($\\psi$)"),
              xlab = 'Cover (c)')
# Plot parturition site selection vs. predator activity
smoothScatter(mprsf$rsf ~ mpp$p, nrpoints = 0, las = 1,
              ylab = TeX("Selection ($\\psi$)"),
              xlab = 'Predator activity (p)')






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's make a plot or two using the results from Model 2
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boxplot(m2$sims.list$beta, names = c(TeX("$\\beta_1$"), TeX("$\\beta_2$"), TeX("$\\beta_3$")))
boxplot(m2$sims.list$gamma, names = c(TeX("$\\gamma_1$")))

n.iter <- length(m1$sims.list$alpha[,1])
res <- 100
xeta <- seq(-3*m2$q50$sigma[4], 3*m2$q50$sigma[4], length.out = res)
pc <- matrix(NA, n.iter, res)
pp <- matrix(NA, n.iter, res)
pd <- matrix(NA, n.iter, res)
prsf <- matrix(NA, n.iter, res)

qc <- matrix(NA, res, 5)
qp <- matrix(NA, res, 5)
qd <- matrix(NA, res, 5)
qrsf <- matrix(NA, res, 5)

for (j in 1:res){
  # expected parameter values
  pp[,j] <- exp(m2$sims.list$alpha[,1] + m2$sims.list$beta[,1] * xeta[j])
  pc[,j] <- m2$sims.list$alpha[,2] + m2$sims.list$beta[,2] * xeta[j]
  pd[,j] <- m2$sims.list$alpha[,3] + m2$sims.list$beta[,3] * xeta[j]  
  prsf[,j] <- plogis(m2$sims.list$alpha[,4] + m2$sims.list$gamma * xeta[j])
  
  
  # quantiles of derived quantities  
  qp[j,] <- quantile(pp[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))
  qc[j,] <- quantile(pc[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))
  qd[j,] <- quantile(pd[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))  
  qrsf[j,] <- quantile(prsf[,j], c(0.025, 0.05, 0.5, 0.95, 0.975))  
}

mpp <- melt(pp); names(mpp) = c('row','j','p')
mpc <- melt(pc); names(mpc) = c('row','j','c')
mpd <- melt(pd); names(mpd) = c('row','j','d')
mprsf <- melt(prsf); names(mprsf) = c('row','j','rsf')

# Plot parturition site selection vs. disturbance
smoothScatter(mprsf$rsf ~ mpd$d, nrpoints = 0, las = 1,
              ylab = TeX("Selection ($\\psi$)"),
              xlab = 'Disturbance (d)')
# Plot parturition site selection vs. cover
smoothScatter(mprsf$rsf ~ mpc$c, nrpoints = 0, las = 1,
              ylab = TeX("Selection ($\\psi$)"),
              xlab = 'Cover (c)')
# Plot parturition site selection vs. predator activity
smoothScatter(mprsf$rsf ~ mpp$p, nrpoints = 0, las = 1,
              ylab = TeX("Selection ($\\psi$)"),
              xlab = 'Predator activity (p)')













