# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we'll simulate a composite covariate (disturbance!)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

set.seed(23) # Jordan
inv.logit=plogis # redefine a function
par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up an empty raster and pull coordinates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grid.res <- 50
n.cells <- grid.res^2
r=rast(ncol=grid.res,nrow=grid.res,extent=c(-grid.res,grid.res,-grid.res,grid.res))                      # make an empty raster
r[]=0                                                                # cero; 0
s=crds(r)                                                          # pull coordinates

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# first, we'll simulate random variation in road density
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rf=grf(1, grid = xyFromCell(r,1:ncell(r)), cov.pars=c(1,10))         # gaussian (auto-correlated) random field
d=r                                                               # make dist a spatial raster
road <- exp(rf$data)
d$road = road                                                      # assign auto-correlated disturbance values
plot(d['road'], main='simulated road density across a ~homogenous national forest')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# note that there is spatial auto-correlation in road density (b/c duh)
# some areas (negative values) have low road densities
# others have higher road densities...
#
# now we need to simulate 5 things that lead to disturbance
# 1) cars: linear with road density
# 2) hunters: 
# 3) mtb-ers:      quadratic with road density
# 4) nature folks: quadratic with road density
# 5) paragliders: just in one weird spot!
#
# Then we'll simulate resource use as a 
# function of a composite 'disturbance' index...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# create hunters as a function of road density and longitude
hunt <- as.numeric(exp(1 + 0.5 * road + -0.0625 * road^2 + rnorm(n.cells, 0, 0.5) + 0.001 * s[,1]))
plot(hunt ~ road, ylab = 'Hunters', xlab = 'Road density', cex.lab = 2)

# create cars
cars <- exp(as.numeric(0.2 * road + rnorm(n.cells, 0, 0.25)))
plot(cars ~ road, ylab = 'Cars', xlab = 'Road density', cex.lab = 2)

# create mtb-ers as a function of road density and longitude
# they're rocking the edge of town...
mtbs <- as.numeric(exp(1 + 0.125 * road + -0.0075 * road^2 + rnorm(n.cells, 0, 0.2) + -0.01 * s[,1]))
plot(mtbs ~ road, ylab = 'MTBers', xlab = 'Road density', cex.lab = 2)


# create birders as a function of road density and latitude
bird <- as.numeric(exp(-1 + 0.125 * road + -0.005 * road^2 + rnorm(n.cells, 0, 0.15) + 0.01 * s[,2]))
plot(bird ~ road, ylab = 'Birders [nerds]', xlab = 'Road density', cex.lab = 2)

# create paragliders... they're just in the corners
para <- exp(0 + 0.0005 * s[,1] + 0.0005 * s[,2] + 0.0005 * s[,1] * s[,2])

# assign to spatial data.frame
d$hunt <- hunt
d$cars <- cars
d$mtbs <- mtbs
d$bird <- bird
d$para <- para


# road density
plot(d['road'], main = 'Road density')

# hunters hunt the edge of town, not the most realistic, change the simulation if you don't like it 
plot(d['hunt'], main = 'Hunter density')
plot(hunt ~ road, ylab = 'Hunters', xlab = 'Road density', cex.lab = 2)

# lots of cars in the towns
plot(d['cars'], main = 'Car density')
plot(cars ~ road, ylab = 'Cars', xlab = 'Road density', cex.lab = 2)

# gosh darn mtbers out west...
plot(d['mtbs'], main = 'MTBer density')
plot(mtbs ~ road, ylab = 'MTBers', xlab = 'Road density', cex.lab = 2)

# more birders in the north
plot(d['bird'], main = 'Birders')
plot(bird ~ road, ylab = 'Birders [nerds]', xlab = 'Road density', cex.lab = 2)

# paragliders are just in the corners... freaking weird
plot(d['para'])

par(mfrow = c(1,1))
plot(para ~ hunt, ylab = 'Paraglider density',
     xlab = "Hunter density", cex.lab = 1.5)


# now we'll make a disturbance composite for ungulate resource selection...
# i just made up the parameter values
dist <- 1 * cars + 0.5 * mtbs + 3 * para + 1 * bird + 2 * hunt
d$dist <- dist

plot(d$dist, main = 'Disturbance')
plot(d$road, main = 'Roads')
hist(dist, xlim = c(0,90))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# well... that's interesting...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cor(cbind(hunt,bird,cars,road,para,mtbs))
plot(dist ~ road, ylab = 'Disturbance index', xlab = 'Road density', cex.lab = 2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll simulate habitat use by samsquantches (?) as a function of these different
# covariates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
psi <- exp(6 + -0.2 * dist)
plot(psi ~ dist)
used <- rpois(n.cells, psi)
d$used <- used
plot(d$used, main = 'Bigfoot locations')
plot(d$road, main = 'Road density')



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's prep our data for JAGS and see if we can recover our parameter estimates
# we'll subsample the total grid
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samp <- sample(1:n.cells, 1000, replace = F)
y <- used[samp]
c <- cars[samp]
h <- hunt[samp]
b <- bird[samp]
p <- para[samp]
m <- mtbs[samp]

plot(y ~ c)
plot(y ~ h)
plot(y ~ b)
plot(y ~ p)
plot(y ~ m)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# data-generating model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("m_comp.jags")
cat("
model{

  alpha[1] ~ dnorm(5, 0.1)
  alpha[2] ~ dnorm(0, 0.1)

  beta[1] = 1
  for (k in 2:5){
    beta[k] ~ dnorm(1, 0.1) T(0,)
  }

  for (i in 1:n){
    d[i] = beta[1] * c[i] + beta[2] * h[i] + beta[3] * p[i] + beta[4] * b[i] + beta[5] * m[i]
    z[i] = d[i]/sd(d[])
    psi[i] = exp(alpha[1] + alpha[2] * z[i])
    y[i] ~ dpois(psi[i])
  }

}
",fill = TRUE)
sink()



dat = list(y = y, h = h, c = c, m = m, b = b, p = p, n = length(samp))
pars = c('alpha','beta')
inits <- function(){list()}  
# this will take a few minutes...
mod <- jags(data = dat, inits=inits, model.file = "m_comp.jags",
           parameters.to.save = pars, n.chains = 1, n.iter = 50000, 
           n.thin = 5, n.burnin = 25000, parallel = T)

print(summary(mod), digits = 3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's look at convergence
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

par(mfrow = c(1,1), mar = c(5,5,2,2))
plot(mod$sims.list$beta[,2])
plot(mod$sims.list$beta[,3])
plot(mod$sims.list$beta[,4])
plot(mod$sims.list$beta[,5])

plot(mod$sims.list$alpha[,1])
plot(mod$sims.list$alpha[,2])

par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(used ~ dist, ylab = 'Bigfoot locations', xlab = 'Disturbance index',
     cex.lab = 1.5, las = 1)
plot(used ~ road, ylab = 'Bigfoot locations', xlab = 'Road density',
     cex.lab = 1.5, las = 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK QUESTIONS
# 
# Make a figure showing the relationships between different types of disturbance
# and composite disturbance.
# 
# Make a map showing expected disturbance and expected bigfoot use
#
# Make a map showing residuals of expected disturbance and expected bigfoot use
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

