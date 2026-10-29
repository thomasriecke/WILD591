library(reshape2)
library(latex2exp)
library(jagsUI)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate size(s) of individuals 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 500
x <- rnorm(n, 35, 2)

par(mfrow = c(1,2))
hist(x, xlab = 'Female tarsus (mm)', main = NULL, breaks = 50, cex.lab = 1.5)
z <- as.numeric(scale(x))
hist(z, xlab = 'Z(Female tarsus)', main = NULL, breaks = 50, cex.lab = 1.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate clutch size(s)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mu is each individual's latent value
# zeta is the regression parameter linking body size to 
zeta <- 2
mu <- zeta * z

# c: cut-off points
# nC: number of clutch types
nC <- 5
cut <- seq(-3,3,length.out = nC-1)


par(mfrow = c(1,1))
plot(mu ~ z, ylab = TeX('$\\mu$'), xlab = 'Z(tarsus)')


par(mfrow = c(1,1))
hist(mu, xlab = TeX('$\\mu$'), breaks = 50, main = NULL)
abline(v = cut, lwd = 3, col = 'red')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# create matrices for relative probabilities of each outcome (p)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
q <- p <- matrix(NA, n, nC)
q[,1] <- plogis(cut[1] - mu)
p[,1] <- q[,1]

for (j in 2:(nC-1)){
  q[,j] <- plogis(cut[j] - mu)
  p[,j] <- q[,j] - q[,j-1]
}
p[,nC] <- 1 - q[,nC-1]
rowSums(p)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# draw clutch sizes (y)
# there are five categories
# note that I'm adding 4 to each category to get back to actual
# clutch sizes (i.e., 5 to 9)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
y <- NULL
for (i in 1:n){
  y[i] <- which(rmultinom(1, 1, p[i,]) == 1) + 4
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# scatterplot of our real data (eggs in nest and tarsus length in mm)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(jitter(y) ~ x, 
     cex.lab = 1.5,
     ylab = 'Eggs', 
     xlab = 'Tarsus (mm)')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model(s)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("poisson.jags")
cat("
model{

  alpha ~ dnorm(log(7), 1)
  beta ~ dnorm(0, 1)

  for (i in 1:n){
    mu[i] = exp(alpha + beta * z[i])
    y[i] ~ dpois(mu[i])
    new[i] ~ dpois(mu[i])
  }

}
",fill = TRUE)
sink()



sink("ordinal.jags")
cat("
model{

  beta ~ dnorm(0, 1)
  for (j in 1:J){
    alpha.star[j] ~ dnorm(0, 0.001)
  }
  alpha[1:J] <- sort(alpha.star[1:J])

  for (i in 1:n){
    mu[i] = beta * z[i]
    y[i] ~ dordered.logit(mu[i], alpha[])
    new[i] ~ dordered.logit(mu[i], alpha[])
  }

}
",fill = TRUE)
sink()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings for the poisson regression model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(y = y, n = n, z = z)
pars = c('alpha','beta','new')
inits <- function(){list()}  
Sys.time()
library(jagsUI)
# this will take at least a few minutes depending on your machine
m1 <- jags(data = dat, inits=inits, model.file = "poisson.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings for the ordinal regression model
# note that we scale eggs (5 - 9) back to categories (1 - 5)
# by subtracting 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(y = y - 4, n = n, z = z, J = 4)
pars = c('alpha','beta','new')
inits <- function(){list(alpha.star = c(-3,-1,1,3))}  
Sys.time()
# this will take at least a few minutes depending on your machine
m2 <- jags(data = dat, inits=inits, model.file = "ordinal.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()
print(summary(m2), digits = 3)
m2$mcmc.info$elapsed.mins

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) we'll make a plot of expected values for the ordinal model
#    a) explore how we do this, e.g., calc prob for each outcome
#       then multiply prob. by outcome
#    b) why do the confidence intervals get tighter near the ends?
#    c) note that we back-transform axis, ensure you understand how
#
# 2) goodness-of-fit
#
# 3) add an intercept to the ordinal model, e.g.,
#    mu[i] = alpha + beta * z[i]
#    what happens? to traceplots of alpha and the cutpoints...
#    why?
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1st, we'll make a plot of expected values for the ordinal model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
iter <- m2$mcmc.info$n.samples
pz <- seq(min(z), max(z), length.out = res)
em1 <- em2 <- matrix(NA, iter, res)
ec1 <- ec2 <- matrix(NA, iter, res)
qc1 <- qc2 <- matrix(NA, 5, res)

for (j in 1:res){
  
  ec1[,j] <- exp(m1$sims.list$alpha + m1$sims.list$beta * pz[j])
  qc1[,j] <- quantile(ec1[,j], c(0.025, 0.25, 0.5, 0.75, 0.975))  
  
  em2[,j] <- m2$sims.list$beta * pz[j]
  # the expected clutch size is equal to the 
  # probability of a five 
  ec2[,j] <- plogis(m2$sims.list$alpha[,1] - em2[,j]) * 5 + 
    (plogis(m2$sims.list$alpha[,2] - em2[,j]) - plogis(m2$sims.list$alpha[,1] - em2[,j])) * 6 +
    (plogis(m2$sims.list$alpha[,3] - em2[,j]) - plogis(m2$sims.list$alpha[,2] - em2[,j])) * 7 +
    (plogis(m2$sims.list$alpha[,4] - em2[,j]) - plogis(m2$sims.list$alpha[,3] - em2[,j])) * 8 +
    (1 - plogis(m2$sims.list$alpha[,4] - em2[,j])) * 9
  qc2[,j] <- quantile(ec2[,j], c(0.025, 0.25, 0.5, 0.75, 0.975))
}

library(reshape2)
mc1 <- melt(ec1)
mc2 <- melt(ec2)
names(mc1) <- names(mc2) <- c('i','x','y')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here we calculate positions on the x-axis to insert our 
# own tickmarks (back-transforming)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(x)
labX <- seq(round(min(x)),round(max(x)))
backZ <- (labX - mean(x))/sd(x)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
smoothScatter(mc1$y ~ pz[mc$x], ylab = 'Expected clutch size (eggs)',
              xlab = 'Tarsus (mm)', xaxt = 'n', nrpoints = 0, ylim = c(3,12))
lines(qc1[3,] ~ pz, lwd = 3, col = 'white')
lines(qc1[1,] ~ pz, lwd = 3, lty = 2, col = 'white')
lines(qc1[5,] ~ pz, lwd = 3, lty = 2, col = 'white')
axis(side = 1, at = backZ, labels = labX) 
points(y ~ z, cex = 0.5, col = 'navy')


smoothScatter(mc2$y ~ pz[mc$x], ylab = 'Expected clutch size (eggs)',
              xlab = 'Tarsus (mm)', xaxt = 'n', nrpoints = 0, ylim = c(3,12))
lines(qc2[3,] ~ pz, lwd = 3, col = 'white')
lines(qc2[1,] ~ pz, lwd = 3, lty = 2, col = 'white')
lines(qc2[5,] ~ pz, lwd = 3, lty = 2, col = 'white')
axis(side = 1, at = backZ, labels = labX) 
points(y ~ z, cex = 0.5, col = 'navy')






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate goodness-of-fit...
# chi.y <- real data
# chi.s <- simulated data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iter <- m1$mcmc.info$n.samples
chi.y <- matrix(NA, iter, n)
chi.s <- matrix(NA, iter, n)

for (i in 1:n){
  chi.y[,i] <- (y[i] - exp(m1$sims.list$alpha + m1$sims.list$beta * z[i]))^2 /
    exp(m1$sims.list$alpha + m1$sims.list$beta * z[i])
  chi.s[,i] <- (m1$sims.list$new[,i] - exp(m1$sims.list$alpha + m1$sims.list$beta * z[i]))^2 /
    exp(m1$sims.list$alpha + m1$sims.list$beta * z[i])  
}
dim(chi.y)

chi.dat <- rowSums(chi.y)
chi.new <- rowSums(chi.s)

plot(chi.new ~ chi.dat, ylim = c(0,600), xlim = c(0,600))
abline(0,1)
# goodness of fit is terrible... Bayes-p = 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate goodness-of-fit...
# there are lot's and lot's of ways to do this,
# many of them feel hacky to me.
#
# for today we'll just plot one 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,3))
plot(c(m2$sims.list$new[1,] + 4) ~ z, ylim = c(0,20),
     ylab = 'Simulated clutches (ordinal)')

plot(c(m1$sims.list$new[305,]) ~ z, ylim = c(0,20),
     ylab = 'Simulated clutches (poisson)')

plot(y ~ z, ylim = c(0,20),
     ylab = 'Actual clutches')

















