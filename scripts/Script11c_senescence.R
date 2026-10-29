# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read the data directly from github
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
url <- "https://raw.githubusercontent.com/thomasriecke/WILD591/refs/heads/main/data/cs_age_data.csv"
data <- read.csv(url, header = T)
str(data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model(s)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("poisson.jags")
cat("
model{

  alpha ~ dnorm(log(7), 1)
  beta[1] ~ dnorm(0, 1)
  beta[2] ~ dnorm(0, 1)
  
  sigma ~ dgamma(1,1)
  for (k in 1:K){
    epsilon.star[k] ~ dnorm(0,1)
    epsilon[k] = sigma * epsilon.star[k]
  }

  for (i in 1:n){
    mu[i] = exp(alpha + beta[1] * z[i] + beta[2] * z[i] * z[i] + epsilon[x[i]])
    y[i] ~ dpois(mu[i])
    new[i] ~ dpois(mu[i])
  }

}
",fill = TRUE)
sink()



sink("ordinal.jags")
cat("
model{

  beta[1] ~ dnorm(0, 1)
  beta[2] ~ dnorm(0, 1)
  for (j in 1:J){
    alpha.star[j] ~ dnorm(0, 0.001)
  }
  alpha[1:J] <- sort(alpha.star[1:J])
  sigma ~ dgamma(1,1)
  for (k in 1:K){
    epsilon.star[k] ~ dnorm(0,1)
    epsilon[k] = sigma * epsilon.star[k]
  }


  for (i in 1:n){
    mu[i] = beta[1] * z[i] + beta[2] * z[i] * z[i] + epsilon[x[i]]
    y[i] ~ dordered.logit(mu[i], alpha[])
    new[i] ~ dordered.logit(mu[i], alpha[])
  }

}
",fill = TRUE)
sink()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings for the poisson regression model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(y = data$CS, n = nrow(data), 
           z = as.numeric(scale(data$AGE)),
           x = data$YEAR,
           J = 4,
           K = max(data$YEAR))
pars = c('alpha','beta','new','sigma','epsilon')
inits <- function(){list()}  
Sys.time()
library(jagsUI)
# this will take at least a few minutes depending on your machine
m1 <- jags(data = dat, inits=inits, model.file = "poisson.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()
m1$mcmc.info$elapsed.mins

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MCMC settings for the ordinal regression model
# note that we scale eggs (5 - 9) back to categories (1 - 5)
# by subtracting 4
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat = list(y = data$CS-1, n = nrow(data), 
           z = as.numeric(scale(data$AGE)),
           x = data$YEAR,
           J = 4,
           K = max(data$YEAR))
pars = c('alpha','beta','new','sigma','epsilon')
inits <- function(){list(alpha.star = c(-3,-1,1,3))}  
Sys.time()
# this will take at least a few minutes depending on your machine
m2 <- jags(data = dat, inits=inits, model.file = "ordinal.jags",
           parameters.to.save = pars, n.chains = 2, n.iter = 10000, 
           n.thin = 5, n.burnin = 5000, parallel = T)
Sys.time()
print(summary(m2), digits = 3)
m2$mcmc.info$elapsed.mins



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make figure(s)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 100
iter <- m2$mcmc.info$n.samples
za <- scale(data$AGE)
pz <- seq(min(za), max(za), length.out = res)
em1 <- em2 <- matrix(NA, iter, res)
ec1 <- ec2 <- matrix(NA, iter, res)
qc1 <- qc2 <- matrix(NA, 5, res)

for (j in 1:res){
  
  ec1[,j] <- exp(m1$sims.list$alpha + m1$sims.list$beta[,1] * pz[j] +
               m1$sims.list$beta[,2] * pz[j] * pz[j])
  qc1[,j] <- quantile(ec1[,j], c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # calculate for ordinal model
  # note that we multiply prob. of each clutch by real clutch sizes
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  em2[,j] <- m2$sims.list$beta[,1] * pz[j] + m2$sims.list$beta[,2] * pz[j] * pz[j]
  ec2[,j] <- plogis(m2$sims.list$alpha[,1] - em2[,j]) * 2 + 
    (plogis(m2$sims.list$alpha[,2] - em2[,j]) - plogis(m2$sims.list$alpha[,1] - em2[,j])) * 3 +
    (plogis(m2$sims.list$alpha[,3] - em2[,j]) - plogis(m2$sims.list$alpha[,2] - em2[,j])) * 4 +
    (plogis(m2$sims.list$alpha[,4] - em2[,j]) - plogis(m2$sims.list$alpha[,3] - em2[,j])) * 5 +
    (1 - plogis(m2$sims.list$alpha[,4] - em2[,j])) * 6
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
summary(za)
labX <- seq(2,26,by = 2)
backZ <- (labX - mean(data$AGE))/sd(data$AGE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a figure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
par(mfrow = c(1,2), mar = c(5,5,2,2))
smoothScatter(mc1$y ~ pz[mc1$x], ylab = 'Expected clutch size (eggs)',
              xlab = 'Age', xaxt = 'n', nrpoints = 0, ylim = c(1.5,6.5))
lines(qc1[3,] ~ pz, lwd = 3, col = 'white')
lines(qc1[1,] ~ pz, lwd = 3, lty = 2, col = 'white')
lines(qc1[5,] ~ pz, lwd = 3, lty = 2, col = 'white')
axis(side = 1, at = backZ, labels = labX) 
points(jitter(data$CS) ~ za, cex = 0.5, col = 'navy')


smoothScatter(mc2$y ~ pz[mc2$x], ylab = 'Expected clutch size (eggs)',
              xlab = 'Age', xaxt = 'n', nrpoints = 0, ylim = c(1.5,6.5))
lines(qc2[3,] ~ pz, lwd = 3, col = 'white')
lines(qc2[1,] ~ pz, lwd = 3, lty = 2, col = 'white')
lines(qc2[5,] ~ pz, lwd = 3, lty = 2, col = 'white')
axis(side = 1, at = backZ, labels = labX) 
points(jitter(data$CS) ~ za, cex = 0.5, col = 'navy')



par(mfrow = c(1,1))
plot(jitter(data$CS) ~ jitter(data$AGE), xlab = 'Female age',
     ylab = 'Clutch Size', col = 'navy', cex = 0.5)





par(mfrow = c(1,2), mar = c(5,5,2,2))
smoothScatter(mc1$y ~ pz[mc1$x], ylab = 'Expected clutch size (eggs)',
              xlab = 'Age', xaxt = 'n', nrpoints = 0, ylim = c(0,15))
lines(qc1[3,] ~ pz, lwd = 3, col = 'white')
lines(qc1[1,] ~ pz, lwd = 3, lty = 2, col = 'white')
lines(qc1[5,] ~ pz, lwd = 3, lty = 2, col = 'white')
axis(side = 1, at = backZ, labels = labX) 
points(jitter(m1$sims.list$new[1,]) ~ jitter(za), cex = 0.5, col = 'navy')


smoothScatter(mc2$y ~ pz[mc2$x], ylab = 'Expected clutch size (eggs)',
              xlab = 'Age', xaxt = 'n', nrpoints = 0, ylim = c(0,15))
lines(qc2[3,] ~ pz, lwd = 3, col = 'white')
lines(qc2[1,] ~ pz, lwd = 3, lty = 2, col = 'white')
lines(qc2[5,] ~ pz, lwd = 3, lty = 2, col = 'white')
axis(side = 1, at = backZ, labels = labX) 
points(jitter(m2$sims.list$new[1,]+1) ~ jitter(za), cex = 0.5, col = 'navy')


