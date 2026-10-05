library(latex2exp)
library(jagsUI)
library(MuMIn)
library(loo)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# simulate some data
# e: number of eggs laid
# b: a z-standardized index of body size...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)

# sample size
n <- 200
b <- rnorm(n, 0, 1)
alpha <- 1
beta <- 0.5
psi <- exp(alpha + beta*b)

e <- round(psi)

par(mfrow = c(1,1))
plot(jitter(e) ~ b, ylab = 'Number of Eggs (e)', xlab = 'Body size (b)', las = 1, cex.lab = 1.5)
mean(e)
var(e)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AIC (frequentist/information theoretic)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m0 <- glm(e ~ 1, family = 'poisson')
m1 <- glm(e ~ b, family = 'poisson')
AICc(m0)
AICc(m1)

summary(m0)
summary(m1)

p0 <- predict(m0)
p1 <- predict(m1)

plot(e ~ b, ylab = 'Number of Eggs (e)', xlab = 'Body size (b)', las = 1, cex.lab = 1.5)
points(exp(p0) ~ b, col = 'purple')

plot(e ~ b, ylab = 'Number of Eggs (e)', xlab = 'Body size (b)', las = 1, cex.lab = 1.5)
points(exp(p1) ~ b, col = 'red')

mean(e)
var(e)

res <- e - exp(p1)
plot(res ~ b)
var(res)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# testing the poisson
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sink("clutch_size.jags")
cat("
      model {

      alpha ~ dnorm(1, 1)
      beta ~ dnorm(0, 0.1)


      for (i in 1:n){
        psi[i] = exp(alpha + beta * b[i])
        e[i] ~ dpois(psi[i])
        new[i] ~ dpois(psi[i])
        
        logL[i] = log(dpois(e[i], psi[i]))
  
      }

      
      }
      ",fill = TRUE)
sink()


sink("clutch_size_no_covariate.jags")
cat("
      model {

      alpha ~ dnorm(1, 1)


      for (i in 1:n){
        psi[i] = exp(alpha)
        e[i] ~ dpois(psi[i])
        new[i] ~ dpois(psi[i])
        
        logL[i] = log(dpois(e[i], psi[i]))
      }

      }
      ",fill = TRUE)
sink()



jags.data <- list(e = e, b = b, n = n)

# Initial values
inits <- function(){list()}  

# Parameters monitored
parameters <- c('alpha','beta','fit','fit.new','new','logL')

nc <- 4
nt <- 10
ni <- 25000
nb <- 10000

library(jagsUI)
Sys.time()

m0 <- jags(jags.data, inits, parameters, "clutch_size_no_covariate.jags", 
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
           parallel = T)


m1 <- jags(jags.data, inits, parameters, "clutch_size.jags", 
          n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
          parallel = T)
Sys.time()


print(m0, digits = 3)
print(m1, digits = 3)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# model selection via DIC
# will return nearly identical results to AIC with simple model types
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m0$DIC
m1$DIC
# the model with the covariate (m1) explains far more of the information in the data

# deviance
m0$q50$deviance
m1$q50$deviance

# calculation for effective number of parameters
m0.pD = var(m0$sims.list$deviance)/2
m0$pD

m1.pD = var(m1$sims.list$deviance)/2
m1$pD





par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bayesian p-values
# Pearson's residuals
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
eps0 <- matrix(NA, length(m0$sims.list$alpha), n)
eps0n <- matrix(NA, length(m0$sims.list$alpha), n)
eps1 <- matrix(NA, length(m1$sims.list$alpha), n)
eps1n <- matrix(NA, length(m1$sims.list$alpha), n)
for (j in 1:n){
  # Pearson's residuals for the null model (eps0) and the new data (eps0n)
  eps0[,j] <- (e[j] - exp(m0$sims.list$alpha))/sd(e)
  eps0n[,j] <- (m0$sims.list$new[,j] - exp(m0$sims.list$alpha))/sd(m0$sims.list$new[,j])  
  
  # Pearson's residuals for the model with a covariate (eps1) and the new data (eps1n)
  eps1[,j] <- (e[j] - exp(m1$sims.list$alpha + m1$sims.list$beta * b[j]))/sd(e)
  eps1n[,j] <- (m1$sims.list$new[,j] - exp(m1$sims.list$alpha + m1$sims.list$beta * b[j]))/sd(m1$sims.list$new[,j])    
}

fit0 <- rowSums(eps0^2)
fit0n <- rowSums(eps0n^2)

fit1 <- rowSums(eps1^2)
fit1n <- rowSums(eps1n^2)

par(mfrow = c(1,1))
plot(fit0n ~ fit0, xlim = c(0,300), ylim = c(0,300), 
     ylab = 'RSS of generated (new) data', xlab = 'RSS of actual data (e)')
abline(0,1)
mean(fit0n > fit0)

plot(fit1n ~ fit1, xlim = c(0,300), ylim = c(0,300), 
     ylab = 'RSS of generated (new) data', xlab = 'RSS of actual data (e)')
abline(0,1)
mean(fit1n > fit1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) So... the model with a covariate was clearly superior to the null 
#    (i.e., intercept only model) when using AICc or DIC...
#    But it seems to have terrible fit to the data...?!
#    Why?! Here's a hint from the third iteration of the model
par(mfrow = c(1,3))
plot(m0$sims.list$new[3,] ~ b, ylab = 'New data from m0')
plot(e ~ b, ylab = 'Actual data')
plot(m1$sims.list$new[3,] ~ b, ylab = 'New data from m1')
#    
# 2) Below we will actually calculate the log-likelihood by hand
#    (across all data points [DIC] and for EACH data point [WAIC])
#    review those equations, ask questions! It's useful to understand
#    what these programs are doing inside the 'black box'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Deviance information criterion
# note JAGS calculates as D + pD
# elsewhere often D + 2*pD
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.iter <- nc * (ni-nb)/nt
lik <- NULL
E <- NULL
for (i in 1:n.iter){
  E <- exp(m1$sims.list$alpha[i] + m1$sims.list$beta[i] * b)
  lik[i] <- sum(e * log(E)) - sum(E) - sum(log(factorial(e)))
}
hist(lik)
D <- -2 * lik
plot(m1$sims.list$deviance ~ D)

pD.hand <- var(D)/2

jags.DIC <- mean(D) + pD.hand
m$DIC
DIC <- mean(D) + 2*pD.hand



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WAIC
# Gelman, Hwang, and Vehtari (2013) Statistics and Computing
# here we will calculate the log-likelihood for EACH DATA POINT (not the sum
# across data points)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n.iter <- nc * (ni-nb)/nt
lik <- matrix(NA, n.iter, n)
E <- NULL

for (i in 1:n.iter){
  E <- exp(m1$sims.list$alpha[i] + m1$sims.list$beta[i] * b)
  for (j in 1:n){
    lik[i,j] <- sum(e[j] * log(E[j])) - sum(E[j]) - sum(log(factorial(e[j])))    
  }
}


waic(m1$sims.list$logL)
waic(lik)

lppd <- sum(log(colMeans(exp(m1$sims.list$logL))))
p.waic <- sum(apply(m1$sims.list$logL, 2, var))
waic <- -2*(lppd - p.waic)


par(mfrow = c(1,2))
plot(e ~ b, ylab = 'Number of Eggs (e)', xlab = 'Body size (b)', 
     ylim = c(0,25), type = 'n')
points(m0$q50$new ~ b, col = 'purple')
arrows(b, m0$q2.5$new, b, m0$q97.5$new, length = 0, 
       lty = 2, lwd = 0.5, col = 'purple')
points(e ~ b)

plot(e ~ b, ylab = 'Number of Eggs (e)', xlab = 'Body size (b)', 
     ylim = c(0,25), type = 'n')
points(m1$q50$new ~ b, col = 'red')
arrows(b, m1$q2.5$new, b, m1$q97.5$new, length = 0, 
       lty = 2, lwd = 0.5, col = 'red')
points(e ~ b)



