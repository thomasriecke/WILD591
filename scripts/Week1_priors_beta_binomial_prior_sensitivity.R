library(vioplot)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# example from slides
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res <- 500
prior <- dbeta(seq(0,1,length.out = res), 1, 1)
data <- dbeta(seq(0,1,length.out = res), 137, 863)
posterior <- dbeta(seq(0,1,length.out = res), 138, 864)

par(mfrow = c(1,1), mar = c(5.1,5.1,2.1,2.1))
plot(prior ~ seq(0,1,length.out = res), type = 'l', ylim = c(0,40), lwd = 2,
     col = 'red', xlab = 'f', ylab = 'Density', cex.lab = 2)
points(data ~ seq(0,1,length.out = res), type = 'l', lwd = 2)
points(posterior ~ seq(0,1,length.out = res), type = 'l', lwd = 2, 
       col = 'blue')
legend('topright', c('Prior','Data','Posterior'),
       lwd = 2, col = c('red','black','blue'))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's imagine we're going to mark and release 1000 ducks (k)...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n <- 1000

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# we'll assume that our 'true' band-recovery probability is 0.05, 
# or that hunters will shoot, recover, and report 5% of the total 
# population (we'll assume that our marked population is representative)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
f <- 0.05

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll generate some data (y), or the number of ducks that are actually
# killed, retrived, and reported by hunters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
?rbinom()
y <- rbinom(1, n, f)

# sample 1000 marked birds with same f for 50 years (z)
z <- rbinom(50, n, f)
hist(z, breaks = 20)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# great, now we have our data!
# let's consider three priors...
# the two numbers are the alpha and beta parameters for a beta distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
weak <- c(1,1)
mod <- c(1,9)
strong <- c(100,900)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's visualize our priors!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samples <- 10000
vioplot(rbeta(samples, weak[1], weak[2]),
        rbeta(samples, mod[1], mod[2]),
        rbeta(samples, strong[1], strong[2]),
        names = c('Beta(1,1)', 'Beta(1,9)','Beta(5,95)'))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lets generate some posteriors (along with an MC maximum likelihood estimate)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
samples <- 10000
p.weak <- rbeta(samples, weak[1]+y, weak[2]+(n-y))
p.mod <- rbeta(samples, mod[1]+y, mod[2]+(n-y))
p.strong <- rbeta(samples, strong[1]+y, strong[2]+(n-y))
ml.est <- rbeta(samples, y, (n-y))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's generate some summary statistics of the posterior distribution
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
q.weak <- quantile(p.weak, c(0.025,0.5,0.975))
q.mod <- quantile(p.mod, c(0.025,0.5,0.975))
q.strong <- quantile(p.strong, c(0.025,0.5,0.975))
q.ml <- quantile(ml.est, c(0.025,0.5,0.975))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's plot our posteriors and our MLE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vioplot(p.weak, p.mod, p.strong, ml.est,
        drawRect = F, ylab = 'Estimates of band-recovery probability',
        xlab = 'Priors (and non-Bayes)', line = 3, cex.lab = 2,
        names = c('Beta(1,1)','Beta(1,9)','Beta(5,95)','ML'), ylim = c(0,0.1))
arrows(1:4, c(q.weak[1],q.mod[1],q.strong[1],q.ml[1]), 
       1:4, c(q.weak[3],q.mod[3],q.strong[3],q.ml[3]),
       length = 0)
points(c(q.weak[2],q.mod[2],q.strong[2],q.ml[2]), pch = 21, bg = 'white', cex = 2)

q.weak[2];q.mod[2];q.strong[2];q.ml[2]



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# using a prior from the literature,
# imagine an estimate of 0.0487 with an sd of 0.01 has been published
# and you want to use that (or something like that) as hyperpriors
# for a beta distribution.... how would you do that?
#
# moment matching! 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mu <- 0.0487
sd <- 0.01

alpha = ((1-mu)/(sd*sd) - (1/mu)) * mu^2
beta = alpha * (1/mu - 1)

pub <- rnorm(samples, 0.0487, 0.01)
pri <- rbeta(samples, alpha, beta)

q.pub <- quantile(pub, c(0.025,0.5,0.975))
q.pri <- quantile(pri, c(0.025,0.5,0.975))

vioplot(pub, pri, drawRect = F,
        names = c('Normal(0.0487, sd = 0.01)','Beta(22.51,439.77)'))
arrows(1:2, c(q.pub[1],q.pri[1]), 1:2, c(q.pub[3],q.pri[3]), length = 0, lty = 1, col = 'white')
points(c(q.pub[2],q.pri[2]) ~ c(1:2), pch = 21, bg = 'white', cex = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# HOMEWORK:
#
# 1) generate a new posterior distribution given the normal distribution from the literature and moment matching!
#
# 2) play with priors! try >5 different priors to generate posteriors... try some really bad ones!
#    assume you think you know everything and you're deeply wrong... see what happens :)
#    simple models like this are the easiest to use to understand the effects of bad priors
#
# 3) Visualize, visualize, visualize! My favorite way to get to know distributions (after
#    a cursory inspection of expected mean and variance) is to simply plot them, e.g.,
hist(rbeta(10000, 22, 78), breaks = 1000)
#    this will help you rapidly understand relationships between values used in the distribution
#    e.g., the mean of a beta is equal to alpha/(alpha + beta), you can see it in the histogram
#    try 50 and 50, 75 and 25, 5 and 5, etc. as the sum of the two numbers decline variance will
#    increase, as alpha gets bigger relative to beta the mean will increase, etc. Make some panel
#    plots, etc.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

