library(reshape2)
library(vioplot)
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# a simple linear model, lm()
# with raw length data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(m1 <- lm(y ~ x))

plot(y ~ x, ylab = 'Mass (lbs)', xlab = 'Length (ft)', 
     xlim = c(-1,10), ylim = c(-100,200),
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
abline(m1$coefficients[1], m1$coefficients[2], lwd = 2, col = 'black')
points(0, m1$coefficients[1], pch = 21, bg = 'red', cex = 3) # the intercept

# a second point at a centered intercept
mean(x)
points(mean(x), m1$coefficients[1] + m1$coefficients[2] * mean(x), 
       pch = 21, bg = 'dodgerblue4', cex = 4)


# this code just shows the best fit if the intercept was fixed to 0...
# plot(y ~ x, ylab = 'Mass (lbs)', xlab = 'Length (ft)', 
#      xlim = c(-1,10), ylim = c(-100,200),
#      pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
# summary(m1b <- lm(y ~ 0 + x))
# abline(0, m1b$coefficients[1], lwd = 2, col = 'black')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll center x...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c <- x - mean(x)

plot(y ~ c, ylab = 'Mass (lbs)', xlab = 'Distance from average length (ft)', 
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)

summary(m2 <- lm(y ~ c))
# the intercept changed, but the slope is the same? why?

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll transform x (length in feet) to length in mm (m)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
m <- c * 25.4 * 12

plot(y ~ m, ylab = 'Mass (lbs)', xlab = 'Distance from average length (mm)', 
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
summary(m3 <- lm(y ~ m))
# note that the intercept is the same, but the slope changed dramatically. Why?

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now we'll z-standardize
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
z <- (x - mean(x))/sd(x)
plot(y ~ z, ylab = 'Mass (lbs)', xlab = 'Z-standardized length', 
     pch = 21, bg = 'forestgreen', cex = 1.5, cex.lab = 2, las = 1)
summary(m4 <- lm(y ~ z))
# The intercept is the same, the slope is slighlty different (it's now 
# on the scale of standard deviations!)















# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Let's write a Stan model (or two...)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# so... now we have the same covariate on four different scales:
# x: alligator measurements in feet
# c: 'centered' alligator measurements in feet
# m: 'centered' alligator measurements in millimeters
# z: 'z-standardized' alligator measurements in feet
#
# We have a measurement of a response variable (y; mass in lbs)
# 
# and we have a sample size (n = 200)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# the first step is writing the model itself, to do this we'll 'sink'
# a text file into our working directory that contains the code necessary
# to run the model. This model code can also be written and 
# stored in a separate file (via Notepad or any text editor),
# but I find it easiest to find my mistakes if it's write here in the code
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Specify model in Stan
library(rstan)
m5.stan <- stan_model(model_code = 
"                        
data {
  int n;                 // in Stan we write comments like this 
  vector[n] y;           // every line has to end with a semi-colon, or you'll get an 'ill-formed' error
  vector[n] x;           
}

parameters {
  // we also 
  real<lower=-100, upper = 100> beta0;     // we can put limits on the potential values of things
  real<lower=-100, upper = 100> beta1;
  real<lower=0> sigma;
}

// in JAGS we can just lump everything willy-nilly into the model file
// Stan requires some careful thought about what a parameter is
// what the data are, what parameters are 'transformed,' or simply functions
// of other parameters
// and what our prior and data models look like
// Stan can also incorporate 'functions' (so can JAGS, we just aren't going 
// to ever get into it), 'transformed data', and 'generated quantities'
// blocks of code

transformed parameters {
  vector[n] mu;
  mu = beta0 + beta1 * x;
  // this could also be written inside a for-loop, e.g., 
  // for (i in 1:n){
  //   mu[i] = beta0 + beta1 * x[i]
  // }
  // but Stan can just vectorize it (and it is faster to do this!)
}

model {
  // priors
  sigma ~ gamma(1,1);
  beta0 ~ normal(0,100);
  beta1 ~ normal(0,100);
  // data model
  y ~ normal(mu, sigma);
}
")
# It will take some time to compile this model in Stan!
# generally longer than it takes to run an entire JAGS model
# however, when running more complex models Stan will sample much, much,
# more efficiently!


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# here is where we provide JAGS with the data it will need to run the model,
# in this instance it needs y (our response variable, mass), x (our covariate, length
# in feet), and n (our sample size)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stan_data <- list(y = y, x = x, n = n)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we provide initial values. This is where JAGS will begin sampling
# to build a posterior distribution, in this case, I won't provide any 
# initial values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list()}  

# sample
nc <- 4
nt <- 5
ni <- 5000
nb <- 2500


# Call JAGS from R 
# 50k iterations takes 1m on an i9
library(jagsUI)
Sys.time()


m5 <- sampling(m5.stan, 
               data = stan_data,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               open_progress = FALSE, refresh = -1, cores = 4)
# here I'm using 4 cores to sample in parallel
Sys.time()

###
### all of our results from this run are now stored in list 'm5'
###

# this command will print a summary of parameter estimates
print(m5)
# Stan will automatically save everything if we don't specify what to save...

m5s <- extract(m5) # this will generate a list similar to the 'sims.list$
# object generated by jags

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make a plot of the effect of length on mass
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# let's visualize our estimate of beta1, or the relationship between mass
# and length in feet, first we'll plot the posterior with vioplot()
vioplot(m5s$beta1, drawRect = F,
        names = c(''))
# add an axis label
mtext(expression(beta[1]), side = 2, line = 2.5, cex = 2)
# draw a 95% credible interval
arrows(1, quantile(m5s$beta1, c(0.025)), 1, quantile(m5s$beta1, c(0.975)), length = 0, col = 'white', lwd = 2)
# add a point at the median of the distribution
points(quantile(m5s$beta1, c(0.5)) ~ 1, cex = 3, pch = 21, bg = 'white')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# make the same plots for beta0 (and sigma)!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# now let's make a plot of the expected value of mass across a range of
# potential covariate values
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define the resolution of your plot
res <- 100
# create a seq of covariate values to predict from
px <- seq(3,10,length.out = res)
# the length of the posterior distribution = (nc * (ni-nb))/nt
n.iter <- length(m5s$beta0)

# empty matrix to store results
py <- matrix(NA, n.iter, res)
# empty matrix to store quantile
qy <- matrix(NA, res, 5)

hist(m5s$beta0)

for (j in 1:res){
  py[,j] <- m5s$beta0 + m5s$beta1 * px[j]
  qy[j,] <- quantile(py[,j], c(0.025,0.05,0.5,0.95,0.975))
}

# melt the matrix to plot 
mpy <- melt(py)
names(mpy) <- c('ni','x','y')

# make a smoothScatter plot (similar alternatives exist in ggplot, etc.)
# note the use of 'nested indexing' to place values on the x axis. In the 
# mpy data.frame, we have a 
smoothScatter(mpy$y ~ px[mpy$x],
              colramp = colorRampPalette(c('white','forestgreen')),
              xlab = 'Length (ft)', ylab = 'Mass (lbs)',
              nrpoints = 0, cex.lab = 2, las = 1)
lines(qy[,3] ~ px, col = 'white', lwd = 2)
lines(qy[,1] ~ px, col = 'white', lwd = 2, lty = 2)
lines(qy[,5] ~ px, col = 'white', lwd = 2, lty = 2)

points(y ~ x, pch = 1, col = 'grey20', cex = 1.5)

boxplot(py[,2], py[,16])
dim(py)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Homework
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# 1) build two additional models, model 6 and model 7.
#    for model 6, try using the centered transformation of the covariate (c)
#    for model 7, try using the z-standardized transformation of the covariate (z)
#    
#    General guidance, remember that you will have to change the data you supply
#    to JAGS in the jags.data list (e.g., c rather than x)
#
#    Make the same plots as were made for m5. Think about the scale of the covariate
#    used to predict
#
# 2) Go back to the original version of model 5 [above]
#    Purposefully mess something up, e.g., remove a key piece of data from
#    jags.data list(). What error message do you receive? Fix the error, and 
#    attempt to create a new error. This may seem a little goofy [don't do this
#    for hours! :)], but intentionally creating problems to see how programs
#    respond can help you learn to effectively error trap.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



hist(c(2,4), breaks = 100, xlim = c(0,10), cex.lab = 2,
     ylim = c(0,3), main = NULL, xlab = 'Population means')












