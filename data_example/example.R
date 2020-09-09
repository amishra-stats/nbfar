###########################################################
## R code for simulation examples and applications       ##
## in "Negative binomial co-sparse factor regression"    ##
###########################################################

## initial value NBRRR: dispesionn parameter: moment estimaate

rm(list = ls())
## Install and load the package "gofar" from GitHub
library(devtools)
library(nbfar)

## Load required packages
library(MASS)
require(ggplot2)
require(reshape2)
library(magrittr)
# library(rrpack)
library(glmnet)
library(tictoc)


# Rcpp::sourceCpp('src/nbfar.cpp')
# source('R/nbfar.R')

## --------------------------------------------------------
## Simulation setting:
## n: sample size
## p: number of predictors
## pz: number of control variable other than intercept parameter
## q1: number of gaussian responses
## q2: number of bernoulli responses
## q3: number of poisson responses
## nrank: true rank of the model
## snr: signal to noise ratio to be used in generating gaussian outcome variables
## nlam: number of lambda values to be fitted
## rank.est: maximum estimated rank to be specified to the model
## s: multiplicative factor of singular values

## The simulation was replicated 100 times under each setting as detailed in the paper.


## ## -----------------------Simulation settings --------------------
# Setting specific to simulated examples:

# Setup I
p <- 100; n <- 200
## The following code is for a single replication.
## ----------------------single replication ------------------------
## Other parameters
SD <- 123
set.seed(SD)
pz <- 0
nrank <- 3                # true rank
rank.est <- 5             # estimated rank
nlam <- 40                # number of tuning parameter
snr <- 0.50               # SNR for variance Gaussian error
s  = 0.5
q <- 30
control <- nbfar_control()  # control parameters




## -------------- Generate low-rank and sparse coefficient matrix ---
## D: singular values
## U: sparse left singular vectors
## V: sparse right singular vectors

D <- rep(0, nrank)
V <- matrix(0, ncol = nrank, nrow = q)
U <- matrix(0, ncol = nrank, nrow = p)
#
U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#

# for similar type response type setting
V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8, replace =
                                TRUE) * runif(8, 0.3, 1), rep(0, q - 16))
V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8, replace =
                                 TRUE) * runif(8, 0.3, 1), rep(0, q - 28))
V[, 3] <- c(
  sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
  sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
)


U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#
D <- s * c(4, 6, 5) # signal strength varries as per the value of s
or <- order(D, decreasing = T)
U <- U[, or]
V <- V[, or]
D <- D[or]
C <- U %*% (D * t(V)) # simulated coefficient matrix
intercept <- rep(0.5, q) # specifying intercept to the model:
C0 <- rbind(intercept, C)


## ----- Simulate data -----
Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = 3, depth = 10)  # Simulated sample
X <- sim.sample$X[1:n, ]                    # simulated predictors (training)
Y <- sim.sample$Y[1:n, ]                    # simulated responses (training)
sum(Y==0)
# 1000 test sample data
sim.sample <- nbfar_sim(U, D, V, 1000, Xsigma, C0, disp = 3, depth = 10)
Xt <- sim.sample$X                    # simulated predictors (test)
Yt <- sim.sample$Y                    # simulated predictors (test)

X0 <- cbind(1, X)                     # 1st column accounting for intercept
# THETA <- X0 %*% C0
# exp(0.5*range(THETA))

# Simulate data with 20% missing entries
miss <- 0.20          # Proportion of entries missing
t.ind <- sample.int(n * q, size = miss * n * q)
y <- as.vector(Y)
y[t.ind] <- NA
Ym <- matrix(y, n, q)   # 20% of entries are missing at random





## ----- Negative binomial reduced rank reggression
# Model with offset
set.seed(1234)
offset <- log(10)*matrix(1,n,q)
control_nbrr <- nbfar_control(initmaxit = 5000, initepsilon = 1e-4)
# Model with known offset
nbrrr_test_of <- nbrrr(Y, X, maxrank = 5, cIndex = NULL, ofset = offset,
                       control = control_nbrr, nfold = 5)

# Model without offset
set.seed(1234)
nbrrr_test <- nbrrr(Y, X, maxrank = 5, control = control_nbrr, nfold = 5)

save(list = ls(), file = '../model_check.rda')


plot(colMeans(nbrrr_test$cv.err))
plot(log(nbrrr_test$diffobj[nbrrr_test$diffobj > 0]))


## ----- Negative binomial co-sparse factor regression
# Model with offset
set.seed(1234)
control_nbfar <- nbfar_control(initmaxit = 5000, gamma0 = 2,
                               spU = 0.5, spV = 0.6,
                               lamMinFac = 1e-10, epsilon = 1e-5)
nbfar_test_of <- nbfar(Y, X, maxrank = 5, nlambda = 20, cIndex = NULL,
                       ofset = offset, control = control_nbfar, nfold = 5,
                       PATH = F)
plot(colMeans(nbfar_test_of$fit[[1]]$dev))



set.seed(1234)
nbfar_test <- nbfar(Y, X, maxrank = 5, nlambda = 20, cIndex = NULL,
                    ofset = NULL, control = control_nbfar, nfold = 5,
                    PATH = FALSE)
plot(colMeans(nbfar_test$fit[[1]]$dev))
save(list = ls(), file = '../model_check.rda')
































# Model fitting begins:
control$epsilon <- 1e-7                       # error tolerance
control$spU <- 20/p;                          # maximum support of U
control$spV <- 20/q;                          # maximum support of V
if (p > 50) control$spU <- 0.5
control$maxit <- 2000                         # maximum number of itreration
control$objI <- 1                             # convergence on the basis of objective function
control$gamma0 <- 1                           # power factor for constructing weights
control$se1 <- 1                              # 1 satndard deviation rule for cross validation


## Initalize NBRRR model for only inttercept in the model
## count distribution model
## Compute null deviance of tthe model


# Model fitting
# GOFAR(S) (full data)

set.seed(SD)
tic()
fit.seq <- gofar_s(Y, X, nrank = rank.est, family = family,
                   nlambda = nlam, familygroup = familygroup,
                   control = control, nfold = 5, PATH = TRUE)
fit.seq$time = toc()


# GOFAR(S) (missing data)
set.seed(SD)
tic()
fit.seq.m <- gofar_s(Ym, X, nrank = rank.est, family = family,
                     nlambda = nlam, familygroup = familygroup,
                     control = control, nfold = 5)
fit.seq.m$time = toc()

