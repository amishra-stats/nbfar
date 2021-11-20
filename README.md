# Negative binomial factor regression models with application to microbiome data analysis

The human microbiome provides essential physiological functions and helps maintain host homeostasis via the 
formation of intricate ecological host-microbiome relationships. The R package **nbfar** implements 
factor regression models that allow the estimation of structured parsimonious associations between 
host-related features X and amplicon-derived microbial taxa Y. 
To account for the overdispersed nature of the amplicon sequencing count data, we propose
**Negative Binomial reduced rank regression (NB-RRR)** and **Negative Binomial co-sparse factor regression (NB-FAR)**.

Both of these models encode the underlying dependency among the microbial abundance data as outcomes and 
the host-associated features as predictors through a structured coefficient matrix. 
NB-RRR does so via a rank constrained coefficient matrix whereas NB-FAR does so via a sparse singular value decomposition (SSVD) 
of the coefficient matrix. The R package **nbfar** provides an efficient implementation of the proposed 
procedure for parameter estimation. Here are the details of the R subroutines:

- **nbfar** : Estimate the parameters of the negative binomial factor regression model
- **nbrrr** : Estimate the parameters of the negative reduced rank regression model 

<img src="https://i.imgur.com/ytq5qZK.jpg" alt="nbfar" height="600" align="left"/>


# Simulated examples 

We showcase the usage of **nbfar** and  **nbrrr** on  simulated data. 


### Data simulation

```
rm(list = ls())


# Install packages 
# install.packages("Rcpp", repos="https://rcppcore.github.io/drat")
# devtools::install_github('amishra-stats/nbfar/nbfar', force = TRUE)
# load library
library(nbfar)
library(RcppParallel)

```

Check check

```
## ## -----------------------Simulation settings --------------------
## Simulation setting:
## n: sample size
## p: number of predictors
## pz: number of control variable other than intercept parameter
## q: number of outcome variables
## nrank: true rank of the model
## snr: signal to noise ratio to be used in generating gaussian outcome variables
## nlam: number of lambda values to be fitted
## rank.est: maximum estimated rank to be specified to the model
## s: multiplicative factor of singular values
## nthread: number of parallel thread can be used in  parallel for cross validation  
## The simulation was replicated 100 times under each setting as detailed in the paper.

#
## Model specification:
p <- 50;
example_seed <- 123
xp = 30
set.seed(example_seed)
n <- 200
nrank <- 3                # true rank
rank.est <- 5             # estimated rank
nlam <- 40                # number of tuning parameter
s  = 0.5; q <- 30
sp  = xp/p
nthread = 1



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
V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8, replace = TRUE) *
              runif(8, 0.3, 1), rep(0, q - 16))
V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8, replace = TRUE) *
              runif(8, 0.3, 1), rep(0, q - 28))
V[, 3] <- c( sample(c(1, -1), 5, replace = TRUE) *
               runif(5, 0.3, 1), rep(0, 23),
             sample(c(1, -1), 2, replace = TRUE) *
               runif(2, 0.3, 1), rep(0, q - 30))
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
sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = 0.75, depth = 10)  # Simulated sample
X <- sim.sample$X[1:n, ]                    # simulated predictors (training)
Y <- sim.sample$Y[1:n, ]                    # simulated responses (training)
# 1000 test sample data
sim.sample <- nbfar_sim(U, D, V, 1000, Xsigma, C0, disp = 0.75, depth = 10)
Xt <- sim.sample$X                    # simulated predictors (test)
Yt <- sim.sample$Y                    # simulated predictors (test)
X0 <- cbind(1, X)                     # 1st column accounting for intercept


# Simulate data with 20% missing entries
miss <- 0.10          # Proportion of entries missing
t.ind <- sample.int(n * q, size = miss * n * q)
y <- as.vector(Y)
y[t.ind] <- NA
Ym <- matrix(y, n, q)   # 20% of entries are missing at random

```

###  Negative binomial reduced rank regression: **nbrrr** 


```

# Model fit: (full data)
set.seed(example_seed)
control_r3 <- nbfar_control(initmaxit = 10000, initepsilon = 1e-5,
                            objI = 1)
nbrrr_test <- nbrrr(Y, X, maxrank = 5, control = control_r3, nfold = 5, trace = F)


# Model fit:  (missing data)
set.seed(example_seed)
control_r3 <- nbfar_control(initmaxit = 10000, initepsilon = 1e-5,
                            objI = 1)
nbrrr_testm <- nbrrr(Ym, X, maxrank = 5, control = control_r3, nfold = 5,trace = F)


```
###  Negative co-sparse factor  regression: **nbfar** 

```
# Model fit: (full data)
RcppParallel::setThreadOptions(numThreads = nthread)
set.seed(example_seed)
control_nbfar <- nbfar_control(gamma0 = 1, spU = sp, spV = 20/q,
                               maxit = 2000, lamMaxFac = 1e-2,
                               lamMinFac = 1e-7, epsilon = 1e-4,
                               objI = 0,
                               initmaxit = 10000, initepsilon = 1e-7)
nbfar_test <- nbfar(Y, X, maxrank = rank.est, nlambda = nlam,
                    cIndex = NULL,
                    ofset = NULL, control = control_nbfar, nfold = 5,
                    PATH = FALSE, nthread = nthread,trace = F)

# Model fit: (missing data)
RcppParallel::setThreadOptions(numThreads = nthread)
set.seed(example_seed)
control_nbfar <- nbfar_control(gamma0 = 1, spU = sp, spV = 20/q,
                               maxit = 2000, lamMaxFac = 1e-2,
                               lamMinFac = 1e-7, epsilon = 1e-4,
                               objI = 0,
                               initmaxit = 10000, initepsilon = 1e-7)
nbfar_testm <- nbfar(Ym, X, maxrank = rank.est, nlambda = nlam,
                    cIndex = NULL,
                    ofset = NULL, control = control_nbfar, nfold = 5,
                    PATH = FALSE, nthread = nthread,trace = F)
```




## Queries
Please contact authors and creators for any queries related to using the analysis 


-   Aditya Mishra: [mailto](mailto:amishra@flatironinstitute.org)
-   Christian Mueller: [mailto](mailto:cmueller@flatironinstitute.org)
