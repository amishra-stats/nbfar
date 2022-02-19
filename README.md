# Negative binomial factor regression models with application to microbiome data analysis

The R package `nbfar` implements Negative Binomial factor regression models that allow the estimation of 
structured (sparse) associations between a feature matrix X and overdispersed count data Y. 

The package has been developed with microbiome count data Y in mind and can be used, e.g., to associate 
host or environmental covariates with microbial abundances.

Currently, two models are available

- **Negative Binomial reduced rank regression (NB-RRR)**
- **Negative Binomial co-sparse factor regression (NB-FAR)**.

The underlying structure of the models are illustrated in the Figure below.

<img src="https://i.imgur.com/ytq5qZK.jpg" alt="nbfar" height="80%" align="center"/>

Both models result in a type of joint biclustering structure linking features to count responses.

Microbiome data example - linking host phenotype data to microbial abundances of the American Gut project
--------------

Using `nbfar` we analyzed the [American Gut Project data](https://journals.asm.org/doi/10.1128/mSystems.00031-18) 
and [Vioscreen](https://www.vioscreen.com/Home) information to identify robust links between diet and life style 
features and broad abundance patterns of microbial families. 

The manually curated data file as 'phyloseq' object is available [here](https://github.com/amishra-stats/nbfar/raw/master/manuscript_file/final_ag_vio.rda). 

Some of our findings are summarized in the Figure below where we show microbial families and host-associated feature 
biclusters automatically identified by `nbfar`.

<img src="https://i.imgur.com/XK37aHm.jpg" alt="nbfar" width=100%/>



Getting started  
--------------

The `nbfar` package is currently available on GitHub and can be installed as follows.
The package `Rcpp` is required for installation.

```
# Install packages
install.packages("Rcpp", repos="https://rcppcore.github.io/drat")
devtools::install_github('amishra-stats/nbfar/nbfar', force = TRUE)
```

To use the library, we also rely on parallel computation. 

```
# load library
library(nbfar)
library(RcppParallel)
```


Simulation examples 
--------------

We showcase the usage of **nbfar** and  **nbrrr** on simulated data. 


### Data simulation


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

###  Negative Binomial reduced rank regression: **nbrrr** 


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
###  Negative Binomial co-sparse factor regression: **nbfar** 

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


Community Guidelines
--------------------

1.  Contributions and suggestions to the software are always welcome.
    Please consult our [contribution guidelines](https://github.com/mingzehuang/latentcor/blob/master/CONTRIBUTING.md) prior
    to submitting a pull request.
2.  Report issues or problems with the software using github’s [issue
    tracker](https://github.com/amishra-stats/nbfar/issues).
3.  Contributors must adhere to the [Code of Conduct](https://github.com/amishra-stats/latentcor/blob/master/CODE_OF_CONDUCT.md).

Acknowledgments
--------------

We thank Dr. Andreas Buja for useful comments on the project.

## Inquiries

You can also contact us via email

- [Aditya Mishra](mailto:amishra@flatironinstitute.org)
- [Christian L. Mueller](mailto:cmueller@flatironinstitute.org)
