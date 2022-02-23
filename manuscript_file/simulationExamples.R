###########################################################
## R code for simulation examples and applications       ##
##     Negative binomial factor regression model         ##
##    with application to microbiome data anaysis        ##
###########################################################


rm(list = ls())
library(devtools)
devtools::install_github("amishra-stats/nbfar/nbfar")
devtools::install_github("amishra-stats/gofar/gofar")


## Load required packages
library(MASS)
require(ggplot2)
library(magrittr)
library(mpath)
library(RcppParallel)
library(tictoc)
library(gofar)
library(nbfar)



## ## -----------------------Simulation settings --------------------
## Simulation setting:
## n: sample size
## p: number of predictors
## q: number of outcome variables
## nrank: true rank of the model
## snr: signal to noise ratio to be used in generating gaussian outcome variables
## nlam: number of lambda values to be fitted
## rank.est: maximum estimated rank to be specified to the model
## s: multiplicative factor of singular values
## nthread: number of thread to execute cross-validation parallel
## The simulation was replicated 100 times under each setting as detailed in the paper.


## ## -----------------------Simulation settings --------------------
# Setting specific to simulated examples:

# Setup I
p <- 50; n <- 200
# Setup II
# p <- 300; n <- 200


## The following code is for a single replication.
## ----------------------single replication ------------------------
## Other parameters
example_seed <- 123
set.seed(example_seed)
nrank <- 3                # true rank
rank.est <- 5             # estimated rank
nlam <- 40                # number of tuning parameter
snr <- 0.50               # SNR for variance Gaussian error
s  = 0.5; q <- 30
sp  = 0.5
nthread = 1
disp = 0.75



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
sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = disp, depth = 10)  # Simulated sample
X <- sim.sample$X[1:n, ]                    # simulated predictors (training)
Y <- sim.sample$Y[1:n, ]                    # simulated responses (training)
# 1000 test sample data
sim.sample <- nbfar_sim(U, D, V, 1000, Xsigma, C0, disp = disp, depth = 10)
Xt <- sim.sample$X                    # simulated predictors (test)
Yt <- sim.sample$Y                    # simulated predictors (test)
X0 <- cbind(1, X)                     # 1st column accounting for intercept


# Simulate data with 20% missing entries
miss <- 0.10          # Proportion of entries missing
t.ind <- sample.int(n * q, size = miss * n * q)
y <- as.vector(Y)
y[t.ind] <- NA
Ym <- matrix(y, n, q)   # 20% of entries are missing at random




# Model fitting: full data
# NBRRR
tic()
set.seed(example_seed)
control_r3 <- nbfar_control(initmaxit = 10000, initepsilon = 1e-5,
                            objI = 1)
offset_v <- matrix(0, nrow = n, ncol = q)
nbrrr_test <- nbrrr(Y, X, maxrank = rank.est, control = control_r3, nfold = 5, trace = F,
                    ofset = offset_v)
nbrrr_test$time = toc()



# NBFAR
tic()
RcppParallel::setThreadOptions(numThreads = nthread)
set.seed(example_seed)
offset_v <- matrix(0, nrow = n, ncol = q)
control_nbfar <- nbfar_control(gamma0 = 1, spU = sp, spV = 20/q,
                               maxit = 2000, lamMaxFac = 1e-2,
                               lamMinFac = 1e-7, epsilon = 1e-4,
                               objI = 0,
                               initmaxit = 10000, initepsilon = 1e-7)
nbfar_test <- nbfar(Y, X, maxrank = rank.est, nlambda = nlam,  ofset = offset_v,
                    cIndex = NULL, control = control_nbfar, nfold = 5,
                    PATH = FALSE, nthread = nthread,trace = F)
nbfar_test$time = toc()


# GOFAR: module for poisson outcomes
tic()
set.seed(example_seed)
control_poisson <- gofar_control(gamma0 = 1, spU = sp, spV = 20/q,
                                 maxit = 2000,
                                 lamMinFac = 1e-6, epsilon = 1e-6, objI = 1,
                                 initmaxit = 2000, initepsilon = 1e-6, alp = 10)
gofar_test <- gofar_s(Y, X, nrank = rank.est, nlambda = nlam,
                      familygroup = rep(3, q),
                      family = list(gaussian(), binomial(), poisson()),
                      control = control_poisson, nfold = 5)
gofar_test$time = toc()




# Model fitting: missing data
# NBRRR
tic()
set.seed(example_seed)
control_r3 <- nbfar_control(initmaxit = 10000, initepsilon = 1e-5,
                            objI = 1)
offset_v <- matrix(0, nrow = n, ncol = q)
nbrrr_testm <- nbrrr(Ym, X, maxrank = 5, control = control_r3, nfold = 5,
                     ofset = offset_v,
                     trace = F)
nbrrr_testm$time = toc()


# NBFAR:
tic()
RcppParallel::setThreadOptions(numThreads = nthread)
set.seed(example_seed)
control_nbfar <- nbfar_control(gamma0 = 1, spU = sp, spV = 20/q,
                               maxit = 2000, lamMaxFac = 1e-2,
                               lamMinFac = 1e-7, epsilon = 1e-4,
                               objI = 0,
                               initmaxit = 10000, initepsilon = 1e-7)
offset_v <- matrix(0, nrow = n, ncol = q)
nbfar_testm <- nbfar(Ym, X, maxrank = rank.est, nlambda = nlam,
                     cIndex = NULL, control = control_nbfar, nfold = 5,
                     ofset = offset_v,
                     PATH = FALSE, nthread = nthread,trace = F)
nbfar_testm$time = toc()




# GOFAR: module for poisson outcomes   [missinng data]
tic()
set.seed(example_seed)
control_poisson <- gofar_control(gamma0 = 1, spU = sp, spV = 20/q,
                                 maxit = 2000,
                                 lamMinFac = 1e-6, epsilon = 1e-6, objI = 1,
                                 initmaxit = 2000, initepsilon = 1e-6, alp = 10)
gofar_testm <- gofar_s(Ym, X, nrank = rank.est, nlambda = nlam,
                       familygroup = rep(3, q),
                       family = list(gaussian(), binomial(), poisson()),
                       control = control_poisson, nfold = 5)
gofar_testm$time = toc()





## -----------------------Model comparisons------------
## A function to compute the error measures for comparing different methods
## Comparison statistics are reported in Table 1
## Error in C, Error in XC, FPR, FNR,  R%, rank measure, Error in dispersion


model_compare = function(X, U, D, V, intercept, n, fit, pHI){
  p <- nrow(U); q <- nrow(V)
  intercept <- matrix(intercept, ncol = q)
  pz <- nrow(intercept)
  pt <- p + pz

  exectime <- with(fit$time, as.numeric(toc -tic))
  if(length(fit) == 0){
    return(NA)
  } else {
    names(fit)[names(fit) == "Phi"] <- 'PHI'
    ErPhi = sqrt(sum(pHI - fit$PHI)^2/q)
    # C <- rbind(intercept, U %*% (D * t(V)))
    # C.est <- rbind(fit$Z, fit$C)
    C <- U %*% (D * t(V))
    C.est <- fit$U %*% (fit$D * t(fit$V))
    nrank <- ncol(U)
    ErC <- 100 * (norm(C - C.est, "f"))^2/(pt * q)     # compute Er(C)

    MU0 <-  X %*% C.est
    Ym <- X %*% C
    ErY <- sqrt(sum((Ym - MU0)^2)/ (n * q))

    # test for rank and sparsity
    beta.sim <- c(as.vector(U), as.vector(V))
    truezero <- sum(beta.sim == 0)
    truenonzero <- sum(beta.sim != 0)
    fit$V <- as.matrix(fit$V)
    fit$U <- as.matrix(fit$U)
    rank.est <- length(fit$D)           # estimated rank
    if (rank.est == 1) {
      U.est <- matrix(fit$U)
      V.est <- matrix(fit$V)
      D.est <- fit$D
    } else {
      ord1 <- order(fit$D,decreasing = T)
      U.est <- fit$U[,ord1]
      V.est <- fit$V[,ord1]
      D.est <- fit$D[ord1]
    }
    rank.est <- NA
    if (is.null(U.est) | is.null(V.est)) {
      fp <- truezero
      fn <- 0
    } else {
      Ue <- matrix(0, p, nrank)
      Ve <- matrix(0, q, nrank)
      rank.est <- ncol(U.est)
      if (nrank < rank.est) {
        Ue <- fit$U[, nrank]
        Ve <- fit$V[, nrank]
      } else {
        Ue[, 1:rank.est] <- U.est
        Ve[, 1:rank.est] <- V.est
      }
      beta.est <- c(as.vector(Ue), as.vector(Ve))
      fp <- sum(beta.sim == 0 & beta.est != 0)
      fn <- sum(beta.sim != 0 & beta.est == 0)
    }
    if (is.na(rank.est)) {
      return(c( 100*ErC, 100*ErY, 100*fp/truezero,
                100*fn/truenonzero, 0, rank.est,ErPhi,exectime))
    }
    else {
      if (nrank >= rank.est) {
        return(c( 100*ErC, 100*ErY, 100*fp/truezero,
                  100*fn/truenonzero, 0, rank.est,ErPhi,exectime))
      }
      else {
        return(c( 100*ErC, 100*ErY,  100*fp/truezero, 100*fn/truenonzero,
                  100 * sum((D.est[-(1:nrank)])^2)/sum(D.est^2),
                  rank.est,ErPhi,exectime))
      }
    }
  }
}


## Lesser value of statistics suggest better model
pHI <- rep(disp,q)
compare.fit.full <- rbind(model_compare(X, U, D, V, intercept, n, nbrrr_test, pHI),
                          model_compare(X, U, D, V, intercept, n, nbfar_test, pHI),
                          model_compare(X, U, D, V, intercept, n, gofar_test, pHI))

compare.fit.miss <- rbind(model_compare(X, U, D, V, intercept, n, nbrrr_testm, pHI),
                          model_compare(X, U, D, V, intercept, n, nbfar_testm, pHI),
                          model_compare(X, U, D, V, intercept, n, gofar_testm, pHI))

colnames(compare.fit.full) <- colnames(compare.fit.miss)  <- c('Er(C)', 'Er(Theta)', 'FPR', 'FNR', 'R%', 'r', 'Er(Phi)','Time(s)')
rownames(compare.fit.full) <- rownames(compare.fit.miss) <- c('NB-RRR','NB-FAR', 'GO-FAR')
compare.fit.full
compare.fit.miss





## ----------------- Figure 2 of the manuscript ------------------
## We have processed the results from model comparison and save it on the
## GitHub page. One can download it and see the plots


# # Load data
data_url = "https://github.com/amishra-stats/nbfar/raw/master/manuscript_file/fig2.rda"
download.file(data_url, destfile='temp.rda')
load("temp.rda")


ggplot(dt, aes(x = Method, y = ErY, fill = M)) +
  facet_wrap( ~ case, scales = "free") +
  scale_y_continuous(trans='log10') +
  geom_boxplot(notch = TRUE,notchwidth = 0.1)+ xlab("") + theme_bw()+ scale_fill_grey() +
  ylab("") + theme(legend.direction = "horizontal",legend.position =  "top" ) +
  theme(strip.text = element_text(size=10, face="bold"),
    axis.text=element_text(size=8),
    legend.text=element_text(size=12),
    axis.text.x = element_text(angle = 45,
                               hjust = 1,size = 14,
                               face = "bold"),
    axis.text.y = element_text( hjust = 1,size = 14,
                                face = "bold"),
    panel.border = element_rect(colour = "black"),
    strip.text.x = element_text(size = 14,face = "bold"))









