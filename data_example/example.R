###########################################################
## R code for simulation examples and applications       ##
## in "Negative binomial co-sparse factor regression"    ##
###########################################################

## initial value NBRRR: dispesionn parameter: moment estimaate
## try the idea of warm start with condition:
## missing data case, and simulation speciefied mean and variance structure
# https://cran.r-project.org/web/packages/mpath/mpath.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4201648/
# https://cran.r-project.org/web/packages/mpath/vignettes/static_german.pdf
# https://arxiv.org/pdf/1712.03412.pdf
# https://www.fs.fed.us/psw/publications/ritchie/psw_2014_ritchie001_crotteau.pdf
# LD setting: epsilon = 1e-7, objI = 0
# control_nbfar <- nbfar_control(gamma0 = 1, spU = 40/p, spV = 20/q,
# lamMinFac = 1e-7, epsilon = 1e-7, objI = 0)
# # c
#
rm(list = ls())
## Install and load the package "gofar" from GitHub
# devtools::install_github("amishra-stats/nbfar/nbfar")
library(nbfar)

## Load required packages
library(MASS)
require(ggplot2)
require(reshape2)
library(magrittr)
library(glmnet)
library(tictoc)
library(RcppParallel)
# Rcpp::sourceCpp('src/nbfar.cpp')
# source('R/nbfar.R')


## --------------------------------------------------------
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
sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = 1, depth = 10)  # Simulated sample
X <- sim.sample$X[1:n, ]                    # simulated predictors (training)
Y <- sim.sample$Y[1:n, ]                    # simulated responses (training)
sum(Y==0)
# 1000 test sample data
sim.sample <- nbfar_sim(U, D, V, 1000, Xsigma, C0, disp = 1, depth = 10)
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
control_nbrr <- nbfar_control(initmaxit = 300, initepsilon = 1e-4)
# Model with known offset
nbrrr_test_of <- nbrrr(Y, X, maxrank = 5, cIndex = NULL, ofset = offset,
                       control = control_nbrr, nfold = 5)

plot(colMeans(nbrrr_test_of$cv.err))

# # Model without offset
# set.seed(1234)
# nbrrr_test <- nbrrr(Y, X, maxrank = 5, control = control_nbrr, nfold = 5)
# save(list = ls(), file = '../model_check.rda')
# plot(colMeans(nbrrr_test$cv.err))
# plot(log(nbrrr_test$diffobj[nbrrr_test$diffobj > 0]))




## ----- Negative binomial co-sparse factor regression
# Model with offset
# defaultNumThreads() / 2
set.seed(1234)
offset <- log(10)*matrix(1,n,q)
control_nbfar <- nbfar_control(gamma0 = 1, spU = 100/p, spV = 20/q, maxit = 2000,
                               lamMinFac = 1e-7, epsilon = 1e-7, objI = 0,
                               initmaxit = 3000, initepsilon = 1e-8)
tic()
nbfar_test_of1 <- nbfar(Y, X, maxrank = 1, nlambda = 20, cIndex = NULL,
                       ofset = offset, control = control_nbfar, nfold = 5,
                       PATH = F, nthread = 1)
toc()
plot(colMeans(nbfar_test_of1$fit[[1]]$dev))
nbfar_test_of1$fit[[1]]$ExecTimekpath



set.seed(1234)
RcppParallel::setThreadOptions(numThreads = 3)
offset <- log(10)*matrix(1,n,q)
control_nbfar <- nbfar_control(gamma0 = 1, spU = 100/p, spV = 20/q, maxit = 2000,
                               lamMinFac = 1e-7, epsilon = 1e-7, objI = 0,
                               initmaxit = 3000, initepsilon = 1e-8)
tic()
nbfar_test_of <- nbfar(Y, X, maxrank = 1, nlambda = 20, cIndex = NULL,
                       ofset = offset, control = control_nbfar, nfold = 5,
                       PATH = F, nthread = 3)
toc()
plot(colMeans(nbfar_test_of$fit[[1]]$dev))
nbfar_test_of$fit[[1]]$ExecTimekpath


nbfar_test_of$fit[[1]]$dev - nbfar_test_of1$fit[[1]]$dev

#
# set.seed(1234)
# nbfar_test <- nbfar(Y, X, maxrank = 5, nlambda = 20, cIndex = NULL,
#                     ofset = NULL, control = control_nbfar, nfold = 5,
#                     PATH = FALSE)
# plot(colMeans(nbfar_test$fit[[1]]$dev))
# save(list = ls(), file = '../model_check.rda')




## ----- Marginal model ------ non-sparse
# nbCol <- function(Y, X0, ofset, naind) {
#   q  <- ncol(Y)
#   PHI <- rep(0, q)
#   Cini <- matrix(0, nrow = ncol(X0), ncol = q)
#   for (i in 1:q) {
#     print(i)
#     qqq <- naind[,i] == 1
#     ft <- MASS::glm.nb(Y[qqq, i]~X0[qqq, ] + offset(ofset[qqq, i]) -1,
#                        control = stats::glm.control(maxit = 10000, epsilon = 1e-5),
#                        init.theta = MASS::glm.nb(Y[qqq, i]~1+
#                                                    offset(ofset[qqq,i]))$theta)
#     tem = coef(ft); tem[is.na(tem)] <- 0
#     Cini[, i] <- ft$coefficients
#     PHI[i] <- ft$theta
#   }
#   return(list(C = Cini, PHI = PHI))
# }


## ----- Marginal model ------ sparse ---






















#
#
#
# # Model fitting begins:
# control$epsilon <- 1e-7                       # error tolerance
# control$spU <- 20/p;                          # maximum support of U
# control$spV <- 20/q;                          # maximum support of V
# if (p > 50) control$spU <- 0.5
# control$maxit <- 2000                         # maximum number of itreration
# control$objI <- 1                             # convergence on the basis of objective function
# control$gamma0 <- 1                           # power factor for constructing weights
# control$se1 <- 1                              # 1 satndard deviation rule for cross validation
#
#
# ## Initalize NBRRR model for only inttercept in the model
# ## count distribution model
# ## Compute null deviance of tthe model
#
#
# # Model fitting
# # GOFAR(S) (full data)
#
# set.seed(SD)
# tic()
# fit.seq <- gofar_s(Y, X, nrank = rank.est, family = family,
#                    nlambda = nlam, familygroup = familygroup,
#                    control = control, nfold = 5, PATH = TRUE)
# fit.seq$time = toc()
#
#
# # GOFAR(S) (missing data)
# set.seed(SD)
# tic()
# fit.seq.m <- gofar_s(Ym, X, nrank = rank.est, family = family,
#                      nlambda = nlam, familygroup = familygroup,
#                      control = control, nfold = 5)
# fit.seq.m$time = toc()
#






# rm(list = ls())
# task_id <- Sys.getenv('DISBATCH_REPEAT_INDEX')
# ij <- as.numeric(task_id)
# print(Sys.getenv())
# cores <- as.numeric(Sys.getenv('CORES'))
# cores <- min(cores, 25)
#
# #ij=4; cores=4   # 4:43
# caseid <- ij%/%4 # 4:43%/%4
# jid <- (ij%%4) #(4:43%%4)+1
# caseSindex <- (25*jid)+(1:25)
# load('cond_siml.rda');scase <- data.frame(scase)
# inpar <- parse(text=paste0(paste(names(scase),scase[caseid,],sep = '='),collapse=';'))
#
#
# data.case <- vector('list',cores)
# for(i in 1:length(data.case)) {
#   data.case[[i]]$inpar <- inpar
#   data.case[[i]]$fid <- caseSindex[i]
#   data.case[[i]]$caseid <- caseid
#   data.case[[i]]$SD <- (2*caseid+1)*(2*caseSindex[i]+1)
#   data.case[[i]]$mmit <- mmiitt[caseid]
# }
#
#
# getKappaC0 = function(X0,familygroup, alpx){
#   family <- list(gaussian(),binomial(),poisson())
#   temp <- rep(0,length(family))
#   svdX0d1 <- svd(X0)$d[1]
#   cfamily <- unique(familygroup)
#   for (j in cfamily)
#     temp[j] <- switch(family[[j]]$family,'gaussian' = svdX0d1,
#                       'binomial' = svdX0d1/2,'poisson' = svdX0d1*alpx)
#   1.0*max(temp)
# }
#
#
# library(snow)
# library(parallel)
# cl <- makeSOCKcluster(c(rep('localhost', cores)))
# clusterExport(cl,list('getKappaC0'))
# caseT <- parLapply(cl,data.case, function(mylist){ # mylist <- data.case[[1]]
#   tryCatch({
#     ## loading package and data
#     ## # mylist <- data.case[[2]]
#     library(glmnet)
#     library(MASS)
#     library(Rcpp)
#     library(RcppArmadillo)
#     library(gofar)
#     library(rrpack)
#     library(tictoc)
#     # source('~/Desktop/googledrive/gsecure/simulation/cal500/miscFun.R')
#     # source('~/Desktop/googledrive/gsecure/simulation/cal500/grrr-all-fun.R')
#     source('miscFun.R')
#     source('grrr-all-fun.R')
#     dt <- mylist
#     for (i in 1:length(mylist)) {
#       tempobj = mylist[[i]]
#       eval(parse(text = paste(names(mylist)[[i]],"= tempobj")))
#     }
#     eval(inpar)
#     ## Other model fitting parameters:
#     nrank <- 3;rank.est <- 5
#     nlam <- 40
#     snr <- 0.5
#     set.seed(SD)
#     # generate dataname and file name;
#     folname <- paste0('out/codes',caseid)
#     ffname <- paste("data",fid,n,p,q1,q2,q3,"mat.RData",sep = '_')
#     ffname2 <- paste("all",ffname,sep = '_')
#     ffname <- file.path(folname,ffname)
#     ffname2 <- file.path(folname,ffname2)
#     #
#     q <- q1+q2+q3
#     respFamily <-c("gaussian","binomial","poisson")
#     family <- list(gaussian(),binomial(),poisson())
#     familygroup <- c(rep(1,q1),rep(2,q2),rep(3,q3))
#     cfamily <- unique(familygroup)
#     nfamily <- length(cfamily)
#     #
#     control <- gofar_control()
#     #
#     # source('~/Desktop/googledrive/gsecure/simulation/testsiml/dataGen.R', echo=TRUE)
#     # source('dataGen.R')
#     ## --------------------------------------------------------
#     ## Generate data
#     D <- rep(0,nrank)
#     V <- matrix(0,ncol=nrank ,nrow=q)
#     U <- matrix(0,ncol=nrank ,nrow=p)
#     #
#     U[,1]<-c(sample(c(1,-1),8,replace=TRUE),rep(0,p-8))
#     U[,2]<-c(rep(0,5),sample(c(1,-1),9,replace=TRUE),rep(0,p-14))
#     U[,3]<-c(rep(0,11),sample(c(1,-1),9,replace=TRUE),rep(0,p-20))
#     #
#     if(nfamily==1){
#       # for similar type response type setting
#       V[,1]<-c(rep(0,8),sample(c(1,-1),8,replace=TRUE)*runif(8,0.3,1),rep(0,q-16))
#       V[,2]<-c(rep(0,20),sample(c(1,-1),8,replace=TRUE)*runif(8,0.3,1),rep(0,q-28))
#       V[,3]<-c(sample(c(1,-1),5,replace=TRUE)*runif(5,0.3,1),rep(0,23),
#                sample(c(1,-1),2,replace=TRUE)*runif(2,0.3,1),rep(0,q-30))
#     } else {
#       # for mixed type response setting
#       # V is generated such that joint learning can be emphasised
#       V1 <- matrix(0,ncol=nrank ,nrow=q/2);
#       V1[,1]<-c(sample(c(1,-1),5,replace=TRUE),rep(0,q/2-5))
#       V1[,2]<-c(rep(0,3),V1[4,1],-1*V1[5,1],sample(c(1,-1),3,replace=TRUE),rep(0,q/2-8))
#       V1[,3]<-c(V1[1,1],-1*V1[2,1],rep(0,4),V1[7,2],-1*V1[8,2],sample(c(1,-1),2,replace=TRUE),rep(0,q/2-10))
#       #
#       V2 <- matrix(0,ncol=nrank ,nrow=q/2);
#       V2[,1]<-c(sample(c(1,-1),5,replace=TRUE),rep(0,q/2-5))
#       V2[,2]<-c(rep(0,3),V2[4,1],-1*V2[5,1],sample(c(1,-1),3,replace=TRUE),rep(0,q/2-8))
#       V2[,3]<-c(V2[1,1],-1*V2[2,1],rep(0,4),V2[7,2],-1*V2[8,2],sample(c(1,-1),2,replace=TRUE),rep(0,q/2-10))
#       #
#       V <- rbind(V1,V2)
#     }
#     U[,1:3]<- apply(U[,1:3],2,function(x)x/sqrt(sum(x^2)))
#     V[,1:3]<- apply(V[,1:3],2,function(x)x/sqrt(sum(x^2)))
#     #
#     D <-   s*c(4,6,5)     # signal strength varries as per the value of s
#     or <- order(D,decreasing = T)
#     U <- U[,or];V<-V[,or];D <- D[or]
#     C <- U %*% (D*t(V))     # simulated coefficient matrix
#     intercept <- rep(0.5,q)   # specifying intercept to the model:
#     C0 <- rbind(intercept,C)
#     #
#     Xsigma <- 0.5^abs(outer(1:p, 1:p,FUN = "-"))
#     sim.sample <- gofar_sim(U,D,V,n,Xsigma,C0,familygroup,snr) ## data simulated
#     pHI <-  c(rep(sim.sample$sigmaG,q1),rep(1,q2),rep(1,q3))     ##  similated stantard deviation
#     X <- sim.sample$X[1:n,]
#     Y <- sim.sample$Y[1:n,]
#     sim.sample <- gofar_sim(U,D,V,1000,Xsigma,C0,familygroup,snr)
#     Xt <- sim.sample$X
#     Yt <- sim.sample$Y
#     #
#     crossprod(X %*% U)
#     apply(X,2,norm,'2')
#     X0 <- cbind(1,X)
#     #----------------------------------
#     #
#     out <- vector('list',2)
#     if ( nfamily == 1) {
#       mname <- c("GSeCURE","GSeCURE.cv","GSeCURE.eea","mRRR.r","mRRR.rs",
#                  "mRRR.n","uGLM")    ## Mpdel Name
#     } else  {
#       mname <- c("GSeCURE","GSeCURE.cv","GSeCURE.eea","GSeCURE.s","GSeCURE.eea.s",
#                  "mRRR.r","mRRR.rs","mRRR.n","uGLM")    ## Mpdel Name
#     }
#     pName <- c('Dev','ErC','ErY','ErYg','ErYng','ErPhi','FPR',
#                'FNR','R%','r','FDR')      # Model accuracy measumrement parameter
#     #
#     ## M% = 20%
#     # Generating required proportion of missing data
#     miss <- 0.20
#     t.ind <- sample.int(n*q, size = miss*n*q)
#     y <- as.vector(Y); y[t.ind] <- NA;  Ym <- matrix(y,n,q)
#     naind <- (!is.na(Ym)) + 0 #matrix(1,n,q)
#     misind <- any(naind == 0) + 0
#
#     # Model fitting begins:
#     control$epsilon <- 1e-7
#     control$spU <- 20/p; control$spV <- 20/q;
#     if (p > 50) control$spU <- 0.5
#     control$maxit <- mmit
#     control$sv.tol <- 1e-2
#     control$objI <- 1
#     control$gamma0 <- 1
#     control$se1 <- 1
#     control$alp <- 10
#     save(list = ls(),file = ffname)
#
#     # GSeCURE
#     tic()
#     set.seed(SD)
#     control$initepsilon = 1e-6;control$initmaxit <- 2000
#     fit.seq  <- gofar_s(Y, X, nrank = rank.est, family = family, nlambda = nlam,
#                         familygroup=familygroup, #cIndex = cidd2,
#                         control = control,nfold=5)
#     fit.seq$time = toc()
#
#
#     # GRRR
#     tic()
#     set.seed(SD)
#     kappaC0 <- getKappaC0(cbind(1,X),familygroup, control$alp)
#     grrtune_hard <- grrr.tune(Y, X, ctrl.id = c(),
#                               family = family,
#                               familygroup = familygroup,
#                               maxrank = 10,
#                               penstr = list(penaltySVD="count",penaltyS="soft",
#                                             lambdaSVD=c(1:10), lambdaS=500),
#                               control = list(epsilon= 1e-6,sv.tol=1e-2,
#                                              maxit=300,trace=FALSE),
#                               init = list(kappaC0=kappaC0), gammaC0 = 1.05,
#                               conv.obj = FALSE, is.pca = NULL,
#                               nfold = 5, foldid=NULL, nlam=nlam,
#                               toplot = FALSE, warm=FALSE)$fit
#     grrtune_hard$Z <- matrix(grrtune_hard$coef[1,],ncol=ncol(Y))
#     grrtune_hard$C <- grrtune_hard$coef[-1,]
#     grrtune_hard$familygroup <- familygroup
#     grrtune_hard$Phi <- grrtune_hard$dispersion
#     nrank <- grrtune_hard$nrank
#     XXtC <- X %*% grrtune_hard$C
#     svdxc <- svd(XXtC)
#     grrtune_hard$V <- as.matrix(svdxc$v[,1:nrank])
#     grrtune_hard$D <- svdxc$d[1:nrank]/sqrt(n)
#     grrtune_hard$U <- grrtune_hard$C%*%grrtune_hard$V%*%diag(1/grrtune_hard$D,nrow=nrank,ncol=nrank)
#     grrtune_hard$time = toc()
#
#
#     # GEEA <-  use estimates from GRRR as initializer for EEA algorithm
#     set.seed(SD)
#     tic()
#     if(all(familygroup==2)){
#       control$initepsilon = 1e-4;control$initmaxit <- 300
#     }
#     fit.eea <- gofar_p(Y, X, nrank = rank.est, nlambda = nlam,
#                        family = family,
#                        familygroup = familygroup, cIndex = NULL,
#                        control = control, nfold = 5)
#     fit.eea$time = toc()
#
#
#     # GSeCURE-sep
#     # Separate model fitting with GSeCURE, when outcome variablesare of mixed types
#
#     if (nfamily > 1 ) {
#       tic()
#       fitT <- vector('list',length(cfamily))
#       C.grrr <- matrix(0,nrow = p + 1, ncol = q)
#       tPhi <- NULL
#       for (i in cfamily) {
#         qq <- familygroup == i
#         set.seed(SD)
#         fitT[[i]]  <- gofar_s(Y[,qq], X, nrank = rank.est, family = family, nlambda = nlam,
#                               familygroup=familygroup[qq], #cIndex = cidd2,
#                               control = control,nfold=5)
#         C.grrr[,qq] <- rbind(fitT[[i]]$Z,fitT[[i]]$C)
#         tPhi <- c(tPhi,fitT[[i]]$Phi)
#       }
#
#       fit.sGRRR <- getUVfromCoef(X,C.grrr)
#       fit.sGRRR$familygroup <- familygroup
#       fit.sGRRR$Phi <- tPhi
#       fit.sGRRR$time = toc()
#     }
#
#
#     # GLMNET
#     tic()
#     CK <- matrix(0,p+1,q)
#     phhh <- rep(1,q)
#     for(i in 1:q){
#       print(i)
#       qq <- !is.na(Y[,i])
#       cv <- cv.glmnet(X[qq,],Y[qq,i],family=family[[familygroup[i]]]$family,
#                       standardize = F,intercept=T, nfold=5,#nlambda = 10,
#                       alpha=0.9,maxit=10000)
#       CK[,i] <- as.vector(coef(cv$glmnet.fit)[,which(cv$lambda == cv$lambda.min)])
#       if(familygroup[i]==1)
#         phhh[i] <- mean((Y[qq,i]-cbind(1,X[qq,])%*%CK[,i])^2)
#     }
#     fit.glmnet <- list()
#     fit.glmnet$Z <- matrix(CK[1,],ncol=q)
#     fit.glmnet$C <- CK[-1,]
#     fit.glmnet$familygroup <- familygroup
#     fit.glmnet$Phi <- phhh  # grrfit_ore2$dispersion # ****************
#     # nrank <- rank.est
#     XXtC <- X %*% fit.glmnet$C
#     svdxc <- svd(XXtC)
#     nrank <- sum(svdxc$d > 1e-07)
#     fit.glmnet$V <- as.matrix(svdxc$v[,1:nrank])
#     fit.glmnet$D <- svdxc$d[1:nrank]/sqrt(n)
#     fit.glmnet$U <- fit.glmnet$C%*%fit.glmnet$V%*%diag(1/fit.glmnet$D,nrow=nrank,ncol=nrank)
#     fit.glmnet$time = toc()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#     ## for missing data case implementation check
#     tic()
#     # GSeCURE
#     set.seed(SD)
#     control$initepsilon = 1e-6;control$initmaxit <- 2000
#     fit.seq3  <- gofar_s(Ym, X, nrank = rank.est, family = family, nlambda = nlam,
#                          familygroup=familygroup, #cIndex = cidd2,
#                          control = control,nfold=5)
#     fit.seq3$time = toc()
#
#     # GRRR
#     set.seed(SD)
#     tic()
#     kappaC0 <- getKappaC0(cbind(1,X),familygroup, control$alp)
#     grrtune_hard3 <- grrr.tune(Ym, X, ctrl.id = c(),
#                                family = family,
#                                familygroup = familygroup,
#                                maxrank = 10,
#                                penstr = list(penaltySVD="count",penaltyS="soft",
#                                              lambdaSVD=c(1:10), lambdaS=50),
#                                control = list(epsilon= 1e-6,sv.tol=1e-2,
#                                               maxit=300,trace=FALSE),
#                                init = list(kappaC0=kappaC0), gammaC0 = 1.05,
#                                conv.obj = FALSE, is.pca = NULL,
#                                nfold = 5, foldid=NULL, nlam=nlam,
#                                toplot = FALSE, warm=FALSE)$fit
#     grrtune_hard3$Z <- matrix(grrtune_hard3$coef[1,],ncol=ncol(Ym))
#     grrtune_hard3$C <- grrtune_hard3$coef[-1,]
#     grrtune_hard3$familygroup <- familygroup
#     grrtune_hard3$Phi <- grrtune_hard3$dispersion
#     nrank <- grrtune_hard3$nrank
#     XXtC <- X %*% grrtune_hard3$C
#     svdxc <- svd(XXtC)
#     grrtune_hard3$V <- as.matrix(svdxc$v[,1:nrank])
#     grrtune_hard3$D <- svdxc$d[1:nrank]/sqrt(n)
#     grrtune_hard3$U <- grrtune_hard3$C%*%grrtune_hard3$V%*%diag(1/grrtune_hard3$D,nrow=nrank,ncol=nrank)
#     grrtune_hard3$time = toc()
#
#     # GEEA <-  use estimates from GRRR as initializer for EEA algorithm
#     set.seed(SD)
#
#     if(all(familygroup==2)){
#       control$initepsilon = 1e-4;control$initmaxit <- 300
#     }
#     tic()
#     fit.eea3 <- gofar_p(Ym, X, nrank = rank.est, nlambda = nlam,
#                         family = family,
#                         familygroup = familygroup, cIndex = NULL,
#                         control = control, nfold = 5)
#     fit.eea3$time = toc()
#
#
#
#     # GSeCURE-sep
#     # Separate model fitting with GSeCURE, when outcome variablesare of mixed types
#
#     if (nfamily > 1 ) {
#       tic()
#       fitT3 <- vector('list',length(cfamily))
#       C.grrr3 <- matrix(0,nrow = p + 1, ncol = q)
#       tPhi3 <- NULL
#       for (i in cfamily) {
#         qq <- familygroup == i
#         set.seed(SD)
#         fitT3[[i]]  <- gofar_s(Ym[,qq], X, nrank = rank.est, family = family, nlambda = nlam,
#                                familygroup=familygroup[qq], #cIndex = cidd2,
#                                control = control,nfold=5)
#         C.grrr3[,qq] <- rbind(fitT3[[i]]$Z,fitT3[[i]]$C)
#         tPhi3 <- c(tPhi3,fitT3[[i]]$Phi)
#       }
#
#       fit.sGRRR3 <- getUVfromCoef(X,C.grrr3)
#       fit.sGRRR3$familygroup <- familygroup
#       fit.sGRRR3$Phi <- tPhi3
#       fit.sGRRR3$time = toc()
#     }
#
#
#
#     # GLMNET
#     tic()
#     CK3 <- matrix(0,p+1,q)
#     phhh3 <- rep(1,q)
#     for(i in 1:q){
#       print(i)
#       qq <- !is.na(Ym[,i])
#       cv <- cv.glmnet(X[qq,],Ym[qq,i],family=family[[familygroup[i]]]$family,
#                       standardize = F,intercept=T, nfold=5,#nlambda = 10,
#                       alpha=0.9,maxit=10000)
#       CK3[,i] <- as.vector(coef(cv$glmnet.fit)[,which(cv$lambda == cv$lambda.min)])
#       if(familygroup[i]==1)
#         phhh3[i] <- mean((Ym[qq,i]-cbind(1,X[qq,])%*%CK3[,i])^2)
#     }
#     fit.glmnet3 <- list()
#     fit.glmnet3$Z <- matrix(CK3[1,],ncol=q)
#     fit.glmnet3$C <- CK3[-1,]
#     fit.glmnet3$familygroup <- familygroup
#     fit.glmnet3$Phi <- phhh3  # grrfit_ore2$dispersion # ****************
#     # nrank <- rank.est
#     XXtC <- X %*% fit.glmnet3$C
#     svdxc <- svd(XXtC)
#     nrank <- sum(svdxc$d > 1e-07)
#     fit.glmnet3$V <- as.matrix(svdxc$v[,1:nrank])
#     fit.glmnet3$D <- svdxc$d[1:nrank]/sqrt(n)
#     fit.glmnet3$U <- fit.glmnet3$C%*%fit.glmnet3$V%*%diag(1/fit.glmnet3$D,nrow=nrank,ncol=nrank)
#     fit.glmnet3$time = toc()
#
#     ## store output locally
#     save(list = ls(),file = ffname2)
#
#     # dt$model <- modT
#     return(1)
#   }, error=function(error_message) {
#     # list(x,message(error_message))
#     return(NA)
#   }
#   )
# })
# stopCluster(cl)



# Here we describe the setting for model simulation setting;
# Simulation task for the project -
# n < p [O% = 5, 10, 20] [n = 200, p =100] [high and moderate]
# n > p [O% = 5, 10, 20] [n = 200, p =300] [high and moderate]
# nsim = 100
# write simulation benchmark with these settings
#DISBATCH REPEAT 100 start 1 module purge ; module load slurm gcc R ; OMP_NUM_THREADS=1 CORES=$(nproc --all) Rscript siml-sp.R > logfile/${DISBATCH_REPEAT_INDEX}.log 2>&1
# module purge ; module load slurm gcc R ; OMP_NUM_THREADS=1 Rscript siml_biostat.R 10 250.0 598.742 0.26703 0.1 49 100 4 > logfile/250.log 2>&1
#
# # Generate simulation stetting code test file and check
# n <- 200
# p <- c(100, 300)
# shFac <- c(6, 8)
# L <- c(0, 1)
# O <- c(0.05, 0.10, 0.20)
# nsim <- 100
# temp <- 'module purge ; module load slurm gcc ; OMP_NUM_THREADS=1 Rscript siml_biostat.R'
# j = 1
# fwname <- "siml_bio"
# zz <- file(fwname, "w")
# for (tp in p) {
#   for (i in 1:nsim) {
#     tem <- paste(temp, n, tp, 0, 6, 0, i, j, paste('> logfile/',j,'.log 2>&1',sep = ''))
#     j <- j + 1
#     # print(tem)
#     writeLines(tem,zz)
#   }
# }
# for (tp in p) {
#   for (tl in L) {
#     for (s in shFac) {
#       for (o in O) {
#         for (i in 1:nsim) {
#           tem <- paste(temp, n, tp, tl, s, o, i, j,  paste('> logfile/',j,'.log 2>&1',sep = ''))
#           j <- j + 1
#           # print(tem)
#           writeLines(tem,zz)
#         }
#       }
#     }
#   }
# }
# close(zz)
# module purge ; module load slurm gcc R ; OMP_NUM_THREADS=1 Rscript siml_biostat.R 200 100 0 6 0 1 1 > logfile/1.log 2>&1
# # sbatch -N 10 -p ccm disBatch.py -t 50 siml_bio

args <- commandArgs(trailingOnly = TRUE)
print(args)
# n, tp, tl, s, o, i, j
n <- as.numeric(args[1])
p <- as.numeric(args[2])
L <- as.numeric(args[3])
shFac <- as.numeric(args[4])
O <- n*as.numeric(args[5])
index <- as.numeric(args[6])
example_seed <- as.numeric(args[7])
print(c(n,p,L,shFac,O,index,example_seed))
temp <- paste('siml',n,p,L,shFac,O,index,example_seed,sep = '_')
fname <- paste0('out/',temp,'.rda')  # out file name
set.seed(example_seed)
# install.packages("../robregcc_1.0.tar.gz")
# rm(list = ls())
# n = 200; p = 300; L = 1; shFac = 6; O = 0.20*n
# index = 12; example_seed = 123



ngrp <- 4                           # number of sub-composition
snr <- 3                            # Signal to noise ratio
## -----------------------------------------------------------------------------
## Load package
library(robregcc)
library(magrittr)

#--------------------------------------------------------------
# ## Define parameters to simulate example
# p <- 300                            # number of predictors
# n <- 200                            # number of sample
# O <- 0.10*n                         # number of outlier, e.g. 15% of observation
# L <- 1                              # indicator variable for outlier type,
# # L = {0,1} => leveraged {no, yes}
#
# # generate outlier by shifting "O"observation by amount equals to shFac times
# # true error variance sigma.
# # shFac = {6,8} corresponds to {moderate, high} outlier
# shFac <- 6
# example_seed <- 2*p+1               # example seed
# set.seed(example_seed)
# ngrp <- 4                           # number of sub-composition
# snr <- 3                            # Signal to noise ratio

#--------------------------------------------------------------
## Simulate true model variables, i.e., y, X, C, beta
## Follow examples from [Pixu Shi 2016]

# Simulate subcomposition matrix
C1 <- matrix(0,ngrp,23)
tind <- c(0,10,16,20,23)
for (ii in 1:ngrp)
  C1[ii,(tind[ii] + 1):tind[ii + 1]] <- 1
C <- matrix(0,ngrp,p)
C[,1:ncol(C1)] <- C1


# model parameter beta
beta0 <- 0.5
beta <- c(1, -0.8, 0.4, 0, 0, -0.6, 0, 0, 0, 0, -1.5, 0, 1.2, 0, 0, 0.3)
beta <- c(beta,rep(0,p - length(beta)))

# Simulate response and predictor, i.e., X, y
Sigma  <- 1:p %>% outer(.,.,'-') %>% abs(); Sigma  <- 0.5^Sigma
data.case <- vector("list",1)
data.case <- robregcc_sim(n,beta,beta0, O = O,Sigma,levg = L,
                          snr,shft = shFac,0,
                          C,out = data.case)


#--------------------------------------------------------------
## Data preprocessing:

X <- data.case$X                          # predictor matrix
y <- data.case$y                          # model response
#

# Predictor transformation due to compositional constraint:
# Equivalent to performing centered log-ratio transform
# Xt <- svd(t(C))$u %>% tcrossprod() %>% subtract(diag(p),.) %>%
#   crossprod(t(X),.)
# Xm <- colMeans(Xt)
# Xt <- scale(Xt,Xm,FALSE)                  # centering of predictors
# #
# mean.y <- mean(y)
# y <- y - mean.y                           # centering of response
#
# Account for intercept in the model
Xt <- cbind(1,X)                          # accounting for intercept in predictor
C <- cbind(0,C)                           # accounting for intercept in constraint
bw <- c(0,rep(1,p))                       # weight matrix to not penalize intercept

#--------------------------------------------------------------
# Prescreening
betaT <- beta
if (FALSE) {
  control <- robregcc_option(spb = 0.5, spy = 0.5)
  # Robust regression using lasso penalty [Huber equivalent]
  fit.soft <- robregcc_sp(Xt,y,C, cindex = 1,
                          control = control, penalty.index = 2,
                          alpha = 0.5)
  supp_index <- with(fit.soft, betapath[,ncol(betapath)]) != 0
  Xt = Xt[,supp_index]; C = C[,supp_index]; bw <- bw[supp_index]
  betaT <- c(beta0,beta)[supp_index]
  betaT <- betaT[-1]
  p <- length(betaT)
}

#--------------------------------------------------------------
## Initialization

# Breakdown point for tukey Bisquare loss function
#b1 = 0.5                    # 50% breakdown point
#cc1 =  1.567                # corresponding model parameter
b1 = 0.25; cc1 =  2.937      # initalization for scale parameter

set.seed(example_seed)      # unique seed
#--------------------------------------------------------------
# Non-robust regression, [Pixu Shi 2016]
control <- robregcc_option(maxiter = 5000, tol = 1e-8, lminfac = 1e-8)
fit.nr <- classo(Xt, y, C, we = bw, type = 1, control = control)
fit.nr$beta
if (O > 0) {
  fit.oracle <- classo(Xt[-(1:O),], data.case$yo[-(1:O)], C,
                       we = bw, type = 1, control = control)
  fit.oracle$beta
} else {
  fit.oracle <- fit.nr
}


# control parameter for intialization method
control <- robregcc_option(maxiter = 1000,tol = 1e-4,lminfac = 1e-7)

# intialization
fit.init <- cpsc_sp(Xt, y,alp = 0.4, cfac = 2, b1 = b1,
                    cc1 = cc1,C,bw,1,control)
# plot(fit.init$residualR)
# plot(fit.init$residuals)
# plot(with(fit.init, residuals/sd ))


# plot(fit.init$residuals/fit.init$sd)



#--------------------------------------------------------------
## Model fitting

# control parameters
control <- robregcc_option()
beta.wt <- fit.init$betaR          # Set weight for model parameter beta
beta.wt[1] <- 0
control$gamma = 1                   # gamma for constructing  weighted penalty
control$spb = 40/p                  # fraction of maximum non-zero model parameter beta
control$outMiter = 1000             # Outer loop iteration
control$inMiter = 3000              # Inner loop iteration
control$nlam = 50                   # Number of tuning parameter lambda to be explored
control$lmaxfac = 1                 # Parameter for constructing sequence of lambda
control$lminfac = 1e-8              # Parameter for constructing sequence of lambda
control$tol = 1e-20;                # tolrence parameter for converging [inner  loop]
control$out.tol = 1e-16             # tolerence parameter for convergence [outer loop]
control$kfold = 10                   # number of fold of crossvalidation
control$sigmafac = 2#1.345
# Robust regression using adaptive lasso penalty
fit.ada <- robregcc_sp(Xt,y,C,
                       beta.init = beta.wt,  cindex = 1,
                       gamma.init = fit.init$residuals,
                       control = control,
                       penalty.index = 1, alpha = 0.95)


# Robust regression using lasso penalty [Huber equivalent]
fit.soft <- robregcc_sp(Xt,y,C, cindex = 1,
                        control = control, penalty.index = 2,
                        alpha = 0.95)



# Robust regression using hard thresholding penalty
control$lmaxfac = 1e2               # Parameter for constructing sequence of lambda
control$lminfac = 1e-3              # Parameter for constructing sequence of lambda
control$sigmafac = 2#1.345
fit.hard <- robregcc_sp(Xt,y,C, beta.init = fit.init$betaf,
                        gamma.init = fit.init$residuals,
                        cindex = 1,
                        control = control, penalty.index = 3,
                        alpha = 0.95)


save(list = ls(), file = fname)


