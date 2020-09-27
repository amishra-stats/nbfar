
#' Default control parameters for Generalized co-sparse factor regresion
#'
#' Control parameters for NBFAR
#'
#' @param lamMaxFac a multiplier of calculated lambda_max
#' @param lamMinFac a multiplier of determing lambda_min as a fraction of lambda_max
#' @param gamma0 power parameter in the adaptive weights
#' @param elnetAlpha elastic net penalty parameter
#' @param spU maximum proportion of nonzero elements in each column of U
#' @param spV maximum proportion of nonzero elements in each column of V
#' @param maxit maximum iteration for each sequential steps
#' @param epsilon tolerence value set for convergene of gcure
#' @param objI 1 or 0 convergence on the basis of objective function or not
#' @param initmaxit maximum iteration for initialization problem
#' @param initepsilon tolerence value for convergene in the initialization problem
#' @return a list of controling parameter.
#' @export
#' @examples
#' control <- nbfar_control()
#' @useDynLib nbfar
nbfar_control <- function(maxit = 5000, epsilon = 1e-6,
                          elnetAlpha = 0.95,
                          gamma0 = 1,
                          spU = 0.5, spV = 0.5,
                          lamMaxFac = 1, lamMinFac = 1e-6,
                          initmaxit = 2000, initepsilon = 1e-6,
                          objI = 1) {
  list(
    lamMaxFac = lamMaxFac,
    lamMinFac = lamMinFac,
    gamma0 = gamma0, elnetAlpha = elnetAlpha, spU = spU, spV = spV,
    maxit = maxit, epsilon = epsilon,
    objI = objI,
    initmaxit = initmaxit, initepsilon = initepsilon
  )
}

# Fit glm.nb Columnwise on the control variable
#
# @param Y Multivariate response matrix
# @param X0 design matrix
# @param ofset offset  matrix
# @param naind indicator vector for missing data in Y
#' @useDynLib nbfar
#' @import magrittr
#' @importFrom MASS glm.nb
#' @importFrom stats glm.control rnorm offset coef
nbCol <- function(Y, X0, ofset, naind) {
  q  <- ncol(Y)
  PHI <- rep(0, q)
  Cini <- matrix(0, nrow = ncol(X0), ncol = q)
  for (i in 1:q) {
    print(i)
    qqq <- naind[,i] == 1
    ft <- MASS::glm.nb(Y[qqq, i]~X0[qqq, ] + offset(ofset[qqq, i]) -1,
                       control = stats::glm.control(maxit = 10000, epsilon = 1e-5),
                       init.theta = MASS::glm.nb(Y[qqq, i]~1+
                                                   offset(ofset[qqq,i]))$theta)
    tem = coef(ft); tem[is.na(tem)] <- 0
    Cini[, i] <- ft$coefficients
    PHI[i] <- ft$theta
  }
  return(list(C = Cini, PHI = PHI))
}




#' Simulate data
#'
#' Genertate random samples from a negative binomial sparse factor regression model
#'
#' @param U specified value of U
#' @param V specified value of V
#' @param D specified value of D
#' @param n sample size
#' @param Xsigma covariance matrix for generating sample of X
#' @param C0 Specified coefficient matrix with first row being intercept
#' @param disp dispersion parameter of the generative model
#' @param depth sequencing depth of the microbiome data
#' @return
#'   \item{Y}{Generated response matrix}
#'   \item{X}{Generated predictor matrix}
#' @export
#' @useDynLib nbfar
#' @import magrittr
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma rpois
#' @examples
#' ## Model specification:
#' SD <- 123
#' set.seed(SD)
#' p <- 100; n <- 200
#' pz <- 0
#' nrank <- 3                # true rank
#' rank.est <- 5             # estimated rank
#' nlam <- 20                # number of tuning parameter
#' s  = 0.5
#' q <- 30
#' control <- nbfar_control()  # control parameters
#' #
#' #
#' ## Generate data
#' D <- rep(0, nrank)
#' V <- matrix(0, ncol = nrank, nrow = q)
#' U <- matrix(0, ncol = nrank, nrow = p)
#' #
#' U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
#' U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
#' U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#' #
#'   # for similar type response type setting
#'   V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 16))
#'   V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 28))
#'   V[, 3] <- c(
#'     sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
#'     sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
#'   )
#' U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' #
#' D <- s * c(4, 6, 5) # signal strength varries as per the value of s
#' or <- order(D, decreasing = TRUE)
#' U <- U[, or]
#' V <- V[, or]
#' D <- D[or]
#' C <- U %*% (D * t(V)) # simulated coefficient matrix
#' intercept <- rep(0.5, q) # specifying intercept to the model:
#' C0 <- rbind(intercept, C)
#' #
#' Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
#' # Simulated data
#' sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = 3, depth = 10)  # Simulated sample
#' # Dispersion parameter
#' X <- sim.sample$X[1:n, ]
#' Y <- sim.sample$Y[1:n, ]
#'
nbfar_sim <- function(U, D, V, n, Xsigma, C0,disp,depth) {
  ## finding basis along more number of columns of data vector
  basis.vec <- function(x) {
    # require(Matrix)
    if (diff(dim(x)) < 0) x <- t(x)
    qd <- qr(x)
    k <- qr.Q(qd) %*% qr.R(qd)[, 1:qd$rank]
    k[abs(k) < 1e-6] <- 0
    b.ind <- vector()
    for (i in 1:qd$rank) {
      b.ind[i] <- which(apply(x, 2, function(x, y) sum(abs(x - y)), k[, i]) < 1e-6)[1]
    }
    return(list(ind = b.ind, vec = x[, b.ind]))
  }

  p <- nrow(U)
  q <- nrow(V)
  nrank <- ncol(U)

  U.t <- diag(max(dim(U)))
  U.t <- U.t[, -basis.vec(U)$ind]
  P <- cbind(U, U.t)
  UtXsUt <- t(U.t) %*% Xsigma %*% U.t
  UtXsU <- t(U.t) %*% Xsigma %*% U
  UXsU <- t(U) %*% Xsigma %*% U
  UXsUinv <- solve(UXsU)
  ## sigma.X2 <- t(U.t)%*%Xsigma%*%U.t - t(U.t)%*%Xsigma%*%U%*%solve(t(U)%*%Xsigma%*%U)%*%t(U)%*%Xsigma%*%U.t
  sigma.X2 <- UtXsUt - UtXsU %*% UXsUinv %*% t(UtXsU)
  sigma.X2 <- (sigma.X2 + t(sigma.X2)) / 2


  X1 <- matrix(nrow = nrank, ncol = n, rnorm(n * nrank))
  ## X1 <- t(mvrnorm(n,rep(0,ncol(U)),diag(ncol(U)) ))
  mean.X2 <- UtXsU %*% UXsUinv %*% X1
  ## mean.X2 <- t(U.t)%*%Xsigma%*%U%*%solve(t(U)%*%Xsigma%*%U)%*%X1
  X2 <- mean.X2 + t(MASS::mvrnorm(ncol(mean.X2), rep(0, nrow(mean.X2)), sigma.X2))
  X <- t(solve(t(P)) %*% rbind(X1, X2)) # /sqrt(n)
  # crossprod(X%*%U)

  X0 <- cbind(1, X)
  MU <- X0 %*% C0
  ## Generation  of Y and ty.lam to find lambda max in the selection process
  Y <- matrix(nrow = n, ncol = q, 0)
  # sigma <- rho
  lam <- depth*exp(MU)
  for(i in 1:nrow(lam)) {
    for(j in 1:ncol(lam)) {
      Y[i,j] <- stats::rgamma(1,shape= disp , scale= lam[i,j]/disp)
      Y[i,j] <- stats::rpois(1,lambda=Y[i,j])
    }
  }
  return(list(Y = Y, X = X))
}


# Yt = Y; maxrank = 5; cIndex = NULL; ofset = NULL
# control = nbfar_control();  nfold = 5
# set.seed(1234)
# nbrrr_test <- nbrrr(Yt, X, maxrank = 7, cIndex = NULL, ofset = NULL,
#       control = gofar_control(initmaxit = 5000, initepsilon = 1e-3), nfold = 5)
# plot(colMeans(nbrrr_test$cv.err))
# plot(log(nbrrr_test$objval[nbrrr_test$objval>0]))
# plot(log(nbrrr_test$diffobj[nbrrr_test$diffobj>0]))



#' Generalize Exclusive factor extraction via co-sparse unit-rank estimation (GOFAR(P)) using k-fold crossvalidation
#'
#' Divide and conquer approach for low-rank and sparse coefficent matrix estimation: Exclusive extraction
#'
#' @param Yt response matrix
#' @param X covariate matrix; when X = NULL, the fucntion performs unsupervised learning
#' @param maxrank an integer specifying the desired rank/number of factors
#' @param cIndex control index, specifying index of control variable in the design matrix X
#' @param ofset offset matrix specified
#' @param control a list of internal parameters controlling the model fitting
#' @param nfold number of fold for cross-validation
#' @return
#'   \item{C}{estimated coefficient matrix; based on GIC}
#'   \item{Z}{estimated control variable coefficient matrix}
#'   \item{PHI}{estimted dispersion parameters}
#'   \item{U}{estimated U matrix (generalize latent factor weights)}
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated V matrix (factor loadings)}
#' @export
#' @import magrittr
#' @useDynLib nbfar
#' @examples
#' \donttest{
#' ## Model specification:
#' SD <- 123
#' set.seed(SD)
#' p <- 100; n <- 200
#' pz <- 0
#' nrank <- 3                # true rank
#' rank.est <- 5             # estimated rank
#' nlam <- 20                # number of tuning parameter
#' s  = 0.5
#' q <- 30
#' control <- nbfar_control()  # control parameters
#' #
#' #
#' ## Generate data
#' D <- rep(0, nrank)
#' V <- matrix(0, ncol = nrank, nrow = q)
#' U <- matrix(0, ncol = nrank, nrow = p)
#' #
#' U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
#' U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
#' U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#' #
#'   # for similar type response type setting
#'   V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 16))
#'   V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 28))
#'   V[, 3] <- c(
#'     sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
#'     sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
#'   )
#' U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' #
#' D <- s * c(4, 6, 5) # signal strength varries as per the value of s
#' or <- order(D, decreasing = TRUE)
#' U <- U[, or]
#' V <- V[, or]
#' D <- D[or]
#' C <- U %*% (D * t(V)) # simulated coefficient matrix
#' intercept <- rep(0.5, q) # specifying intercept to the model:
#' C0 <- rbind(intercept, C)
#' #
#' Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
#' # Simulated data
#' sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = 3, depth = 10)  # Simulated sample
#' # Dispersion parameter
#' X <- sim.sample$X[1:n, ]
#' Y <- sim.sample$Y[1:n, ]
#' X0 <- cbind(1, X)                     # 1st column accounting for intercept
#'
#' # Model with known offset
#' set.seed(1234)
#' offset <- log(10)*matrix(1,n,q)
#' control_nbrr <- nbfar_control(initmaxit = 5000, initepsilon = 1e-4)
#' # nbrrr_test <- nbrrr(Y, X, maxrank = 5, cIndex = NULL, ofset = offset,
#' #                       control = control_nbrr, nfold = 5)
#' }
nbrrr <- function(Yt, X, maxrank = 10,
                  cIndex = NULL, ofset = NULL,
                  control = list(), nfold = 5) {
  cat("Initializing...", "\n")
  n <- nrow(Yt)
  p <- ncol(X)
  q <- ncol(Yt)
  X0 <- cbind(1, X)
  Y <- Yt

  if (is.null(cIndex)) {
    cIndex <- 1
    pz <- 0
    p <- ncol(X)
  } else {
    pz <- length(cIndex)
    p <- ncol(X) - pz
    cIndex <- c(1, cIndex + 1)
  }

  ## Initialization
  if (is.null(ofset)) ofset <- matrix(0, nrow = n, ncol = q)
  Yin <- Y
  naind <- (!is.na(Y)) + 0 # matrix(1,n,q)
  misind <- any(naind == 0) + 0
  if (misind == 1) Yin[is.na(Y)] <- 0
  # init_model <- nbCol(Yin, X0, ofset, naind)
  init_model <- nbZeroSol(Yin, X0, cIndex, ofset, naind)
  # Z0 <- init_model$C[cIndex,, drop = FALSE]
  Z0 <- init_model$Z
  PHI0 <- init_model$PHI

  ## Begin exclusive extraction procedure
  N.nna <- sum(!is.na(Y))
  naind <- (!is.na(Y)) + 0
  ind.nna <- which(!is.na(Y))
  fit.nfold <- vector("list", maxrank)
  # generate kfold index
  ID <- rep(1:nfold, len = N.nna)
  ID <- sample(ID, N.nna, replace = FALSE)
  ## store the deviance of the test data
  dev <- matrix(NA, nfold, maxrank)
  sdcal <- matrix(NA, nfold, maxrank)
  tec <- rep(0, nfold)

  # rank selection via cross validation
  for (k in 1:maxrank) { # desired rank extraction
    cat('Cross-validation for rank r: ', k, '\n')
    fitT <- vector("list", nfold)
    for (ifold in 1:nfold) { # ifold=3; k = 3
      cat('Fold ', ifold, '\n')
      ind.test <- ind.nna[which(ID == ifold)]
      Yte <- Y
      Yte[-ind.test] <- NA
      Ytr <- Y
      Ytr[ind.test] <- NA

      naind2 <- (!is.na(Ytr)) + 0 # matrix(1,n,q)
      misind <- any(naind2 == 0) + 0
      Ytr[is.na(Ytr)] <- 0
      # tryCatch({
      fitT[[ifold]] <- nbrrr_cpp(Ytr, X0, k, cIndex, ofset,
                                 Z0,  PHI0, control,
                                 misind, naind2)
      # compute test sample error
      tem  <- nbrrr_likelihood(Yte, fitT[[ifold]]$mu, fitT[[ifold]]$eta,
                               fitT[[ifold]]$PHI, (!is.na(Yte)) + 0)
      dev[ifold, k] = tem[1];
      sdcal[ifold, k] = tem[2];
      tec[ifold] <- tem[3];
      # }, error=function(e) print("nbrrr execution issue!") )

    }
    fit.nfold[[k]] <- fitT
  }
  # select rank with lowest mean;
  dev.mean <- colMeans(dev, na.rm = T)
  rank_sel <- which.max(dev.mean)

  # compute model estimate for the selected rank
  misind <- any(naind == 0) + 0
  control_nbrr <- nbfar_control(initmaxit = 5000, initepsilon = 1e-7,
                                objI = 1)
  Y[is.na(Y)] <- 0
  out <- nbrrr_cpp(Y, X0, rank_sel, cIndex, ofset, Z0,  PHI0,
                   control_nbrr, misind, naind)
  out$cv.err <- dev
  XC <- X0[, -cIndex] %*% out$C[-cIndex, ]
  out$Z <- matrix(out$C[cIndex, ], ncol = q)
  svdxc <- svd(XC)
  nkran <- sum(svdxc$d > 1e-07)
  out$V <- svdxc$v[, 1:nkran, drop = FALSE]
  out$D <- svdxc$d[1:nkran] / sqrt(n)
  out$U <- out$C[-cIndex, ] %*% out$V %*%
    diag(1 / out$D, nrow = nkran, ncol = nkran)
  return(out)
}




#
# Get lambda max for the negattive binomial disttribution
#
# @param Y Multivariate response matrix
# @param X design matrix
# @param offset offset  matrix
#' @useDynLib nbfar
#' @import magrittr
#' @importFrom MASS glm.nb
lmax_nb = function(Y, X,offset){
  Y[is.na(Y)] <- 0
  mu <- exp(offset)
  ttt <- abs(crossprod(X,Y - mu)); maxtt <- max(ttt)
  out <- 2*maxtt
  out
}

# Get zero solution for the negative binomial regression
#
# @param Y Multivariate response matrix
# @param X0 design matrix
# @param c_index control index in the model
# @param ofset offset  matrix
# @param naind indicator vector for missing data in Y
#' @useDynLib nbfar
#' @import magrittr
#' @importFrom MASS glm.nb
#' @importFrom stats glm.control rnorm
nbZeroSol <- function(Y, X0, c_index, ofset, naind) {
  q  <- ncol(Y)
  PHI <- rep(1, q)
  Cini <- matrix(0, nrow = length(c_index), ncol = q)
  for (i in 1:q) {
    tryCatch({
      qqq <- naind[,i] == 1
      xx <- MASS::glm.nb(Y[qqq, i]~1)
      ft <- MASS::glm.nb(Y[qqq, i]~X0[qqq, c_index,drop = FALSE] +
                           offset(ofset[qqq, i]) -1,
                         control = glm.control(maxit = 1000, epsilon = 1e-3,
                                               trace = FALSE),
                         init.theta = xx$theta)
      Cini[, i] <- ft$coefficients
      PHI[i] <- ft$theta
    }, error = function(error_message) {
      return(NA)
    })
  }
  return(list(Z = Cini, PHI = PHI))
}
# nbZeroSol <- function(Y, X0, c_index, ofset, naind) {
#   q  <- ncol(Y)
#   PHI <- rep(0, q)
#   Cini <- matrix(0, nrow = length(c_index), ncol = q)
#   for (i in 1:q) {
#     qqq <- naind[,i] == 1
#     xx <- MASS::glm.nb(Y[qqq, i]~1,
#                        control = glm.control(maxit = 10000, epsilon = 1e-3,
#                                              trace = FALSE))
#     ft <- MASS::glm.nb(Y[qqq, i]~X0[qqq, c_index,drop = FALSE] +
#                          offset(ofset[qqq, i]) -1,
#                        control = glm.control(maxit = 10000, epsilon = 1e-3,
#                                              trace = FALSE),
#                        init.theta = xx$theta)
#     Cini[, i] <- ft$coefficients
#     PHI[i] <- ft$theta
#   }
#   return(list(Z = Cini, PHI = PHI))
# }







# load('aditya.rda')
# Yt = Y; X = X; maxrank = 5; cIndex = NULL; ofset = NULL
# control = gofar_control();  nfold = 5; nlambda = 40;
# PATH = FALSE

# set.seed(1234)
# nbrrr_test <- nbrrr(Yt, X, maxrank = 7, cIndex = NULL, ofset = NULL,
#       control = gofar_control(initmaxit = 5000, initepsilon = 1e-3), nfold = 5)
# plot(colMeans(nbrrr_test$cv.err))
# plot(log(nbrrr_test$objval[nbrrr_test$objval>0]))
# plot(log(nbrrr_test$diffobj[nbrrr_test$diffobj>0]))


#' Generalize Sequential factor extraction via co-sparse unit-rank estimation (GOFAR(S)) using k-fold crossvalidation
#'
#' Divide and conquer approach for low-rank and sparse coefficent matrix estimation: Sequential
#'
#' @param Yt response matrix
#' @param X covariate matrix; when X = NULL, the fucntion performs unsupervised learning
#' @param maxrank an integer specifying the desired rank/number of factors
#' @param nlambda number of lambda values to be used along each path
#' @param cIndex control index, specifying index of control variable in the design matrix X
#' @param ofset offset matrix specified
#' @param control a list of internal parameters controlling the model fitting
#' @param nfold number of folds in k-fold crossvalidation
#' @param PATH TRUE/FALSE for generating solution path of sequential estimate after cross-validation step
#' @param nthread number of thread to be used for parallelization
#' @return
#'   \item{C}{estimated coefficient matrix; based on GIC}
#'   \item{Z}{estimated control variable coefficient matrix}
#'   \item{Phi}{estimted dispersion parameters}
#'   \item{U}{estimated U matrix (generalize latent factor weights)}
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated V matrix (factor loadings)}
#' @export
#' @useDynLib nbfar
#' @importFrom RcppParallel RcppParallelLibs setThreadOptions
#' @examples
#' \donttest{
#' ## Model specification:
#' SD <- 123
#' set.seed(SD)
#' p <- 100; n <- 200
#' pz <- 0
#' nrank <- 3                # true rank
#' rank.est <- 5             # estimated rank
#' nlam <- 20                # number of tuning parameter
#' s  = 0.5
#' q <- 30
#' control <- nbfar_control()  # control parameters
#' #
#' #
#' ## Generate data
#' D <- rep(0, nrank)
#' V <- matrix(0, ncol = nrank, nrow = q)
#' U <- matrix(0, ncol = nrank, nrow = p)
#' #
#' U[, 1] <- c(sample(c(1, -1), 8, replace = TRUE), rep(0, p - 8))
#' U[, 2] <- c(rep(0, 5), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 14))
#' U[, 3] <- c(rep(0, 11), sample(c(1, -1), 9, replace = TRUE), rep(0, p - 20))
#' #
#'   # for similar type response type setting
#'   V[, 1] <- c(rep(0, 8), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 16))
#'   V[, 2] <- c(rep(0, 20), sample(c(1, -1), 8,
#'     replace =
#'       TRUE
#'   ) * runif(8, 0.3, 1), rep(0, q - 28))
#'   V[, 3] <- c(
#'     sample(c(1, -1), 5, replace = TRUE) * runif(5, 0.3, 1), rep(0, 23),
#'     sample(c(1, -1), 2, replace = TRUE) * runif(2, 0.3, 1), rep(0, q - 30)
#'   )
#' U[, 1:3] <- apply(U[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' V[, 1:3] <- apply(V[, 1:3], 2, function(x) x / sqrt(sum(x^2)))
#' #
#' D <- s * c(4, 6, 5) # signal strength varries as per the value of s
#' or <- order(D, decreasing = TRUE)
#' U <- U[, or]
#' V <- V[, or]
#' D <- D[or]
#' C <- U %*% (D * t(V)) # simulated coefficient matrix
#' intercept <- rep(0.5, q) # specifying intercept to the model:
#' C0 <- rbind(intercept, C)
#' #
#' Xsigma <- 0.5^abs(outer(1:p, 1:p, FUN = "-"))
#' # Simulated data
#' sim.sample <- nbfar_sim(U, D, V, n, Xsigma, C0,disp = 3, depth = 10)  # Simulated sample
#' # Dispersion parameter
#' X <- sim.sample$X[1:n, ]
#' Y <- sim.sample$Y[1:n, ]
#' X0 <- cbind(1, X)                     # 1st column accounting for intercept
#'
#' # Model with known offset
#' set.seed(1234)
#' offset <- log(10)*matrix(1,n,q)
#' control_nbfar <- nbfar_control(initmaxit = 5000, gamma0 = 2, spU = 0.5,
#' spV = 0.6, lamMinFac = 1e-10, epsilon = 1e-5)
#' # nbfar_test <- nbfar(Y, X, maxrank = 5, nlambda = 20, cIndex = NULL,
#' # ofset = offset, control = control_nbfar, nfold = 5, PATH = F)
#' }
nbfar <- function(Yt, X, maxrank = 3, nlambda = 40, cIndex = NULL,
                  ofset = NULL, control = list(), nfold = 5,
                  PATH = FALSE, nthread = 1) {
  cat("Initializing...", "\n")
  n <- nrow(Yt)
  p <- ncol(X)
  q <- ncol(Yt)
  X0 <- cbind(1, X)
  Y <- Yt

  if (is.null(cIndex)) {
    cIndex <- 1
    pz <- 0
    p <- ncol(X)
  } else {
    pz <- length(cIndex)
    p <- ncol(X) - pz
    cIndex <- c(1, cIndex + 1)
  }

  # define
  U <- matrix(0, nrow = p, ncol = maxrank)
  V <- matrix(0, nrow = q, ncol = maxrank)
  D <- rep(0, maxrank)
  Z <- matrix(0, nrow = pz + 1, ncol = q)
  PHI <- rep(0, q)
  lamSel <- rep(0, maxrank)


  if (is.null(ofset)) ofset <- matrix(0, nrow = n, ncol = q)
  naind <- (!is.na(Y)) + 0
  aft <- nbZeroSol(Y, X0, c_index = cIndex, ofset, naind)
  Z <- aft$Z
  PHI <- aft$PHI


  N.nna <- sum(!is.na(Y))
  ind.nna <- which(!is.na(Y))
  Yf2 <- Y
  naind22 <- (!is.na(Yf2)) + 0
  misind22 <- any(naind22 == 0) + 0
  Yf2[is.na(Yf2)] <- 0
  fit.nlayer <- fit.nfold <- vector("list", maxrank)




  # Implementation of k-fold cross validation:
  for (k in 1:maxrank) { # desired rank extraction # k = 1
    cat("Initializing unit-rank unit", k, "\n")
    control_nbrrr <- control # nbfar_control()
    xx <- nbrrr_cpp(Yf2, X0, 1, cIndex, ofset,
                    Z,  PHI, control_nbrrr,
                    misind22, naind22)
    # print(c(xx$maxit,xx$converge,ndev))
    XC <- X0[, -cIndex] %*% xx$C[-cIndex, ]
    Z0 <- matrix(xx$C[cIndex, ], ncol = q)
    Phi0 <- xx$PHI
    svdxc <- svd(XC)
    nkran <- sum(svdxc$d > 1e-07)
    xx$V <- svdxc$v[, 1:nkran, drop = FALSE]
    xx$D <- svdxc$d[1:nkran] / sqrt(n)
    xx$U <- xx$C[-cIndex, ] %*% xx$V %*% diag(1 / xx$D, nrow = nkran, ncol = nkran)
    initW <- list(
      wu = abs(xx$U)^-control$gamma0,
      wd = abs(xx$D)^-control$gamma0,
      wv = abs(xx$V)^-control$gamma0)
    lambda.max <- lmax_nb(Y, X,ofset)
    zerosol <- nbZeroSol(Yf2, X0, c_index = cIndex, ofset, naind22)
    cat(xx$D, '\n')
    cat("Cross validation:", k, "\n")


    # tic()
    ## store the deviance of the test data
    if (nthread > 1){
      outcv = cv_nbfar_par(Y, X0, nlam = nlambda, cindex = cIndex,
                           ofset = ofset, initw = initW,
                           Zini = Z0, PhiIni = Phi0,
                           xx = xx,
                           lmax = lambda.max, control = control,
                           zerosol,
                           control$maxit, epsilon = control$epsilon,
                           nfold = 5);
    } else {
      outcv = cv_nbfar_cpp(Y, X0, nlam = nlambda, cindex = cIndex,
                           ofset = ofset, initw = initW,
                           Zini = Z0, PhiIni = Phi0,
                           xx = xx,
                           lmax = lambda.max, control = control,
                           zerosol,
                           control$maxit, epsilon = control$epsilon,
                           nfold = 5);
    }
    # toc()
    dev.mean <- colMeans(outcv$dev, na.rm = FALSE)
    l.mean <- which.max(dev.mean)
    lamS <- outcv$lamseq[l.mean]


    Yf <- Y
    naind <- (!is.na(Yf)) + 0 # matrix(1,n,q)
    misind <- any(naind == 0) + 0
    Yf[is.na(Yf)] <- 0
    # save(list=ls(),file= 'aditya.rda')
    # zerosol <- nbZeroSol(Yf, X0, c_index= cIndex, ofset, naind)
    if ( PATH ) {
      control1 <- nbfar_control(lamMinFac = 1e-10, spU  = 0.9, spV = 0.9)
      # control1$spU  = 0.9; control1$spV = 0.9
      fit.nlayer[[k]] <- fit.layer <- nbfar_cpp(Yf,X0, nlam = nlambda,
                                                cindex = cIndex,
                                                ofset = ofset,
                                                initw = initW,
                                                Dini = xx$D,
                                                Zini = Z0 ,
                                                PhiIni = Phi0,
                                                Uini = xx$U,
                                                Vini = xx$V,
                                                lmax = lambda.max,
                                                control1,misind,
                                                naind,zerosol,
                                                control$maxit,
                                                control$epsilon)
      l.mean <- which(fit.layer$lamKpath == lamS)
      if ( length(l.mean) == 0) {
        l.mean <- 1;
      }
    } else {
      control1 <- control
      control1$lamMinFac <- 1; control1$epsilon = 1e-12
      fit.nlayer[[k]] <- fit.layer <- nbfar_cpp(Yf,X0, nlam = 1,
                                                cindex = cIndex,
                                                ofset = ofset,
                                                initw = initW,
                                                Dini = xx$D,
                                                Zini = Z0 ,
                                                PhiIni = Phi0,
                                                Uini = xx$U,
                                                Vini = xx$V,
                                                lmax = lamS,
                                                control1,misind,
                                                naind,zerosol,
                                                control$maxit,
                                                control$epsilon)
      l.mean <- 1
    }

    ## Same as before
    fit.nlayer[[k]]$initw <- initW
    fit.nlayer[[k]]$lamseq <- outcv$lamseq   ## save sequence of lambda path
    fit.nlayer[[k]]$dev <- outcv$dev
    fit.nlayer[[k]]$lamS <- lamS
    U[, k] <- fit.layer$ukpath[, l.mean]
    V[, k] <- fit.layer$vkpath[, l.mean]
    D[k] <- fit.layer$dkpath[l.mean, 1]
    cat(D[k], '\n')
    lamSel[k] <- fit.layer$lamKpath[l.mean, 1]

    if (D[k] == 0) {
      U <- matrix(U[, 1:(k - 1)], ncol = k - 1)
      V <- matrix(V[, 1:(k - 1)], ncol = k - 1)
      D <- D[1:(k - 1)]
      lamSel <- lamSel[1:(k - 1)]
      break
    }

    Ck <- D[k] * tcrossprod(U[, k], V[, k])
    PHI <- fit.layer$phipath[, l.mean]
    if (pz == 0) {
      Z[1, ] <-  drop(fit.layer$zpath[, ,l.mean])
    } else {
      Z <- fit.layer$zpath[, , l.mean]
    }
    ofset <- ofset + crossprod(t(X0[, -cIndex]), Ck) # crossprod(t(X), Ck)
  }

  cat("Estimated rank =", sum(D != 0), "\n")
  if (sum(D != 0) == maxrank) {
    cat("Increase maxrank value!")
  }
  ft1 <- list(fit = fit.nlayer, C = U %*% (D * t(V)), Z = Z, Phi = PHI,
              U = U, V = V, D = D,
              lam = lamSel)
  return(ft1)
}





# ## store the deviance of the test data
# dev <- matrix(NA, nfold, nlambda)
# sdcal <- matrix(NA, nfold, nlambda)
# tec <- rep(0, nfold)
# ID <- rep(1:nfold, len = N.nna)
# ID <- sample(ID, N.nna, replace = FALSE)
# fitT <- vector("list", nfold)
#
# for (ifold in 1:nfold) { # ifold=5
#   ind.test <- ind.nna[which(ID == ifold)]
#   Yte <- Y
#   Yte[-ind.test] <- NA
#   Ytr <- Y
#   Ytr[ind.test] <- NA
#
#   naind <- (!is.na(Ytr)) + 0 # matrix(1,n,q)
#   misind <- any(naind == 0) + 0
#   Ytr[is.na(Ytr)] <- 0
#   # zerosol <- nbZeroSol(Ytr, X0, c_index= cIndex, ofset, naind)
#   # compromised convergence criteria for enhanced speed
#   control_cv <- control; #control_cv$lamMinFac = 1e-10;
#   # control_cv$spU  = 0.8; control_cv$spV = 0.9
#   cat('Fold ', ifold, '\n')
#   tic()
#   fitT[[ifold]] <- nbfar_cpp(Ytr, X0,
#                              nlam = nlambda,
#                              cindex = cIndex,
#                              ofset = ofset, initw = initW,
#                              Dini = xx$D, Zini = Z0, PhiIni = Phi0,
#                              Uini = xx$U, Vini = xx$V,
#                              lmax = lambda.max,
#                              control_cv, misind,
#                              naind, zerosol,
#                              control$maxit, epsilon = control$epsilon)
#   fitF <- fitT[[ifold]]
#   a = toc()
#   cat('Fold ', ifold, 'done ', a$toc -a$tic,'\n')
#
#   for (im in 1:fitF$nkpath) {
#     lam <- fitF$lamKpath[im, 1]
#     insel <- which(fitF$lamseq == lam)
#     mu.test <- fitF$mukpath[, , im]; mu.test[-ind.test] <- NA
#     eta.test <- fitF$etapath[, , im]; eta.test[-ind.test] <- NA
#     tttval <- nbrrr_likelihood(Yte, mu.test, eta.test, fitF$phipath[, im],
#                                (!is.na(Yte)) + 0)
#     dev[ifold, insel] <- tttval[1]
#     sdcal[ifold, insel] <- tttval[2]
#     tec[ifold] <- tttval[3]
#   }
# }
# fit.nfold[[k]] <- fitT


