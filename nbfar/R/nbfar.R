#' Control parameters for NBFAR and NBRRR
#'
#' Default value for a list of control parameters that are used to estimate the parameters of negative binomial co-sparse factor regression (NBFAR) and negative binomial reduced rank regression (NBRRR).
#'
#'
#' @param lamMaxFac a multiplier of the computed maximum value (lambda_max) of the tuning parameter
#' @param lamMinFac a multiplier to determine lambda_min as a fraction of lambda_max
#' @param gamma0 power parameter for generating the adaptive weights
#' @param elnetAlpha elastic net penalty parameter
#' @param spU maximum proportion of nonzero elements in each column of U
#' @param spV maximum proportion of nonzero elements in each column of V
#' @param maxit maximum iteration for each sequential steps
#' @param epsilon tolerance value required for convergence of inner loop in GCURE
#' @param objI 1 or 0  to indicate that the convergence will be on the basis of objective function or not
#' @param initmaxit maximum iteration for minimizing the objective function while computing the initial estimates of the model parameter
#' @param initepsilon tolerance value required for the convergence of the objective function while computing the initial estimates of the model parameter
#' @return a list of controlling parameter.
#' @export
#' @examples
#' control <- nbfar_control()
#' @useDynLib nbfar
#' @references
#' Mishra, A., Müller, C. (2022) \emph{Negative binomial factor regression models with application to microbiome data analysis.  https://doi.org/10.1101/2021.11.29.470304}
nbfar_control <- function(maxit = 5000, epsilon = 1e-7,
                          elnetAlpha = 0.95,
                          gamma0 = 1,
                          spU = 0.5, spV = 0.5,
                          lamMaxFac = 1, lamMinFac = 1e-6,
                          initmaxit = 10000, initepsilon = 1e-8,
                          objI = 0) {
  list(
    lamMaxFac = lamMaxFac,
    lamMinFac = lamMinFac,
    gamma0 = gamma0, elnetAlpha = elnetAlpha, spU = spU, spV = spV,
    maxit = maxit, epsilon = epsilon,
    objI = objI,
    initmaxit = initmaxit, initepsilon = initepsilon
  )
}







#' Simulated data for testing NBFAR and NBRRR model
#'
#' Simulate response and covariates for multivariate negative binomial regression with a low-rank and sparse coefficient matrix. Coefficient matrix is expressed in terms of U (left singular vector), D (singular values) and V (right singular vector).
#'
#' @param U specified value of U
#' @param V specified value of V
#' @param D specified value of D
#' @param n sample size
#' @param Xsigma covariance matrix used to generate predictors in X
#' @param C0 intercept value in the coefficient matrix
#' @param disp dispersion parameter of the generative model
#' @param depth log of the sequencing depth of the microbiome data (used as an offset in the simulated multivariate negative binomial regression model)
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
#' # disp = 3; depth = 10;
#' # simulate_nbfar <- list(Y = Y,X = X, U = U, D = D, V = V, n=n,
#' # Xsigma = Xsigma, C0 = C0,disp =disp, depth =depth)
#' # save(simulate_nbfar, file = 'data/simulate_nbfar.RData')
#' @references
#' Mishra, A., Müller, C. (2022) \emph{Negative binomial factor regression models with application to microbiome data analysis.  https://doi.org/10.1101/2021.11.29.470304}
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




#' Negative binomial reduced rank regression (NBRRR)
#'
#' In the range of 1 to maxrank, the estimation procedure selects the rank r of the coefficient matrix using a cross-validation approach. For the selected rank, a rank r coefficient matrix is estimated that best fits the observations.
#'
#' @param Yt response matrix
#' @param X design matrix; when X = NULL, we set X as identity matrix and perform generalized PCA.
#' @param maxrank an integer specifying the maximum possible rank of the coefficient matrix or the number of factors
#' @param cIndex specify index of control variable in the design matrix X
#' @param ofset offset matrix or microbiome data analysis specific scaling: common sum scaling = CSS (default), total sum scaling = TSS, median-ratio scaling = MRS, centered-log-ratio scaling  = CLR
#' @param control a list of internal parameters controlling the model fitting
#' @param nfold number of folds in k-fold crossvalidation
#' @param trace TRUE/FALSE checking progress of cross validation error
#' @param verbose TRUE/FALSE checking progress of estimation procedure
#' @return
#'   \item{C}{estimated coefficient matrix}
#'   \item{Z}{estimated control variable coefficient matrix}
#'   \item{PHI}{estimted dispersion parameters}
#'   \item{U}{estimated U matrix (generalize latent factor weights)}
#'   \item{D}{estimated singular values}
#'   \item{V}{estimated V matrix (factor loadings)}
#' @export
#' @import magrittr
#' @importFrom Rcpp evalCpp
#' @importFrom graphics title
#' @useDynLib nbfar
#' @examples
#' \donttest{
#' ## Load simulated data set:
#' data('simulate_nbfar')
#' attach(simulate_nbfar)
#'
#' # Model with known offset
#' set.seed(1234)
#' offset <- log(10)*matrix(1,n,ncol(Y))
#' control_nbrr <- nbfar_control(initmaxit = 5000, initepsilon = 1e-4)
#' # nbrrr_test <- nbrrr(Y, X, maxrank = 5, cIndex = NULL, ofset = offset,
#' #                       control = control_nbrr, nfold = 5)
#' }
#' @references
#' Mishra, A., Müller, C. (2022) \emph{Negative binomial factor regression models with application to microbiome data analysis.  https://doi.org/10.1101/2021.11.29.470304}
nbrrr <- function(Yt, X, maxrank = 10,
                  cIndex = NULL, ofset = 'CSS',
                  control = list(), nfold = 5, trace = FALSE, verbose = TRUE) {
  if(verbose) cat("Initializing...", "\n")
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
  ofset <- offset_sacling(Y, ofset)
  # if (is.null(ofset)) ofset <- matrix(0, nrow = n, ncol = q)
  Yin <- Y
  naind <- (!is.na(Y)) + 0 # matrix(1,n,q)
  misind <- any(naind == 0) + 0
  if (misind == 1) Yin[is.na(Y)] <- 0
  # init_model <- nbCol(Yin, X0, ofset, naind)
  # init_model <- nbColSp(Y, X0, ofset, cIndex,  naind)
  # init_model <- nbZeroSol(Yin, X0, cIndex, ofset, naind)
  # Z0 <- init_model$C[cIndex,, drop = FALSE]
  # init_model <- nbColSp(Y, X0, ofset, cIndex,  naind)
  init_model <- nbzerosol_cpp(Yin, X0[,cIndex,drop = FALSE], ofset, control,
                              misind, naind)
  Z0 <- init_model$Z
  PHI0 <- init_model$PHI
  ETA0 <- ofset + X0[, cIndex] %*% Z0
  MU0 <- exp(ETA0)
  # C0 <- init_model$C
  # XC <- X0[, -cIndex] %*% init_model$C; svdxc <- svd(XC)
  # Vk <- svdxc$v[, 1:maxrank, drop = FALSE]
  # Dk <- svdxc$d[1:maxrank] / sqrt(n)
  # Uk <- init_model$C %*% Vk %*%
  #   diag(1 / Dk, nrow = maxrank, ncol = maxrank)
  # cat('Singular values inittial:', Dk, '\n')


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
  dev0 <- rep(0, nfold)

  # rank selection via cross validation
  for (k in 1:maxrank) { # desired rank extraction
    if(verbose) cat('Cross-validation for rank r: ', k, '\n')
    fitT <- vector("list", nfold)
    # C0 <- Uk[, 1:k] %*% (Dk[1:k]*t(Vk[, 1:k]))
    for (ifold in 1:nfold) { # ifold=3; k = 3
      ind.test <- ind.nna[which(ID == ifold)]
      Yte <- Y
      Yte[-ind.test] <- NA
      Ytr <- Y
      Ytr[ind.test] <- NA

      naind2 <- (!is.na(Ytr)) + 0 # matrix(1,n,q)
      misind <- any(naind2 == 0) + 0
      Ytr[is.na(Ytr)] <- 0

      fitT[[ifold]] <- nbrrr_cpp(Ytr, X0, k, cIndex, ofset,
                                 Z0,  PHI0, matrix(0,p,q), control,
                                 misind, naind2)

      cat('Fold ', ifold, ': [Error,iteration] = [',
          with(fitT[[ifold]],diffobj[maxit]),
          with(fitT[[ifold]],maxit), ']', '\n')

      # compute test sample error
      tem  <- nbrrr_likelihood(Yte, fitT[[ifold]]$mu, fitT[[ifold]]$eta,
                               fitT[[ifold]]$PHI, (!is.na(Yte)) + 0)
      dev0[ifold]  <- nbrrr_likelihood(Yte, MU0, ETA0, PHI0, (!is.na(Yte)) + 0)[1]


      dev[ifold, k] = tem[1];
      sdcal[ifold, k] = tem[2];
      tec[ifold] <- tem[3];

    }
    if (trace) {
      plot(colMeans(dev, na.rm = T),
           xlab = 'Rank ', ylab = 'Log-Likelihood')
      graphics::title(paste('Component ', k))
      Sys.sleep(0.5)
      }

    fit.nfold[[k]] <- fitT
  }
  # select rank with lowest mean;
  dev.mean <- colMeans(cbind(dev0,dev), na.rm = T)
  rank_sel <- which.max(dev.mean)-1

  # compute model estimate for the selected rank
  misind <- any(naind == 0) + 0
  control_nbrr <- control
  control_nbrr$initepsilon <- 0.01*control_nbrr$initepsilon
  Y[is.na(Y)] <- 0
  # C0 <- Uk[, 1:rank_sel] %*% (Dk[1:rank_sel]*t(Vk[, 1:rank_sel]))
  if (rank_sel > 0){
    out <- nbrrr_cpp(Y, X0, rank_sel, cIndex, ofset, Z0,  PHI0, matrix(0,p,q),
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
    out$Y <- Yt;  out$X = X
  } else {
    C <- matrix(rep(0,ncol(X0)*q), ncol(X0),q)
    C[cIndex,] <- Z0
    out = list(C = C, PHI = PHI0, sc = NA, sb =NA, mu = MU0,
               objval = NA, diffobj =NA, converged = NA, ExecTimekpath = NA,
               maxit = NA, converge = NA, Z = Z0, D = 0,
               U = matrix(rep(0,ncol(X0) -length(cIndex))),
               V = matrix(rep(0,q)), Y = Yt, X =X)
  }



  return(out)
}





#' Negative binomial co-sparse factor regression  (NBFAR)
#'
#' To estimate a low-rank and sparse coefficient matrix in large/high dimensional setting, the approach extracts unit-rank components of required matrix in sequential order. The algorithm automatically stops after extracting sufficient unit rank components.
#'
#' @param Yt response matrix
#' @param X design matrix; when X = NULL, we set X as identity matrix and perform generalized sparse PCA.
#' @param maxrank an integer specifying the maximum possible rank of the coefficient matrix or the number of factors
#' @param nlambda number of lambda values to be used along each path
#' @param cIndex specify index of control variables in the design matrix X
#' @param ofset offset matrix or microbiome data analysis specific scaling: common sum scaling = CSS (default), total sum scaling = TSS, median-ratio scaling = MRS, centered-log-ratio scaling  = CLR
#' @param control a list of internal parameters controlling the model fitting
#' @param nfold number of folds in k-fold crossvalidation
#' @param PATH TRUE/FALSE for generating solution path of sequential estimate after cross-validation step
#' @param nthread number of thread to be used for parallelizing the crossvalidation procedure
#' @param trace TRUE/FALSE checking progress of cross validation error
#' @param verbose TRUE/FALSE checking progress of estimation procedure
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
#' @importFrom Rcpp evalCpp
#' @importFrom graphics title
#' @examples
#' \donttest{
#' ## Load simulated data set:
#' data('simulate_nbfar')
#' attach(simulate_nbfar)
#'
#' # Model with known offset
#' set.seed(1234)
#' offset <- log(10)*matrix(1,n,ncol(Y))
#' control_nbfar <- nbfar_control(initmaxit = 5000, gamma0 = 2, spU = 0.5,
#' spV = 0.6, lamMinFac = 1e-10, epsilon = 1e-5)
#' # nbfar_test <- nbfar(Y, X, maxrank = 5, nlambda = 20, cIndex = NULL,
#' # ofset = offset, control = control_nbfar, nfold = 5, PATH = F)
#' }
#' @references
#' Mishra, A., Müller, C. (2022) \emph{Negative binomial factor regression models with application to microbiome data analysis.  https://doi.org/10.1101/2021.11.29.470304}
nbfar <- function(Yt, X, maxrank = 3, nlambda = 40, cIndex = NULL,
                  ofset = 'CSS', control = list(), nfold = 5,
                  PATH = FALSE, nthread = 1, trace = FALSE, verbose = TRUE) {
  if(verbose) cat("Initializing...", "\n")
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


  ofset <- offset_sacling(Y, ofset)
  # if (is.null(ofset)) ofset <- matrix(0, nrow = n, ncol = q)
  # naind <- (!is.na(Y)) + 0
  # aft <- nbZeroSol(Y, X0, c_index = cIndex, ofset, naind)
  # Z <- aft$Z
  # PHI <- aft$PHI

  # aft <- nbzerosol_cpp(Y, X0[,cIndex,drop = FALSE], ofset, control,
  #                      any(naind == 0) + 0, naind)
  # Z <- aft$Z
  # PHI <- aft$PHI
  # aft <- nbColSp(Y , X0, ofset, cIndex,  naind)
  # Z <- aft$Z
  # PHI <- aft$PHI
  # XC <- X0[, -cIndex] %*% aft$C; svdxc <- svd(XC)
  # Vk <- svdxc$v[, 1:maxrank, drop = FALSE]
  # Dk <- svdxc$d[1:maxrank] / sqrt(n)
  # Uk <- aft$C %*% Vk %*% diag(1 / Dk, nrow = maxrank, ncol = maxrank)
  # cat('Singular values initial:', Dk, '\n')
  # save(list = ls(), file = 'aditya.rda')

  Yf2 <- Y
  naind22 <- (!is.na(Yf2)) + 0
  misind22 <- any(naind22 == 0) + 0
  Yf2[is.na(Yf2)] <- 0
  aft <- nbzerosol_cpp(Yf2, X0[,cIndex,drop = FALSE], ofset, control,
                       misind22, naind22)
  Z <- aft$Z
  PHI <- aft$PHI
  fit.nlayer <- vector("list", maxrank)

  # Implementation of k-fold cross validation:
  for (k in 1:maxrank) { # desired rank extraction # k = 1
    control_nbrrr <- control # nbfar_control()
    control_nbrrr$objI <- 1
    xx <- nbrrr_cpp(Yf2, X0, 1, cIndex, ofset,
                    Z,  PHI, matrix(0,p,q),#0*aft$C, # -  0*U[, 1:k] %*% (D[1:k]*t( V[, 1:k])),
                    control_nbrrr,
                    misind22, naind22)

    XC <- X0[, -cIndex] %*% xx$C[-cIndex, ]
    Z0 <- matrix(xx$C[cIndex, ], ncol = q)
    Phi0 <- xx$PHI
    svdxc <- svd(XC)
    nkran <- sum(svdxc$d > 1e-07)
    xx$V <- svdxc$v[, 1:nkran, drop = FALSE]
    xx$D <- svdxc$d[1:nkran] / sqrt(n)
    xx$U <- xx$C[-cIndex, ] %*% xx$V %*%
      diag(1 / xx$D, nrow = nkran, ncol = nkran)
    initW <- list(
      wu = abs(xx$U)^-control$gamma0,
      wd = abs(xx$D)^-control$gamma0,
      wv = abs(xx$V)^-control$gamma0)
    lambda.max <- lmax_nb(Y, X,ofset)
    if(verbose) cat("Initializing unit-rank unit ", k, ': [Error,iteration,D] = [',
        with(xx,diffobj[maxit]), xx$maxit, xx$D, ']', '\n')
    if(verbose) cat("Cross validation for component:", k, "\n")

    ## store the deviance of the test data
    if (nthread > 1) {
      outcv = cv_nbfar_par(Y, X0, nlam = nlambda, cindex = cIndex,
                           ofset = ofset, initw = initW,
                           Zini = Z0, PhiIni = Phi0,
                           xx = xx,
                           lmax = lambda.max, control = control,
                           control$maxit, epsilon = control$epsilon,
                           nfold = nfold);
    } else {

      outcv = cv_nbfar_cpp(Y, X0, nlam = nlambda, cindex = cIndex,
                           ofset = ofset, initw = initW,
                           Zini = Z0, PhiIni = Phi0,
                           xx = xx,
                           lmax = lambda.max, control = control,
                           control$maxit, epsilon = control$epsilon,
                           nfold = nfold);
    }
    # toc()
    dev.mean <- colMeans(outcv$dev, na.rm = FALSE)
    l.mean <- max( which.max(dev.mean) - 1 , 1 )
    lamS <- outcv$lamseq[l.mean]
    if (trace) {
      plot(colMeans(outcv$dev, na.rm = FALSE),
           xlab = 'Lambda index', ylab = 'Log-likelihood')
      graphics::title(paste('Component ', k))
      Sys.sleep(0.5)
      }

    Yf <- Y
    naind <- (!is.na(Yf)) + 0 # matrix(1,n,q)
    misind <- any(naind == 0) + 0
    Yf[is.na(Yf)] <- 0

    # zerosol <- nbZeroSol(Yf, X0, c_index= cIndex, ofset, naind)
    if ( PATH ) {
      control1 <- control
      control1$epsilon = 0.001*control$epsilon
      control1$spU  = 0.9; control1$spV = 0.9
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
                                                naind,
                                                control1$maxit,
                                                control1$epsilon)
      l.mean <- which(fit.layer$lamKpath == lamS)
      if ( length(l.mean) == 0) {
        l.mean <- 1;
      }
    } else {
      control1 <- control
      control1$lamMinFac <- 1; control1$epsilon = 0.001*control$epsilon

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
                                                naind,
                                                control1$maxit,
                                                control1$epsilon)
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

    if ((D[k] == 0) ) {
      if(k > 1){
        U <- matrix(U[, 1:(k - 1)], ncol = k - 1)
        V <- matrix(V[, 1:(k - 1)], ncol = k - 1)
        D <- D[1:(k - 1)]
        lamSel <- lamSel[1:(k - 1)]
      }
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
              lam = lamSel, Y = Yt, X = X)
  return(ft1)
}

