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



#' Suitably generates offset matrix for the multivariate regression problem
#'
#' @param Y outcome matrix
#' @param ofset offset matrix or microbiome data analysis specific scaling: common sum scaling = CSS (default), total sum scaling = TSS, median-ratio scaling = MRS, centered-log-ratio scaling  = CLR
#' @importFrom stats median
offset_sacling = function(Y, ofset){
  n <- nrow(Y)
  q <- ncol(Y)
  if (is.matrix(ofset)){
    return(ofset)
  } else if (tolower(ofset) == 'css'){
    # ofset <- matrix(0, nrow = n, ncol = q)
    ofset <- matrix(log(rowSums(Y, na.rm = T)),n,q)
    return(ofset)
  } else if (tolower(ofset) == 'tss'){
    ofset <- matrix(log(rowSums(Y,na.rm = T)),n,q)
    return(ofset)
  } else if (tolower(ofset) == 'mrs'){
    tem <- exp(rowSums(log(Y + (Y == 0)), na.rm = T)/q)
    ofset <- matrix(log(apply(Y/tem,1,median, na.rm=T)),n,q)
    return(ofset)
  } else if (tolower(ofset) == 'clr'){
    ofset <- matrix(rowSums(log(Y + (Y == 0)), na.rm = T)/q,n,q)
    return(ofset)
  }
}






#' @importFrom mpath cv.glmregNB
nbColSp <- function(Y, X0, ofset, cindex,  naind) {
  q  <- ncol(Y)
  PHI <- rep(1, q)
  Cini <- matrix(0, nrow = ncol(X0), ncol = q)
  pf = rep(1, ncol(X0)); pf[cindex] <- 0
  for (i in 1:q) {
    if ((i %% 5) == 0 ) cat('Initialization upto column index: ',i,'\n')
    tryCatch({
      qqq <- naind[,i] == 1
      dt = data.frame(y = Y[qqq, i], X = X0[qqq, -1])
      ft <- mpath::glmregNB(y~., data = dt, maxit = 100, nlambda = 2,
                            offset = ofset[qqq, i], thresh = 1e-2,
                            penalty.factor  = pf[-1],maxit.theta = 10,
                            # lambda = seq(0.001, 1, length.out = 5),
                            # lambda.min.ratio = 0.0001,
                            # init.theta = MASS::glm.nb(Y[qqq, i]~1)$theta,
                            alpha = 0.5,
                            rescale = FALSE, standardize = FALSE)
      tem = coef(ft)[,2]; tem[is.na(tem)] <- 0
      Cini[, i] <- tem
      PHI[i] <- ft$theta[2]
    }, error = function(error_message) {
      return(NA)
    })
  }
  return(list(Z = Cini[cindex,,drop = FALSE],
              C = Cini[-cindex,,drop = FALSE],
              PHI = PHI))
}







#' @importFrom mpath cv.glmregNB
nbColSpx <- function(Y, X0, ofset, cindex,  naind) {
  q  <- ncol(Y)
  PHI <- rep(1, q)
  Cini <- matrix(0, nrow = ncol(X0), ncol = q)
  pf = rep(1, ncol(X0)); pf[cindex] <- 0
  for (i in 1:q) {
    if ((i %% 5) == 0 ) cat('Initialization upto column index: ',i,'\n')
    tryCatch({
      qqq <- naind[,i] == 1
      dt = data.frame(y = Y[qqq, i], X = X0[qqq, -1])
      ft <- mpath::cv.glmregNB(y~., data = dt, maxit = 100, nfolds = 4,
                               offset = ofset[qqq, i], thresh = 1e-2,
                               penalty.factor  = pf[-1],maxit.theta = 10,
                               # nlambda = 5,
                               lambda = seq(0.001, 1, length.out = 5),
                               # init.theta = MASS::glm.nb(Y[qqq, i]~1)$theta,
                               alpha = 0.95,
                               rescale = FALSE, standardize = FALSE,
                               plot.it = FALSE)
      tem = coef(ft); tem[is.na(tem)] <- 0
      Cini[, i] <- tem
      PHI[i] <- ft$fit$theta[ft$lambda.which]
    }, error = function(error_message) {
      return(NA)
    })
  }
  return(list(Z = Cini[cindex,,drop = FALSE],
              C = Cini[-cindex,,drop = FALSE],
              PHI = PHI))
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

