bridge_select <- function(type1 = "trunc", type2 = "continuous") {
  if (type1 == "binary" & type2 == "binary") { bridge_select <- bridgeF_bb
  } else if (type1 == "trunc" & type2 == "trunc") { bridge_select <- bridgeF_tt
  } else if (type1 == "trunc" & type2 == "continuous") { bridge_select <- bridgeF_tc
  } else if (type1 == "continuous" & type2 == "trunc") { bridge_select <- bridgeF_ct
  } else if (type1 == "binary" & type2 == "continuous") { bridge_select <- bridgeF_bc
  } else if (type1 == "continuous" & type2 == "binary") { bridge_select <- bridgeF_cb
  } else if (type1 == "trunc" & type2 == "binary") { bridge_select <- bridgeF_tb
  } else if (type1 == "binary" & type2 == "trunc") { bridge_select <- bridgeF_bt
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}
bridgeF_bc <- function(r, zratio1, zratio2 = NULL){
  # binary and continuous
  de1 <- stats::qnorm(zratio1)
  res <- as.numeric( 4*fMultivar::pnorm2d(de1, 0, rho = r/sqrt(2)) - 2*zratio1 )
  return(res)
}
bridgeF_cb <- function(r, zratio1 = NULL, zratio2){
  # continuous and binary
  de2 <- stats::qnorm(zratio2)
  res <- as.numeric( 4*fMultivar::pnorm2d(0, de2, rho = r/sqrt(2)) - 2*zratio2 )
  return(res)
}
bridgeF_tb <- function(r, zratio1, zratio2){
  # truncated and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  mat1 <- matrix(c(1, -r, 1/sqrt(2),
                   -r, 1, -r/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 <- matrix(c(1, 0, -1/sqrt(2),
                   0, 1, -r/sqrt(2),
                   -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(
    2*(1-zratio1)*(zratio2)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2)
  )
  return(res)
}
bridgeF_bt <- function(r, zratio1, zratio2){
  # binary and truncated
  de1 <- stats::qnorm(zratio2)
  de2 <- stats::qnorm(zratio1)
  mat1 <- matrix(c(1, -r, 1/sqrt(2),
                   -r, 1, -r/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 <- matrix(c(1, 0, -1/sqrt(2),
                   0, 1, -r/sqrt(2),
                   -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(
    2*(1-zratio2)*(zratio1)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2)
  )
  return(res)
}
bridgeF_tc <- function(r, zratio1, zratio2 = NULL){
  # truncated and continuous
  de1 <- stats::qnorm(zratio1)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric( -2*fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                       4*mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2) )
  return(res)
}
bridgeF_ct <- function(r, zratio1 = NULL, zratio2){
  # continuous and truncated
  de1 <- stats::qnorm(zratio2)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric( -2*fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                       4*mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2) )
  return(res)
}
bridgeF_bb <- function(r, zratio1, zratio2){
  # binary and binary
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  res <- as.numeric(2*(fMultivar::pnorm2d(de1, de2, rho = r) - zratio1*zratio2))
  return(res)
}
bridgeF_tt <- function(r, zratio1, zratio2){
  # truncated and truncated
  de1 <- stats::qnorm(zratio1)
  de2 <- stats::qnorm(zratio2)
  
  mat1 <- matrix(c(1, 0, 1/sqrt(2), -r/sqrt(2),
                   0, 1, -r/sqrt(2), 1/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1, -r,
                   -r/sqrt(2), 1/sqrt(2), -r, 1), nrow = 4)
  mat2 <- matrix(c(1, r, 1/sqrt(2), r/sqrt(2),
                   r, 1, r/sqrt(2), 1/sqrt(2),
                   1/sqrt(2), r/sqrt(2), 1, r,
                   r/sqrt(2), 1/sqrt(2), r, 1), nrow = 4)
  
  res <- as.numeric( -2*mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat1) +
                       2*mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat2)
  )
  return(res)
}

# Calculate norm of canonical vector w. It is calculated as t(w) %*% Sigma %*% w.
# If w is all zero vector, then 0 will be returned.
w_norm <- function(W, Sigma){
  W <- as.vector(W)
  # Check first if the inputs are valid.
  if(!isSymmetric(Sigma)){
    stop("Sigma should be symmetric.")
  }
  if(length(W) != nrow(Sigma)){
    stop("non-conformable arguments. Check W and Sigma.")
  }
  wnorm <- as.numeric(crossprod(W, Sigma %*% W))
  return(wnorm)
}



# Function to normalize the canonical vector W1. The returned vector should satisfy t(W) %*% Sigma %*% W = 1.
normalizedW <- function(W, Sigma){
  W <- as.vector(W)
  # Check first if the inputs are valid.
  if(!isSymmetric(Sigma)){
    stop("Sigma should be symmetric.")
  }
  if(length(W) != nrow(Sigma)){
    stop("non-conformable arguments. Check W and Sigma.")
  }
  normW <- w_norm(W, Sigma)
  if(normW == 0){
    normalizedW <- 0
  } else {
    normalizedW <- W/sqrt(normW)
  }
  return(normalizedW)
}


# Estimate canonical correlation based on two estimates w1 and w2
# Sigma1, Sigma2, Sigma12 are all true latent correlation matrices or from out-of-sample correlation matrices.
# w1 is a vector of length p1. w2 is a vector of length p2.
# Sigma1 is a matrix of size p1 by p1. Sigma2 is a matrix of size p2 by p2. Sigma12 is a matrix of size p1 by p2.
cancorhat <- function(w1, w2, Sigma1, Sigma2, Sigma12){
  w1 <- as.vector(w1)
  w2 <- as.vector(w2)
  
  # Check first if the inputs are valid.
  if(length(w1) != ncol(Sigma1)){
    stop("non-conformable arguments: xcoef and Sigma1")
  }
  if(length(w2) != ncol(Sigma2)){
    stop("non-conformable arguments: ycoef and Sigma2")
  }
  if(ncol(Sigma1) != nrow(Sigma12)){
    stop("non-conformable arguments: Sigma1 and Sigma12")
  }
  if(ncol(Sigma2) != ncol(Sigma12)){
    stop("non-conformable arguments: Sigma2 and Sigma12")
  }
  # when one of the estimates w1 or w2 is all zero vector, the output will be zero.
  norm1 <- w_norm(w1, Sigma1)
  norm2 <- w_norm(w2, Sigma2)
  if (norm1 == 0 | norm2 == 0){
    output <- 0
  } else {
    output <- as.numeric(crossprod(w1, Sigma12 %*% w2))/sqrt(norm1)/sqrt(norm2)
  }
  return(output)
}
#' @title Estimate latent correlation matrix
#'
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#'
#' @aliases estimateR_mixed
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables.
#' @param type A type of variables in \code{X}, must be one of "continuous", "binary" or "trunc".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param rho Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @return \code{estimateR} returns
#' \itemize{
#'       \item{type: }{Type of the data matrix \code{X}}
#'       \item{R: }{Estimated p by p latent correlation matrix of \code{X}}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12168}{"High dimensional semiparametric latent graphicalmodel for mixed data"}, \emph{J. R. Statist. Soc. B}, 79: 405-421.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2018+) \href{http://arxiv.org/abs/1807.05274}{"Sparse semiparametric canonical correlation analysis for data of mixed types"}, \emph{arXiv 1807.05274}.
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/estimateR_ex.R
estimateR <- function(X, type = "trunc", use.nearPD = TRUE, rho = 0.01, tol = 1e-3){
  X <- as.matrix(X)
  
  n <- nrow(X)
  p1 <- ncol(X)
  
  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }
  
  if (type == "trunc"){
    if(sum(X < 0) > 0) {stop("The data contains negative values.")}
    if(sum(colSums(X == 0)) == 0){
      message("The data does not contain zeros. Consider changing the type to \"continuous\".")
    }
  }
  if (type == "binary"){
    if(sum(!(X %in% c(0, 1))) > 0) {stop("The data is not \"binary\".")}
  }
  if (type == "continuous"){
    R1 <- sin(pi/2 * pcaPP::cor.fk(X))
  } else {
    zratio1 <- colMeans(X == 0)
    R1 <- fromKtoR(Kendall_matrix(X), zratio = zratio1, type = type, tol = tol)
  }
  
  if ( use.nearPD == TRUE & min(eigen(R1)$values) < 0 ) {
    message(" minimum eigenvalue of correlation estimator is ", min(eigen(R1)$values), "\n nearPD is used")
    R1 <- as.matrix(Matrix::nearPD(R1, corr = TRUE)$mat)
  }
  # shrinkage method
  if(rho < 0 | rho > 1){ stop("rho must be be between 0 and 1.") }
  R1 <- (1 - rho)*R1 + rho*diag(p1)
  
  ### Want to keep the correct column names for each matrices
  colnames(R1) <- rownames(R1) <- c(colnames(X))
  
  return(list(type = type, R = R1))
}
#'
#'
#'
#' @title Estimate latent correlation matrix
#' @rdname estimateR
#' @aliases estimateR_mixed
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of variables in \code{X1}, must be one of "continuous", "binary" or "trunc".
#' @param type2 A type of variables in \code{X2}, must be one of "continuous", "binary" or "trunc".
#' @inheritParams use.nearPD
#' @inheritParams rho
#' @inheritParams tol
#'
#' @return \code{estimateR_mixed} returns
#' \itemize{
#'       \item{type1: }{Type of the data matrix \code{X1}}
#'       \item{type2: }{Type of the data matrix \code{X2}}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} = (\code{X1}, \code{X2}) (p1+p2 by p1+p2)}
#'       \item{R1: }{Estimated latent correlation matrix of \code{X1} (p1 by p1)}
#'       \item{R2: }{Estimated latent correlation matrix of \code{X2} (p2 by p2)}
#'       \item{R12: }{Estimated latent correlation matrix between \code{X1} and \code{X2} (p1 by p2)}
#' }
#'
#' @export
#' @importFrom Matrix nearPD
estimateR_mixed <- function(X1, X2, type1 = "trunc", type2 = "continuous", use.nearPD = TRUE, rho = 0.01, tol = 1e-3){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
  }
  
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  if (sum(c(type1, type2) %in% c("continuous", "binary", "trunc")) != 2){
    stop("Unrecognised type of variables. Should be one of continuous, binary or trunc.")
  }
  
  if (type1 == "trunc"){
    if(sum(X1 < 0) > 0) {stop("The data X1 contains negative values.")}
    if(sum(colSums(X1 == 0)) == 0){
      message("The data X1 does not contain zeros. Consider changing the type to \"continuous\".")
    }
  }
  if (type1 == "binary"){
    if(sum(!(X1 %in% c(0, 1))) > 0) {stop("The data X1 is not \"binary\".")}
  }
  
  if (type2 == "trunc"){
    if(sum(X2 < 0) > 0) {stop("The data X2 contains negative values.")}
    if(sum(colSums(X2 == 0)) == 0){
      message("The data X2 does not contain zeros. Consider changing the type to \"continuous\".")
    }
  }
  if (type2 == "binary"){
    if(sum(!(X2 %in% c(0, 1))) > 0) {stop("The data X2 is not \"binary\".")}
  }
  
  if (type1 == type2) {
    # both datasets are of the same type. CC, TT or BB case.
    
    Xcomb <- cbind(X1, X2)
    Rall <- estimateR(Xcomb, type = type1, use.nearPD = use.nearPD, rho = rho, tol = tol)$R
    R1 <- Rall[1:p1, 1:p1]
    R2 <- Rall[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- Rall[1:p1, (p1 + 1):(p1 + p2)]
    
  } else {
    # datasets are of different type
    if (type1 == "continuous"){ # These are CT or CB case.
      zratio2 <- colMeans(X2 == 0)
      R1 <- sin(pi/2 * pcaPP::cor.fk(X1))
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type2, tol = tol)
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    } else if (type2 == "continuous"){ # These are TC or BC case.
      zratio1 <- colMeans(X1 == 0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type1, tol = tol)
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
      R2 <- sin(pi/2 * pcaPP::cor.fk(X2))
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    } else { # These are TB or BT case.
      zratio1 <- colMeans(X1 == 0)
      zratio2 <- colMeans(X2 == 0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type1, tol = tol)
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type2, tol = tol)
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    }
    
    if ( use.nearPD == TRUE & min(eigen(Rall)$values) < 0 ) {
      message(" minimum eigenvalue of correlation estimator is ", min(eigen(Rall)$values), "\n nearPD is used")
      Rall <- as.matrix(Matrix::nearPD(Rall, corr = TRUE, posd.tol = 1e-04)$mat)
    }
    # shrinkage method
    if(rho < 0 | rho > 1){ stop("rho must be be between 0 and 1.") }
    Rall <- (1 - rho)*Rall + rho*diag(p1 + p2)
    
    ### Want to keep the correct column names for each matrices
    colnames(Rall) <- rownames(Rall) <- c(colnames(X1), colnames(X2))
    
    # For conveninence, split the R matrices
    R1 <- Rall[1:p1, 1:p1]
    R2 <- Rall[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- Rall[1:p1, (p1 + 1):(p1 + p2)]
  }
  
  return(list(type = c(type1, type2), R1 = R1, R2 = R2, R12 = R12, R = Rall))
}
fromKtoR <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  de1 <- NULL
  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else {
    if (is.null(zratio)){ stop ("The input for zratio is required for \"trunc\" and \"binary\" types.") }
    if (length(zratio)!=d1){ stop ("The length of zratio must match with the number of columns in K.") }
    bridge <- bridge_select(type1 = type, type2 = type)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d1)
    
    if (d1 <= 1){ # for scalar K: one Kendall tau value
      i <- j <- 1
      f1 <- function(r)(bridge(r, zratio = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      if(op == 100) {
        hatR[i, j] <- hatR[j, i] <-0
      } else {
        hatR[i, j] <- hatR[j, i] <- unlist(op)
      }
    } else { # for matrix K
      for(i in 1:(d1-1)) {
        for(j in (i+1):d1){
          # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
          f1 <- function(r)(bridge(r, zratio1 = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            hatR[i, j] <- hatR[j, i] <-0
          }else {
            hatR[i, j] <- hatR[j, i] <- unlist(op)
          }
        }
      }
    }
  }
  return(hatR)
}

fromKtoR_mixed <- function(K12, zratio1 = NULL, zratio2 = NULL, type1 = "trunc", type2 = "continuous", tol = 1e-3) {
  
  K12 <- as.matrix(K12)
  d1 <- nrow(K12)
  d2 <- ncol(K12)
  ###################################################################
  de1 <- de2 <- NULL
  
  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    if (type1 != "continuous"){
      if (is.null(zratio1)){ stop ("The input for zratio1 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio1)!=d1){ stop ("The length of zratio1 must match with the number of rows in K12.") }
    }
    if (type2 != "continuous"){
      if (is.null(zratio2)){ stop ("The input for zratio2 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio2)!=d2){ stop ("The length of zratio2 must match with the number of columns in K12.") }
    }
    
    bridge <- bridge_select(type1 = type1, type2 = type2)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d2)
    
    if ( d1 <= 1 & d2 <= 1 ){
      i <- j <- 1
      f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      if(op == 100) {
        hatR[i, j] <- hatR[j, i] <- 0
      } else {
        hatR[i, j] <- hatR[j, i] <- unlist(op)
      }
    } else {
      for(i in 1:d1) {
        for(j in 1:d2){
          # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
          f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            hatR[i,j] <- 0
          } else {
            hatR[i,j] <- unlist(op)
          }
        }
      }
    }
  }
  return(hatR)
}

#' Kendall's tau correlation
#'
#' Calculate Kendall's tau correlation.
#' \deqn{ \hat{\tau}_{jk} = \frac{2}{n(n-1)}\sum_{1\le i<i'\le n} sign(X_{ji}-X_{ji'}) sign(X_{ki}-X_{ki'}) }
#' The function \code{KendallTau} calculates Kendall's tau correlation between two variables, returning a single correlation value. The function \code{Kendall_matrix} returns a correlation matrix (p1 by p2).
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @rdname KendallTau
#' @examples
#'
#' n <- 100 # sample size
#' r <- 0.8 # true correlation
#'
#' ### vector input
#' # Data generation (X1: truncated continuous, X2: continuous)
#' Z <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, r, r, 1), nrow = 2))
#' X1 <- Z[,1]
#' X1[Z[,1] < 1] <- 0
#' X2 <- Z[,2]
#'
#' KendallTau(X1, X2)
#' Kendall_matrix(X1, X2)
#'
#' ### matrix data input
#' p1 <- 3; p2 <- 4 # dimension of X1 and X2
#' JSigma <- matrix(r, nrow = p1+p2, ncol = p1+p2); diag(JSigma) <- 1
#' Z <- mvrnorm(n, mu = rep(0, p1+p2), Sigma = JSigma)
#' X1 <- Z[,1:p1]
#' X1[Z[,1:p1] < 0] <- 0
#' X2 <- Z[,(p1+1):(p1+p2)]
#'
#' Kendall_matrix(X1, X2)
#'
#' @importFrom pcaPP cor.fk
#' @export
KendallTau <- function(x, y){ # both x and y are vectors, not matrix.
  # Based on cor.fk function from pcaPP package to make the computation faster.
  # It can handle ties.
  if (length(x) != length(y)){ # Check of they have the same length.
    stop ("x and y must have same length.")
  }
  n <- length(x)
  n0 <- n*(n-1)/2
  if (length(unique(x)) != n) {
    x <- as.vector(x) # sometimes input x is a matrix n by 1, which gives errors for rle function below.
    x.info <- rle(sort(x))
    t1 <- x.info$lengths[x.info$lengths>1]
    n1 <- sum(t1*(t1-1)/2)
  } else {
    n1 <- 0
  }
  if (length(unique(y)) != n) {
    y <- as.vector(y) # sometimes input y is a matrix n by 1, which gives errors for rle function below.
    y.info <- rle(sort(y))
    u1 <- y.info$lengths[y.info$lengths>1]
    n2 <- sum(u1*(u1-1)/2)
  } else {
    n2 <- 0
  }
  
  tau <- pcaPP::cor.fk(x, y)*sqrt(n0-n1)*sqrt(n0-n2)/n0
  
  return(tau)
}

#' @param X A numeric matrix (n by p1).
#' @param Y A numeric matrix (n by p2).
#'
#' @rdname KendallTau
#' @export
#' @importFrom pcaPP cor.fk
Kendall_matrix <- function(X, Y = NULL){ # X and Y are matrix.
  if (is.null(Y)) {
    X <- as.matrix(X) # In case that X is vector, this line gives you one column matrix.
    n <- nrow(X)
    p1 <- ncol(X)
    if (p1 <= 1){
      tau <- KendallTau(X, X)
    } else {
      tau <- matrix(1, p1, p1)
      for (i in 1:(p1-1)){
        for (j in (i+1):p1){ # calculate KendallTau elementwise.
          tau[i, j] <- tau[j, i] <- KendallTau(X[, i], X[, j])
        }
      }
    }
  } else {
    X <- as.matrix(X); Y <- as.matrix(Y) # In case that X and Y are vectors, this line gives you one column matrix.
    n <- nrow(X)
    p1 <- ncol(X); p2 <- ncol(Y)
    if ( p1 <= 1 & p2 <= 1){
      tau <- KendallTau(X, Y)
    } else {
      tau <- matrix(0, p1, p2)
      for (i in 1:p1){
        for (j in 1:p2){ # calculate KendallTau elementwise.
          tau[i, j] <- KendallTau(X[, i], Y[, j])
        }
      }
    }
  }
  return(tau)
}
