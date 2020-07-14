require(igraph)
require(tmvtnorm)
require(gdata)
require(RANN)
require(ppcor)
require(glasso)
require(huge)

# ----- Auxiliary functions -----

pdParCor <- function(A, zeta, eps=1e-04){
  p <- nrow(A)
  Ip <- diag(1, p)
  B <- Ip - A
  while (min(eigen(B)$val)< eps){
    B <- B + diag(zeta, p)
  }
  A <- Ip - cov2cor(B)
  return(A)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

gm_mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

#' The geometric mean for positive components only
gm_mean_pos <- function(x){
  y <- x[x > 0 & !is.na(x)]
  if (length(y)==0){
    stop("The input vector is all zero!")
  } else {
    exp(sum(log(y)) / length(y))
  }
}

clr <- function(x, base=exp(1)){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

#' x is a matrix
clr_epsilon <- function(mat, base=exp(1), const=0.1){
  index <- which(mat>0 & !is.na(mat))
  a <- apply(mat,1,function(x){
    y <- x[x > 0 & !is.na(x)]
    y <- log(y / gm_mean_pos(y), base)
    x[x > 0 & !is.na(x)] <- y
    x
  })
  a <- t(a)
  mat[index] <- a[index] + abs(min(a)) + const
  
  return(mat)
}

#' Generate a random graph with p nodes
#' For more details, check the function sample_gnm in the igraph package.
model_gnm <- function(p, m, const=0.01){
  g <- sample_gnm(n=p,m=m,directed = FALSE)
  A <- as.matrix(as_adjacency_matrix(g,type = "both"))
  weights <- matrix(0, p, p)
  upperTriangle(weights, diag = F) <- runif((p*(p - 1))/2, 0.5, 1)*(2*rbinom((p*(p - 1))/2, 1, 0.5) - 1)
  weights <- weights + t(weights)
  Omega <- A*weights
  diag(Omega) <- ifelse(min(eigen(Omega)$val) > 0, 0, - min(eigen(Omega)$val)) + const
  Sigma <- solve(Omega)
  Sigma <- cov2cor(Sigma)
  Omega <- solve(Sigma)
  class(A) <- "graph"
  return(list(Omega = Omega, Sigma = Sigma, Adj=A))
}

#' @details This function is designed for simulating from log-normal or multivariate Gamma distribution. 
#' To simulate count data, we need a different function. 
#' @param pars model parameters 
data.simulator <- function(pars, sid=NULL) {
  p <- pars$p

  if (!is.null(sid)){
    set.seed(sid)
  }
  
  model <- pars$model

  if (identical(model,"gamma")){
    S.svd <- svd(pars$S)
    FM <- S.svd$u %*% diag(sqrt(S.svd$d))
    U <- matrix(rgamma(pars$n * (pars$p), shape=pars$shape, scale = pars$scale), ncol=pars$n)
    Y <- t(matrix(rep(pars$mu,pars$n),ncol=pars$n) + FM %*% U/sqrt(pars$shape))
  } else if (identical(model,"normal")){
    # Simulate from MVN; with some contamination
    Y <- mvrnorm(n=pars$n, mu=pars$mu, Sigma=pars$S) 
  }
  
  # Get the basis and compositions
  basis <- exp(Y)
  theta <- sweep(basis,1,rowSums(basis),FUN='/')

  return(list(basis = basis, composition=theta))
}

#' @details This function simulates zero-inflated negative binomial count data from a given data set. 
#' @param pars model parameters 
data.simulator.zinegbin <- function(pars, sid=NULL, pseudo = 1) {
  if (!is.null(sid)){
    set.seed(sid)
  }
  
  Y <- SpiecEasi::synth_comm_from_counts(pars$dat, mar=2, distr='zinegbin',Sigma = pars$S)
  # Apply CLR transformation to observed counts
  obs_clr <- t(apply(Y + pseudo, 1, clr)) # as done in SPIEC-EASI
  obs_mclr <- clr_epsilon(Y)
  return(list(clr = obs_clr, mclr = obs_mclr,censCounts=Y,
              percent.zero=sum(Y==0)/prod(dim(Y))))
  
}

#' @details Function to compute the trace of a matrix
matTr <- function(z) sum(diag(z))

#' @details Function to compute the nearest positive definite matrix
nearPDMat <- function(A){
  # replace -ve eigen values with small +ve number
  newEig <- eigen(A)
  newEig2 <- ifelse(newEig$values < 0, 0, newEig$values)
  
  # create modified matrix eqn 5 from Brissette et al 2007, inv = transp foreig vectors
  A <- newEig$vectors %*% diag(newEig2) %*% t(newEig$vectors)
  
  # normalize modified matrix eqn 6 from Brissette et al 2007
  A <- A/sqrt(diag(A) %*% t(diag(A)))
  return(A)
}


#---- Regularization paths ----
#' @details This function returns the path of the graphical lasso estimator. 
#' @param X The n by p data matrix
#' @param lambda A sequence of tuning parameters
gm.glasso <- function(X,lambda, use.nearPD=TRUE){
  lambda.sorted <- sort(lambda,decreasing = TRUE)
  X <- scale(X, center = T, scale = T)
  n <- nrow(X)
  empcov <- (1/n) * t(X) %*% X
  if ( use.nearPD == TRUE & min(eigen(empcov)$values) < 0 ) {
    message(" minimum eigenvalue of correlation estimator is ", min(eigen(empcov)$values), "\n nearPD is used")
    empcov <- as.matrix(Matrix::nearPD(empcov, corr = TRUE, posd.tol=1e-03)$mat)
  }
  
  pathList <- lapply(lambda.sorted, function(m) glasso(empcov, rho=m, penalize.diagonal = FALSE, approx=FALSE)$wi)
  pathList <- lapply(pathList, function(a) (a+t(a))/2)
  adjList <- lapply(pathList, function(a) 1*(abs(a)>1e-08) - diag(1,ncol(X)))
  return(list(path=adjList, icov=pathList, cov=empcov, lambda=lambda.sorted,data=X, method='glasso'))
}

#' @details This function returns the path of the SPRING estimator. 
gm.spring <- function(X,lambda,method=c('clime','glasso'),approx=FALSE){
  lambda.sorted <- sort(lambda,decreasing = TRUE)
  
  method <- match.arg(method)
  hat_sigma <- spring.corr(X)
  
  if (method=='clime'){
    ## Fit clime
    pathList <- clime(hat_sigma, sigma=TRUE, lambda=lambda.sorted)
  }
  if (method=='glasso'){
    ## Fit the graphical lasso
    pathList <- lapply(lambda.sorted, function(m) glasso(hat_sigma, rho=m, penalize.diagonal = FALSE, approx=approx)$wi)
  }
  
  pathList <- lapply(pathList, function(a) (a+t(a))/2)
  adjList <- lapply(pathList, function(a) 1*(abs(a)>1e-08) - diag(1,ncol(X)))
  return(list(path=adjList, icov=pathList,cov=hat_sigma,lambda=lambda.sorted,data=X, method='spring'))
}

#' @details This function returns the path of the gCoDa estimator. 
gm.gcoda <- function(X, lambda,counts=T){
  lambda.sorted <- sort(lambda,decreasing = TRUE)
  res_gcoda_count <- gcoda(x = X, counts = counts, lambda=lambda.sorted);
  return(list(path=res_gcoda_count$path,icov=res_gcoda_count$icov,lambda=res_gcoda_count$lambda, data=X, method='gcoda'))
}

#' Function to compare the AUC between SpiecEasi and ProbitDirect
#' @param dd A list with count and CLR-transformed data. The first p columns of the count data are microbiome counts
#' @param nlambda Number of tuning parameters.
#' @param lambda The vector of tuning parameters. Either nlambda or lambda needs to be specified.
#' @param graph The conditional dependence graph

compare.roc <- function(dd, lambda, graph, OTU.only=TRUE){
  ncond <- 6
  out <- vector("list",ncond)
  paths <- vector("list",ncond)
  names(out) <- c('gcoda', "oracle", "spiec", "spring")
  names(paths) <- names(out)
  paths$gcoda <- NULL
  out$gcoda <- NULL
  if (OTU.only){
    paths$gcoda <- gcoda(x = dd$noisyCounts, counts = T, lambda=lambda)$path;
    out$gcoda <- huge.roc(paths$gcoda, graph, verbose = FALSE)
  }
  
  paths$oracle  <- gm.glasso(log(dd$basis), lambda)$path
  paths$spiec <- gm.glasso(dd$clr, lambda)$path
  paths$spring <- gm.spring(dd$mclr, lambda, method='glasso')$path
  for (nn in names(paths)){
    # cat('...',nn,'...\n')
    out[[nn]] <- huge.roc(paths[[nn]], graph, verbose = FALSE)
  }
  return(out)
}
