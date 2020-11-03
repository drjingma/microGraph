##-------------------------------------------------------------
## collection of functions for generating simulation data set
##
## Kun Yue, yuek@uw.edu
## Date: 2020/11/02
## 
## (modify the file paths in this script to use on your computer)
##-------------------------------------------------------------

library(entropy)
library(gtools)
library(SpiecEasi)
library(seqgroup)





##---------------------------------------
## Null model data set generation
##---------------------------------------


#---------------
# option 1: generate each OUT independently from marginal distributions, 
# with parameter values referencing some empirical data
#---------------

## directly use copula model generation function with null graph matrix

#---------------
# option 1.1: shuffle a reference data set to obtain null model data
#---------------

null_data_generate_1.1 = function(x){ # x is n by p real data matrix
  # shuffle within taxa, so keep its marginal distribution
  n = nrow(x)
  p = ncol(x)
  cat('Shuffle reference dataset, nsub ', n, ', ntaxa ', p, '\n')
  index = replicate(p, sample(n, replace=F))
  x_shuffle = sapply(1:p, function(i)x[index[,i], i])
  return(x_shuffle)
}


#---------------
# option 2: generate composition from Dirichlet distribution, 
# then generate counts from multinomial. Equivalent to independent Gamma distributions
#---------------

null_data_generate_2 = function(n, p, library_scale, alpha, mu){ 
  # library_scale is total count
  # generate Dirichelt(alpha*X_i) with a fixed alpha; where X is compositional from log-normal(mu, sd=1.5)
  
  Phi = mvtnorm::rmvnorm(n,mean=mu,sigma = diag(1.5,p)) 
  X = sweep(exp(Phi),1,STATS = rowSums(exp(Phi)), FUN='/') 
  W <- t(sapply(1:n, function(i) {
    a <- dirmult::rdirichlet(1,alpha*X[i,])
    rmultinom(n=1,size=library_scale[i],prob=a) # generate counts from multinomial-dirichlet
  }))
  
  return(W)
}

#---------------
# option 3: based on option 2 model, subtract all entries by the data matrix mean, 
# set negative values to zero -> increase sparsity level 
#---------------
null_data_generate_3 = function(n, p, library_scale){
  data = null_data_generate_2(n, p, library_scale)
  data = data-mean(data)
  data[data<0]<-0
  return(data)
}


##---------------------------------------
## Informative network model data set generation
##---------------------------------------

#-----------------
# option 1: logistic normal multinomial:
#-----------------
para_data_generate_1 = function(n, p, library_scale){
  mu_vec = rnorm(p-1, sd=1) 
  Sigma_mat_decomp = matrix(rnorm((p-1)^2, sd=0.1), p-1, p-1)
  Sigma_mat = Sigma_mat_decomp%*%t(Sigma_mat_decomp)
  z = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*(p-1), mean = 0), n, p-1, byrow=T)) + mu_vec
  z = cbind(t(z), rep(0, n))
  x = exp(z)/rowSums(exp(z))
  data = t(sapply(1:n, function(i)rmultinom(n=1, size=rpois(n=1,library_scale), prob=x[i,])))
  return(data)
}

#-----------------
# option 2: log normal distribution for generating abundance data; takes in mu and Sigma, output n by p
#-----------------
para_data_generate_2 = function(n, p, mu, Sigma){
  Sigma_mat_decomp = t(chol(Sigma))
  # Sigma_mat = Sigma_mat_decomp%*%t(Sigma_mat_decomp)
  x = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*p, mean = 0), n, p, byrow=T)) + mu
  y = t(exp(x))
  data=y  
  return(data)
}

#-----------------
# option3: multinomial log normal
#-----------------
para_data_generate_3 = function(n, p, mu, Sigma, library_scale){
  Sigma_mat_decomp = t(chol(Sigma))
  # Sigma_mat = Sigma_mat_decomp%*%t(Sigma_mat_decomp)
  x = Sigma_mat_decomp %*% t(matrix(rnorm(n=n*p, mean = 0), n, p, byrow=T)) + mu
  y = t(exp(x))
  comp_mat = sweep(y,1,STATS = rowSums(y), FUN='/')
  data = t(do.call(cbind,lapply(1:n, function(i)rmultinom(n=1,size=library_scale[i], prob = comp_mat[i,]))))
  return(data)
}

#----------------
# option 4: Copula model 
# generate OTU values with specified marginal distributions; 
# distribution parameters estimated based on reference data set; 
# relationship among OTUs based on specified graph matrix
#----------------

setwd('..../SpiecEasi-master/SpiecEasi-master/R')
import_files = list.files()
sapply(import_files, source)


SpiecEasi_graph_Sigma = function(p,            # number of OTUs
                                 e,            # number of edges in the graph (if used)
                                 type,         # c('band', 'cluster', 'erdos_renyi', 'hub', 'scale_free', 'block', 'chain_small', 'chain_large', 'cov_erdos_renyi')
                                 graph = NULL, # if specify this graph matrix, will override other arguments
                                 network_condition_number # the targeted precision matrix condition number
                                 ){
  
  if(!is.null(graph)){
    Prec  <- graph2prec(graph, posThetaLims = c(2, 3), # this is precision strength range
                        targetCondition = network_condition_number, # this is condition number kappa
                        epsBin = 0.01, numBinSearch = 100) # here can set different generating parameters
    Omega = diag(1/sqrt(diag(Prec))) %*% Prec %*% diag(1/sqrt(diag(Prec)))
    Cor = cov2cor(prec2cov(Prec))
    A_cov = (abs(Cor)>1e-11)*1
    diag(A_cov) = 0
  }
  
  if(is.null(graph)){
    if(type %in% c('band', 'cluster', 'erdos_renyi', 'hub', 'scale_free', 'block')){
      d <- p
      graph <- make_graph(type, d, e) # choose from band, cluster, erdos_renyi, hub, scale_free, block; can also write my own graph (essentially an adjacency matrix)
      Prec  <- graph2prec(graph, posThetaLims = c(2, 3), # this is precision strength range
                          targetCondition = network_condition_number, # this is condition number kappa
                          epsBin = 0.01, numBinSearch = 100) # here can set different generating parameters
      Omega = diag(1/sqrt(diag(Prec))) %*% Prec %*% diag(1/sqrt(diag(Prec)))
      Cor = cov2cor(prec2cov(Prec))
      A_cov = (abs(Cor)>1e-11)*1
      diag(A_cov) = 0
      
    }else if(type %in% c('chain_small', 'chain_large')){ # chain graph for correlation
      library(stats)
      
      if(type=='chain_small'){
        # create a toeplitz matrix
        Cor<- toeplitz(0.5^(0:(p-1)))
      }else{
        Cor<-  toeplitz(0.8^(0:(p-1)))
      }
      
      Omega = solve(Cor)
      graph = 1*(abs(Omega)>1e-11)
      A_cov = 1*(abs(Cor)>1e-11)
      diag(graph)<-diag(A_cov)<-0
      
    }else if(type %in% c('cov_erdos_renyi')){ # erdos_renyi graph for correlation 
      d <- p
      graph <- make_graph('erdos_renyi', d, e) # choose from band, cluster, erdos_renyi, hub, scale_free, block; can also write my own graph (essentially an adjacency matrix)
      Cov  <- graph2prec(graph, posThetaLims = c(2, 3), # this is precision strength range
                         targetCondition = network_condition_number, # this is condition number kappa
                         epsBin = 0.01, numBinSearch = 100) # here can set different generating parameters
      Cor = diag(1/sqrt(diag(Cov))) %*% Cov %*% diag(1/sqrt(diag(Cov)))  
      Omega = solve(Cor)
      A_cov = (abs(Cor)>1e-11)*1
      diag(A_cov) = 0
      graph = (abs(Omega)>1e-11)*1
      diag(graph) <-0
      
    }else{
      stop('error: type of graph do not match pool of available graphs')
    }
    
  }
  
  
  return(list(Sigma = Cor, A_inv = graph, Omega = Omega, A_cov = A_cov))
  
}


SpiecEasi_generate = function(data,         # reference data set
                              graph_Sigma,  # list(Sigma) that contains the Covariance corresponding to network specified
                              distr,        # marginal distribution
                              n             # desired sample size
                              ){
  depths <- rowSums(data) # raw counts data, rows are obs, n by p
  data.n  <- t(apply(data, 1, norm_to_total))
  data.cs <- round(data.n * min(depths))
  
  d <- ncol(data.cs)
  e <- d
  
  
  X <- synth_comm_from_counts(data.cs, mar=2, 
                              distr=distr, # choose from zipois, zinegbin, pois, lognorm, negbin
                              Sigma=graph_Sigma$Sigma, n=n) # output data n by p
  
  return(list(data=X, Sigma_list = graph_Sigma))
}

#------------------
# wrap up generation functions for simulation study
#------------------

data_generation = function(n, p, option){
  
  model_data = NA
  other_data = list()
  
  if(option$hypothesis == 'null'){
    
    if(option$model =='shuffle'){
      #input real data set of n by p, shuffle real data
      model_data = null_data_generate_1.1(option$reference_data) 
      n = nrow(model_data)
      p = ncol(model_data)
      
    }
    
    if(option$model == 'copula'){
      data_list = SpiecEasi_generate(data = option$reference_data, 
                                     graph_Sigma = option$Sigma_list,
                                     distr = option$distr,
                                     n=n)
      model_data = data_list$data
      n = nrow(model_data)
      p = ncol(model_data)
    }
    
    if(option$model == 'Dirichlet'){
      model_data = null_data_generate_2(n, p, library_scale = option$library_scale, alpha = option$alpha, mu=option$mu) 
    }
    
  }
  
  if(option$hypothesis == 'alternative'){
    if(option$model == 'copula'){
      # generate with SpiecEasi, with a chosen graph type and zero-inflated negative binomial marginal distribution
      data_list = SpiecEasi_generate(data = option$reference_data, 
                                     graph_Sigma = option$Sigma_list,
                                     distr = option$distr,
                                     n=n)
      model_data = data_list$data
      
    }
    if(option$model == 'log-normal'){
      model_data = para_data_generate_2(n, p, mu = option$mu, Sigma = option$Sigma_list$Sigma)
    }
    
    if(option$model == 'multinomial-log-normal'){
      model_data = para_data_generate_3(n, p, mu = option$mu, Sigma = option$Sigma_list$Sigma, 
                                        library_scale = option$library_scale)
    }
  }
  
  
  # check proportion of zeros
  W_zeros <- range(apply(model_data,1,function(a) sum(a==0)/p))
  cat('range of zero proportions across cells: ', W_zeros, '\n')
  
  return(list(data = model_data))
}






##-----------------
## Extra generation functions from literature review (not used)
##-----------------

#-----------------
# Conet naive generative model
#-----------------
CoNet_data_generate = function(n, p, library_scale){
  noise = sample(c(-1, 1), replace=T, size=n*p)*rpois(n*p, library_scale/p/20)
  t1 = round(seq(library_scale/4/3/3, library_scale/4/3*3, length=n))
  t2 = rep(library_scale/4/2, n)
  t3 = round(seq(library_scale/4/3*3, library_scale/4/3/3, length=n))
  t4 = rep(library_scale/4/2,n)
  t_rest = round(library_scale/p/(p-4))
  data = cbind(t1, t2, t3, t4, matrix(t_rest, nrow=n, ncol=p-4)) + matrix(noise, nrow=n, ncol=p)
  data[data<0] = 0
  return(data)
}


#-----------------
# SparCC generative model 
#-----------------
library(entropy)

SparCC_data_generate = function(n, p, mu0, mu1, prob){
  mu = c(mu1, rep(mu0, p-1))
  
  Sigma = matrix(0, p, p)
  rand_corr = as.numeric(runif(p*(p-1)/2)<prob)
  rand_corr[rand_corr!=0] = sample(c(-1, 1), replace=T, size=length(rand_corr[rand_corr!=0]))*rand_corr[rand_corr!=0]
  Sigma[lower.tri(Sigma, diag=F)] <- rand_corr*0.1*0.1 # the sd of each taxon is set to be 0.1
  Sigma = Sigma+ t(Sigma)
  diag(Sigma)<-0.01
  # make sure it is psd
  tmp = eigen(Sigma)
  eigenvalue = tmp$values
  eigenvalue[eigenvalue<=0]<-.Machine$double.eps
  Sigma = tmp$vector %*% diag(eigenvalue)%*%t(tmp$vectors)
  corr_mat = diag(1/sqrt(diag(Sigma)))%*%Sigma%*%diag(1/sqrt(diag(Sigma)))
  
  data = exp(tmp$vector%*%diag(sqrt(eigenvalue))%*%matrix(rnorm(n*p),nrow=p, ncol=n) + mu)
  n_eff = mean(sapply(1:n, function(i)entropy.ChaoShen(data[,i])))
  return(list(n_eff=n_eff, prob=prob, mu=mu, Sigma=Sigma, corr_mat = corr_mat, data=data))
}


#-----------------------
# CClasso generative model (parametric, log-normal for abundance)
#-----------------------
CClasso_mu = function(n,p){
  runif(n=p, -0.5, 0.5)
}
# option: random, neighbor, AR, hub, block
CClasso_Sigma = function(n, p, option){
  Sigma = matrix(0, p, p)
  if(option=='random'){
    Sigma[lower.tri(Sigma, diag=F)] <- sample(c(-0.15, 0.15), replace=T, size=p*(p-1)/2) * as.numeric(runif(p*(p-1)/2)<0.3)
    Sigma = Sigma+t(Sigma) + diag(1, p)
  }
  
  if(option=='neighbor'){
    library(FNN)
    pos = matrix(runif(p*2), p, 2)
    for (i in 1:p){
      neighbor=((1:p)[-i])[knnx.index(pos[-i,], pos[i,,drop=F], k=10)] # gives index of 10 nearest neighbours
      Sigma[i, neighbor] = 0.5
    }
    
    Sigma = Sigma + t(Sigma) + diag(1, p)
  }
  
  if(option == 'AR'){
    library(stats)
    # create a toeplitz matrix
    Sigma = toeplitz(c(1, 0.4, 0.2, 0.2, 0.1, rep(0, p-5)))
  }
  
  if(option == 'hub'){
    hub = sample(1:p, size=3, replace = F)
    Sigma[hub, 1:p] = 0.2 * as.numeric(runif(3*p)<0.7)
    Sigma[(1:p)[-hub], 1:p] = 0.2 * as.numeric(runif((p-3)*p)<0.2)
    Sigma[lower.tri(Sigma, diag=T)] <-0
    Sigma = Sigma + t(Sigma) + diag(1, p)
  }
  
  if(option=='block'){ # assume we just have p that can be divided by 5
    groupsize = p/5
    flds = lapply(1:5, function(i)((i-1)*groupsize+1):(i*groupsize))
    for(i in 1:5){
      Sigma[as.matrix(expand.grid(flds[[i]], flds[[i]]))] <- 0.4*as.numeric(runif(length(flds[[i]])^2)<0.6)
      Sigma[as.matrix(expand.grid(flds[[i]], (1:p)[-flds[[i]]] ))] <- 0.2*as.numeric(runif(length(flds[[i]])*(p-length(flds[[i]]))   )<0.2)
    }
    Sigma[lower.tri(Sigma, diag=T)] <-0
    Sigma = Sigma + t(Sigma) + diag(1, p)
  }
  
  # make Sigma psd
  if(!all(eigen(Sigma)$value>0))
    stop('error, not PSD, not sure how to correct it for now')
  
  return(Sigma)
}


#-----------------------
# COAT generative model (parametric, log-normal / log-gamma)
#-----------------------
source('../COAT-master/COAT-master/simulation.R') # this contains all different data generating models

# data_list = generateDataCell(n, modelCov, p1 = 50, p2 = 100, p3 = 200, nRep = 100) #modelCov = 1: hub cov; modelCov = 2: block cov; modelCov = 3: sparse cov











