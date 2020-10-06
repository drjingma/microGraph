#######################
##
##
##  Data Examples for book chapter
##
##
########################

######
## updated on 9/24/2020: 
## 1. adding generating partial correlation network from AR(1)/chain graph; use rho = 0.5 (small) and rho=0.8 (large) cases, use small n and large p
## 2. try to change the mean parameters under null2 and alt2, to target overall mean sparsity level 65%-70%
 
# filepath = 'C:\\Users\\yuek\\Dropbox\\Microbial_Networks\\microGraph' #BOX
# 
# filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph' # MAC
# 
# filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph' #PC

filepath = '~/Desktop/micro_net' #bayes


setwd(filepath)
source('lib/func_libs.R')
source('Kun_code/generation_function - Copy.R')
source('Kun_code/method_function - Copy.R')
library(SpiecEasi)
library(MASS)

data(amgut1.filt)
reference_data = amgut1.filt
# dim(reference_data) #data format n by p
# n = c(100, 200, 500)[1]
# p = 200

args = commandArgs(trailingOnly = T) # (n p null1 2),  n and p override by reference data set, e.g. 
print(args)
# args = c(50, 127, 'alt1', 20, 'zinegbin', 'chain_small', 0)
# args = c(50, 40, 'alt2', 20, 'none', 'chain_small', 0)
# args = c(20, 30, 'null2', 2)
# args = c(100, 127, 'null1.1', 100)
n = as.integer(args[1])
p = as.integer(args[2])
choose_model = as.character(args[3]) # choose among 'null1' (shuffle reference data), 
                                    # 'null1.1' (use copula model with NB marginal distr, based on reference data), 
                                    # 'null2' (generate from Dirichelt distribution with varying total count), 
                                    # 'alt1' (generate with copula model and reference data), 
                                    # 'alt2' (generate from log-normal distribution directly)
                                    # 'alt3' (generate from logisitic normal multinomial with fixed total count)
nreps = as.integer(args[4]) # just run 200 for all methods
distr = as.character(args[5]) # specify the distribution for copula model, if used 
network_option = as.character(args[6]) # 'erdos_renyi', or 'chain_small', 'chain_large' for AR(1) (small: rho=0.5, large:rho=0.8), 'cov_erdos_renyi' for correlation based graph
network_condition_number = as.numeric(args[7]) # specify the condition number of the inverse covariance matrix

set.seed(102)

#---------------
# part I: null hypothesis
#---------------

if(choose_model == 'null1'){
  ## Null 1
  # shuffle reference data set, empty graph for both cov and inv-cov
  reference_data = amgut1.filt[sample(1:nrow(amgut1.filt), size=n, replace = F), ]
  p = ncol(reference_data)
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='shuffle', reference_data = reference_data, Sigma_list = Sigma_list)
}

if(choose_model == 'null1.1'){
  ## Null1 1.1
  # Use reference data set, marginal distributions are ZINB, and the network is set to empty
  reference_data = amgut1.filt[sample(1:nrow(amgut1.filt), size=n, replace = F), ]
  p = ncol(reference_data)
  # generate a null graph object
  option = list(hypothesis = 'null', model = 'copula', 
                reference_data = reference_data, Sigma_list = NULL, distr = distr)
  option$Sigma_list$A_inv <-  option$Sigma_list$A_cov <- matrix(0, p, p)
  option$Sigma_list$Sigma <- option$Sigma_list$Omega <- diag(1, p)
}

if(choose_model == 'null2'){
  ## Null 2
  # generate from Dirichlet distribution, so will assume both inv and cov based graph being empty?
  mu = runif(p, 0, 8) # the mu is for underlying log-normal of Dirichlet parameters
  alpha=100
  library_scale = rnegbin(n,mu=3e4, theta=9e8/(3e6-3e4)) # use library_scale to directly specify the total counts per sample
  # library_scale = rep(3e4, n)
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='Dirichlet', library_scale = library_scale, alpha = alpha, mu=mu, Sigma_list = Sigma_list, distr = distr)
}



#---------------
# part II: Alternative hypothesis
#---------------


if(choose_model=='alt1'){
  ## Alternative 1
  # use Copula generative model, graph is based on inverse covariance matrix.
  reference_data = amgut1.filt[1:n, ]
  p = ncol(reference_data)
  Sigma_list = SpiecEasi_graph_Sigma(p,e=p, type = network_option, graph=NULL, network_condition_number = network_condition_number)
  option = list(hypothesis = 'alternative', model = 'copula', reference_data = reference_data, 
                Sigma_list = Sigma_list, distr = distr, network_option = network_option, network_condition_number = network_condition_number)
}


if(choose_model == 'alt2'){
  ## Alternative 2: log normal
  # generate graph, p nodes and e edges; 
  # this graph corresponds to inverse covariance
  Sigma_list  = SpiecEasi_graph_Sigma(p,e=p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  mu = runif(p, 0, 4)
  option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma_list = Sigma_list, 
                distr = distr, network_option = network_option, network_condition_number = network_condition_number)
  
}

if(choose_model == 'alt3'){
  ## Alternative 3: logistic normal multinomial
  # generate graph, p nodes and e edges; 
  # this graph corresponds to inverse covariance
  Sigma_list  = SpiecEasi_graph_Sigma(p,e=p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  mu = runif(p, 0, 4)
  library_scale = rep(3e4,n)
  option = list(hypothesis = 'alternative', model = 'multinomial-log-normal', 
                mu = mu, Sigma_list = Sigma_list, distr = distr, 
                library_scale = library_scale, network_option = network_option, network_condition_number = network_condition_number)
  
}
# data_list = data_generation(n, p, option)
# dim(data_list$data)
# 
# 
# data = data_list$data
# colnames(data) <- paste('taxa', 1:ncol(data))
# rownames(data) <- paste('cell', 1:nrow(data))


## wrap nrep replications for running the methods; data_rep[[1,k]] is composition matrix; data_rep[[2,k]] is count matrix
## the data matrix is not zero-corrected. 

data_rep = matrix(list(), 2, nreps)
for(k in 1:nreps){
  set.seed(k*100)
  data_list = data_generation(n, p, option)
  dim(data_list$data)
  
  data = data_list$data
  colnames(data) <- paste('taxa', 1:ncol(data))
  rownames(data) <- paste('cell', 1:nrow(data))
  
  data_rep[[1, k]] = count_to_comp(data) # compositions
  data_rep[[2, k]] = data # counts
}

save(data_rep, file = 
       paste0('dist_data/', distr,'/', network_option, '/cond_', network_condition_number, '/n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))
save(list=c('data_rep', 'option', 'n', 'p'), file = 
       paste0('dist_data/', distr,'/', network_option, '/cond_', network_condition_number,  '/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))






