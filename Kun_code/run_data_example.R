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
 
## update on 10/10/2020:
## 1. change the reference dataset as from reference_data.Rdata. Can have p=200 and n=100, 200, 500 for all settings
## 2. use #edge=3p for erdos_renyi graphs, under alt1 and alt2
## 3. change null2 library_scale (previously based on NB distribution) to be based on real data: sample from 1*min(seq_depth) to 10*min(seq_depth)
## 4. change the mean parameters under null2 with the new library_scale generation, still target 65% sparsity
## 5. change the penalize lambda sequence, with lambda.min.ratio=1e-5 (may help extend the ROC towards (1,1))

## update on 10/15/2020:
## 1. for only null1.1, set the sampled reference data set as nested; still keep random sampling for alt1
## 2. for null1.1, try zinegbin, and try pois
## 3. for alt1, try posi, for both inv and cov, erdos_renyi cond_1000

## update on 10/19/2020
## 1. for null1.1 and alt1, use american gut project data as previously did
## 2. for null2 and alt2, generate total count based on the minimum sample depth we see in the data provided by Jing
  ## notes: for alt2 and null2, run p=200 but the available reference data only have 127 available. But still just use its minimum seq depth
## 3. use zinegbin for copular model; use 3*p number of edges for alternative models

# filepath = '\\\\fs2-vip\\students\\yuek\\Desktop\\micro_net' #BOX
# 
# filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph' # MAC
# 
# filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph' #PC



filepath = '~/Desktop/micro_net' #bayes


setwd(filepath)
source('lib/func_libs.R')
source('Kun_code/generation_function.R')
source('Kun_code/method_function.R')
library(SpiecEasi)
library(MASS)

# use this as reference count data set
data(amgut1.filt)
reference_data = amgut1.filt



# use this for minimum sampling depth reference
# load('Kun_code/reference_data.RData') 
# reference_data = t(sub_counts)



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
save_folder_name = as.character(args[8]) # subfolder name to save files



set.seed(102)
## rarify the data
library(phyloseq)
tmp = otu_table(reference_data, taxa_are_rows = F)
out = rarefy_even_depth(tmp, sample.size = 1000)
mean(as.matrix(out)==0)
reference_data = out


# do not need to subsample the samples of reference data; only subsample the taxa
if(choose_model %in% c('null1', 'null1.1', 'alt1')) reference_data = reference_data[,sample(1:ncol(reference_data), size=p, replace = F)] # pick out p taxa

min_seq_depth = round(min(rowSums(reference_data))) # determine the min_seq_depth


#---------------
# part I: null hypothesis
#---------------

if(choose_model == 'null1'){
  ## Null 1
  # shuffle reference data set, empty graph for both cov and inv-cov

  p = ncol(reference_data)
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='shuffle', reference_data = reference_data, Sigma_list = Sigma_list)
}

if(choose_model == 'null1.1'){
  ## Null1 1.1
  # Use reference data set, marginal distributions are ZINB, and the network is set to empty
  
  p = ncol(reference_data)
  # generate a null graph object
  option = list(hypothesis = 'null', model = 'copula', 
                reference_data = reference_data, Sigma_list = NULL, distr = distr)
  option$Sigma_list$A_inv <-  option$Sigma_list$A_cov <- matrix(0, p, p)
  option$Sigma_list$Sigma <- option$Sigma_list$Omega <- diag(1, p)
}

if(choose_model == 'null2'){
  ## Null 2
  ## generate from Dirichlet distribution, so will assume both inv and cov based graph being empty?
  
  # mu = runif(p, 0, 8) # the mu is for underlying log-normal of Dirichlet parameters
  mu = runif(p, 0, 4)
  alpha=100
  # library_scale = rnegbin(n,mu=3e4, theta=9e8/(3e6-3e4)) # use library_scale to directly specify the total counts per sample
  # library_scale = rep(3e4, n)
  library_scale = runif(n, min=1, max=10) * min_seq_depth # generate total abundance from real data
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='Dirichlet', library_scale = library_scale, alpha = alpha, mu=mu, Sigma_list = Sigma_list, distr = distr)
}



#---------------
# part II: Alternative hypothesis
#---------------


if(choose_model=='alt1'){
  ## Alternative 1
  ## use Copula generative model, graph is based on inverse covariance matrix.
  p = ncol(reference_data)
  # Sigma_list = SpiecEasi_graph_Sigma(p,e=p, type = network_option, graph=NULL, network_condition_number = network_condition_number)
  Sigma_list = SpiecEasi_graph_Sigma(p,e=3*p, type = network_option, graph=NULL, network_condition_number = network_condition_number)
  
  
  option = list(hypothesis = 'alternative', model = 'copula', reference_data = reference_data, 
                Sigma_list = Sigma_list, distr = distr, network_option = network_option, network_condition_number = network_condition_number)
}


if(choose_model == 'alt2'){
  ## Alternative 2: log normal
  ## generate graph, p nodes and e edges; 
  ## this graph corresponds to inverse covariance
  
  # Sigma_list  = SpiecEasi_graph_Sigma(p,e=p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  Sigma_list  = SpiecEasi_graph_Sigma(p,e=3*p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  
  mu = runif(p, 0, 4)
  option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma_list = Sigma_list, 
                distr = distr, network_option = network_option, network_condition_number = network_condition_number)
  
}

if(choose_model == 'alt3'){
  ## Alternative 3: logistic normal multinomial
  ## generate graph, p nodes and e edges; 

  # Sigma_list  = SpiecEasi_graph_Sigma(p,e=p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  Sigma_list  = SpiecEasi_graph_Sigma(p,e=3*p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
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
  ## add this additional step of total count variation
  for(i in 1:nrow(data))data[i,] = data[i,]*round(runif(1,min=1, max=10))
  
  colnames(data) <- paste('taxa', 1:ncol(data))
  rownames(data) <- paste('cell', 1:nrow(data))
  
  data_rep[[1, k]] = count_to_comp(data) # compositions
  data_rep[[2, k]] = data # counts
}

save(data_rep, file = 
       paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number, '/vary_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))
save(list=c('data_rep', 'option', 'n', 'p'), file = 
       paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number,  '/vary_image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))






