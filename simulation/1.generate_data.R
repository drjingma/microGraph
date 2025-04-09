#######################
##
##
##  Data Examples for book chapter (Part 1: generate data sets)
##
## Kun Yue, yuek@uw.edu
## Date: 2020/11/02
##
## (modify the file path before running on your computer)
########################


#### This script can be run with command line arguments.
#### command format: 
# Rscript 1.generate_data.R n p choose_model nreps marginal_distribution network_structure network_condition_number folder_name


#### Examples to run the script: 
# Rscript 1.generate_data.R 100 127 alt1 200 zinegbin erdos_renyi 100 simulation1
# Rscript 1.generate_data.R 100 127 null1.1 200 zinegbin none 0 simulation1
# Rscript 1.generate_data.R 100 200 null2 200 none none 0 simulation1
# Rscript ../microGraph/simulation/1.generate_data.R 100 127 alt1 10 zinegbin erdos_renyi 1000 simulation1


######
## updated on 9/24/2020: 
## 1. adding generating partial correlation network from AR(1)/chain graph; use rho = 0.5 (small) and rho=0.8 (large) cases, use small n and large p

## update on 10/10/2020:
## 1. use #edge=3p for erdos_renyi graphs, under alt1 and alt2 models
## 2. change null2 library_scale to be based on real data: sample from 1*min(total_count) to 10*min(total_count)


# filepath = '......' 
# setwd(filepath)

source('~/Dropbox/Projects/Microbiome/microGraph/simulation/generation_function.R')
source('~/Dropbox/Projects/Microbiome/microGraph/simulation/method_function.R')

library(SpiecEasi)
library(MASS)

# use this as reference count data set
data(amgut1.filt)
reference_data = amgut1.filt



args = commandArgs(trailingOnly = T)
print(args)

n = as.integer(args[1])
p = as.integer(args[2]) # n and p override by reference data set, e.g.
choose_model = as.character(args[3]) # choose among 'null1' (shuffle reference data),
                                    # 'null1.1' (use copula model with NB marginal distr, based on reference data),
                                    # 'null2' (generate from Dirichelt distribution with varying total count),
                                    # 'alt1' (generate with copula model and reference data),
                                    # 'alt2' (generate from log-normal distribution directly)
                                    # 'alt3' (generate from logisitic normal multinomial with fixed total count)
nreps = as.integer(args[4]) # run 200 for all methods
distr = as.character(args[5]) # specify the distribution for copula model, if used
network_option = as.character(args[6]) # 'erdos_renyi', or 'chain_small', 'chain_large' for AR(1) (small: rho=0.5, large:rho=0.8), 'cov_erdos_renyi' for correlation based graph
network_condition_number = as.numeric(args[7]) # specify the condition number of the inverse covariance matrix
save_folder_name = as.character(args[8]) # subfolder name to save files
# n <- 100
# p <- 127
# choose_model <- "alt1"
# nreps <- 10
# distr <- 'zinegbin'
# network_option <- 'erdos_renyi'
# network_condition_number <- 100
# save_folder_name <- 'simulation1'

set.seed(102)
## rarify the data (saved to a separate file)
library(phyloseq)
tmp = otu_table(reference_data, taxa_are_rows = F)
out = rarefy_even_depth(tmp, sample.size = 1000)
mean(as.matrix(out)==0)
reference_data = out


# only subsample the OTUs if have a small p
if(choose_model %in% c('null1', 'null1.1', 'alt1')) reference_data = reference_data[,sample(1:ncol(reference_data), size=p, replace = F)] # pick out p OTUs

min_seq_depth = round(min(rowSums(reference_data))) # determine the min of total count

##--------------
## create appropriate option variable for each setting
##--------------

#---------------
# part I: null model
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
  ## generate from Dirichlet distribution, so will assume both inv and cov based graph being empty
  mu = runif(p, 0, 4)
  alpha=100
  library_scale = runif(n, min=1, max=10) * min_seq_depth # generate total abundance from real data
  Sigma_list = list(Sigma=NULL, Omega=NULL, A_inv = diag(0, p), A_cov = diag(0, p))
  option = list(hypothesis = 'null', model='Dirichlet', library_scale = library_scale, alpha = alpha, mu=mu, Sigma_list = Sigma_list, distr = distr)
}


#---------------
# part II: Alternative models (non-emtpy network)
#---------------

if(choose_model=='alt1'){
  ## Alternative 1
  ## use Copula generative model, graph is based on inverse covariance matrix.
  p = ncol(reference_data)
  Sigma_list = SpiecEasi_graph_Sigma(p,e=3*p, type = network_option, graph=NULL, network_condition_number = network_condition_number)
  option = list(hypothesis = 'alternative', model = 'copula', reference_data = reference_data, 
                Sigma_list = Sigma_list, distr = distr, network_option = network_option, network_condition_number = network_condition_number)
}


if(choose_model == 'alt2'){
  ## Alternative 2: log normal
  Sigma_list  = SpiecEasi_graph_Sigma(p,e=3*p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  mu = runif(p, 0, 4)
  option = list(hypothesis = 'alternative', model = 'log-normal', mu = mu, Sigma_list = Sigma_list, 
                distr = distr, network_option = network_option, network_condition_number = network_condition_number)
}

if(choose_model == 'alt3'){
  ## Alternative 3: logistic normal multinomial
  Sigma_list  = SpiecEasi_graph_Sigma(p,e=3*p, type=network_option, graph = NULL, network_condition_number = network_condition_number)  # output=list(Sigma for covariance/Correlation, Omega for inv-cov, A for adjacency matrix)
  mu = runif(p, 0, 4)
  library_scale = rep(3e4,n)
  option = list(hypothesis = 'alternative', model = 'multinomial-log-normal', 
                mu = mu, Sigma_list = Sigma_list, distr = distr, 
                library_scale = library_scale, network_option = network_option, network_condition_number = network_condition_number)
  
}

##-----------------------------------------------------------
##
## Generate nrep replications for analysis; data_rep[[1,k]] is composition matrix; data_rep[[2,k]] is corresponding count matrix
##
## (the data matrix is not zero-corrected)
##-------------------------------------------------------------


data_rep = matrix(list(), 2, nreps)
for(k in 1:nreps){
  set.seed(k*100)
  data_list = data_generation(n, p, option)
  dim(data_list$data)
  data = data_list$data
  
  ## add this additional step of total count variation, for null1.1 and alt1 models only
  if(choose_model %in% c('alt1', 'null1.1')){
    for(i in 1:nrow(data))data[i,] = data[i,]*round(runif(1,min=1, max=10))
  }
  colnames(data) <- paste('taxa', 1:ncol(data))
  rownames(data) <- paste('cell', 1:nrow(data))
  
  data_rep[[1, k]] = count_to_comp(data) # compositions
  data_rep[[2, k]] = data # counts
}

output_folder = paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number)
if (!dir.exists(output_folder)){
  dir.create(output_folder, showWarnings = T, recursive = TRUE)
  print(output_folder)
}

saveRDS(c('data_rep', 'option', 'n', 'p'),
        file =
       paste0(output_folder, '/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.rds'))






