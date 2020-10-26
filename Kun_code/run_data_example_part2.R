# filepath = '//fs2-vip/students/yuek/Desktop/micro_net' #BOX
# 
# filepath = '/Users/Kun/Desktop/Dropbox/Microbial_Networks/microGraph' # MAC
# 
# filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph' #PC

filepath = '~/Desktop/micro_net' #bayes

setwd(filepath)



args = commandArgs(trailingOnly = T) # (n p null1 2 1),  n and p override by reference data set, e.g. args = c(100, 200, 'null2', 100, 8, 1)
# args = c(100, 200, 'alt2', 200, 'none', 'chain_large', 0, 6, 6)
n = as.integer(args[2])
p = as.integer(args[3])
choose_model = as.character(args[4])
nreps = as.integer(args[5])
distr = as.character(args[6])
network_option = as.character(args[7])
network_condition_number = as.numeric(args[8]) 
run_rep = as.integer(args[9])
part = as.character(args[10])
save_folder_name = args[1]

print(args)

# check if file exist
check_file=
       paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number, '/', choose_model,'/res_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData')

if_have_file = tryCatch(
  {
    load(check_file)
  }
  ,error =  function(e) e
)

if('simpleError' %in% class(if_have_file) ){
  

load(paste0('data/',save_folder_name,'/',  distr,'/', network_option, '/cond_', network_condition_number, '/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))


source('lib/func_libs.R')
source('Kun_code/generation_function.R')
source('Kun_code/method_function.R')
library(SpiecEasi)
library(MASS)

set.seed(2020)

# use gcoda implemented lambda path computation, based on the first rep data
x = data_rep[[2,1]]
if (any(x==0)){
  x = x+1
}
x = sweep(x,1,STATS = rowSums(x), FUN='/')

# Log transformation for compositional data
S <- var(log(x) - rowMeans(log(x)));
# Generate lambda via lambda.min.ratio and nlambda
lambda.max <- max(max(S - diag(p)), -min(S - diag(p)));
lambda.min.ratio = 1e-5
lambda.min <- lambda.min.ratio * lambda.max;


lambda_seq = exp(c(seq(log(lambda.min*0.9),log(lambda.max*1.1),length = 40)))
lambda_seq = sort(lambda_seq, decreasing = T)


# compare methods for their ROC
switch(
  part,
  '1' = {
    cat('CoNet/ReBoot \n')
    roc = compare_methods(data_rep[,run_rep, drop=F],
                                    est_mat = c('covariance', 'precision')[1],
                                    method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                               'SpiecEasi', 'gCoDa','Spring')[1],
                                    target_graph_cov = option$Sigma_list$A_cov,
                                    target_graph_inv = option$Sigma_list$A_inv,
                                    option = option,
                                    lambda_seq = lambda_seq)
  }, 
  '2' = {
    cat('SparCC \n')
    
    roc = compare_methods(data_rep[,run_rep, drop=F],  
                                     est_mat = c('covariance', 'precision')[1],
                                     method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                                'SpiecEasi', 'gCoDa','Spring')[2],
                                     target_graph_cov = option$Sigma_list$A_cov,
                                     target_graph_inv = option$Sigma_list$A_inv,
                                     option = option,
                                     lambda_seq = lambda_seq)
  },
  '3' = {
    cat('CCLasso \n')
    
    roc = compare_methods(data_rep[,run_rep, drop=F], 
                                      est_mat = c('covariance', 'precision')[1],
                                      method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                                 'SpiecEasi', 'gCoDa','Spring')[3],
                                      target_graph_cov = option$Sigma_list$A_cov,
                                      target_graph_inv = option$Sigma_list$A_inv,
                                      option = option,
                                      lambda_seq = lambda_seq)
  },
  '4' = {
    cat('coat \n')
    
    roc = compare_methods(data_rep[,run_rep, drop=F],  
                                   est_mat = c('covariance', 'precision')[1],
                                   method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                              'SpiecEasi', 'gCoDa','Spring')[4],
                                   target_graph_cov = option$Sigma_list$A_cov,
                                   target_graph_inv = option$Sigma_list$A_inv,
                                   option = option,
                                   lambda_seq = lambda_seq)
  },
  '5' = {
    cat('Spieceasi \n')
    
    roc = compare_methods(data_rep[,run_rep, drop=F],  
                                        est_mat = c('covariance', 'precision')[2],
                                        method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                                   'SpiecEasi', 'gCoDa','Spring')[5],
                                        target_graph_cov = option$Sigma_list$A_cov,
                                        target_graph_inv = option$Sigma_list$A_inv,
                                        option = option,
                                        lambda_seq = lambda_seq)
  },
  '6' = {
    cat('gCoDa \n')
    
    roc = compare_methods(data_rep[,run_rep, drop=F],  
                                    est_mat = c('covariance', 'precision')[2],
                                    method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                               'SpiecEasi', 'gCoDa','Spring')[6],
                                    target_graph_cov = option$Sigma_list$A_cov,
                                    target_graph_inv = option$Sigma_list$A_inv,
                                    option = option,
                                    lambda_seq = lambda_seq)
  },
  '7' = {
    cat('Spring \n')
    
    roc = compare_methods(data_rep[,run_rep, drop=F],  
                                     est_mat = c('covariance', 'precision')[2],
                                     method = c('CoNet', 'SparCC','CCLasso', 'COAT',
                                                'SpiecEasi', 'gCoDa','Spring')[7],
                                     target_graph_cov = option$Sigma_list$A_cov,
                                     target_graph_inv = option$Sigma_list$A_inv,
                                     option = option,
                                     lambda_seq = lambda_seq)
  }
)


cat('Done \n')


save(list = c('roc'), file=
       paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number, '/', choose_model,'/res_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData'))


}

