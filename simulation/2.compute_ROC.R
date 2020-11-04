#######################
##
##
##  Data Examples for book chapter (Part 2: compute ROC for each method)
##
## Kun Yue, yuek@uw.edu
## Date: 2020/11/02
##
## (modify the file path before running on your computer)
########################


#### This script can be run with command line arguments. 
#### command format: 
# Rscript 2.compute_ROC.R n p choose_model nreps marginal_distribution network_structure network_condition_number which_replication which_method folder_name


#### Examples to run the script: 
# Rscript 2.compute_ROC.R 100 127 alt1 200 zinegbin erdos_renyi 100 1 2 simulation1
# Rscript 2.compute_ROC.R 100 127 null1.1 200 zinegbin none 0 1 2 simulation1
# Rscript 2.compute_ROC.R 100 200 null2 200 none none 0 1 2 simulation1




filepath = '....' 

setwd(filepath)



args = commandArgs(trailingOnly = T) 
n = as.integer(args[1])
p = as.integer(args[2])
choose_model = as.character(args[3])
nreps = as.integer(args[4])
distr = as.character(args[5])
network_option = as.character(args[6])
network_condition_number = as.numeric(args[7]) 
run_rep = as.integer(args[8])
part = as.character(args[9])
save_folder_name = args[10]

print(args)

# check if file exist
check_file=
       paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number, '/', choose_model,'/fix_lambda_res_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData')

if_have_file = tryCatch(
  {
    load(check_file)
  }
  ,error =  function(e) e
)

if('simpleError' %in% class(if_have_file) ){
  
  load(paste0('data/',save_folder_name,'/',  distr,'/', network_option, '/cond_', network_condition_number, '/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))

  source('generation_function.R')
  source('method_function.R')
  library(SpiecEasi)
  library(MASS)
  
  set.seed(2020)
  
  # use a fixed lambda sequence
 
  lambda_seq = c(seq(0.01,0.1,0.015),seq(0.1,6,0.25))*sqrt(log(p)/n) 
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
  output_folder = paste0('data/',save_folder_name,'/', distr,'/', network_option, '/cond_', network_condition_number,'/', choose_model)
  dir.create(output_folder, showWarnings = F)
  
  save(list = c('roc'), file=
         paste0(output_folder, '/res_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData'))
  
  
}

