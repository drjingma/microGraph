filepath = 'E:\\Dropbox\\Microbial_Networks\\microGraph' #PC

setwd(filepath)

args = c('newref2p', 100, 200, 'alt2', 200, 'none', 'erdos_renyi', 1000, 1, 3) # this is the updated data we are using
args = c('newref2p', 300, 200, 'alt2', 200, 'none', 'erdos_renyi', 1000, 1, 3) # this is the updated data we are using
args = c('newref2p', 100, 127, 'alt1', 200, 'zinegbin', 'erdos_renyi', 100, 1, 3) # this is the updated data we are using
args = c('newref2p', 500, 127, 'null1.1', 200, 'zinegbin', 'none', 0, 1, 3) # this is the updated data we are using

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


  
  
  load(paste0('data/',save_folder_name,'/',  distr,'/', network_option, '/cond_', network_condition_number, '/image_n_', n, '_p_', p, '_', choose_model, '_nreps_', nreps, '_data_rep.RData'))

  source('lib/func_libs.R')
  source('Kun_code/generation_function.R')
  source('Kun_code/method_function.R')
  library(SpiecEasi)
  library(MASS)
  
  set.seed(2020)
  
  
  target_graph_cov = option$Sigma_list$A_cov
  target_graph_inv = option$Sigma_list$A_inv
  
  i=4
  
  
  
  
  
  # time = sum(system.time({
  #   # use glasso for alternative model estimation and evaluation
  #   sparcc_res =  sparcc(data = date_rep[[2,i]],iter = 20, inner_iter = 10, th = 0.1)
  # })[2:3])
  # 
  # 
  # sparpath = huge(sparcc_res$Cor, method='glasso', lambda = lambda_seq2)
  # tmp2 = huge.roc(sparpath$path, target_graph_inv)
  
  
  source(paste0(filepath, "/Kun_code/CCLasso-master/R/cclasso.R"));
  
  # for CCLasso, the input need to be zero-corrected. so add peudocount 1 if there are zero counts
  for( k in 1:ncol(data_rep)){
    if (any(data_rep[[2,k]]==0)){
      print('zero counts corrected by adding pseudocount=1')
      data_rep[[2,k]] = data_rep[[2,k]]+1
      data_rep[[1,k]] = sweep(data_rep[[2,k]],1,STATS = rowSums(data_rep[[2,k]]), FUN='/')
    }
  }
  
  
  
  source(paste0(filepath, "/Kun_code/CCLasso-master/R/cclasso.R"));
  keep = list()
  t1 = NULL
  check_n = c(c(200,250,300,350,400,450,500),305,310, 320, 325, 330, 340, 345)
  check_n = sort(check_n)
  for(i in 1:length(check_n)){
    set.seed(i)
    print(check_n[i])
    time= sum(system.time({
      res_ccl_count <- cclasso(x = data_rep[[1,10]][1:(check_n[i]),],# input format n by p
                               counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                               k_cv = 3, # cv folds, min=2
                               lam_int = c(1e-6, 3), #tuning parameter value range; need to vary this for ROC curve
                               k_max=20, 
                               n_boot =2)  # n_boot for estimateing the p_value
      
    })[2:3])
    # res_ccl_count$cor_w[1:5, 1:5]
    keep[[i]] = res_ccl_count
    t1 = c(t1, mean(abs(res_ccl_count$cor_w-diag(1, p))>1e-11))
  }
    t1
  
    modify = list()
  set.seed(2)
  modify[[1]] <- cclasso(x = data_rep[[1,10]][1:250,],# input format n by p
                           counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                           k_cv = 3, # cv folds, min=2
                           lam_int = c(1e-6, 8), #tuning parameter value range; need to vary this for ROC curve
                           k_max=20, 
                           n_boot =2)  # n_boot for estimateing the p_value
  
  set.seed(3)
  modify[[2]] <- cclasso(x = data_rep[[1,10]][1:300,],# input format n by p
                         counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                         k_cv = 3, # cv folds, min=2
                         lam_int = c(1e-6, 8), #tuning parameter value range; need to vary this for ROC curve
                         k_max=20, 
                         n_boot =2)
  
  set.seed(4)
  modify[[3]] <- cclasso(x = data_rep[[1,10]][1:305,],# input format n by p
                         counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                         k_cv = 3, # cv folds, min=2
                         lam_int = c(1e-6, 8), #tuning parameter value range; need to vary this for ROC curve
                         k_max=20, 
                         n_boot =2)
  
  mean(abs(modify[[1]]$cor_w-diag(1, p))>1e-11)
  mean(abs(modify[[2]]$cor_w-diag(1, p))>1e-11)
  mean(abs(modify[[3]]$cor_w-diag(1, p))>1e-11)
  
  t2 = t1
  t2[2:4]<-c(  mean(abs(modify[[1]]$cor_w-diag(1, p))>1e-11),
               mean(abs(modify[[2]]$cor_w-diag(1, p))>1e-11),
               mean(abs(modify[[3]]$cor_w-diag(1, p))>1e-11))
  plot(check_n, t2)
  
  save.image('inspect_CCLasso.RData')
  
  
  tmp = cclasso(x = data_rep[[1,10]][1:305,],# input format n by p
          counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
          k_cv = 3, # cv folds, min=2
          lam_int = c(1e-6, 8), #tuning parameter value range; need to vary this for ROC curve
          k_max=20, 
          n_boot =500)
  
  
  i=3
  check_n[i] # corresponding 
  cor = keep[[i]]
  diag(cor)<-0
  ggplot(data.frame(x=as.vector(row(cor)), y=as.vector(col(cor)), z=as.vector(cor)),
         aes(x=x, y=y, fill=z))+
    geom_tile()+
    scale_fill_gradient2()
  
  mean(abs(cor)^2) # F2 norm
  
  
  time= sum(system.time({
    res_ccl_count <- cclasso(x = data_rep[[1,10]][1:(check_n[i]),],# input format n by p
                             counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                             k_cv = 8, # cv folds, min=2
                             lam_int = c(1e-6, 3), #tuning parameter value range; need to vary this for ROC curve
                             k_max=20, 
                             n_boot =2)  # n_boot for estimateing the p_value
    
  })[2:3])
  res_ccl_count$cor_w[1:5, 1:5]
  mean(abs(res_ccl_count$cor_w-diag(1, p))>1e-11)
  
  
  time= sum(system.time({
    res_ccl_count <- cclasso(x = data_rep[[1,10]][1:(check_n[i]),],# input format n by p
                             counts = F, pseudo = 1, # for correction of count matrix; ignored if we direclty supply composition matrix
                             k_cv =4, # cv folds, min=2
                             lam_int = c(1e-6, 3), #tuning parameter value range; need to vary this for ROC curve
                             k_max=50, 
                             n_boot =2)  # n_boot for estimateing the p_value
    
  })[2:3])
  mean(abs(res_ccl_count$cor_w-diag(1, p))>1e-11)
  
  
  plot( res_ccl_count$info_cv$lams,res_ccl_count$info_cv$fvals)
  
  
  
  
  
  
  plot(c(200,250,300,350,400,450,500),t1)
  plot(c(305,310, 320, 325, 330, 340, 345),t2)
  
    ccpath = huge(res_ccl_count$cor_w, method='glasso', 
                  lambda = sort(lambda_seq2, decreasing = T))
    sum(res_ccl_count$cor_w-diag(1, p))
    tmp1 = huge.roc(ccpath$path, target_graph_inv)
    tmp1$tp
    # ccpath = huge(res_ccl_count$cor_w, method='glasso', 
    #               lambda = sort(lambda_seq2, decreasing = F))
    # huge.roc(ccpath$path, target_graph_inv)
    
    tmp = huge.roc(ccpath$path, target_graph_inv)
    tmp
    keep[[i]] = list(fp =tmp$fp, tp=tmp$tp)
  }
  plot(tmp$fp, tmp$tp)
  
  plot(c(0, 1), c(0, 1))
  for(i in 1:20){
    points(keep[[i]]$fp, keep[[i]]$tp, color='gray', type='l')
  }
  
  
  # ROC_inv = find_glasso_ROC(res_ccl_count$cor_w, lambda_seq, target_graph_inv)
  
  ccpath = huge(res_ccl_count$cor_w, method='glasso', lambda = lambda_seq2)
  tmp = huge.roc(ccpath$path, target_graph_inv)
  abline(a=0, b=1)
  plot(tmp$fp, tmp$tp, type='l')
  

  points(tmp2$fp, tmp2$tp, type='l', col=2)
  
  
  lambda_seq2 = c(seq(0.01,0.1,0.015),seq(0.1,6,0.25))*sqrt(log(p)/n) 
  
  
  1
  
  
  
  
  
  
  
  
  
  
  
  
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
  
  target_graph_cov = option$Sigma_list$A_cov
  target_graph_inv = option$Sigma_list$A_inv
  
  use_data = data_rep[[2, 3]]
  
  Spiec_network <- spiec.easi(use_data, method='glasso', # choose from 'mb' for neighbourhood selection, or 'glasso'
                               # if to perform model selection
                               pulsar.select = T,
                               pulsar.params = list(
                                 thresh=0.05,# Threshold for StARS criterion.
                                 subsample.ratio=0.8, # Subsample size for StARS.
                                 rep.num = 20), # Number of subsamples for StARS.
                               lambda.min.ratio= min(lambda_seq)/max(lambda_seq), lambda.max = max(lambda_seq), nlambda=length(lambda_seq)
  )
                               
                               
  Spiec_ROC_res = huge::huge.roc(path = Spiec_network$est$path, 
                      theta = target_graph_inv, verbose=FALSE)
  
  
  
  #####
  ## another way: run glasso on clr transformed data
  ?clr
  
  dim(use_data)
  use_data[1:4, 1:4]
  tmp = clr(use_data, mar=1)
  tmp[1:4, 1:4]
  colSums(tmp) # colSums are centered to 0, attention to the rotation of the matrix
  
  # for huge, need x is n by d
  set.seed(1)
  use_clr = huge(t(tmp), method='glasso', lambda = lambda_seq)
  huge.roc(use_clr$path, theta = target_graph_inv)  # not much difference
  
  lambda_seq2 = c(seq(0.01,0.1,0.015),seq(0.1,6,0.25))*sqrt(log(p)/n) 
  
  ## one problem might be the reason: if I do not force the lambda sequence, seems to be god
  set.seed(1)
  use_clr_180 = huge(t(tmp)[1:180,], method='glasso')
  huge.roc(use_clr_180$path, theta = target_graph_inv)  # good
  
  
  set.seed(1)
  use_clr_180_2 = huge(t(tmp)[1:180,], method='glasso', lambda = use_clr_180$lambda)
  huge.roc(use_clr_180_2$path, theta = target_graph_inv)  # good
  
  set.seed(1)
  use_clr_180_3 = huge(t(tmp)[1:180,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_180_3$path, theta = target_graph_inv)  # not good
  
  set.seed(1)
  use_clr_180_3_2 = huge(t(tmp)[1:180,], method='glasso', lambda = lambda_seq2[15:36])
  huge.roc(use_clr_180_3_2$path, theta = target_graph_inv)  # not good
  
  set.seed(1)
  tmp_lambda = c(lambda_seq,use_clr_180$lambda)
  tmp_lambda = sort(tmp_lambda, decreasing = T)
  use_clr_180_4 = huge(t(tmp)[1:180,], method='glasso', lambda =  tmp_lambda)
  huge.roc(use_clr_180_4$path, theta = target_graph_inv)  # not good?? why
  

  
  # how about we just use the first x samples, will it get better?
  set.seed(1)
  use_clr_100 = huge(t(tmp)[1:100,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_100$path, theta = target_graph_inv)  # really good                     
  
  set.seed(1)
  use_clr_150 = huge(t(tmp)[1:150,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_150$path, theta = target_graph_inv)  # really good
  
  set.seed(1)
  use_clr_last_150 = huge(t(tmp)[50:200,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_last_150$path, theta = target_graph_inv)  # really good!!!
  
  set.seed(1)
  use_clr_151 = huge(t(tmp)[1:151,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_151$path, theta = target_graph_inv)  # bad??!!
  
  set.seed(1)
  use_clr_150_rand = huge(t(tmp)[sample(1:200, 150, replace=F),], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_150_rand$path, theta = target_graph_inv)  # ??
  
  
  set.seed(2)
  use_clr_151_rand = huge(t(tmp)[sample(1:200, 151, replace=F),], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_151_rand$path, theta = target_graph_inv)  # ??
  
  
  set.seed(1)
  use_clr_155 = huge(t(tmp)[1:155,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_155$path, theta = target_graph_inv)  #  not good
  
  set.seed(1)
  use_clr_160 = huge(t(tmp)[1:160,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_160$path, theta = target_graph_inv)  # not good
  
  set.seed(1)
  use_clr_180 = huge(t(tmp)[1:180,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_180$path, theta = target_graph_inv)  # not good
  
  set.seed(1)
  use_clr_190 = huge(t(tmp)[1:190,], method='glasso', lambda = lambda_seq)
  huge.roc(use_clr_190$path, theta = target_graph_inv)  # not good 

  
  

  