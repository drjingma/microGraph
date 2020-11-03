###-----------------------
### this file contains utility functions for reading files and generating plots
###-----------------------

##------------
## for getting null model results and plots
##------------
# non-existing result files will be printed
get_summary_fp = function(add_prex,folder_name, choose_model, n, p, nreps, part, part_name, copula_distr, network_option, network_condition_number, include){
  # conet and Sparcc have p values, others only fp_null
  res_for_200 = do.call(rbind, lapply(1:include, function(run_rep){
    load_file = paste0(filepath,
                       '/', folder_name, '/', copula_distr,'/', network_option, '/cond_', network_condition_number,'/', choose_model, '/',add_prex ,'res_n_', n, '_p_', p, '_', choose_model, 
                       '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData')
    # print(load_file)
    iferror = tryCatch(load(load_file), error = function(e){
      print(load_file)
      return(e)
    })
    if('simpleError' %in% class(iferror)) return(NULL)
    if(part %in% c(1, 2)){
      rbind(data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                       fp=roc$fp_null$fp_null, adjust = 'none', rep = run_rep, time = roc$time),
            data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                       fp=roc$fp_null$fp_null_FDR, adjust = 'FDR', rep=run_rep, time = roc$time))
    }else if(part %in% c(3, 4)){
      data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                 fp=roc$fp_null, adjust = 'none', rep=run_rep, time = roc$time)
    }else{
      data.frame(n=n, p=p, model = choose_model, method = part_name[[part]], 
                 fp=roc$fp_inv_null, adjust = 'none', rep=run_rep, time = roc$time)
    }
  }))
  if(!is.null(res_for_200)){
    res_for_200$method = factor(res_for_200$method, levels=do.call(c,part_name))
    
    output = with(res_for_200, aggregate(fp, by = list(n=n, p=p, model= model, method=method, adjust=adjust), 
                                         function(x)mean(x,na.rm=T) ))
    
    time_out = with(res_for_200, aggregate(time, by = list(n=n, p=p, model= model, method=method, adjust=adjust), 
                                           function(x)mean(x,na.rm=T) ))
    
    output$time = time_out$x
    output
  }else{
    NULL
  }
  
  
  
}

##------------
## for getting network recovery ROC results
##------------
get_summary_roc = function(add_prex, folder_name, choose_model, n, p, nreps, part, part_name, copula_distr, network_option, network_condition_number, include){
  # all models we compare the inverse cov graph ROC
  # to have an 'averaged' ROC curve, we take average of fp and tp at each grid,(since we are using the same sequence of lambda)
  res_for_20 = lapply(1:include, function(run_rep){
    load_file = paste0(folder_name, '/', copula_distr,'/', network_option, '/cond_', network_condition_number,'/', choose_model, '/',add_prex,'res_n_', n, '_p_', p, '_', choose_model, 
                       '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData')
    #print(load_file)
    iferror = tryCatch(load(load_file), error=function(e){
      print(load_file)
      e}
    )
    if('simpleError' %in% class(iferror)) return(NULL)
    list(n=n, p=p, model = choose_model, method = part_name[[part]],rep=run_rep,
         tp = roc$ROC_inv$tp, fp = roc$ROC_inv$fp, auc = roc$ROC_inv$AUC,
         tp_cov = roc$ROC_cov$tp, fp_cov = roc$ROC_cov$fp, auc_cov = roc$ROC_cov$AUC, time = roc$time)
    
  })
  # average time information
  
  res_for_20 = res_for_20[!sapply(res_for_20, is.null)]
  if(length(res_for_20)==0) return(NULL)
  
  ave_time = mean(do.call(c, sapply(res_for_20, `[`, 'time')))
  
  # inv ROC results
  tps = sapply(res_for_20, `[`, 'tp')
  if(length(unique(sapply(tps, FUN = length)))>1) stop(paste('check lambda sequence, tp not same length across replications: n',n,'p',p, 'part', part))
  fps = sapply(res_for_20, `[`, 'fp')
  if(length(unique(sapply(fps, FUN = length)))>1) stop(paste('check lambda sequence, fp not same length across replications: n',n,'p',p, 'part', part))
  
  if(is.null(do.call(cbind, tps))){ave_tp<-NA}else{ave_tp = rowMeans(do.call(cbind, tps))}
  if(is.null(do.call(cbind, fps))){ave_fp =  NA}else{ave_fp =  rowMeans(do.call(cbind, fps))}
  if(is.null(do.call(cbind, sapply(res_for_20, `[`, 'auc')))){ave_auc = NA}else{ave_auc =  mean(do.call(cbind, sapply(res_for_20, `[`, 'auc')))}
  
  ave_res = list(n=n, p=p, model = choose_model, method = part_name[[part]],
                 ave_tp = ave_tp, ave_fp = ave_fp, ave_auc = ave_auc, ave_time = ave_time)
  
  # cov ROC results (not really used, discard)
  ave_res_cov = NULL
  if(F){
    tps_cov = sapply(res_for_20, `[`, 'tp_cov')
    if(length(unique(sapply(tps_cov, FUN = length)))>1) stop(paste('check lambda sequence, tp_cov not same length across replications: n',n,'p',p, 'part', part))
    fps_cov = sapply(res_for_20, `[`, 'fp_cov')
    if(length(unique(sapply(fps_cov, FUN = length)))>1) stop(paste('check lambda sequence, fp_cov not same length across replications: n',n,'p',p, 'part', part))
    
    ave_tp_cov = rowMeans(do.call(cbind, tps_cov))
    ave_fp_cov = rowMeans(do.call(cbind, fps_cov))
    if(is.null(do.call(cbind, sapply(res_for_20, `[`, 'auc_cov')))){ave_auc_cov = NA}else{ave_auc_cov =  mean(do.call(cbind, sapply(res_for_20, `[`, 'auc_cov')))}
    
    ave_res_cov = list(n=n, p=p, model = choose_model, method = part_name[[part]],
                       ave_tp = ave_tp_cov, ave_fp = ave_fp_cov, ave_auc = ave_auc_cov, ave_time = ave_time)
  }

  
  
  return(list(ave_res = ave_res, ave_res_cov = ave_res_cov, res_for_20 = res_for_20))
}

plot_roc = function(roc_res, target){
  if(target=='inv'){
    tmp1=data.frame(tp = roc_res$ave_res$ave_tp, fp = roc_res$ave_res$ave_fp, group=rep('mean', length(roc_res$ave_res$ave_fp)))
    tmp2= do.call(rbind, lapply(1:length(roc_res$res_for_20), function(i)
      data.frame(tp = roc_res$res_for_20[[i]]$tp, fp = roc_res$res_for_20[[i]]$fp, 
                 group=rep(as.character(i), length(roc_res$res_for_20[[i]]$fp)))))
    
  }
  if(target=='cov'){
    tmp1=data.frame(tp = roc_res$ave_res_cov$ave_tp, fp = roc_res$ave_res_cov$ave_fp, 
                    group=rep('mean', length(roc_res$ave_res_cov$ave_fp)))
    tmp2= do.call(rbind, lapply(1:length(roc_res$res_for_20), function(i)
      data.frame(tp = roc_res$res_for_20[[i]]$tp_cov, fp = roc_res$res_for_20[[i]]$fp_cov, 
                 group=rep(as.character(i), length(roc_res$res_for_20[[i]]$fp_cov)))))
    
    
  }
  
  tmp = rbind(tmp1, tmp2)  
  tmp$n = roc_res$ave_res$n[1]
  tmp$p = roc_res$ave_res$p[1]
  tmp$model = roc_res$ave_res$model[1]
  tmp$method = roc_res$ave_res$method[1]
  tmp$time = roc_res$ave_res$ave_time
  
  
  gplot = ggplot(tmp, aes(x=fp, y=tp, group=group, color=group))+
    geom_line()+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    ggtitle(paste('n',tmp$n[1], 'p', tmp$p[1], tmp$model[1], tmp$method[1]))+
    xlab('False Positive')+
    ylab('True Positive')+
    scale_color_manual(values = c('black', rep('grey',10)))+
    scale_linetype_manual(values = c(1, rep(4, 10)))+
    geom_line(data = tmp[tmp$group=='mean',], aes(x=fp, y=tp))+
    theme(legend.position = 'none')
  
  return(list(plot = gplot, data = tmp))
}

output_plot_all = function(add_prex,folder_name, choose_model, n, p, nreps, part_name, copula_distr, network_option, network_condition_number, include, target){
  data = do.call(rbind, lapply(c(2:5,6, 7), function(part){
    roc_res = get_summary_roc(add_prex,folder_name, choose_model, n, p, nreps, part, part_name, copula_distr, network_option, network_condition_number, include)
    if(is.null(roc_res)){NULL}else{plot_roc(roc_res, target)$data}
  }))
  data$group = factor(data$group, levels = c(1:nreps, 'mean'))
  data$method = factor(data$method, levels=do.call(c,part_name))
  
  plot_all = ggplot(data, aes(x=fp, y=tp, group=group, color = group))+
    geom_line()+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    ggtitle(paste('dist', copula_distr, 'n',data$n[1], 'p', data$p[1], data$model[1]))+
    theme(legend.position = 'none')+
    scale_color_manual(values=c(rep('grey', 200), 2))+
    facet_wrap(~method)
  
  # ggsave(plot_all, filename = paste0('plot/', copula_distr, '_n',data$n[1], 'p', data$p[1], data$model[1],'.png'))
  return(list(data= data, plots = plot_all))
}


###---------------------
### for each of the four models in the manuscript, generate summary plots
###---------------------

library(ggplot2)
library(kableExtra)
library(cowplot)
theme_set(theme_cowplot())

# specify a color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7","#F0E442")
part_name = list('1'='ReBoot', '2' = 'SparCC', '3' = 'CCLasso', '4' = 'COAT','5' = 'SPIEC-EASI', '6' = 'gCoda', '7' = 'SPRING')


##------------
## null1.1 model
##------------
plt_null1.1 = function(
  n_list,
  p,
  folder_name, 
  nreps,
  include, # number of reps to summarize
  copula_distr, 
  choose_model, 
  network_option, 
  network_condition_number,
  fix = NULL){
  
  
  null1.1= lapply(n_list, 
                  function(n) cbind(do.call(rbind, lapply(1:7,  
                                                          function(part)
                                                            get_summary_fp(fix, folder_name, choose_model,n=n, p=p, nreps=nreps, part=part, part_name, copula_distr=copula_distr,
                                                                           network_option = network_option, network_condition_number = network_condition_number, include= include)
                  )
                  ),
                  dist = copula_distr, 
                  network_option = network_option, 
                  network_condition_number = network_condition_number)
  )
  
  null1.1 = do.call(rbind, null1.1)
  
  
  null1.1_none = null1.1[null1.1$adjust=='none', ]
  pp = ggplot(data = null1.1_none, aes(x=n, y=x, group=method, color=method, 
                                       linetype = method))+
    coord_cartesian(ylim = c(0, min(max(null1.1_none$x), 0.4)))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    ylab('Type I Error')+
    scale_shape_manual(values=1:7)+
    scale_color_manual(values=cbPalette[1:7])+
    geom_hline(yintercept = 0.05, color=1, linetype=2)+
    theme(legend.title = element_blank())+
    scale_linetype_manual(values=1:7)+
    # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
    scale_x_continuous(breaks=c(100,  200, 300,  400,  500))
  
  
  # summarize time information for each method under varying sample size
  pp_time = ggplot(data = null1.1_none, aes(x=n, y=log10(time), group=method, color=method, 
                                            linetype = method))+
    # coord_cartesian(ylim = c(0, min(max(null1.1_none$time), 0.2)))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    ylab('log10(time) (seconds)')+
    scale_shape_manual(values=1:7)+
    scale_color_manual(values=cbPalette[1:7])+
    theme(legend.title = element_blank())+
    scale_linetype_manual(values=1:7)+
    # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
    scale_x_continuous(breaks=c(100, 200,300,  400, 500))
  
  return(list(res = pp, time=pp_time))
}

##------------
## null2 model
##------------

plt_null2 = function(
  n_list,
  p,
  folder_name, 
  nreps,
  include, # number of reps to summarize
  copula_distr, 
  choose_model, 
  network_option, 
  network_condition_number,
  fix = NULL){
  
  null2 = lapply(p, function(p){
    null2= lapply(n_list, function(n)
      cbind(
        do.call(rbind, lapply(1:7,  function(part)
          get_summary_fp(fix, folder_name, choose_model,n=n, p=p, nreps=nreps, part=part, part_name, 
                       copula_distr=copula_distr,network_option = network_option, 
                       network_condition_number = network_condition_number, include= include)
                              )
                ),
        dist = copula_distr, 
        network_option = network_option, 
        network_condition_number = network_condition_number)
      )
    
    null2 = do.call(rbind, null2)
  })
  
  null2 = do.call(rbind, null2)
  
  null2_none = null2[null2$adjust=='none' & null2$p==200,]
  pp1 =ggplot(data = null2_none, aes(x=n, y=x, group=method, color=method, 
                                     linetype=method))+
    coord_cartesian(ylim = c(0, min(max(null2_none$x), 0.2)))+
    geom_line(size=1)+
    scale_linetype_manual(values=1:7)+
    geom_point(size=2, shape=18)+
    geom_hline(yintercept = 0.05, color=1, linetype=2)+
    theme(legend.title = element_blank())+
    ylab('Type I Error')+
    scale_shape_manual(values=1:7)+
    scale_color_manual(values=cbPalette[c(1:7)])+
    scale_x_continuous(breaks=c(100,  200, 300,400,500))+
    scale_y_continuous(breaks=c(0,0.01,0.02,0.03,0.04,0.05,0.06))
  
  
  pp_time = ggplot(data = null2_none, aes(x=n, y=log10(time), group=method, color=method, 
                                          linetype = method))+
    # coord_cartesian(ylim = c(0, min(max(null1.1_none$time), 0.2)))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    ylab('log10(time) (seconds)')+
    scale_shape_manual(values=1:7)+
    scale_color_manual(values=cbPalette[1:7])+
    theme(legend.title = element_blank())+
    scale_linetype_manual(values=1:7)+
    # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
    scale_x_continuous(breaks=c(100,  200, 300,400,500))
  
  
  return(list(res = pp1, time = pp_time))
}

##------------
## alt1 / alt2 model plot
##------------

plt_alt = function(
  n_list,
  p,
  folder_name, 
  nreps,
  include, # number of reps to summarize
  copula_distr, 
  choose_model, 
  network_option, 
  network_condition_number,
  fix = NULL,
  target = 'inv')
{
  
  data = lapply(n_list, function(n){
    output_plot_all(fix, folder_name, choose_model, n=n, p=p, nreps=nreps, 
                    part_name, copula_distr, network_option, network_condition_number,include= include, target)$data
  })
  data = do.call(rbind, data)
  data$n = as.factor(data$n)
  
  tmpdata = data[data$group=='mean' & data$method != 'ReBoot',] 
  plot_all = ggplot(tmpdata, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    #ggtitle(paste('dist', copula_distr, 'n',nn, data$model[1]))+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    facet_wrap(~n, labeller = label_both, nrow = 2)+
    scale_color_manual(values = cbPalette[2:7])+
    scale_linetype_manual(values = c(2:7))
  
  tmpdata2 = with(tmpdata, aggregate(time, by = list(n=n, p=p, method=method), FUN=mean, ))
  for(i in c(1))tmpdata2[,i] = as.numeric(levels(tmpdata2[,i]))[tmpdata2[,i]]
  plot_time = ggplot(tmpdata2, aes(x=n, y=log10(x), group=method, color = method,linetype=method))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    scale_shape_manual(values=2:7)+
    xlab('n')+
    ylab('log10(time) (seconds)')+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    scale_color_manual(values = cbPalette[2:7])+
    scale_linetype_manual(values = c(2:7))+
    scale_x_continuous(breaks=c(100, 200, 300, 400, 500))
  
  return(list(res = plot_all, time = plot_time))
}

