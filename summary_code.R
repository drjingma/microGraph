filepath = '\\\\fs2-vip\\students\\yuek\\Desktop\\micro_net' #BOX
setwd(filepath)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",  "#0072B2", "#D55E00", "#CC79A7","#F0E442")

# choose_model = c('null1', 'null2', 'alt1', 'alt2')[2]
# n = c(100, 200, 500, 289)[1]
# p = c(200, 127)[1]
# nreps=100
# run_rep = 1:10 # just compute for 10 repetitions
# part = (1:7)[1] # for the 7 methods
part_name = list('1'='ReBoot', '2' = 'SparCC', '3' = 'CCLasso', '4' = 'COAT','5' = 'SPIEC-EASI', '6' = 'gCoda', '7' = 'SPRING')


copula_distr = 'zinegbin'         # for null1.1 and alt1
copula_distr = 'none'             # for alt2; also for null2


# compute fp summary for null models, dist_data files
get_summary_fp = function(choose_model, n, p, nreps, part, part_name, copula_distr, network_option, network_condition_number){
  # conet and Sparcc have p values, others only fp_null
  res_for_200 = do.call(rbind, lapply(1:200, function(run_rep){
    load_file = paste0(filepath,
                       '/dist_data/', copula_distr,'/', network_option, '/cond_', network_condition_number,'/', choose_model, '/res_n_', n, '_p_', p, '_', choose_model, 
                       '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData')
    # print(load_file)
    tryCatch(load(load_file), error = function(e){
      print(e)
      print(load_file)})
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
  res_for_200$method = factor(res_for_200$method, levels=do.call(c,part_name))
  
  output = with(res_for_200, aggregate(fp, by = list(n=n, p=p, model= model, method=method, adjust=adjust), 
                              function(x)mean(x,na.rm=T) ))
  
  time_out = with(res_for_200, aggregate(time, by = list(n=n, p=p, model= model, method=method, adjust=adjust), 
                              function(x)mean(x,na.rm=T) ))
  
  output$time = time_out$x
  output
  
  
}

library(ggplot2)
library(kableExtra)
library(cowplot)
theme_set(theme_cowplot())


# copula_distr = 'none'
# 
# choose_model = 'null1'
# 
# null1= lapply(c(100, 289), function(n)cbind(do.call(rbind, lapply(1:7,  function(part)
#   get_summary_fp(choose_model,n=n, p=127, nreps=100, part=part, part_name, copula_distr=copula_distr))),
#   dist = copula_distr))
# null1 = do.call(rbind, null1)
# head(null1)
# 
# 
# pp = ggplot(data = null1, aes(x=n, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
#   coord_cartesian(ylim = c(0, min(max(null1$x), 0.2)))+
#   geom_line(data=null1, aes(linetype = interaction(method, adjust)))+
#   geom_hline(yintercept = 0.05, color=1, linetype=2)+
#   theme(legend.title = element_blank())+
#   ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))
# ggsave(pp, filename = paste0('plot/', copula_distr, '_', choose_model,'.png'))



# null1.1
copula_distr = 'zinegbin'

choose_model = 'null1.1'

network_option='none'

network_condition_number = 0

null1.1= lapply(c(100, 150, 200, 289), function(n)cbind(do.call(rbind, lapply(1:7,  function(part)
  get_summary_fp(choose_model,n=n, p=127, nreps=200, part=part, part_name, copula_distr=copula_distr,
                 network_option = network_option, network_condition_number = network_condition_number))),
  dist = copula_distr, network_option = network_option, network_condition_number = network_condition_number))
null1.1 = do.call(rbind, null1.1)
head(null1.1)



# pp = ggplot(data = null1.1, aes(x=n, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
#   coord_cartesian(ylim = c(0, min(max(null1.1$x), 0.2)))+
#   geom_line(data=null1.1, aes(linetype = interaction(method, adjust)))+
#   geom_hline(yintercept = 0.05, color=1, linetype=2)+
#   theme(legend.title = element_blank())+
#   ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))
# ggsave(pp, filename = paste0('plot/', copula_distr, '_', choose_model,'.png'))


null1.1_none = null1.1[null1.1$adjust=='none', ]
pp = ggplot(data = null1.1_none, aes(x=n, y=x, group=method, color=method, 
                                linetype = method))+
  coord_cartesian(ylim = c(0, min(max(null1.1_none$x), 0.2)))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  ylab('Type I Error')+
  scale_shape_manual(values=1:7)+
  scale_color_manual(values=cbPalette[1:7])+
  geom_hline(yintercept = 0.05, color=1, linetype=2)+
  theme(legend.title = element_blank())+
  scale_linetype_manual(values=1:7)+
  # ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
  scale_x_continuous(breaks=c(100, 150, 200, 250))
ggsave(pp, filename = paste0('plot/type_I_null1.1.png'),width =11.7*1.7, height= 11.7*1.7, units = 'cm')

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
  scale_x_continuous(breaks=c(100, 150, 200, 250))
ggsave(pp_time, filename = paste0('plot/time_null1.1.png'),width =11.7*1.7, height= 11.7*1.7, units = 'cm')



# # null1 models
# null1.2 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null1',n=289, p=127, nreps=100, part=part, part_name)))
# null1.2$x = round(null1.2$x, 4)
# kable(null1.2, 'latex', booktab=T, col.names = c(names(null1.2)[1:5], 'False Positive Rate'))%>%
#   kable_styling()
# 
# 
# null1.2 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null1',n=100, p=127, nreps=100, part=part, part_name)))
# null1.2$x = round(null1.2$x, 4)
# kable(null1.2, 'latex', booktab=T, col.names = c(names(null1.2)[1:5], 'False Positive Rate'))%>%
#   kable_styling()



# null2 model
copula_distr = 'none'

choose_model = 'null2'
network_option='none'

network_condition_number = 0

# null2 = lapply(c(100, 200), function(p){
#                   null2= lapply(c(100, 200, 500), function(n)
#                     cbind(do.call(rbind, lapply(1:7,  function(part)
#                                                       get_summary_fp(choose_model,n=n, p=p, nreps=200, part=part, part_name, 
#                                                                      copula_distr=copula_distr,network_option = network_option, 
#                                                                      network_condition_number = network_condition_number))),
#                           dist = copula_distr, 
#                           network_option = network_option, 
#                           network_condition_number = network_condition_number))
# 
#                   null2 = do.call(rbind, null2)
#                   null2
#                 })
# 
# null2 = do.call(rbind, null2)
# head(null2)
# 
# pp1 = ggplot(data = null2, aes(x=n, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
#   coord_cartesian(ylim = c(0, min(max(null2$x), 0.2)))+
#   geom_line(data=null2, aes(linetype = interaction(method, adjust)))+
#   geom_hline(yintercept = 0.05, color=1, linetype=2)+
#   theme(legend.title = element_blank())+
#   ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
#   facet_wrap(~p)
# pp2 = ggplot(data = null2, aes(x=p, y=x, group=interaction(method,adjust), color=interaction(method, adjust)))+
#   coord_cartesian(ylim = c(0, min(max(null2$x), 0.2)))+
#   geom_line(data=null2, aes(linetype = interaction(method, adjust)))+
#   geom_hline(yintercept = 0.05, color=1, linetype=2)+
#   theme(legend.title = element_blank())+
#   ggtitle(paste0('dist: ', copula_distr, ', ', choose_model))+
#   facet_wrap(~n)
# 
# ggsave(pp1, filename = paste0('plot/', copula_distr, '_', choose_model,'n.png'))
# ggsave(pp2, filename = paste0('plot/', copula_distr, '_', choose_model,'p.png'))

## null2 p=200

null2 = lapply(c(200), function(p){
  null2= lapply(c(100, 200, 500), function(n)
    cbind(do.call(rbind, lapply(1:7,  function(part)
      get_summary_fp(choose_model,n=n, p=p, nreps=200, part=part, part_name, 
                     copula_distr=copula_distr,network_option = network_option, 
                     network_condition_number = network_condition_number))),
      dist = copula_distr, 
      network_option = network_option, 
      network_condition_number = network_condition_number))
  
  null2 = do.call(rbind, null2)
  null2
})

null2 = do.call(rbind, null2)
head(null2)

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

ggsave(pp1, filename = paste0('plot/type_I_null2.png'), width = 11.7*1.7, height = 11.7*1.7, units = 'cm')

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
ggsave(pp_time, filename = paste0('plot/time_null2.png'),width =11.7*1.7, height= 11.7*1.7, units = 'cm')



library(ggpubr)
null_plot = ggarrange(pp,pp1,labels = c('A', 'B'), common.legend = T, legend = 'right' )
ggsave(null_plot, filename = paste0('plot/type1error_null1.1_null2.png'), width = 11.7*1.7, height = 7*1.7, units = 'cm')


# null2.1 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null2',n=100, p=200, nreps=100, part=part, part_name)))
# null2.1$x = round(null2.1$x, 4)
# kable(null2.1, 'latex', booktab=T, col.names = c(names(null2.1)[1:5], 'False Positive Rate'))%>%
#   kable_styling()
# 
# null2.2 = do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null2',n=200, p=200, nreps=100, part=part, part_name)))
# null2.2$x = round(null2.2$x, 4)
# kable(null2.2, 'latex', booktab=T, col.names = c(names(null2.2)[1:5], 'False Positive Rate'))%>%
#   kable_styling()
# 
# 
# null2.3= do.call(rbind, lapply(1:7,  function(part)get_summary_fp('null2',n=500, p=200, nreps=100, part=part, part_name)))
# null2.3$x = round(null2.3$x, 4)
# kable(null2.3, 'latex', booktab=T, col.names = c(names(null2.3)[1:5], 'False Positive Rate'))%>%
#   kable_styling()

#######################################################################
# get ROC summary for alternative models, dist_data files
get_summary_roc = function(choose_model, n, p, nreps, part, part_name, copula_distr, network_option, network_condition_number){
  # all models we compare the inverse cov graph ROC
  # to have an 'averaged' ROC curve, we may take average of fp and tp at each grid? (since we are using the same sequence of lambda)
  res_for_20 = lapply(1:50, function(run_rep){
    load_file = paste0(filepath,
                       '/dist_data/', copula_distr,'/', network_option, '/cond_', network_condition_number,'/', choose_model, '/res_n_', n, '_p_', p, '_', choose_model, 
                       '_nreps_', nreps, '_run_rep', run_rep,'_part_', part, '.RData')
    #print(load_file)
    tryCatch(load(load_file), error=function(e)print(load_file))
    list(n=n, p=p, model = choose_model, method = part_name[[part]],rep=run_rep,
         tp = roc$ROC_inv$tp, fp = roc$ROC_inv$fp, auc = roc$ROC_inv$AUC,
         tp_cov = roc$ROC_cov$tp, fp_cov = roc$ROC_cov$fp, auc_cov = roc$ROC_cov$AUC, time = roc$time)

  })
  # average time information
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
  
  # cov ROC results
  tps_cov = sapply(res_for_20, `[`, 'tp_cov')
  if(length(unique(sapply(tps_cov, FUN = length)))>1) stop(paste('check lambda sequence, tp_cov not same length across replications: n',n,'p',p, 'part', part))
  fps_cov = sapply(res_for_20, `[`, 'fp_cov')
  if(length(unique(sapply(fps_cov, FUN = length)))>1) stop(paste('check lambda sequence, fp_cov not same length across replications: n',n,'p',p, 'part', part))
  
  ave_tp_cov = rowMeans(do.call(cbind, tps_cov))
  ave_fp_cov = rowMeans(do.call(cbind, fps_cov))
  if(is.null(do.call(cbind, sapply(res_for_20, `[`, 'auc_cov')))){ave_auc_cov = NA}else{ave_auc_cov =  mean(do.call(cbind, sapply(res_for_20, `[`, 'auc_cov')))}
  
  ave_res_cov = list(n=n, p=p, model = choose_model, method = part_name[[part]],
                 ave_tp = ave_tp_cov, ave_fp = ave_fp_cov, ave_auc = ave_auc_cov, ave_time = ave_time)
  
  
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

# if to plot all methods together:
output_plot_all = function(choose_model, n, p, nreps, part_name, copula_distr, network_option, network_condition_number, target){
  data = do.call(rbind, lapply(1:7, function(part){
    roc_res = get_summary_roc(choose_model, n, p, nreps, part, part_name, copula_distr, network_option, network_condition_number)
    plot_roc(roc_res, target)$data
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


### first get results for 50 replicates
#alt1

{
target = 'inv' # or 'cov'
copula_distr = 'zinegbin'
choose_model = 'alt1'
network_option='erdos_renyi'
network_condition_number = 100
}


{
 target = 'inv' # or 'cov'
copula_distr = 'zinegbin'
choose_model = 'alt1'
network_option='erdos_renyi'
network_condition_number = 1000 
}

{
  target = 'inv' # or 'cov'
  copula_distr = 'zinegbin'
  choose_model = 'alt1'
  network_option='chain_small'
  network_condition_number = 0 
}

{
  target = 'inv' # or 'cov'
  copula_distr = 'zinegbin'
  choose_model = 'alt1'
  network_option='chain_large'
  network_condition_number = 0 
}

if(T){
  data = lapply(c(100, 150,200), function(n){
    output_plot_all(choose_model,n=n, p=127, nreps=200, part_name, copula_distr,network_option, network_condition_number, target)$data
  })
  data = do.call(rbind, data)
  data$n = as.factor(data$n)
  
  tmpdata = data[data$group=='mean' & data$method != 'ReBoot',] # for
  plot_all = ggplot(tmpdata[tmpdata$n %in% c(100,150,200,289),], aes(x=fp, y=tp, group=method, color = method,linetype=method))+
    geom_line(size=1.2)+
    geom_abline(slope=1, intercept=0,linetype = 3)+
    xlab('False Positive')+
    ylab('True Positive')+
    #ggtitle(paste('dist', copula_distr, 'n',nn, data$model[1]))+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    facet_wrap(~n, labeller = label_both, nrow = 2, ncol=2)+
    scale_color_manual(values = cbPalette[2:7])+
    scale_linetype_manual(values = c(2:7))
  ggsave(plot_all, filename = paste0('plot/(200 reps)ROC_alt1_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')
  
  tmpdata2 = with(tmpdata, aggregate(time, by = list(n=n, p=p, method=method), FUN=mean, ))
  for(i in c(1))tmpdata2[,i] = as.numeric(levels(tmpdata2[,i]))[tmpdata2[,i]]
  plot_time = ggplot(tmpdata2[tmpdata2$n %in% c(100,150,200,289),], aes(x=n, y=log10(x), group=method, color = method,linetype=method))+
    geom_line(size=1)+
    geom_point(size=2, shape=18)+
    scale_shape_manual(values=2:7)+
    xlab('n')+
    ylab('log10(time) (seconds)')+
    theme(legend.position = 'right')+
    theme(legend.title = element_blank())+
    scale_color_manual(values = cbPalette[2:7])+
    scale_linetype_manual(values = c(2:7))+
    scale_x_continuous(breaks=c(100, 150, 200, 250))
  ggsave(plot_time, filename = paste0('plot/time_alt1_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')
  
} # generate the plots


# alt2
{
  target = 'inv' # or 'cov'
  copula_distr = 'none'
  choose_model = 'alt2'
  network_option='erdos_renyi'
  network_condition_number = 100
}


{
  target = 'inv' # or 'cov'
  copula_distr = 'none'
  choose_model = 'alt2'
  network_option='erdos_renyi'
  network_condition_number = 1000 
}

{
  target = 'inv' # or 'cov'
  copula_distr = 'none'
  choose_model = 'alt2'
  network_option='chain_small'
  network_condition_number = 0 
}

{
  target = 'inv' # or 'cov'
  copula_distr = 'none'
  choose_model = 'alt2'
  network_option='chain_large'
  network_condition_number = 0 
}


if(T){
  data = lapply(c(100,200, 500), function(n){
  output_plot_all(choose_model,n=n, p=200, nreps=200, part_name, copula_distr,network_option, network_condition_number, target)$data
})
data = do.call(rbind, data)
data$n = as.factor(data$n)

tmpdata = data[data$group=='mean' & data$method != 'ReBoot',]
plot_all = ggplot(tmpdata[tmpdata$n %in% c(100,200,500),], aes(x=fp, y=tp, group=method, color = method,linetype=method))+
  geom_line(size=1.2)+
  geom_abline(slope=1, intercept=0,linetype = 3)+
  xlab('False Positive')+
  ylab('True Positive')+
  #ggtitle(paste('dist', copula_distr, 'n',nn, data$model[1]))+
  theme(legend.position = 'right')+
  theme(legend.title = element_blank())+
  facet_wrap(~n, labeller = label_both, nrow = 2, ncol=2)+
  scale_color_manual(values = cbPalette[2:7])+
  scale_linetype_manual(values = c(2:7))
ggsave(plot_all, filename = paste0('plot/(50 reps)ROC_alt2_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')

tmpdata2 = with(tmpdata, aggregate(time, by = list(n=n, p=p, method=method), FUN=mean, ))
for(i in c(1))tmpdata2[,i] = as.numeric(levels(tmpdata2[,i]))[tmpdata2[,i]]
plot_time = ggplot(tmpdata2[tmpdata2$n %in% c(100,200,500),], aes(x=n, y=log10(x), group=method, color = method,linetype=method))+
  geom_line(size=1)+
  geom_point(size=2, shape=18)+
  scale_shape_manual(values=2:7)+
  xlab('n')+
  ylab('log10(time) (seconds)')+
  theme(legend.position = 'right')+
  theme(legend.title = element_blank())+
  scale_color_manual(values = cbPalette[2:7])+
  scale_linetype_manual(values = c(2:7))+
  scale_x_continuous(breaks=c(100, 150, 200, 250))
ggsave(plot_time, filename = paste0('plot/time_alt2_',target, '_',network_option, network_condition_number,'.png'), width = 11.7*1.7, height = 8*1.7, units = 'cm')

}


# 
# data3 = lapply(c(100, 200, 500), function(n){
#   tmp = lapply(c(20, 50, 100, 150, 200), function(p){
#     output_plot_all('alt3',n=n, p=p, nreps=100, part_name, copula_distr)$data
#   })
#   do.call(rbind, tmp)
# })
# data3 = do.call(rbind, data3)
# data3$p = as.factor(data3$p)
# 
# for(nn in c(100, 200, 500)){
#   tmpdata3 = data3[data3$group=='mean'& data3$n==nn,]
#   plot_all = ggplot(tmpdata3, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
#     geom_line(size=1.2)+
#     geom_abline(slope=1, intercept=0,linetype = 3)+
#     xlab('False Positive')+
#     ylab('True Positive')+
# #    ggtitle(paste('sample size n =', nn))+
#     theme(legend.position = 'right', axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+
#     scale_color_manual(values=cbPalette[1:7])+
#     scale_linetype_manual(values=1:7)+
#     
#     facet_wrap(~p, labeller = label_both)
#   ggsave(plot_all, filename = paste0('plot/ROC_alt3_n_',nn,'.png'), width = 11.7*1.9, 
#          height = 8*1.9, units = 'cm')
#   
# }
# 
# for(pp in c(20,50,100,150,200)){
#   tmpdata3 = data3[data3$group=='mean'& data3$p==pp,]
#   plot_all = ggplot(tmpdata3, aes(x=fp, y=tp, group=method, color = method,linetype=method))+
#     geom_line(size=1.2)+
#     geom_abline(slope=1, intercept=0,linetype = 3)+
#     xlab('False Positive')+
#     ylab('True Positive')+
#    # ggtitle(paste('number of genes p =', pp))+
#     theme(legend.position = 'right', axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+
#     scale_color_manual(values=cbPalette[1:7])+
#     scale_linetype_manual(values=1:7)+
#     facet_wrap(~n, labeller = label_both)
#   ggsave(plot_all, filename = paste0('plot/ROC_alt3_p_',pp,'.png'), width = 11.7*1.9,
#          height = 4*1.9, units = 'cm')
#   
# }
# 

# output_plot_all('alt2',n=200, p=200, nreps=100, part_name, copula_distr)
# output_plot_all('alt2',n=500, p=200, nreps=100, part_name, copula_distr)

# inspect CClasso problem
# nn=100
# for(pp in c(20, 50, 100, 150, 200)){
#   plot_all =output_plot_all('alt2',n=nn, p=pp, nreps=100, part_name, copula_distr)$plots
#   ggsave(plot_all, filename = paste0('plot/', copula_distr, '_n_',nn,'_p_',pp,'alt2.png'))
# }


# output_plot_all('alt2',n=200, p=100, nreps=100, part_name, copula_distr)$plots
# output_plot_all('alt2',n=200, p=150, nreps=100, part_name, copula_distr)$plots
# output_plot_all('alt2',n=200, p=200, nreps=100, part_name, copula_distr)$plots



# summarize dataset sparsity level, output mean sparsity level of the whole dataset (across all cells and taxa)
get_data_sparsity = function(choose_model, n, p, nreps, copula_distr, network_option, network_condition_number){
  load_file = paste0(filepath,
                     '/dist_data/', copula_distr,'/', network_option, '/cond_', network_condition_number,'/n_', n, '_p_', p, '_', choose_model, 
                     '_nreps_', nreps, '_data_rep.RData')

  tryCatch(load(load_file), error = function(e)print(load_file))
  c(sparsity = mean(data_rep[[1, 1]]==0), n=n, p=p, model = choose_model, copula_distr = copula_distr, network_option = network_option, network_condition_number = network_condition_number)
}

tmp1 = do.call(rbind, lapply(c(100, 150, 200, 289), function(n){
  get_data_sparsity('null1.1', n, 127, 200, 'zinegbin', 'none', 0)
}))

tmp2 = do.call(rbind, lapply(c(100, 200, 500), function(n){
  rbind(
    get_data_sparsity('null2', n, 100, 200, 'none', 'none', 0),
    get_data_sparsity('null2', n, 200, 200, 'none', 'none', 0)
  )
}))

tmp3 = do.call(rbind, lapply(c(100, 150, 200, 289), function(n){
  do.call(rbind, lapply(c(127), function(p){
    rbind(get_data_sparsity('alt1', n, p, 200, 'zinegbin', 'chain_small', 0),
          get_data_sparsity('alt1', n, p, 200, 'zinegbin', 'chain_large', 0),
          get_data_sparsity('alt1', n, p, 200, 'zinegbin', 'erdos_renyi', 100),
          get_data_sparsity('alt1', n, p, 200, 'zinegbin', 'erdos_renyi', 1000)
    )
  }))
}))

tmp4 = do.call(rbind, lapply(c(100, 200, 500), function(n){
  do.call(rbind, lapply(c(100, 200), function(p){
    rbind(get_data_sparsity('alt2', n, p, 200, 'none', 'chain_small', 0),
          get_data_sparsity('alt2', n, p, 200, 'none', 'chain_large', 0),
          get_data_sparsity('alt2', n, p, 200, 'none', 'erdos_renyi', 100),
          get_data_sparsity('alt2', n, p, 200, 'none', 'erdos_renyi', 1000)
    )
  }))
}))

sparsity_table = data.frame(rbind(tmp1, tmp2, tmp3, tmp4), stringsAsFactors = F)
for(i in c(1, 2, 3, 7)) sparsity_table[,i] = as.numeric(sparsity_table[,i])
head(sparsity_table) # for alt2 model, is log-normal so sparsity=0
write.table(sparsity_table, row.names = F, file=paste0(filepath, '/generated_data_mean_sparsity.txt'))
