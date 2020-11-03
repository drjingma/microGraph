
## Generate parameters
jid=commandArgs(trailingOnly=T)[1]
jid=as.numeric(jid)

CodeFolder <- "../RCode/"
OutFolder <- "../Output/"
source(paste0(CodeFolder,"lib/func_libs.R"))
source(paste0(CodeFolder,"lib/mixedCorr.R"))
source(paste0(CodeFolder,"lib/gcoda_20191005.R"))

filename <- "data_20191022_gnp_p60_q0_n100_nreps20_normal_edge60"
load(paste0(OutFolder,'synData/',filename, '.rda'))
n <- pars$n # sample size
p <- pars$p  # number of microbes
q <- pars$q  # number of metabolites
Omega <- pars$graph 
lambda_vec <- c(seq(0.01,0.1,0.01),seq(0.1,3,0.1)) * sqrt(log(p+q)/n)
lambda_vec <- sort(lambda_vec, decreasing = TRUE)

my.iterations <- function(iters){
  dd <- dataList[[iters]]
  res <- compare.roc(dd, graph=pars$graph,lambda=lambda_vec,OTU.only=T)
  return(res)
}

cat("... Run the main function ... \n")
        
res <- tryCatch({my.iterations(jid)}, error=function(e) NULL)
filename <- gsub('data', 'roc', filename)

if (!is.null(res)) {save(res, file=paste0(OutFolder, filename, '_',jid,'.rda'))}


