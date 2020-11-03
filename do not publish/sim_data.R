rm(list=ls())

## load packages
library(glasso)
library(mvtnorm)
library(MASS)
library(gdata)
library(fields)
library(igraph)
library(RANN)
library(ggplot2)
library(censReg)

CodeFolder <- "../RCode/"
OutFolder <- "../Output/synData/"
source(paste0(CodeFolder,"lib/censoredGGM.R"))
source(paste0(CodeFolder,"lib/mixedCorr.R"))
source(paste0(CodeFolder,"lib/gcoda_20191005.R"))
source(paste0(CodeFolder,"lib/fitdistr.R"))

set.seed(102)

today <- "20191022"
DEBUG <- "sf"
q <- 0
nreps <- 20  # number of replications
mdistr <- 'normal'

p <- 60  #nrow(abundance_even)
n <- 100  #ncol(abundance_even)

if (DEBUG=='sf'){
  myModel <- model_pa(p+q, 1, const=0.1)
}
if (DEBUG=='gnp'){
  myModel <- model_gnm(p+q, 1*(p+q), const=0.1)
}
if (DEBUG=='hubnet'){
  myModel <- hubNet(p+q, 4)
}
if (DEBUG=='nn'){
  myModel <- nn.net(p+q, 5)
}
if (DEBUG=='sw'){
  myModel <- model_smallworld(p+q,2)
}
today <- paste0(today,"_", DEBUG,"_p",p,"_q",q,'_n',n,"_nreps", nreps, '_',mdistr,'_edge',sum(myModel$Adj)/2)

simul.mu <- runif(p+q,-0.5,2)
Sigma <- myModel$Sigma
pars <- list(n=n, p=p, q=q, mu = simul.mu, S=Sigma, Omega=solve(Sigma), graph=myModel$Adj, model=mdistr, lb=1, ub=10, nreps=nreps)

dataList <- lapply(1:nreps, function(j) data.simulator(pars, sid=j))
# save(pars, dataList, file=paste0(OutFolder, 'data_', today, '.rda'))

