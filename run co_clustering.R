rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(1)

setwd('U:\\GIT_models\\co_clustering')
sourceCpp('rcpp_func.cpp')
source('gibbs functions.R')
source('co_clustering_main.R')
dat=data.matrix(read.csv('fake data.csv',as.is=T))
n=read.csv('fake data n.csv',as.is=T)$x

ngibbs=1000
burnin=ngibbs/2
ngroup.loc=20
ngroup.spp=10
gamma.v=0.1
gamma.u=0.1
a.prior=b.prior=0.1

res=co_clustering(dat=dat,ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp,ngibbs=ngibbs,burnin=burnin,n=n,
                  gamma.v=gamma.v,gamma.u=gamma.u,a.prior=a.prior,b.prior=b.prior)
