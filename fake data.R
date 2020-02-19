rm(list=ls(all=TRUE))
set.seed(30)

ngroup.loc=6
ngroup.spp=rpois(ngroup.loc,lambda=3)+1

#get parameters
tmp=runif(ngroup.loc)
theta.true=theta=tmp/sum(tmp)
phi=psi=matrix(NA,ngroup.loc,max(ngroup.spp))
for (i in 1:ngroup.loc){
  tmp=runif(ngroup.spp[i])
  ind=1:ngroup.spp[i]
  phi[i,ind]=tmp/sum(tmp)
  psi[i,ind]=rbeta(ngroup.spp[i],0.1,0.1)
}
phi.true=phi
psi.true=psi

#get latent variables
nloc=1000
tmp=rmultinom(1,size=nloc,prob=theta)
tmp1=rep(1:ngroup.loc,times=tmp)
z.true=z=sample(tmp1,nloc); 

nspp=50
w=matrix(NA,ngroup.loc,nspp)
for (i in 1:ngroup.loc){
  phi.tmp=phi[i,]
  tmp=rmultinom(1,size=nspp,prob=phi.tmp[!is.na(phi.tmp)])  
  tmp1=rep(1:ngroup.spp[i],times=tmp)
  w[i,]=sample(tmp1,nspp)
}
w.true=w

#generate data
n=rbinom(nloc,size=3,prob=0.6)+1; table(n)
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  for (j in 1:nspp){
    zsel=z[i]
    wsel=w[zsel,j]
    y[i,j]=rbinom(1,size=n[i],prob=psi[zsel,wsel])    
  }
}
image(y)

setwd('U:\\GIT_models\\co_clustering')
write.csv(y,'fake data.csv',row.names=F)
write.csv(n,'fake data n.csv',row.names=F)