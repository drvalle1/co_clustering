plot(res$llk,type='l',ylim=quantile(res$llk,c(0.01,1)))

zestim=res$z[nrow(res$z),]
k=data.frame(ztrue=z.true,zestim=zestim)
k1=table(k); k1
ind.loc=numeric()
for (i in 1:nrow(k1)){
  tmp=which(k1[i,]==max(k1[i,]))
  ind.loc=c(ind.loc,colnames(k1)[tmp])
}
ind.loc=as.numeric(ind.loc)
k1[,ind.loc]

theta.estim=res$theta[nrow(res$theta),]
plot(theta.estim,type='h')
theta.estim[ind.loc];theta.true

rango=range(c(theta.true,theta.estim[ind.loc]))
plot(theta.true,theta.estim[ind.loc],xlim=rango,ylim=rango)
lines(rango,rango)
#---------------------------------------
westim=res$w[nrow(res$w),]
k=data.frame(wtrue=w.true,westim=westim)
k1=table(k); k1
ind.spp=numeric()
for (i in 1:nrow(k1)){
  tmp=which(k1[i,]==max(k1[i,]))
  ind.spp=c(ind.spp,colnames(k1)[tmp])
}
k1[,ind.spp]
ind.spp=as.numeric(ind.spp)

phi.estim=res$phi[nrow(res$phi),]
round(phi.estim[ind.spp],3); round(phi.true,3)
plot(phi.estim,type='h')

rango=range(c(phi.estim[ind.spp],phi.true))
plot(phi.true,phi.estim[ind.spp],xlim=rango,ylim=rango)
lines(rango,rango)
#----------------------------------------
tmp=res$psi[nrow(res$psi),]
psi0=matrix(tmp,20,10)
psi=psi0[,ind.spp]
psi1=psi[ind.loc,]
rango=c(0,1)
plot(psi1,data.matrix(psi.true),xlim=rango,ylim=rango)
lines(rango,rango)
