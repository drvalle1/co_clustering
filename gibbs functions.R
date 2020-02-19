#' Sample psi parameters
#' 
#' Sample the presence probability for each location and species group (psi) from 
#' it's full conditional distribution

#' @param dat matrix has L rows (e.g., locations) and S columns (e.g., species),
#'            containing the number of observation occasions that a particular species was observed 
#'            in a particular location 
#' @param nmat  matrix has L rows (e.g., locations) and S columns (e.g., species), 
#'              containing the number of observation occasions for each species and location
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp maximum number of species groups (KS)
#' @param z vector of length L, storing the current membership of each location
#' @param w vector of length S, storing the current membership of each species
#' @return a KL x KS matrix of psi parameters
#' @export

sample.psi=function(z,w,dat,ngroup.loc,ngroup.spp,nmat,a.prior,b.prior,lo,hi){
  #summarize the data by calculating how many observations were assigned to each group
  #for which dat[i,j]=1 and for which dat[i,j]=0
  tmp=getql(z=z-1,w=w-1,dat=dat,ngrloc=ngroup.loc,ngrspp = ngroup.spp,nmat=nmat)
  
  #generate psi from beta distribution
  tmp1=rbeta(ngroup.loc*ngroup.spp,tmp$nqly+a.prior,tmp$nqlny+b.prior)
  tmp1[tmp1>hi]=hi
  tmp1[tmp1<lo]=lo
  matrix(tmp1,ngroup.loc,ngroup.spp)
}

#' Sample theta and vk parameters
#' 
#' Sample vk parameters from their full conditional distributions and convert 
#' these vk parameters into theta parameters
#'
#' @param ngroup.loc maximum number of groups for locations (KL)
#' @param gamma.v truncated stick-breaking prior parameter for the 
#'                location groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param burnin number of MCMC samples that are going to be thrown out as 
#'               part of the burn-in phrase
#' @param gibbs.step current iteration of the gibbs sampler
#' @param theta vector of size KL containing the current estimate of theta (i.e., probability of each location group)                  
#' @param psi KL x KS matrix  containing the current estimate of psi
#' @param z this is a vector of length L, storing the current membership of each location
#' @return this function returns a list containing 4 items (theta, vk, psi, and z)
#' @export

sample.theta=function(ngroup.loc,gamma.v,burnin,gibbs.step,theta,psi,z,lo,hi,w,phi){
  #re-order thetas in decreasing order if we are still in the burn-in phrase. 
  #based on this re-ordering, re-order z's and psi's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    psi=psi[ind,]
    w=w[ind,]
    phi=phi[ind,]
    
    #get z.new
    z.new=z; z.new[]=NA
    for (i in 1:ngroup.loc){
      cond=z==ind[i]
      z.new[cond]=i
    }
    z=z.new
  }
  
  #get the number of locations in each group
  nk=rep(0,ngroup.loc)
  tmp=table(z)
  nk[as.numeric(names(tmp))]=tmp
  
  #sample vk from a beta distribution
  ind=ngroup.loc:1
  invcumsum=cumsum(nk[ind])[ind]
  vk=rbeta(ngroup.loc,nk+1,invcumsum-nk+gamma.v)
  # vk[vk>hi]=hi #to avoid numerical issues
  # vk[vk<lo]=lo
  vk[ngroup.loc]=1
  
  #convert from vk to theta
  theta=convertSBtoNormal(vk)
  
  #output vk, theta, z, and psi
  list(vk=vk,theta=theta,z=z,psi=psi,w=w,phi=phi)
}

#' Sample phi and uk parameters
#' 
#' Sample uk parameters from their full conditional distributions and 
#' then calculate the implied phi parameters
#'
#' @param ngroup.spp maximum number of groups for species (KS)
#' @param gamma.u the truncated stick-breaking prior parameter for 
#'                species groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param burnin number of MCMC samples that are going to be thrown out as 
#'               part of the burn-in phrase
#' @param gibbs.step current iteration of the gibbs sampler
#' @param phi vector of size KS containing the current estimate of phi (i.e., probability of each species group)                  
#' @param psi matrix of size KL x KS containing the current estimate  of psi
#' @param w vector of length S, storing the current membership of each species
#' @return this function returns a list with 4 items (phi, uk, w, psi)
#' @export
#' 

sample.phi=function(ngroup.loc,ngroup.spp,gamma.u,burnin,gibbs.step,phi,psi,w,lo,hi){
  #re-order phi in decreasing order if we are still in burn-in phase 
  #Based on this re-ordering, re-order w's and psi's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    w.new=w;
    for (i in 1:ngroup.loc){
      ind=order(phi[i,],decreasing=T)  
      phi[i,]=phi[i,ind]
      psi[i,]=psi[i,ind]
      
      #get w.new
      for (j in 1:ngroup.spp){
        cond=w[i,]==ind[j]
        w.new[i,cond]=j
      }
    }
    w=w.new
  }
  
  #calculate the number of species in each group and location
  uk.store=phi=mk=matrix(0,ngroup.loc,ngroup.spp)
  for (i in 1:ngroup.loc){
    tmp=table(w[i,])
    mk[i,as.numeric(names(tmp))]=tmp
    
    #sample uk from a beta distribution
    ind=ngroup.spp:1
    invcumsum=cumsum(mk[i,ind])[ind]
    uk=rbeta(ngroup.spp,mk[i,]+1,invcumsum-mk[i,]+gamma.u)
    uk[uk>hi]=hi #for numerical issues
    # uk[uk<lo]=lo #for numerical issues
    uk[ngroup.spp]=1
    uk.store[i,]=uk
    
    #convert uk to phi
    phi[i,]=convertSBtoNormal(uk)
  }

  #return uk, phi, w, and psi
  list(uk=uk.store,phi=phi,w=w,psi=psi)
}
#--------------------------
#' Sample gamma.u
#' 
#' Sample the TSB prior parameter for species groups (gamma.u) from its full conditional distribution
#'
#' @param uk vector of size KS with probabilities
#' @param ngroup.spp maximum number of species groups (KS)
#' @param gamma.possib this is a vector containing the possible values that gamma.u can take
#' @return this function returns a real number (gamma.u) 
#' @export
#' 
sample.gamma.u=function(uk,gamma.possib,ngroup.spp,ngroup.loc){
  #calculate the stick-breaking probabilities for different values of gamma.u
  ngamma=length(gamma.possib)
  soma=sum(log(1-uk[,-ngroup.spp]))
  k=ngroup.loc*(ngroup.spp-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #to check this code: sum(dbeta(v[-ngroup.spp],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize these probabilities to draw from categorical distribution
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
#--------------------------
#' Sample gamma.v
#' 
#' Sample the TSB prior parameter for location groups (gamma.v) from its full conditional distribution
#'
#' @param vk vector of size KL with probabilities
#' @param ngroup.loc maximum number of location groups (KL)
#' @param gamma.possib vector containing the possible values that gamma.u can take
#' @return this function returns a real number (gamma.v) 
#' @export
#' 

sample.gamma.v=function(vk,gamma.possib,ngroup.loc){
  #calculate the stick-breaking probabilities for different values of gamma.v
  ngamma=length(gamma.possib)
  soma=sum(log(1-vk[-ngroup.loc]))
  k=(ngroup.loc-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #to check this code: sum(dbeta(v[-ngroups],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize these probabilities to draw from categorical distribution
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
#--------------------------
#' Sample z
#' 
#' Sample the vector z containing the membership of each location
#'
#' @param ltheta equal to log(theta)
#' @param dat matrix has L rows (e.g., locations) and S columns (e.g., species),
#'            containing the number of observation occasions that a particular species was observed 
#'            in a particular location 
#' @param nmat  matrix has L rows (e.g., locations) and S columns (e.g., species), 
#'              containing the number of observation occasions for each species and location
#' @param lpsi equal to log(psi)
#' @param l1mpsi equal to log(1-psi)
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp maximum number of species groups (KS)
#' @param nloc total number of locations (L)
#' @param nspp total number of species (S)
#' @param w vector of length S, storing the current membership of each species
#' @param z vector of length L, storing the current membership of each locations
#' @return this function returns a vector of length L containing z 
#' @export
#' 

sample.z=function(ltheta,dat,nmat,lpsi,l1mpsi,ngroup.loc,ngroup.spp,nloc,nspp,w,z,
                  a.prior,b.prior,constant){
  #calculation of log probability for groups that already exist
  lprob.exist=matrix(NA,nloc,ngroup.loc)
  for (i in 1:ngroup.loc){
    lpsi1=matrix(lpsi[i,w[i,]],nloc,nspp,byrow=T)
    l1mpsi1=matrix(l1mpsi[i,w[i,]],nloc,nspp,byrow=T)
    lprob.exist[,i]=ltheta[i]+rowSums(dat*lpsi1+(nmat-dat)*l1mpsi1)
  }

  #sample z
  for (i in 1:nloc){
    max.z=max(z)
    lprob=lprob.exist[i,1:max.z]
    if (max.z<ngroup.loc) {
      tmp=ltheta[max.z+1]+sum(lgamma(dat[i,]+a.prior)+lgamma(nmat[i,]-dat[i,]+b.prior)-
                              lgamma(a.prior+b.prior+nmat[i,]))+constant
      lprob=c(lprob,tmp)
    }
    
    #get normalized probs
    tmp1=lprob-max(lprob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1
    
    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    z[i]=ind1
  }
  z
}
#--------------------------
#' Sample w
#' 
#' Sample the vector w containing the membership of each species
#'
#' @param lphi equal to log(phi)
#' @param dat matrix has L rows (e.g., locations) and S columns (e.g., species),
#'            containing the number of observation occasions that a particular species was observed 
#'            in a particular location 
#' @param nmat  matrix has L rows (e.g., locations) and S columns (e.g., species), 
#'              containing the number of observation occasions for each species and location
#' @param lpsi equal to log(psi)
#' @param l1mpsi equal to log(1-psi)
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp thismaximum number of species groups (KS)
#' @param nloc total number of locations (L)
#' @param nspp total number of species (S)
#' @param w vector of length S, storing the current membership of each species
#' @param z vector of length L, storing the current membership of each locations
#' @return this function returns a vector of length S containing w 
#' @export
#' 
sample.w=function(lphi,dat,nmat,lpsi,l1mpsi,ngroup.spp,ngroup.loc,nloc,nspp,w,z,a.prior,b.prior){

  for (i in 1:ngroup.loc){
    cond=z==i
    soma.cond=sum(cond)
    
    if (soma.cond>0){ #if this location group exists
      #calculate log probability for groups that already exist
      lprob.exist=matrix(NA,nspp,ngroup.spp)
      for (j in 1:ngroup.spp){
        lpsi1=matrix(lpsi[i,j],soma.cond,nspp)
        l1mpsi1=matrix(l1mpsi[i,j],soma.cond,nspp)
        lprob.exist[,j]=lphi[i,j]+colSums(dat[cond,]*lpsi1+(nmat[cond,]-dat[cond,])*l1mpsi1)
      }
      
      #sample w 
      for (j in 1:nspp){
        max.w=max(w[i,])  
        lprob=lprob.exist[j,1:max.w]
        if (max.w<ngroup.spp){
          soma.y=sum(dat[cond,j])
          soma.n=sum(nmat[cond,j])
          tmp=lphi[i,max.w+1]+lgamma(soma.y+a.prior)+lgamma(soma.n+soma.y+b.prior)-
              lgamma(soma.n+a.prior+b.prior)
          lprob=c(lprob,tmp)
        }
        
        #get normalized probs
        tmp1=lprob-max(lprob) #for numerical stability
        tmp2=exp(tmp1) #exponentiate log probability
        prob=tmp2/sum(tmp2) #normalize to sum to 1
        
        #draw from multinomial distrib
        ind=rmultinom(1,size=1,prob=prob)
        ind1=which(ind==1)
        w[i,j]=ind1
      }
    }
  }
 w
}