###########################################################################################
## Title: Statistical Inference for Environmental Efficiency and Productivity Dynamics
## Authors: Valentin Zelenyuk and Shirong Zhao
## Date: Aug 29, 2026
##############
## This function estimates the MPI under weak disposability
## x1, x2, is n times p matrix for inputs 
## g1, g2 is n times qg matrix for desirable outputs
## b1, b2 is n times qb matrix for undesirable outputs
###########################################################################################
coverage.mpi.without.CO2.weighted<-function(x1=x1,y1=y1,x2=x2,y2=y2,
                                     XDIREC1=XDIREC1,YDIREC1=YDIREC1,
                                     XDIREC2=XDIREC2,YDIREC2=YDIREC2,
                                     H1=H1,H2=H2,L=50) {
  
  p=ncol(x1)
  q=ncol(y1)
  n=nrow(x1)
  na=floor(n/2)
  nb=n-na
  kappa=2/(p+q+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  # evaluate efficiency using CRS-DEA
  
  c11= 1-Benchmarking::dea.direct(X=x1,Y=y1,
                                  DIREC=cbind(XDIREC1,YDIREC1),
                                  XREF=x1,YREF=y1,
                                  ORIENTATION="in-out",
                                  RTS="crs")$eff[,1]
  c12= 1-Benchmarking::dea.direct(X=x1,Y=y1,
                                  DIREC=cbind(XDIREC1,YDIREC1),
                                  XREF=x2,YREF=y2,
                                  ORIENTATION="in-out",
                                  RTS="crs")$eff[,1]
  
  c21= 1-Benchmarking::dea.direct(X=x2,Y=y2,
                                  DIREC=cbind(XDIREC2,YDIREC2),
                                  XREF=x1,YREF=y1,
                                  ORIENTATION="in-out",
                                  RTS="crs")$eff[,1]
  c22= 1-Benchmarking::dea.direct(X=x2,Y=y2,
                                  DIREC=cbind(XDIREC2,YDIREC2),
                                  XREF=x2,YREF=y2,
                                  ORIENTATION="in-out",
                                  RTS="crs")$eff[,1]
  #
  delta=-0.5*(log(1+c21)-log(1+c11)+log(1+c22)-log(1+c12))
  #
  U1=(1+c21)*as.vector(H2)
  U2=(1+c22)*as.vector(H2)
  U3=(1+c11)*as.vector(H1)
  U4=(1+c12)*as.vector(H1)
  U5=as.vector(H2)
  U6=as.vector(H1)
  mu1=mean(U1)
  mu2=mean(U2)
  mu3=mean(U3)
  mu4=mean(U4)
  mu5=mean(U5)
  mu6=mean(U6)
  Delta.logs=-0.5*(log(mu1)-log(mu3)+log(mu2)-log(mu4))+log(mu5)-log(mu6)
  Delta=exp(Delta.logs)
  #
  Ti=cbind(U1,U2,U3,U4,U5,U6)
  Sigma=cov(Ti)
  gg1=-1/(2*mu1)
  gg2=-1/(2*mu2)
  gg3= 1/(2*mu3)
  gg4= 1/(2*mu4)
  gg5= 1/mu5
  gg6=-1/mu6
  gg=c(gg1,gg2,gg3,gg4,gg5,gg6)
  V.Sigma=t(gg)%*%Sigma%*%gg
  sig=sqrt(V.Sigma)
  sig=as.vector(sig)
  sig=sig*Delta
  if (is.nan(sig)){
    sig=0
  } 
  #####################################################
  # compute bias corrections via generalized jackknife:
  tbar0=0
  tbar=rep(0,n)
  ind=c(1:n)
  for (j in 1:L) {
    #
    ind1=sample(ind,size=n)
    x1.b=matrix(x1[ind1,],ncol=p)
    y1.b=matrix(y1[ind1,],ncol=q) 
    XDIREC1.b=matrix(XDIREC1[ind1,],ncol=p)
    YDIREC1.b=matrix(YDIREC1[ind1,],ncol=q) 
    H1.b=H1[ind1]
    #
    x2.b=matrix(x2[ind1,],ncol=p)
    y2.b=matrix(y2[ind1,],ncol=q) 
    XDIREC2.b=matrix(XDIREC2[ind1,],ncol=p)
    YDIREC2.b=matrix(YDIREC2[ind1,],ncol=q) 
    H2.b=H2[ind1]
    #
    x1a=matrix(x1.b[1:na,],ncol=p)
    y1a=matrix(y1.b[1:na,],ncol=q)   
    XDIREC1a=matrix(XDIREC1.b[1:na,],ncol=p)
    YDIREC1a=matrix(YDIREC1.b[1:na,],ncol=q)   
    H1a=H1.b[1:na]
    
    x1b=matrix(x1.b[(na+1):n,],ncol=p)
    y1b=matrix(y1.b[(na+1):n,],ncol=q)
    XDIREC1b=matrix(XDIREC1.b[(na+1):n,],ncol=p)
    YDIREC1b=matrix(YDIREC1.b[(na+1):n,],ncol=q)
    H1b=H1.b[(na+1):n]
    #
    x2a=matrix(x2.b[1:na,],ncol=p)
    y2a=matrix(y2.b[1:na,],ncol=q)   
    XDIREC2a=matrix(XDIREC2.b[1:na,],ncol=p)
    YDIREC2a=matrix(YDIREC2.b[1:na,],ncol=q)   
    H2a=H2.b[1:na]
    
    x2b=matrix(x2.b[(na+1):n,],ncol=p)
    y2b=matrix(y2.b[(na+1):n,],ncol=q)
    XDIREC2b=matrix(XDIREC2.b[(na+1):n,],ncol=p)
    YDIREC2b=matrix(YDIREC2.b[(na+1):n,],ncol=q)
    H2b=H2.b[(na+1):n]
    #
    c11a= 1-Benchmarking::dea.direct(X=x1a,Y=y1a,
                                     DIREC=cbind(XDIREC1a,YDIREC1a),
                                     XREF=x1a,YREF=y1a,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    c12a= 1-Benchmarking::dea.direct(X=x1a,Y=y1a,
                                     DIREC=cbind(XDIREC1a,YDIREC1a),
                                     XREF=x2a,YREF=y2a,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    
    c21a= 1-Benchmarking::dea.direct(X=x2a,Y=y2a,
                                     DIREC=cbind(XDIREC2a,YDIREC2a),
                                     XREF=x1a,YREF=y1a,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    c22a= 1-Benchmarking::dea.direct(X=x2a,Y=y2a,
                                     DIREC=cbind(XDIREC2a,YDIREC2a),
                                     XREF=x2a,YREF=y2a,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    #
    delta.a=-0.5*(log(1+c21a)-log(1+c11a)+log(1+c22a)-log(1+c12a))
    #
    U1a=(1+c21a)*as.vector(H2a)
    U2a=(1+c22a)*as.vector(H2a)
    U3a=(1+c11a)*as.vector(H1a)
    U4a=(1+c12a)*as.vector(H1a)
    U5a=as.vector(H2a)
    U6a=as.vector(H1a)
    mu1a=mean(U1a)
    mu2a=mean(U2a)
    mu3a=mean(U3a)
    mu4a=mean(U4a)
    mu5a=mean(U5a)
    mu6a=mean(U6a)
    Delta.a.logs=-0.5*(log(mu1a)-log(mu3a)+log(mu2a)-log(mu4a))+log(mu5a)-log(mu6a)
    Delta.a=exp(Delta.a.logs)
    #
    c11b= 1-Benchmarking::dea.direct(X=x1b,Y=y1b,
                                     DIREC=cbind(XDIREC1b,YDIREC1b),
                                     XREF=x1b,YREF=y1b,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    c12b= 1-Benchmarking::dea.direct(X=x1b,Y=y1b,
                                     DIREC=cbind(XDIREC1b,YDIREC1b),
                                     XREF=x2b,YREF=y2b,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    
    c21b= 1-Benchmarking::dea.direct(X=x2b,Y=y2b,
                                     DIREC=cbind(XDIREC2b,YDIREC2b),
                                     XREF=x1b,YREF=y1b,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    c22b= 1-Benchmarking::dea.direct(X=x2b,Y=y2b,
                                     DIREC=cbind(XDIREC2b,YDIREC2b),
                                     XREF=x2b,YREF=y2b,
                                     ORIENTATION="in-out",
                                     RTS="crs")$eff[,1]
    #
    delta.b=-0.5*(log(1+c21b)-log(1+c11b)+log(1+c22b)-log(1+c12b))
    #
    U1b=(1+c21b)*as.vector(H2b)
    U2b=(1+c22b)*as.vector(H2b)
    U3b=(1+c11b)*as.vector(H1b)
    U4b=(1+c12b)*as.vector(H1b)
    U5b=as.vector(H2b)
    U6b=as.vector(H1b)
    mu1b=mean(U1b)
    mu2b=mean(U2b)
    mu3b=mean(U3b)
    mu4b=mean(U4b)
    mu5b=mean(U5b)
    mu6b=mean(U6b)
    Delta.b.logs=-0.5*(log(mu1b)-log(mu3b)+log(mu2b)-log(mu4b))+log(mu5b)-log(mu6b)
    Delta.b=exp(Delta.b.logs)
    #
    tbar0=tbar0+(Delta.a+Delta.b)/2-Delta
    tbar[ind1[1:na]]=tbar[ind1[1:na]] +
      delta.a - delta[ind1[1:na]]
    tbar[ind1[(na+1):n]]=tbar[ind1[(na+1):n]] +
      delta.b - delta[ind1[(na+1):n]]
  }
  # tbar contains the bias for eff of each obs
  tbar=(1/L)*bc.fac*tbar
  tbar0=(1/L)*bc.fac*tbar0
  #
  Delta.bias=tbar0
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  ts=sig/sqrt(n)
  bounds0=matrix((Delta-ts*crit),nrow=3,ncol=2)
  bounds1=matrix((Delta-Delta.bias-ts*crit),nrow=3,ncol=2)

  if (p+q<4) {
    # make a list of results to return to calling routine and then quit:
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             sig=sig,
             estimate=c(Delta,Delta.bias,Delta-Delta.bias),
             estimate_individual=c(exp(delta-tbar))
             )
  } else {
    #Delta.nk=mean(delta[1:nk])
    U1.nk=(1+c21[1:nk])*as.vector(H2[1:nk])
    U2.nk=(1+c22[1:nk])*as.vector(H2[1:nk])
    U3.nk=(1+c11[1:nk])*as.vector(H1[1:nk])
    U4.nk=(1+c12[1:nk])*as.vector(H1[1:nk])
    U5.nk=as.vector(H2[1:nk])
    U6.nk=as.vector(H1[1:nk])
    mu1.nk=mean(U1.nk)
    mu2.nk=mean(U2.nk)
    mu3.nk=mean(U3.nk)
    mu4.nk=mean(U4.nk)
    mu5.nk=mean(U5.nk)
    mu6.nk=mean(U6.nk)
    Delta.nk.logs=-0.5*(log(mu1.nk)-log(mu3.nk)+log(mu2.nk)-log(mu4.nk))+log(mu5.nk)-log(mu6.nk)
    Delta.nk=exp(Delta.nk.logs)
    #
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    #
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             bounds2=bounds2,
             sig=sig,
             estimate=c(Delta.nk,Delta.bias,Delta.nk-Delta.bias,Delta,Delta.bias,Delta-Delta.bias),
             estimate_individual=c(exp(delta-tbar))
            )
  }
  return(res)
}