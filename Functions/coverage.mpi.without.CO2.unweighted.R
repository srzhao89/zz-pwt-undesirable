###########################################################################################
## Title: Statistical Inference for Environmental Efficiency and Productivity Dynamics
## Authors: Valentin Zelenyuk and Shirong Zhao
## Date: Aug 29, 2026
##############
## This function estimates the MPI under weak disposability
## x1, x2, is n times p matrix for inputs 
## y1, g2 is n times qg matrix for desirable outputs
## b1, b2 is n times qb matrix for undesirable outputs
###########################################################################################
coverage.mpi.without.CO2.unweighted<-function(x1=x1,y1=y1,x2=x2,y2=y2,
                                       XDIREC1=XDIREC1,YDIREC1=YDIREC1,
                                       XDIREC2=XDIREC2,YDIREC2=YDIREC2,L=50) {
  
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
  delta.exp=exp(delta)
  Delta=exp(mean(delta)) # Geometric Means of MLPIs
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
    #
    x2.b=matrix(x2[ind1,],ncol=p)
    y2.b=matrix(y2[ind1,],ncol=q) 
    XDIREC2.b=matrix(XDIREC2[ind1,],ncol=p)
    YDIREC2.b=matrix(YDIREC2[ind1,],ncol=q) 
    #
    x1a=matrix(x1.b[1:na,],ncol=p)
    y1a=matrix(y1.b[1:na,],ncol=q)   
    XDIREC1a=matrix(XDIREC1.b[1:na,],ncol=p)
    YDIREC1a=matrix(YDIREC1.b[1:na,],ncol=q)   
    
    x1b=matrix(x1.b[(na+1):n,],ncol=p)
    y1b=matrix(y1.b[(na+1):n,],ncol=q)
    XDIREC1b=matrix(XDIREC1.b[(na+1):n,],ncol=p)
    YDIREC1b=matrix(YDIREC1.b[(na+1):n,],ncol=q)
    #
    x2a=matrix(x2.b[1:na,],ncol=p)
    y2a=matrix(y2.b[1:na,],ncol=q)   
    XDIREC2a=matrix(XDIREC2.b[1:na,],ncol=p)
    YDIREC2a=matrix(YDIREC2.b[1:na,],ncol=q)   
    
    x2b=matrix(x2.b[(na+1):n,],ncol=p)
    y2b=matrix(y2.b[(na+1):n,],ncol=q)
    XDIREC2b=matrix(XDIREC2.b[(na+1):n,],ncol=p)
    YDIREC2b=matrix(YDIREC2.b[(na+1):n,],ncol=q)
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
    delta.exp.a=exp(delta.a)
    Delta.a=exp(mean(delta.a)) # Geometric Means of MLPIs
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
    delta.exp.b=exp(delta.b)
    Delta.b=exp(mean(delta.b)) # Geometric Means of MLPIs
    #
    tbar0=tbar0+(Delta.a+Delta.b)/2-Delta
    #
    tbar[ind1[1:na]]=tbar[ind1[1:na]] +
      delta.exp.a - delta.exp[ind1[1:na]]
    tbar[ind1[(na+1):n]]=tbar[ind1[(na+1):n]] +
      delta.exp.b - delta.exp[ind1[(na+1):n]]
  }
  # tbar contains the bias for eff of each obs
  tbar=(1/L)*bc.fac*tbar
  tbar0=(1/L)*bc.fac*tbar0
  #
  Delta.bias=tbar0
  sig=sd(delta)*Delta
  if (is.nan(sig)){
    sig=0
  } 
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
    Delta.nk=exp(mean(delta[1:nk]))
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    
    ts.nk=sig/sqrt(nk)
    
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