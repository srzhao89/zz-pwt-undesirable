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
coverage.mpi.with.CO2.unweighted<-function(x1,g1,b1,x2,g2,b2,
                                                 XDIREC1=XDIREC1,GDIREC1=GDIREC1,BDIREC1=BDIREC1,
                                                 XDIREC2=XDIREC2,GDIREC2=GDIREC2,BDIREC2=BDIREC2,L=10) {
  
  p=ncol(x1)
  qg=ncol(g1)
  qb=ncol(b1)
  n=nrow(x1)
  na=floor(n/2)
  nb=n-na
  kappa=2/(p+qg+qb+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  # evaluate efficiency using CRS-DEA
  c11=dea.direc.with.CO2.crs(XOBS=x1,GOBS=g1,BOBS=b1,
                              XDIREC=XDIREC1,GDIREC=GDIREC1,BDIREC=BDIREC1,
                              XREF=x1,GREF=g1,BREF=b1)
  c22=dea.direc.with.CO2.crs(XOBS=x2,GOBS=g2,BOBS=b2,
                              XDIREC=XDIREC2,GDIREC=GDIREC2,BDIREC=BDIREC2,
                              XREF=x2,GREF=g2,BREF=b2)
  c12=dea.direc.with.CO2.crs(XOBS=x1,GOBS=g1,BOBS=b1,
                              XDIREC=XDIREC1,GDIREC=GDIREC1,BDIREC=BDIREC1,
                              XREF=x2,GREF=g2,BREF=b2)
  c21=dea.direc.with.CO2.crs(XOBS=x2,GOBS=g2,BOBS=b2,
                              XDIREC=XDIREC2,GDIREC=GDIREC2,BDIREC=BDIREC2,
                              XREF=x1,GREF=g1,BREF=b1)
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
    x1.b=x1
    g1.b=g1
    b1.b=b1
    XDIREC1.b=XDIREC1
    GDIREC1.b=GDIREC1
    BDIREC1.b=BDIREC1
    #
    x2.b=x2
    g2.b=g2
    b2.b=b2
    XDIREC2.b=XDIREC2
    GDIREC2.b=GDIREC2
    BDIREC2.b=BDIREC2
    #
    
    ind1=sample(ind,size=n)
    x1.b[1:n,]=x1[ind1,]
    g1.b[1:n,]=g1[ind1,]
    b1.b[1:n,]=b1[ind1,]
    XDIREC1.b[1:n,]=XDIREC1[ind1,]
    GDIREC1.b[1:n,]=GDIREC1[ind1,]
    BDIREC1.b[1:n,]=BDIREC1[ind1,]
    #
    x2.b[1:n,]=x2[ind1,]
    g2.b[1:n,]=g2[ind1,]
    b2.b[1:n,]=b2[ind1,]
    XDIREC2.b[1:n,]=XDIREC2[ind1,]
    GDIREC2.b[1:n,]=GDIREC2[ind1,]
    BDIREC2.b[1:n,]=BDIREC2[ind1,]
    #
    x1a=matrix(x1.b[1:na,],ncol=p)
    g1a=matrix(g1.b[1:na,],ncol=qg)   
    b1a=matrix(b1.b[1:na,],ncol=qb)
    XDIREC1a=matrix(XDIREC1.b[1:na,],ncol=p)
    GDIREC1a=matrix(GDIREC1.b[1:na,],ncol=qg)   
    BDIREC1a=matrix(BDIREC1.b[1:na,],ncol=qb)
    
    x1b=matrix(x1.b[(na+1):n,],ncol=p)
    g1b=matrix(g1.b[(na+1):n,],ncol=qg)
    b1b=matrix(b1.b[(na+1):n,],ncol=qb)
    XDIREC1b=matrix(XDIREC1.b[(na+1):n,],ncol=p)
    GDIREC1b=matrix(GDIREC1.b[(na+1):n,],ncol=qg)
    BDIREC1b=matrix(BDIREC1.b[(na+1):n,],ncol=qb)
    #
    x2a=matrix(x2.b[1:na,],ncol=p)
    g2a=matrix(g2.b[1:na,],ncol=qg)   
    b2a=matrix(b2.b[1:na,],ncol=qb)
    XDIREC2a=matrix(XDIREC2.b[1:na,],ncol=p)
    GDIREC2a=matrix(GDIREC2.b[1:na,],ncol=qg)   
    BDIREC2a=matrix(BDIREC2.b[1:na,],ncol=qb)
    
    x2b=matrix(x2.b[(na+1):n,],ncol=p)
    g2b=matrix(g2.b[(na+1):n,],ncol=qg)
    b2b=matrix(b2.b[(na+1):n,],ncol=qb)
    XDIREC2b=matrix(XDIREC2.b[(na+1):n,],ncol=p)
    GDIREC2b=matrix(GDIREC2.b[(na+1):n,],ncol=qg)
    BDIREC2b=matrix(BDIREC2.b[(na+1):n,],ncol=qb)
    #
    c11a=dea.direc.with.CO2.crs(XOBS=x1a,GOBS=g1a,BOBS=b1a,
                                 XDIREC=XDIREC1a,GDIREC=GDIREC1a,BDIREC=BDIREC1a,
                                 XREF=x1a,GREF=g1a,BREF=b1a)
    c22a=dea.direc.with.CO2.crs(XOBS=x2a,GOBS=g2a,BOBS=b2a,
                                 XDIREC=XDIREC2a,GDIREC=GDIREC2a,BDIREC=BDIREC2a,
                                 XREF=x2a,GREF=g2a,BREF=b2a)
    c12a=dea.direc.with.CO2.crs(XOBS=x1a,GOBS=g1a,BOBS=b1a,
                                 XDIREC=XDIREC1a,GDIREC=GDIREC1a,BDIREC=BDIREC1a,
                                 XREF=x2a,GREF=g2a,BREF=b2a)
    c21a=dea.direc.with.CO2.crs(XOBS=x2a,GOBS=g2a,BOBS=b2a,
                                 XDIREC=XDIREC2a,GDIREC=GDIREC2a,BDIREC=BDIREC2a,
                                 XREF=x1a,GREF=g1a,BREF=b1a)
    #
    delta.a=-0.5*(log(1+c21a)-log(1+c11a)+log(1+c22a)-log(1+c12a))
    delta.exp.a=exp(delta.a)
    Delta.a=exp(mean(delta.a)) # Geometric Means of MLPIs
    #
    c11b=dea.direc.with.CO2.crs(XOBS=x1b,GOBS=g1b,BOBS=b1b,
                                 XDIREC=XDIREC1b,GDIREC=GDIREC1b,BDIREC=BDIREC1b,
                                 XREF=x1b,GREF=g1b,BREF=b1b)
    c22b=dea.direc.with.CO2.crs(XOBS=x2b,GOBS=g2b,BOBS=b2b,
                                 XDIREC=XDIREC2b,GDIREC=GDIREC2b,BDIREC=BDIREC2b,
                                 XREF=x2b,GREF=g2b,BREF=b2b)
    c12b=dea.direc.with.CO2.crs(XOBS=x1b,GOBS=g1b,BOBS=b1b,
                                 XDIREC=XDIREC1b,GDIREC=GDIREC1b,BDIREC=BDIREC1b,
                                 XREF=x2b,GREF=g2b,BREF=b2b)
    c21b=dea.direc.with.CO2.crs(XOBS=x2b,GOBS=g2b,BOBS=b2b,
                                 XDIREC=XDIREC2b,GDIREC=GDIREC2b,BDIREC=BDIREC2b,
                                 XREF=x1b,GREF=g1b,BREF=b1b)
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
  ##
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  ts=sig/sqrt(n)
  
  bounds0=matrix((Delta-ts*crit),nrow=3,ncol=2)
  bounds1=matrix((Delta-Delta.bias-ts*crit),nrow=3,ncol=2)

  if (p+qg+qb<4) {
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
    #
    bounds2=matrix((Delta-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
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