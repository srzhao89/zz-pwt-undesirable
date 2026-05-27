###########################################################################################
## Title: ???
## Authors:  Valentin Zelenyuk, and Shirong Zhao
## Date: July 17, 2025
##############
## This function estimates the MPI under weak disposability
## x1, x2, is n times p matrix for inputs 
## g1, g2 is n times qg matrix for desirable outputs
## b1, b2 is n times qb matrix for undesirable outputs
###########################################################################################
coverage.mpi.kuosmanen.illu.unweighted<-function(x1,g1,b1,x2,g2,b2,
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
  c11=dea.direc.kuosmanen.crs(XOBS=x1,GOBS=g1,BOBS=b1,
                              XDIREC=XDIREC1,GDIREC=GDIREC1,BDIREC=BDIREC1,
                              XREF=x1,GREF=g1,BREF=b1)
  c22=dea.direc.kuosmanen.crs(XOBS=x2,GOBS=g2,BOBS=b2,
                              XDIREC=XDIREC2,GDIREC=GDIREC2,BDIREC=BDIREC2,
                              XREF=x2,GREF=g2,BREF=b2)
  c12=dea.direc.kuosmanen.crs(XOBS=x1,GOBS=g1,BOBS=b1,
                              XDIREC=XDIREC1,GDIREC=GDIREC1,BDIREC=BDIREC1,
                              XREF=x2,GREF=g2,BREF=b2)
  c21=dea.direc.kuosmanen.crs(XOBS=x2,GOBS=g2,BOBS=b2,
                              XDIREC=XDIREC2,GDIREC=GDIREC2,BDIREC=BDIREC2,
                              XREF=x1,GREF=g1,BREF=b1)
  #
  delta=-0.5*(log(1+c21)-log(1+c11)+log(1+c22)-log(1+c12))
  delta.exp=exp(delta)
  Delta=exp(mean(delta)) # Geometric Means of MLPIs
  ## For OECD
  Delta.OECD=exp(mean(delta[OECD]))
  ## For Non-OECD
  Delta.NOECD=exp(mean(delta[-OECD]))
  #####################################################
  # compute bias corrections via generalized jackknife:
  tbar0=0
  tbar0.OECD=0
  tbar0.NOECD=0
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
    c11a=dea.direc.kuosmanen.crs(XOBS=x1a,GOBS=g1a,BOBS=b1a,
                                 XDIREC=XDIREC1a,GDIREC=GDIREC1a,BDIREC=BDIREC1a,
                                 XREF=x1a,GREF=g1a,BREF=b1a)
    c22a=dea.direc.kuosmanen.crs(XOBS=x2a,GOBS=g2a,BOBS=b2a,
                                 XDIREC=XDIREC2a,GDIREC=GDIREC2a,BDIREC=BDIREC2a,
                                 XREF=x2a,GREF=g2a,BREF=b2a)
    c12a=dea.direc.kuosmanen.crs(XOBS=x1a,GOBS=g1a,BOBS=b1a,
                                 XDIREC=XDIREC1a,GDIREC=GDIREC1a,BDIREC=BDIREC1a,
                                 XREF=x2a,GREF=g2a,BREF=b2a)
    c21a=dea.direc.kuosmanen.crs(XOBS=x2a,GOBS=g2a,BOBS=b2a,
                                 XDIREC=XDIREC2a,GDIREC=GDIREC2a,BDIREC=BDIREC2a,
                                 XREF=x1a,GREF=g1a,BREF=b1a)
    #
    delta.a=-0.5*(log(1+c21a)-log(1+c11a)+log(1+c22a)-log(1+c12a))
    delta.exp.a=exp(delta.a)
    Delta.a=exp(mean(delta.a)) # Geometric Means of MLPIs
    ## For OECD and Non-OECD
    OECD.a=which(ind1[1:na] %in% OECD)
    Delta.OECD.a=exp(mean(delta.a[OECD.a]))
    Delta.NOECD.a=exp(mean(delta.a[-OECD.a]))
    #
    c11b=dea.direc.kuosmanen.crs(XOBS=x1b,GOBS=g1b,BOBS=b1b,
                                 XDIREC=XDIREC1b,GDIREC=GDIREC1b,BDIREC=BDIREC1b,
                                 XREF=x1b,GREF=g1b,BREF=b1b)
    c22b=dea.direc.kuosmanen.crs(XOBS=x2b,GOBS=g2b,BOBS=b2b,
                                 XDIREC=XDIREC2b,GDIREC=GDIREC2b,BDIREC=BDIREC2b,
                                 XREF=x2b,GREF=g2b,BREF=b2b)
    c12b=dea.direc.kuosmanen.crs(XOBS=x1b,GOBS=g1b,BOBS=b1b,
                                 XDIREC=XDIREC1b,GDIREC=GDIREC1b,BDIREC=BDIREC1b,
                                 XREF=x2b,GREF=g2b,BREF=b2b)
    c21b=dea.direc.kuosmanen.crs(XOBS=x2b,GOBS=g2b,BOBS=b2b,
                                 XDIREC=XDIREC2b,GDIREC=GDIREC2b,BDIREC=BDIREC2b,
                                 XREF=x1b,GREF=g1b,BREF=b1b)
    #
    delta.b=-0.5*(log(1+c21b)-log(1+c11b)+log(1+c22b)-log(1+c12b))
    delta.exp.b=exp(delta.b)
    Delta.b=exp(mean(delta.b)) # Geometric Means of MLPIs
    ## For OECD and Non-OECD
    OECD.b=which(ind1[(na+1):n] %in% OECD)
    Delta.OECD.b=exp(mean(delta.b[OECD.b]))
    Delta.NOECD.b=exp(mean(delta.b[-OECD.b]))
    #
    tbar0=tbar0+(Delta.a+Delta.b)/2-Delta
    tbar0.OECD=tbar0.OECD+(Delta.OECD.a+Delta.OECD.b)/2-Delta.OECD
    tbar0.NOECD=tbar0.NOECD+(Delta.NOECD.a+Delta.NOECD.b)/2-Delta.NOECD
    #
    tbar[ind1[1:na]]=tbar[ind1[1:na]] +
      delta.exp.a - delta.exp[ind1[1:na]]
    tbar[ind1[(na+1):n]]=tbar[ind1[(na+1):n]] +
      delta.exp.b - delta.exp[ind1[(na+1):n]]
  }
  # tbar contains the bias for eff of each obs
  tbar=(1/L)*bc.fac*tbar
  tbar0=(1/L)*bc.fac*tbar0
  tbar0.OECD=(1/L)*bc.fac*tbar0.OECD
  tbar0.NOECD=(1/L)*bc.fac*tbar0.NOECD
  #
  Delta.bias=tbar0
  sig=sd(delta)*Delta
  if (is.nan(sig)){
    sig=0
  } 
  ## For OECD
  Delta.OECD.bias=tbar0.OECD
  sig.OECD=sd(delta[OECD])*Delta.OECD
  ## For Non-OECD
  Delta.NOECD.bias=tbar0.NOECD
  sig.NOECD=sd(delta[-OECD])*Delta.NOECD
  ##
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  n.OECD=length(OECD)
  n.NOECD=n-length(OECD)
  ts=sig/sqrt(n)
  ts.OECD=sig.OECD/sqrt(n.OECD)
  ts.NOECD=sig.NOECD/sqrt(n.NOECD)
  bounds0=matrix((Delta-ts*crit),nrow=3,ncol=2)
  bounds1=matrix((Delta-Delta.bias-ts*crit),nrow=3,ncol=2)
  bounds1.OECD=matrix((Delta.OECD-Delta.OECD.bias-ts.OECD*crit),nrow=3,ncol=2)
  bounds1.NOECD=matrix((Delta.NOECD-Delta.NOECD.bias-ts.NOECD*crit),nrow=3,ncol=2)
  if (p+qg+qb<4) {
    that=((Delta.OECD-Delta.OECD.bias)-(Delta.NOECD-Delta.NOECD.bias))/sqrt(sig.OECD^2/n.OECD+sig.NOECD^2/n.NOECD)
    pval=ifelse(that<0,2*pnorm(that,lower.tail=TRUE),2*pnorm(that,lower.tail=FALSE))
    # make a list of results to return to calling routine and then quit:
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             sig=sig,
             estimate=c(Delta,Delta.bias,Delta-Delta.bias),
             estimate_individual=c(exp(delta-tbar)),
             bounds1.OECD=bounds1.OECD,
             estimate.OECD=c(Delta.OECD,Delta.OECD.bias,Delta.OECD-Delta.OECD.bias),
             bounds1.NOECD=bounds1.NOECD,
             estimate.NOECD=c(Delta.NOECD,Delta.NOECD.bias,Delta.NOECD-Delta.NOECD.bias),
             that=that,
             pval=pval
             )
  } else {
    Delta.nk=exp(mean(delta[1:nk]))

    ## For OECD
    nk.OECD=floor(n.OECD**(2*kappa))
    Delta.OECD.nk=exp(mean((delta[OECD])[1:nk.OECD]))
    ## For Non-OECD
    nk.NOECD=floor(n.NOECD**(2*kappa))
    Delta.NOECD.nk=exp(mean((delta[-OECD])[1:nk.NOECD]))
    
    ts.nk=sig/sqrt(nk)
    ts.OECD.nk=sig.OECD/sqrt(nk.OECD)
    ts.NOECD.nk=sig.NOECD/sqrt(nk.NOECD)
    
    #bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    #bounds2.OECD=matrix((Delta.OECD.nk-Delta.OECD.bias-ts.OECD.nk*crit),nrow=3,ncol=2)
    #bounds2.NOECD=matrix((Delta.NOECD.nk-Delta.NOECD.bias-ts.NOECD.nk*crit),nrow=3,ncol=2)
    #
    bounds2=matrix((Delta-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    bounds2.OECD=matrix((Delta.OECD-Delta.OECD.bias-ts.OECD.nk*crit),nrow=3,ncol=2)
    bounds2.NOECD=matrix((Delta.NOECD-Delta.NOECD.bias-ts.NOECD.nk*crit),nrow=3,ncol=2)
    # test equality between OECD and Non-OECD
    #that=((Delta.OECD.nk-Delta.OECD.bias)-(Delta.NOECD.nk-Delta.NOECD.bias))/sqrt(sig.OECD^2/nk.OECD+sig.NOECD^2/nk.NOECD)
    #that=((Delta.OECD-Delta.OECD.bias)-(Delta.NOECD-Delta.NOECD.bias))/sqrt(sig.OECD^2/n.OECD+sig.NOECD^2/n.NOECD)
    that=((Delta.OECD-Delta.OECD.bias)-(Delta.NOECD-Delta.NOECD.bias))/sqrt(sig.OECD^2/nk.OECD+sig.NOECD^2/nk.NOECD)
    pval=ifelse(that<0,2*pnorm(that,lower.tail=TRUE),2*pnorm(that,lower.tail=FALSE))
    #
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             bounds2=bounds2,
             sig=sig,
             estimate=c(Delta.nk,Delta.bias,Delta.nk-Delta.bias,Delta,Delta.bias,Delta-Delta.bias),
             estimate_individual=c(delta.exp-tbar),
             bounds2.OECD=bounds2.OECD,
             sig.OECD=sig.OECD,
             estimate.OECD=c(Delta.OECD.nk,Delta.OECD.bias,Delta.OECD.nk-Delta.OECD.bias,Delta.OECD,Delta.OECD.bias,Delta.OECD-Delta.OECD.bias),
             bounds2.NOECD=bounds2.NOECD,
             sig.NOECD=sig.NOECD,
             estimate.NOECD=c(Delta.NOECD.nk,Delta.NOECD.bias,Delta.NOECD.nk-Delta.NOECD.bias,Delta.NOECD,Delta.NOECD.bias,Delta.NOECD-Delta.NOECD.bias),
             that=that,
             pval=pval)
  }
  return(res)
}