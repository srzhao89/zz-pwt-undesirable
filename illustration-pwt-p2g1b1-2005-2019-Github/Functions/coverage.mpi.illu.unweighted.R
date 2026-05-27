###########################################################################################
## Title: ???
## Authors:  Valentin Zelenyuk, and Shirong Zhao
## Date: July 17, 2025
##############
## This function estimates the MPI under weak disposability
## x1, x2, is n times p matrix for inputs 
## y1, g2 is n times qg matrix for desirable outputs
## b1, b2 is n times qb matrix for undesirable outputs
###########################################################################################
coverage.mpi.illu.unweighted<-function(x1=x1,y1=y1,x2=x2,y2=y2,
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
  #c11=FEAR::dea.direc(XOBS=t(x1),YOBS=t(y1),XDIREC=t(XDIREC1),YDIREC=t(YDIREC1),
  #                    XREF=t(x1),YREF=t(y1),RTS=3)
  #c22=FEAR::dea.direc(XOBS=t(x2),YOBS=t(y2),XDIREC=t(XDIREC2),YDIREC=t(YDIREC2),
  #                    XREF=t(x2),YREF=t(y2),RTS=3)
  #c12=FEAR::dea.direc(XOBS=t(x1),YOBS=t(y1),XDIREC=t(XDIREC1),YDIREC=t(YDIREC1),
  #                    XREF=t(x2),YREF=t(y2),RTS=3)
  #c21=FEAR::dea.direc(XOBS=t(x2),YOBS=t(y2),XDIREC=t(XDIREC2),YDIREC=t(YDIREC2),
  #                    XREF=t(x1),YREF=t(y1),RTS=3)
  
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
    #x1.b=x1
    #y1.b=y1
    #XDIREC1.b=XDIREC1
    #YDIREC1.b=YDIREC1
    #
    #x2.b=x2
    #y2.b=y2
    #XDIREC2.b=XDIREC2
    #YDIREC2.b=YDIREC2
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
    #c11a=FEAR::dea.direc(XOBS=t(x1a),YOBS=t(y1a),XDIREC=t(XDIREC1a),YDIREC=t(YDIREC1a),
    #                    XREF=t(x1a),YREF=t(y1a),RTS=3)
    #c22a=FEAR::dea.direc(XOBS=t(x2a),YOBS=t(y2a),XDIREC=t(XDIREC2a),YDIREC=t(YDIREC2a),
    #                    XREF=t(x2a),YREF=t(y2a),RTS=3)
    #c12a=FEAR::dea.direc(XOBS=t(x1a),YOBS=t(y1a),XDIREC=t(XDIREC1a),YDIREC=t(YDIREC1a),
    #                    XREF=t(x2a),YREF=t(y2a),RTS=3)
    #c21a=FEAR::dea.direc(XOBS=t(x2a),YOBS=t(y2a),XDIREC=t(XDIREC2a),YDIREC=t(YDIREC2a),
    #                    XREF=t(x1a),YREF=t(y1a),RTS=3)
    
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
    ## For OECD and Non-OECD
    OECD.a=which(ind1[1:na] %in% OECD)
    Delta.OECD.a=exp(mean(delta.a[OECD.a]))
    Delta.NOECD.a=exp(mean(delta.a[-OECD.a]))
    #
    #c11b=FEAR::dea.direc(XOBS=t(x1b),YOBS=t(y1b),XDIREC=t(XDIREC1b),YDIREC=t(YDIREC1b),
    #                    XREF=t(x1b),YREF=t(y1b),RTS=3)
    #c22b=FEAR::dea.direc(XOBS=t(x2b),YOBS=t(y2b),XDIREC=t(XDIREC2b),YDIREC=t(YDIREC2b),
    #                    XREF=t(x2b),YREF=t(y2b),RTS=3)
    #c12b=FEAR::dea.direc(XOBS=t(x1b),YOBS=t(y1b),XDIREC=t(XDIREC1b),YDIREC=t(YDIREC1b),
    #                    XREF=t(x2b),YREF=t(y2b),RTS=3)
    #c21b=FEAR::dea.direc(XOBS=t(x2b),YOBS=t(y2b),XDIREC=t(XDIREC2b),YDIREC=t(YDIREC2b),
    #                    XREF=t(x1b),YREF=t(y1b),RTS=3)
    
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
  if (p+q<4) {
    # test equality between OECD and Non-OECD
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
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    
    ## For OECD
    nk.OECD=floor(n.OECD**(2*kappa))
    Delta.OECD.nk=exp(mean((delta[OECD])[1:nk.OECD]))
    ## For Non-OECD
    nk.NOECD=floor(n.NOECD**(2*kappa))
    Delta.NOECD.nk=exp(mean((delta[-OECD])[1:nk.NOECD]))
    
    ts.nk=sig/sqrt(nk)
    ts.OECD.nk=sig.OECD/sqrt(nk.OECD)
    ts.NOECD.nk=sig.NOECD/sqrt(nk.NOECD)
    
    bounds2.OECD=matrix((Delta.OECD.nk-Delta.OECD.bias-ts.OECD.nk*crit),nrow=3,ncol=2)
    bounds2.NOECD=matrix((Delta.NOECD.nk-Delta.NOECD.bias-ts.NOECD.nk*crit),nrow=3,ncol=2)
    # test equality between OECD and Non-OECD
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