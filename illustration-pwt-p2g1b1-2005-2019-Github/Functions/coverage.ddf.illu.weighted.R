###########################################################################################
## Title: Central Limit Theorems for Aggregates of Directional Distance Functions
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: January 22, 2024
##############
## This function estimates the aggregate DDFs under weak disposability, along its
## confidence intervals using the theoretical results in the above mentioned paper.
## x is n times p matrix for inputs 
## g is n times qg matrix for desirable outputs
## b is n times qb matrix for undesirable outputs
###########################################################################################
coverage.ddf.illu.weighted<-function(x,y,XDIREC=XDIREC,YDIREC=YDIREC,H2=H2,L=10) {
  
  p=ncol(x)
  q=ncol(y)
  n=nrow(x)
  na=floor(n/2)
  nb=n-na
  kappa=2/(p+q+1)
  bc.fac=1/(2**kappa - 1)
  nk=floor(n**(2*kappa))
  # evaluate efficiency using VRS-DEA
  #delta=FEAR::dea.direc(XOBS=t(x),YOBS=t(y),XDIREC=t(XDIREC),YDIREC=t(YDIREC),
  #                      RTS=1)
  
  delta= 1-Benchmarking::dea.direct(X=x,Y=y,
                                    DIREC=cbind(XDIREC,YDIREC),
                                    XREF=x,YREF=y,
                                    ORIENTATION="in-out",
                                    RTS="vrs")$eff[,1]
  #
  H1=delta*H2
  mu1=mean(H1)
  mu2=mean(H2)
  Delta=mean(H1)/mean(H2)
  # OECD
  mu1.OECD=mean(H1[OECD])
  mu2.OECD=mean(H2[OECD])
  Delta.OECD=mean(H1[OECD])/mean(H2[OECD])
  # Non-OECD
  mu1.NOECD=mean(H1[-OECD])
  mu2.NOECD=mean(H2[-OECD])
  Delta.NOECD=mean(H1[-OECD])/mean(H2[-OECD])
  #####################################################
  # compute bias corrections via generalized jackknife:
  tbar=rep(0,n)
  ind=c(1:n)
  for (j in 1:L) {
    
    ind1=sample(ind,size=n)
    x.b=matrix(x[ind1,],ncol=p)
    y.b=matrix(y[ind1,],ncol=q)
    XDIREC.b=matrix(XDIREC[ind1,],ncol=p)
    YDIREC.b=matrix(YDIREC[ind1,],ncol=q)
    #
    xa=matrix(x.b[1:na,],ncol=p)
    ya=matrix(y.b[1:na,],ncol=q)
    XDIRECa=matrix(XDIREC.b[1:na,],ncol=p)
    YDIRECa=matrix(YDIREC.b[1:na,],ncol=q)
    #
    xb=matrix(x.b[(na+1):n,],ncol=p)
    yb=matrix(y.b[(na+1):n,],ncol=q)
    XDIRECb=matrix(XDIREC.b[(na+1):n,],ncol=p)
    YDIRECb=matrix(YDIREC.b[(na+1):n,],ncol=q)
    #
    #delta.a=FEAR::dea.direc(XOBS=t(xa),
    #                        YOBS=t(ya),
    #                        XDIREC=t(XDIRECa),
    #                        YDIREC=t(YDIRECa),
    #                        RTS=1
    #                        )
    #delta.b=FEAR::dea.direc(XOBS=t(xb),
    #                        YOBS=t(yb),
    #                        XDIREC=t(XDIRECb),
    #                        YDIREC=t(YDIRECb),
    #                        RTS=1
    #                        )
    delta.a= 1-Benchmarking::dea.direct(X=xa,Y=ya,
                                        DIREC=cbind(XDIRECa,YDIRECa),
                                        XREF=xa,YREF=ya,
                                        ORIENTATION="in-out",
                                        RTS="vrs")$eff[,1]
    
    delta.b= 1-Benchmarking::dea.direct(X=xb,Y=yb,
                                        DIREC=cbind(XDIRECb,YDIRECb),
                                        XREF=xb,YREF=yb,
                                        ORIENTATION="in-out",
                                        RTS="vrs")$eff[,1]
    #
    tbar[ind1[1:na]]=tbar[ind1[1:na]] +
      delta.a - delta[ind1[1:na]]
    tbar[ind1[(na+1):n]]=tbar[ind1[(na+1):n]] +
      delta.b - delta[ind1[(na+1):n]]
  }
  # tbar contains the bias for eff of each obs
  tbar=(1/L)*bc.fac*tbar
  #
  H1=delta*H2
  #
  mu1=mean(H1)
  mu2=mean(H2)
  Delta=mean(H1)/mean(H2)
  Delta.bias=mean(tbar*H2)/mean(H2) 
  var1=var(H1)
  var2=var(H2)
  cov12=cov(H1,H2)
  var=Delta^2*(var1/(mu1^2)+var2/(mu2^2)-2*cov12/(mu1*mu2)) 
  sig=sqrt(var)
  if (is.nan(sig)){
    sig=0
  } 
  #
  # OECD
  Delta.OECD.bias=mean(tbar[OECD]*H2[OECD])/mean(H2[OECD]) 
  var1.OECD=var(H1[OECD])
  var2.OECD=var(H2[OECD])
  cov12.OECD=cov(H1[OECD],H2[OECD])
  var.OECD=Delta.OECD^2*(var1.OECD/(mu1.OECD^2)+var2/(mu2.OECD^2)-2*cov12.OECD/(mu1.OECD*mu2.OECD)) 
  sig.OECD=sqrt(var.OECD)
  # Non-OECD
  Delta.NOECD.bias=mean(tbar[-OECD]*H2[-OECD])/mean(H2[-OECD]) 
  var1.NOECD=var(H1[-OECD])
  var2.NOECD=var(H2[-OECD])
  cov12.NOECD=cov(H1[-OECD],H2[-OECD])
  var.NOECD=Delta.NOECD^2*(var1.NOECD/(mu1.NOECD^2)+var2/(mu2.NOECD^2)-2*cov12.NOECD/(mu1.NOECD*mu2.NOECD)) 
  sig.NOECD=sqrt(var.NOECD)
  ##################
  crit=qnorm(p=c(0.95,0.975,0.995,0.05,0.025,0.005))
  n.OECD=length(OECD)
  n.NOECD=n-length(OECD)
  ts=sig/sqrt(n)
  bounds0=matrix((Delta-ts*crit),nrow=3,ncol=2)
  bounds1=matrix((Delta-Delta.bias-ts*crit),nrow=3,ncol=2)
  if (p+q<4) {
    # make a list of results to return to calling routine and then quit:
    that=((Delta.OECD-Delta.OECD.bias)-(Delta.NOECD-Delta.NOECD.bias))/sqrt(sig.OECD^2/n.OECD+sig.NOECD^2/n.NOECD)
    pval=ifelse(that<0,2*pnorm(that,lower.tail=TRUE),2*pnorm(that,lower.tail=FALSE))
    #
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             sig=sig,
             estimate=c(Delta,Delta.bias,Delta-Delta.bias),
             estimate_individual=c(delta-tbar),
             estimate.OECD=c(Delta.OECD,Delta.OECD.bias,Delta.OECD-Delta.OECD.bias),
             estimate.NOECD=c(Delta.NOECD,Delta.NOECD.bias,Delta.NOECD-Delta.NOECD.bias),
             that=that,
             pval=pval
             )
  } else {
    Delta.nk=mean(H1[1:nk])/mean(H2[1:nk])
    ts.nk=sig/sqrt(nk)
    bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    #
    nk.OECD=floor(n.OECD**(2*kappa))
    Delta.OECD.nk=mean((H1[OECD])[1:nk.OECD])/mean((H2[OECD])[1:nk.OECD])
    #
    nk.NOECD=floor(n.NOECD**(2*kappa))
    Delta.NOECD.nk=mean((H1[-OECD])[1:nk.NOECD])/mean((H2[-OECD])[1:nk.NOECD])
    # test equality between OECD and Non-OECD
    #that=((Delta.OECD.nk-Delta.OECD.bias)-(Delta.NOECD.nk-Delta.NOECD.bias))/sqrt(sig.OECD^2/nk.OECD+sig.NOECD^2/nk.NOECD)
    that=((Delta.OECD-Delta.OECD.bias)-(Delta.NOECD-Delta.NOECD.bias))/sqrt(sig.OECD^2/nk.OECD+sig.NOECD^2/nk.NOECD)
    pval=ifelse(that<0,2*pnorm(that,lower.tail=TRUE),2*pnorm(that,lower.tail=FALSE))
    #
    res=list(bounds0=bounds0,
             bounds1=bounds1,
             bounds2=bounds2,
             sig=sig,
             estimate=c(Delta.nk,Delta.bias,Delta.nk-Delta.bias,Delta,Delta.bias,Delta-Delta.bias),
             estimate_individual=c(delta-tbar),             
             estimate.OECD=c(Delta.OECD.nk,Delta.OECD.bias,Delta.OECD.nk-Delta.OECD.bias,Delta.OECD,Delta.OECD.bias,Delta.OECD-Delta.OECD.bias),
             estimate.NOECD=c(Delta.NOECD.nk,Delta.NOECD.bias,Delta.NOECD.nk-Delta.NOECD.bias,Delta.NOECD,Delta.NOECD.bias,Delta.NOECD-Delta.NOECD.bias),
             that=that,
             pval=pval
             )
  }
  return(res)
}