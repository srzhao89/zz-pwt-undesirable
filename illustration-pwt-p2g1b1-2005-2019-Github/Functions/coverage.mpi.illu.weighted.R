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
coverage.mpi.illu.weighted<-function(x1=x1,y1=y1,x2=x2,y2=y2,
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
  gg1=1/(2*mu1)
  gg2=1/(2*mu2)
  gg3=-1/(2*mu3)
  gg4=-1/(2*mu4)
  gg5=-1/mu5
  gg6=1/mu6
  gg=c(gg1,gg2,gg3,gg4,gg5,gg6)
  V.Sigma=t(gg)%*%Sigma%*%gg
  sig=sqrt(V.Sigma)
  sig=as.vector(sig)
  sig=sig*Delta
  if (is.nan(sig)){
    sig=0
  } 
  ## For OECD
  c11.OECD=c11[OECD]
  c22.OECD=c22[OECD]
  c12.OECD=c12[OECD]
  c21.OECD=c21[OECD]
  H1.OECD=H1[OECD]
  H2.OECD=H2[OECD]
  U1.OECD=(1+c21.OECD)*as.vector(H2.OECD)
  U2.OECD=(1+c22.OECD)*as.vector(H2.OECD)
  U3.OECD=(1+c11.OECD)*as.vector(H1.OECD)
  U4.OECD=(1+c12.OECD)*as.vector(H1.OECD)
  U5.OECD=as.vector(H2.OECD)
  U6.OECD=as.vector(H1.OECD)
  mu1.OECD=mean(U1.OECD)
  mu2.OECD=mean(U2.OECD)
  mu3.OECD=mean(U3.OECD)
  mu4.OECD=mean(U4.OECD)
  mu5.OECD=mean(U5.OECD)
  mu6.OECD=mean(U6.OECD)
  Delta.OECD.logs=-0.5*(log(mu1.OECD)-log(mu3.OECD)+log(mu2.OECD)-log(mu4.OECD))+log(mu5.OECD)-log(mu6.OECD)
  Delta.OECD=exp(Delta.OECD.logs)
  #
  Ti.OECD=cbind(U1.OECD,U2.OECD,U3.OECD,U4.OECD,U5.OECD,U6.OECD)
  Sigma.OECD=cov(Ti.OECD)
  gg1.OECD=1/(2*mu1.OECD)
  gg2.OECD=1/(2*mu2.OECD)
  gg3.OECD=-1/(2*mu3.OECD)
  gg4.OECD=-1/(2*mu4.OECD)
  gg5.OECD=-1/mu5.OECD
  gg6.OECD=1/mu6.OECD
  gg.OECD=c(gg1.OECD,gg2.OECD,gg3.OECD,gg4.OECD,gg5.OECD,gg6.OECD)
  V.Sigma.OECD=t(gg.OECD)%*%Sigma.OECD%*%gg.OECD
  sig.OECD=sqrt(V.Sigma.OECD)
  sig.OECD=as.vector(sig.OECD)
  sig.OECD=sig.OECD*Delta.OECD
  ## For Non-OECD
  c11.NOECD=c11[-OECD]
  c22.NOECD=c22[-OECD]
  c12.NOECD=c12[-OECD]
  c21.NOECD=c21[-OECD]
  H1.NOECD=H1[-OECD]
  H2.NOECD=H2[-OECD]
  U1.NOECD=(1+c21.NOECD)*as.vector(H2.NOECD)
  U2.NOECD=(1+c22.NOECD)*as.vector(H2.NOECD)
  U3.NOECD=(1+c11.NOECD)*as.vector(H1.NOECD)
  U4.NOECD=(1+c12.NOECD)*as.vector(H1.NOECD)
  U5.NOECD=as.vector(H2.NOECD)
  U6.NOECD=as.vector(H1.NOECD)
  mu1.NOECD=mean(U1.NOECD)
  mu2.NOECD=mean(U2.NOECD)
  mu3.NOECD=mean(U3.NOECD)
  mu4.NOECD=mean(U4.NOECD)
  mu5.NOECD=mean(U5.NOECD)
  mu6.NOECD=mean(U6.NOECD)
  Delta.NOECD.logs=-0.5*(log(mu1.NOECD)-log(mu3.NOECD)+log(mu2.NOECD)-log(mu4.NOECD))+log(mu5.NOECD)-log(mu6.NOECD)
  Delta.NOECD=exp(Delta.NOECD.logs)
  #
  Ti.NOECD=cbind(U1.NOECD,U2.NOECD,U3.NOECD,U4.NOECD,U5.NOECD,U6.NOECD)
  Sigma.NOECD=cov(Ti.NOECD)
  gg1.NOECD=1/(2*mu1.NOECD)
  gg2.NOECD=1/(2*mu2.NOECD)
  gg3.NOECD=-1/(2*mu3.NOECD)
  gg4.NOECD=-1/(2*mu4.NOECD)
  gg5.NOECD=-1/mu5.NOECD
  gg6.NOECD=1/mu6.NOECD
  gg.NOECD=c(gg1.NOECD,gg2.NOECD,gg3.NOECD,gg4.NOECD,gg5.NOECD,gg6.NOECD)
  V.Sigma.NOECD=t(gg.NOECD)%*%Sigma.NOECD%*%gg.NOECD
  sig.NOECD=sqrt(V.Sigma.NOECD)
  sig.NOECD=as.vector(sig.NOECD)
  sig.NOECD=sig.NOECD*Delta.NOECD
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
    #H1.b=H1
    #
    #x2.b=x2
    #y2.b=y2
    #XDIREC2.b=XDIREC2
    #YDIREC2.b=YDIREC2
    #H2.b=H2
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
    Delta.a.logs=-0.5*(log(mu1a)-log(mu3a)+log(mu2a)-log(mu4a))+log(mu5)-log(mu6)
    Delta.a=exp(Delta.a.logs)
    ## For OECD and Non-OECD
    OECD.a=which(ind1[1:na] %in% OECD)
    c11.OECD.a=c11a[OECD.a]
    c22.OECD.a=c22a[OECD.a]
    c12.OECD.a=c12a[OECD.a]
    c21.OECD.a=c21a[OECD.a]
    H1.OECD.a=H1a[OECD.a]
    H2.OECD.a=H2a[OECD.a]
    U1.OECD.a=(1+c21.OECD.a)*as.vector(H2.OECD.a)
    U2.OECD.a=(1+c22.OECD.a)*as.vector(H2.OECD.a)
    U3.OECD.a=(1+c11.OECD.a)*as.vector(H1.OECD.a)
    U4.OECD.a=(1+c12.OECD.a)*as.vector(H1.OECD.a)
    U5.OECD.a=as.vector(H2.OECD.a)
    U6.OECD.a=as.vector(H1.OECD.a)
    mu1.OECD.a=mean(U1.OECD.a)
    mu2.OECD.a=mean(U2.OECD.a)
    mu3.OECD.a=mean(U3.OECD.a)
    mu4.OECD.a=mean(U4.OECD.a)
    mu5.OECD.a=mean(U5.OECD.a)
    mu6.OECD.a=mean(U6.OECD.a)
    Delta.OECD.logs.a=-0.5*(log(mu1.OECD.a)-log(mu3.OECD.a)+log(mu2.OECD.a)-log(mu4.OECD.a))+log(mu5.OECD)-log(mu6.OECD)
    Delta.OECD.a=exp(Delta.OECD.logs.a)
    #
    c11.NOECD.a=c11a[-OECD.a]
    c22.NOECD.a=c22a[-OECD.a]
    c12.NOECD.a=c12a[-OECD.a]
    c21.NOECD.a=c21a[-OECD.a]
    H1.NOECD.a=H1a[-OECD.a]
    H2.NOECD.a=H2a[-OECD.a]
    U1.NOECD.a=(1+c21.NOECD.a)*as.vector(H2.NOECD.a)
    U2.NOECD.a=(1+c22.NOECD.a)*as.vector(H2.NOECD.a)
    U3.NOECD.a=(1+c11.NOECD.a)*as.vector(H1.NOECD.a)
    U4.NOECD.a=(1+c12.NOECD.a)*as.vector(H1.NOECD.a)
    U5.NOECD.a=as.vector(H2.NOECD.a)
    U6.NOECD.a=as.vector(H1.NOECD.a)
    mu1.NOECD.a=mean(U1.NOECD.a)
    mu2.NOECD.a=mean(U2.NOECD.a)
    mu3.NOECD.a=mean(U3.NOECD.a)
    mu4.NOECD.a=mean(U4.NOECD.a)
    mu5.NOECD.a=mean(U5.NOECD.a)
    mu6.NOECD.a=mean(U6.NOECD.a)
    Delta.NOECD.logs.a=-0.5*(log(mu1.NOECD.a)-log(mu3.NOECD.a)+log(mu2.NOECD.a)-log(mu4.NOECD.a))+log(mu5.NOECD)-log(mu6.NOECD)
    Delta.NOECD.a=exp(Delta.NOECD.logs.a)
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
    Delta.b.logs=-0.5*(log(mu1b)-log(mu3b)+log(mu2b)-log(mu4b))+log(mu5)-log(mu6)
    Delta.b=exp(Delta.b.logs)
    #
    ## For OECD and Non-OECD
    OECD.b=which(ind1[(na+1):n] %in% OECD)
    c11.OECD.b=c11b[OECD.b]
    c22.OECD.b=c22b[OECD.b]
    c12.OECD.b=c12b[OECD.b]
    c21.OECD.b=c21b[OECD.b]
    H1.OECD.b=H1b[OECD.b]
    H2.OECD.b=H2b[OECD.b]
    U1.OECD.b=(1+c21.OECD.b)*as.vector(H2.OECD.b)
    U2.OECD.b=(1+c22.OECD.b)*as.vector(H2.OECD.b)
    U3.OECD.b=(1+c11.OECD.b)*as.vector(H1.OECD.b)
    U4.OECD.b=(1+c12.OECD.b)*as.vector(H1.OECD.b)
    U5.OECD.b=as.vector(H2.OECD.b)
    U6.OECD.b=as.vector(H1.OECD.b)
    mu1.OECD.b=mean(U1.OECD.b)
    mu2.OECD.b=mean(U2.OECD.b)
    mu3.OECD.b=mean(U3.OECD.b)
    mu4.OECD.b=mean(U4.OECD.b)
    mu5.OECD.b=mean(U5.OECD.b)
    mu6.OECD.b=mean(U6.OECD.b)
    Delta.OECD.logs.b=-0.5*(log(mu1.OECD.b)-log(mu3.OECD.b)+log(mu2.OECD.b)-log(mu4.OECD.b))+log(mu5.OECD)-log(mu6.OECD)
    Delta.OECD.b=exp(Delta.OECD.logs.b)
    #
    c11.NOECD.b=c11b[-OECD.b]
    c22.NOECD.b=c22b[-OECD.b]
    c12.NOECD.b=c12b[-OECD.b]
    c21.NOECD.b=c21b[-OECD.b]
    H1.NOECD.b=H1b[-OECD.b]
    H2.NOECD.b=H2b[-OECD.b]
    U1.NOECD.b=(1+c21.NOECD.b)*as.vector(H2.NOECD.b)
    U2.NOECD.b=(1+c22.NOECD.b)*as.vector(H2.NOECD.b)
    U3.NOECD.b=(1+c11.NOECD.b)*as.vector(H1.NOECD.b)
    U4.NOECD.b=(1+c12.NOECD.b)*as.vector(H1.NOECD.b)
    U5.NOECD.b=as.vector(H2.NOECD.b)
    U6.NOECD.b=as.vector(H1.NOECD.b)
    mu1.NOECD.b=mean(U1.NOECD.b)
    mu2.NOECD.b=mean(U2.NOECD.b)
    mu3.NOECD.b=mean(U3.NOECD.b)
    mu4.NOECD.b=mean(U4.NOECD.b)
    mu5.NOECD.b=mean(U5.NOECD.b)
    mu6.NOECD.b=mean(U6.NOECD.b)
    Delta.NOECD.logs.b=-0.5*(log(mu1.NOECD.b)-log(mu3.NOECD.b)+log(mu2.NOECD.b)-log(mu4.NOECD.b))+log(mu5.NOECD)-log(mu6.NOECD)
    Delta.NOECD.b=exp(Delta.NOECD.logs.b)
    #
    tbar0=tbar0+(Delta.a+Delta.b)/2-Delta
    tbar0.OECD=tbar0.OECD+(Delta.OECD.a+Delta.OECD.b)/2-Delta.OECD
    tbar0.NOECD=tbar0.NOECD+(Delta.NOECD.a+Delta.NOECD.b)/2-Delta.NOECD
    #
    tbar[ind1[1:na]]=tbar[ind1[1:na]] +
      delta.a - delta[ind1[1:na]]
    tbar[ind1[(na+1):n]]=tbar[ind1[(na+1):n]] +
      delta.b - delta[ind1[(na+1):n]]
  }
  # tbar contains the bias for eff of each obs
  tbar=(1/L)*bc.fac*tbar
  tbar0=(1/L)*bc.fac*tbar0
  tbar0.OECD=(1/L)*bc.fac*tbar0.OECD
  tbar0.NOECD=(1/L)*bc.fac*tbar0.NOECD
  #
  Delta.bias=tbar0
  ## For OECD
  Delta.OECD.bias=tbar0.OECD
  ## For Non-OECD
  Delta.NOECD.bias=tbar0.NOECD
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
    Delta.nk.logs=0.5*(log(mu1.nk)-log(mu3.nk)+log(mu2.nk)-log(mu4.nk))+log(mu6.nk)-log(mu5.nk)
    Delta.nk=exp(Delta.nk.logs)
    # OECD
    nk.OECD=floor(n.OECD**(2*kappa))
    U1.OECD.nk=(1+c21.OECD[1:nk.OECD])*as.vector(H2.OECD[1:nk.OECD])
    U2.OECD.nk=(1+c22.OECD[1:nk.OECD])*as.vector(H2.OECD[1:nk.OECD])
    U3.OECD.nk=(1+c11.OECD[1:nk.OECD])*as.vector(H1.OECD[1:nk.OECD])
    U4.OECD.nk=(1+c12.OECD[1:nk.OECD])*as.vector(H1.OECD[1:nk.OECD])
    U5.OECD.nk=as.vector(H2.OECD[1:nk.OECD])
    U6.OECD.nk=as.vector(H1.OECD[1:nk.OECD])
    mu1.OECD.nk=mean(U1.OECD.nk)
    mu2.OECD.nk=mean(U2.OECD.nk)
    mu3.OECD.nk=mean(U3.OECD.nk)
    mu4.OECD.nk=mean(U4.OECD.nk)
    mu5.OECD.nk=mean(U5.OECD.nk)
    mu6.OECD.nk=mean(U6.OECD.nk)
    Delta.OECD.nk.logs=0.5*(log(mu1.OECD.nk)-log(mu3.OECD.nk)+log(mu2.OECD.nk)-log(mu4.OECD.nk))+log(mu6.OECD.nk)-log(mu5.OECD.nk)
    Delta.OECD.nk=exp(Delta.OECD.nk.logs)
    # NOECD
    # OECD
    nk.NOECD=floor(n.NOECD**(2*kappa))
    U1.NOECD.nk=(1+c21.NOECD[1:nk.NOECD])*as.vector(H2.NOECD[1:nk.NOECD])
    U2.NOECD.nk=(1+c22.NOECD[1:nk.NOECD])*as.vector(H2.NOECD[1:nk.NOECD])
    U3.NOECD.nk=(1+c11.NOECD[1:nk.NOECD])*as.vector(H1.NOECD[1:nk.NOECD])
    U4.NOECD.nk=(1+c12.NOECD[1:nk.NOECD])*as.vector(H1.NOECD[1:nk.NOECD])
    U5.NOECD.nk=as.vector(H2.NOECD[1:nk.NOECD])
    U6.NOECD.nk=as.vector(H1.NOECD[1:nk.NOECD])
    mu1.NOECD.nk=mean(U1.NOECD.nk)
    mu2.NOECD.nk=mean(U2.NOECD.nk)
    mu3.NOECD.nk=mean(U3.NOECD.nk)
    mu4.NOECD.nk=mean(U4.NOECD.nk)
    mu5.NOECD.nk=mean(U5.NOECD.nk)
    mu6.NOECD.nk=mean(U6.NOECD.nk)
    Delta.NOECD.nk.logs=0.5*(log(mu1.NOECD.nk)-log(mu3.NOECD.nk)+log(mu2.NOECD.nk)-log(mu4.NOECD.nk))+log(mu6.NOECD.nk)-log(mu5.NOECD.nk)
    Delta.NOECD.nk=exp(Delta.NOECD.nk.logs)
    #
    ts.nk=sig/sqrt(nk)
    ts.OECD.nk=sig.OECD/sqrt(nk.OECD)
    ts.NOECD.nk=sig.NOECD/sqrt(nk.NOECD)
    bounds2=matrix((Delta.nk-Delta.bias-ts.nk*crit),nrow=3,ncol=2)
    bounds2.OECD=matrix((Delta.OECD.nk-Delta.OECD.bias-ts.OECD.nk*crit),nrow=3,ncol=2)
    bounds2.NOECD=matrix((Delta.NOECD.nk-Delta.NOECD.bias-ts.NOECD.nk*crit),nrow=3,ncol=2)
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
             estimate_individual=c(exp(delta-tbar)),
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