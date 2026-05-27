###########################################################################################
## Title: Central Limit Theorems for Aggregates of Directional Distance Functions
## Authors: Leopold Simar, Valentin Zelenyuk, and Shirong Zhao
## Date: January 22, 2024
##############
## This function performs VRS-DEA estimates for the DDF under weak disposability
## i.e., obtaining the estimates of delta using equation (3.6) in the paper
## XOBS is n times p matrix for inputs 
## GOBS is n times qg matrix for desirable outputs
## BOBS is n times qb matrix for undesirable outputs
###########################################################################################

require(Rglpk)

dea.direc.kuosmanen<- function(XOBS,GOBS,BOBS,XDIREC,GDIREC,BDIREC,XREF=NULL,GREF=NULL,BREF=NULL) {
  
  if (is.null(XREF)){
    XREF=XOBS
    GREF=GOBS
    BREF=BOBS
  }
  
  n=nrow(XREF)
  p=ncol(XREF)
  qg=ncol(GREF)
  qb=ncol(BREF)
  
  n1=nrow(XOBS)
  
  res=rep(NA,n1)
  
  for (i in 1:n1) {
    
    xgbi=c(XOBS[i,],GOBS[i,],BOBS[i,])
    XDRECi=XDIREC[i,]
    GDRECi=GDIREC[i,]
    BDRECi=BDIREC[i,]
    
    f.obj=c(rep(0,2*n),1)
    
    gblock=rbind(GREF,matrix(0,nrow=nrow(GREF),ncol=ncol(GREF)),-1*GDRECi)
    bblock=rbind(BREF,matrix(0,nrow=nrow(BREF),ncol=ncol(BREF)),BDRECi)
    xblock=rbind(XREF,XREF,XDRECi)
    xgbblock=rbind(t(xblock),t(gblock),t(bblock))
    f.con=rbind(xgbblock,c(rep(1,2*n),0))
    
    # Set unequality signs
    f.dir=c(rep("<=",p),rep(">=",qg),rep("==",qb),"==")
    
    f.rhs=c(xgbi,1)
    
    f.bounds <- list(lower = list(ind = 2*n+1, val = -Inf))
    
    results=Rglpk::Rglpk_solve_LP(f.obj,f.con,f.dir,f.rhs,f.bounds,max=TRUE)
    
    res[i]=results$optimum
    
  }
  return(res)
}