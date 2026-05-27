###########################################################################################
## When computing MPI with undesirable outputs, we use CRS-DEA estimator
## with Kappa=2/(p+qg+qb+1)
###########################################################################################

require(Rglpk)
require(readxl)

source("./Functions/coverage.mpi.kuosmanen.illu.weighted.R")
source("./Functions/dea.direc.kuosmanen.crs.R")

if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)
######################################################
df<- read_excel("Data/Combined_data_with_price0518.xlsx")
######################################################


Year_min=min(df$Year)
Year_max=max(df$Year)

for (ii in Year_min:Year_max) {
  jj=length(which(df$Year==ii))
  cat(ii,":",jj,"\n")
}

index=df$index
Year=df$Year


#y = df$rgdpo_per_emp
#b = df$CO2_per_emp
#x1 = df$Energy_per_emp
#x2 = df$Capital_per_emp

#w_y=1
#w_b=df$CO2_price_2017
#w_x1=df$Energy_Price
#w_x2=df$pl_n

#y=matrix(y,ncol=1)
#b=matrix(b,ncol=1)
#x1=matrix(x1,ncol=1)
#x2=matrix(x2,ncol=1)

#w_y=matrix(w_y,ncol=1,nrow=nrow(y))
#w_b=matrix(w_b,ncol=1,nrow=nrow(b))
#w_x1=matrix(w_x1,ncol=1,nrow=nrow(x1))
#w_x2=matrix(w_x2,ncol=1,nrow=nrow(x2))

#
#Year_start = 2005
#Year_end   = 2019
#Year_pre = Year_start-1

year1=c(seq(2005,2018,1),2005,2010,2015,2005)
year2=c(seq(2006,2019,1),2010,2015,2019,2019)

#res=matrix(0,nrow=Year_end-Year_start+1,ncol=11)
res=matrix(0,nrow=length(year1),ncol=10)
res.star=rep(0,length=length(year1))
res.OECD.star=rep(0,length=length(year1))
res.NOECD.star=rep(0,length=length(year1))
res.that=rep(0,length=length(year1))

for (i in 1:length(year1)) {
  
  cat (i,"\n")
  
  i1=which(df$Year==year1[i])
  df1=df[i1,]
  x1=cbind(df1$Energy_per_emp,df1$Capital_per_emp)
  g1=matrix(df1$rgdpo_per_emp,ncol=1)
  b1=matrix(df1$CO2_per_emp,ncol=1)
  # computing H1
  wx1=cbind(df1$Energy_Price,df1$pl_n)
  wg1=matrix(1,nrow=nrow(g1),ncol=1)
  wb1=matrix(df1$CO2_price_2017,ncol=1)
  #
  H1=apply(x1*wx1,1,sum)+apply(g1*wg1+b1*wb1,1,sum)
  ##
  OECD=which(df1$Group=="OECD")
  #NOECD=which(df1$Group=="Non_OECD")
  ##
    
  i2=which(df$Year==year2[i])
  df2=df[i2,]
  x2=cbind(df2$Energy_per_emp,df2$Capital_per_emp)
  g2=matrix(df2$rgdpo_per_emp,ncol=1)
  b2=matrix(df2$CO2_per_emp,ncol=1)
  # computing H2
  wx2=cbind(df2$Energy_Price,df2$pl_n)
  wg2=matrix(1,nrow=nrow(g2),ncol=1)
  wb2=matrix(df2$CO2_price_2017,ncol=1)
  #
  H2=apply(x2*wx2,1,sum)+apply(g2*wg2+b2*wb2,1,sum)
    
  #
  XDIREC1=x1
  GDIREC1=g1
  BDIREC1=b1
  XDIREC2=x2
  GDIREC2=g2
  BDIREC2=b2
  #
  tt=coverage.mpi.kuosmanen.illu.weighted(x1=x1,g1=g1,b1=b1,x2=x2,g2=g2,b2=b2,
                                            XDIREC1=XDIREC1,GDIREC1=GDIREC1,BDIREC1=BDIREC1,
                                            XDIREC2=XDIREC2,GDIREC2=GDIREC2,BDIREC2=BDIREC2,
                                            H1=H1,H2=H2,L=100)
  # 
  #res[i,1:3]=tt$estimate[1:3]
  #res[i,4:6]=tt$estimate.OECD[1:3]
  #res[i,7:9]=tt$estimate.NOECD[1:3]
  res[i,1:3]=tt$estimate[4:6]
  res[i,4:6]=tt$estimate.OECD[4:6]
  res[i,7:9]=tt$estimate.NOECD[4:6]
  res[i,10]=tt$estimate.OECD[6]-tt$estimate.NOECD[6]
  res.that[i]=tt$pval
  #res[i,7]=tt$sig
  #res[i,8:9]=tt$bounds0[2,]
  #res[i,4:5]=tt$bounds2[2,]
  
  # obtain the individual estimates
  # note that these are bias-corrected estimates
  print(tt$estimate_individual)
  
  #
  if (tt$bounds2[3,1]>1 | 1>tt$bounds2[3,2]){
    res.star[i]=3
  } else if (tt$bounds2[2,1]>1 | 1>tt$bounds2[2,2]){
    res.star[i]=2
  } else if  (tt$bounds2[1,1]>1 | 1>tt$bounds2[1,2]){
    res.star[i]=1
  } else {
    res.star[i]=0
  }
  
  if (tt$bounds2.OECD[3,1]>1 | 1>tt$bounds2.OECD[3,2]){
    res.OECD.star[i]=3
  } else if (tt$bounds2.OECD[2,1]>1 | 1>tt$bounds2.OECD[2,2]){
    res.OECD.star[i]=2
  } else if  (tt$bounds2.OECD[1,1]>1 | 1>tt$bounds2.OECD[1,2]){
    res.OECD.star[i]=1
  } else {
    res.OECD.star[i]=0
  }
  
  if (tt$bounds2.NOECD[3,1]>1 | 1>tt$bounds2.NOECD[3,2]){
    res.NOECD.star[i]=3
  } else if (tt$bounds2.NOECD[2,1]>1 | 1>tt$bounds2.NOECD[2,2]){
    res.NOECD.star[i]=2
  } else if  (tt$bounds2.NOECD[1,1]>1 | 1>tt$bounds2.NOECD[1,2]){
    res.NOECD.star[i]=1
  } else {
    res.NOECD.star[i]=0
  }
  
}

round(res,3)

res.star

# make latex code for table:

#tab=matrix(nrow=nrow(res),ncol=ncol(res))
tab=formatC(res,width=6,digits=3,format="f")


w1=formatC(res[,3],width=6,digits=3,format="f")
for (i in 1:nrow(res)) {
  
  if (res.star[i]==3) {
    tab[i,3]=paste(w1[i],"$^{***}$",sep="")
  } else if (res.star[i]==2) {
    tab[i,3]=paste(w1[i],"$^{** }$",sep="")
  } else if (res.star[i]==1) {
    tab[i,3]=paste(w1[i],"$^{*  }$",sep="")
  } else {
    tab[i,3]=paste(w1[i],"$^{   }$",sep="")
  }
}


w1.OECD=formatC(res[,6],width=6,digits=3,format="f")
for (i in 1:nrow(res)) {
  
  if (res.OECD.star[i]==3) {
    tab[i,6]=paste(w1.OECD[i],"$^{***}$",sep="")
  } else if (res.OECD.star[i]==2) {
    tab[i,6]=paste(w1.OECD[i],"$^{** }$",sep="")
  } else if (res.OECD.star[i]==1) {
    tab[i,6]=paste(w1.OECD[i],"$^{*  }$",sep="")
  } else {
    tab[i,6]=paste(w1.OECD[i],"$^{   }$",sep="")
  }
}


w1.NOECD=formatC(res[,9],width=6,digits=3,format="f")
for (i in 1:nrow(res)) {
  
  if (res.NOECD.star[i]==3) {
    tab[i,9]=paste(w1.NOECD[i],"$^{***}$",sep="")
  } else if (res.NOECD.star[i]==2) {
    tab[i,9]=paste(w1.NOECD[i],"$^{** }$",sep="")
  } else if (res.NOECD.star[i]==1) {
    tab[i,9]=paste(w1.NOECD[i],"$^{*  }$",sep="")
  } else {
    tab[i,9]=paste(w1.NOECD[i],"$^{   }$",sep="")
  }
}


ii=c(10)
for (k in ii) {
  for (j in 1:nrow(tab)) {
    if (substr(tab[j,k],1,1)==" ") {
      tab[j,k]=paste("\\phs",substr(tab[j,k],2,6),sep="")
    } else {
      tab[j,k]=paste("   ",tab[j,k],sep="")
    }
  }
}


for (i in 1:nrow(res)) {
  
  if (res.that[i]<0.01) {
    tab[i,10]=paste(tab[i,10],"$^{***}$","        ",sep="")
  } else if (res.that[i]<0.05) {
    tab[i,10]=paste(tab[i,10],"$^{** }$","\\phast  ",sep="")
  } else if (res.that[i]<0.10) {
    tab[i,10]=paste(tab[i,10],"$^{*  }$","\\phastt ",sep="")
  } else {
    tab[i,10]=paste(tab[i,10],"$^{   }$","\\phasttt",sep="")
  }
}


### construct the Table ###
#year=Year_start:Year_end
tex1=formatC(year1,width=7,digits=0,format="f")
tex2=formatC(year2,width=7,digits=0,format="f")
tex=paste(tex1,"&",tex2)

for (k in 1:ncol(res)) {
  tex = paste(tex,"&",tab[,k])
}


tex = paste(tex,"\\\\")
write(tex,file="./Output/coverage-mpi-weak-with-CO2-weighted.tex")



