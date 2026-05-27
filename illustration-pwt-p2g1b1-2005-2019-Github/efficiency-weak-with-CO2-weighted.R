###########################################################################################
###########################################################################################

remove(list=ls())
require(Rglpk)
require(readxl)

source("./Functions/coverage.ddf.kuosmanen.illu.weighted.R")
source("./Functions/dea.direc.kuosmanen.R")

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
#Year_start= 2005
#Year_end  = 2019
#Year_pre = Year_start-1

year=seq(2005,2019,1)
#res=matrix(0,nrow=Year_end-Year_start+1,ncol=11)
res=matrix(0,nrow=length(year),ncol=10)
res.that=rep(0,length=length(year))

for (i in 1:length(year)) {
  
  cat(i,"\n")
  ii=which(df$Year==year[i])
  df1=df[ii,]
  
  g=matrix(df1$rgdpo_per_emp,ncol=1)
  b=matrix(df1$CO2_per_emp,ncol=1)
  #x1=matrix(df1$Energy_per_emp,ncol=1)
  #x2=matrix(df1$Capital_per_emp,ncol=1)
  x=cbind(df1$Energy_per_emp,df1$Capital_per_emp)
  
  # computing H1
  wx=cbind(df1$Energy_Price,df1$pl_n)
  wg=matrix(1,nrow=nrow(g),ncol=1)
  wb=matrix(df1$CO2_price_2017,ncol=1)
  #
  H2=apply(x*wx,1,sum)+apply(g*wg+b*wb,1,sum)
  
  ##
  OECD=which(df1$Group=="OECD")
  #
  XDIREC=x
  GDIREC=g
  BDIREC=b
  #
  tt=coverage.ddf.kuosmanen.illu.weighted(x=x,g=g,b=b,XDIREC=XDIREC,GDIREC=GDIREC,BDIREC=BDIREC,H2=H2,L=100)
  # 
  res[i,1:3]=tt$estimate[4:6]
  res[i,4:6]=tt$estimate.OECD[4:6]
  res[i,7:9]=tt$estimate.NOECD[6]
  res[i,10]=tt$estimate.OECD[6]-tt$estimate.NOECD[6]
  res.that[i]=tt$pval
  #res[i,10:11]=c(tt$that,tt$pval)
  #res[i-Year_pre,4:5]=tt$bounds2[2,]
  
  # obtain the individual estimates
  # note that these are bias-corrected estimates
  print(tt$estimate_individual)
  
}

round(res,3)

tab=formatC(res,width=6,digits=3,format="f")

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

tex=formatC(year,width=7,digits=0,format="f")

for (k in 1:ncol(res)) {
  tex = paste(tex,"&",tab[,k])
}

tex = paste(tex,"\\\\")
write(tex,file="./Output/coverage-eff-weak-with-CO2-weighted.tex")



