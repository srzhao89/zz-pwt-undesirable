###########################################################################################
###########################################################################################

require(Rglpk)
require(readxl)

source("./Functions/coverage.mpi.without.CO2.weighted.R")

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



year1=c(seq(2005,2018,1),2005,2010,2015,2005)
year2=c(seq(2006,2019,1),2010,2015,2019,2019)

#res=matrix(0,nrow=Year_end-Year_start+1,ncol=11)
res=matrix(0,nrow=length(year1),ncol=3)
res.star=rep(0,length=length(year1))


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
  H1=apply(x1*wx1,1,sum)+apply(g1*wg1,1,sum)

    
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
  H2=apply(x2*wx2,1,sum)+apply(g2*wg2,1,sum)
    
  #
  XDIREC1=x1
  GDIREC1=g1
  BDIREC1=b1
  XDIREC2=x2
  GDIREC2=g2
  BDIREC2=b2
  #
  y1=g1
  y2=g2
  YDIREC1=g1
  YDIREC2=g2
  #
  tt=coverage.mpi.without.CO2.weighted(x1=x1,y1=y1,x2=x2,y2=y2,
                                XDIREC1=XDIREC1,YDIREC1=YDIREC1,
                                XDIREC2=XDIREC2,YDIREC2=YDIREC2,
                                H1=H1,H2=H2,L=100)
  # 
  res[i,1:3]=tt$estimate[1:3]
  
  #
  if (tt$bounds1[3,1]>1 | 1>tt$bounds1[3,2]){
    res.star[i]=3
  } else if (tt$bounds1[2,1]>1 | 1>tt$bounds1[2,2]){
    res.star[i]=2
  } else if  (tt$bounds1[1,1]>1 | 1>tt$bounds1[1,2]){
    res.star[i]=1
  } else {
    res.star[i]=0
  }

  
}

round(res,3)

res.star

# make latex code for table:

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




### construct the Table ###
#year=Year_start:Year_end
tex1=formatC(year1,width=7,digits=0,format="f")
tex2=formatC(year2,width=7,digits=0,format="f")
tex=paste(tex1,"&",tex2)

for (k in 1:ncol(res)) {
  tex = paste(tex,"&",tab[,k])
}


tex = paste(tex,"\\\\")
write(tex,file="./Output/coverage-mpi-without-CO2-weighted.tex")



