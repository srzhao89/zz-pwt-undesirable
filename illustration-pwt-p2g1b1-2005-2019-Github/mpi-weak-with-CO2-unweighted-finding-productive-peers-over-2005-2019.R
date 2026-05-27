###########################################################################################
## When computing MPI with undesirable outputs, we use CRS-DEA estimator
## with Kappa=2/(p+qg+qb+1)
###########################################################################################

require(Rglpk)
require(readxl)

source("./Functions/coverage.mpi.kuosmanen.illu.unweighted.R")
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

year1=c(2005)
year2=c(2019)

#res=matrix(0,nrow=Year_end-Year_start+1,ncol=11)
res=matrix(0,nrow=1,ncol=10)
res.star=rep(0,length=1)
res.OECD.star=rep(0,length=1)
res.NOECD.star=rep(0,length=1)
res.that=rep(0,length=1)

for (i in 1:length(year1)) {
  
  cat (i,"\n")
  
  i1=which(df$Year==2005)
  df1=df[i1,]
  x1=cbind(df1$Energy_per_emp,df1$Capital_per_emp)
  g1=matrix(df1$rgdpo_per_emp,ncol=1)
  b1=matrix(df1$CO2_per_emp,ncol=1)
  ##
  OECD=which(df1$Group=="OECD")
  #NOECD=which(df1$Group=="Non_OECD")
  ##
    
  i2=which(df$Year==2019)
  df2=df[i2,]
  x2=cbind(df2$Energy_per_emp,df2$Capital_per_emp)
  g2=matrix(df2$rgdpo_per_emp,ncol=1)
  b2=matrix(df2$CO2_per_emp,ncol=1)
    
  #
  XDIREC1=x1
  GDIREC1=g1
  BDIREC1=b1
  XDIREC2=x2
  GDIREC2=g2
  BDIREC2=b2
  #
  tt=coverage.mpi.kuosmanen.illu.unweighted(x1=x1,g1=g1,b1=b1,x2=x2,g2=g2,b2=b2,
                                            XDIREC1=XDIREC1,GDIREC1=GDIREC1,BDIREC1=BDIREC1,
                                            XDIREC2=XDIREC2,GDIREC2=GDIREC2,BDIREC2=BDIREC2,L=100)
  # 
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
  
  
  
}


ii1=order(tt$estimate_individual,decreasing = TRUE)
estimate_individual=tt$estimate_individual[ii1]
Countries=df1$Country_x[ii1]

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

library(stringr)
library(dplyr)


Countries= Countries %>%  str_to_lower()
Countries=sapply(Countries, simpleCap)
Countries=as.vector(Countries)

Countries

tex=formatC(Countries,width=max(nchar(Countries)))
tab=formatC(estimate_individual, digits = 3, format = "f")


tex=paste(tex,"&",tab)
tex


tex=paste(tex,"\\\\")
write(tex,file="./Output/prod-2005-2019.tex")




