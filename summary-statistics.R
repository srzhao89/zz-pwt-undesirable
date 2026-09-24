
require(Rglpk)
require(readxl)
require(dplyr)


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


i1=which(df$Year==2019 & df$Group=="OECD")
df1=df[i1,]

i2=which(df$Year==2019 & df$Group=="Non_OECD")
df2=df[i2,]

df1$CO2_per_emp
df1$Energy_per_emp
df1$Country_x

## the amount of per worker for OECD
ii1=order(df1$Energy_per_emp,decreasing = TRUE)
cbind(df1$Country_x[ii1],df1$Energy_per_emp[ii1])

ii2=order(df1$CO2_per_emp,decreasing = TRUE)
cbind(df1$Country_x[ii2],df1$CO2_per_emp[ii2])

## the whole amount for OECD
ii1=order(df1$`Primary energy consumption (TWh)`,decreasing = TRUE)
cbind(df1$Country_x[ii1],df1$`Primary energy consumption (TWh)`[ii1])

ii2=order(df1$CO2,decreasing = TRUE)
cbind(df1$Country_x[ii2],df2$CO1[ii2])


## the amount of per worker for Non-OECD
jj1=order(df2$Energy_per_emp,decreasing = TRUE)
cbind(df2$Country_x[jj1],df2$Energy_per_emp[jj1])

jj2=order(df2$CO2_per_emp,decreasing = TRUE)
cbind(df2$Country_x[jj2],df2$CO2_per_emp[jj2])

## the whole amount for Non-OECD
jj1=order(df2$`Primary energy consumption (TWh)`,decreasing = TRUE)
cbind(df2$Country_x[jj1],df2$`Primary energy consumption (TWh)`[jj1])

jj2=order(df2$CO2,decreasing = TRUE)
cbind(df2$Country_x[jj2],df2$CO2[jj2])

######## construct the table

###############  for OECD  ########################
ii1=order(df1$Energy_per_emp,decreasing = TRUE)

res.OECD=cbind(df1$Energy_per_emp[ii1],df1$CO2_per_emp[ii1],df1$`Primary energy consumption (TWh)`[ii1],df1$CO2[ii1])

tab=formatC(res.OECD,width=10,digits=3,format="f")



rank1=rank(desc(df1$Energy_per_emp[ii1]))
rank2=rank(desc(df1$CO2_per_emp[ii1]))
rank3=rank(desc(df1$`Primary energy consumption (TWh)`[ii1]))
rank4=rank(desc(df1$CO2[ii1]))
rank.all=cbind(rank1,rank2,rank3,rank4)
tab.rank=formatC(rank.all,width=3)
tab.rank


###### Change countries name to a proper name, like CHINA to China
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

library(stringr)
library(dplyr)
OECD.countries= df1$Country_x[ii1] %>%  str_to_lower()
OECD.countries=sapply(OECD.countries, simpleCap)
OECD.countries=as.vector(OECD.countries)
################

tex=formatC(OECD.countries,width=max(nchar(df1$Country_x[ii1])))

  
for (k in 1:ncol(tab)) {
  tex=paste(tex,"&",tab[,k], "&",tab.rank[,k])
}

tex=paste(tex,"\\\\")
write(tex,file="./Output/summary-stats-OECD.tex")





###############  for Non-OECD  ########################
jj1=order(df2$Energy_per_emp,decreasing = TRUE)
res.NonOECD=cbind(df2$Energy_per_emp[jj1],df2$CO2_per_emp[jj1],df2$`Primary energy consumption (TWh)`[jj1],df2$CO2[jj1])

tab=formatC(res.NonOECD,width=10,digits=3,format="f")

rank1=rank(desc(df2$Energy_per_emp[jj1]))
rank2=rank(desc(df2$CO2_per_emp[jj1]))
rank3=rank(desc(df2$`Primary energy consumption (TWh)`[jj1]))
rank4=rank(desc(df2$CO2[jj1]))
rank.all=cbind(rank1,rank2,rank3,rank4)
tab.rank=formatC(rank.all,width=3)
tab.rank


simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

library(stringr)
library(dplyr)
NonOECD.countries= df2$Country_x[jj1] %>%  str_to_lower()
NonOECD.countries=sapply(NonOECD.countries, simpleCap)
NonOECD.countries=as.vector(NonOECD.countries)

tex=formatC(NonOECD.countries,width=max(nchar(df2$Country_x[jj1])))

for (k in 1:ncol(tab)) {
  tex=paste(tex,"&",tab[,k], "&",tab.rank[,k])
}


tex=paste(tex,"\\\\")
write(tex,file="./Output/summary-stats-NonOECD.tex")
