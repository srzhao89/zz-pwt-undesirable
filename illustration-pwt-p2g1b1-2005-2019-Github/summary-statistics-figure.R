###########################################################################################
## When computing MPI with undesirable outputs, we use CRS-DEA estimator
## with Kappa=2/(p+qg+qb+1)
###########################################################################################

require(Rglpk)
require(readxl)
require(dplyr)

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
df$Group[df$Group=="Non_OECD"]="Non-OECD"


#df$emp*df$Energy_per_emp/df$`Primary energy consumption (TWh)`


ii=which(df$Year==2019)
df2019=df[ii,]

df2019$CO2_per_emp
df2019$Energy_per_emp
df2019$Country_x

df2019$Group





###### Change countries name to a proper name, like CHINA to China
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

library(stringr)
library(dplyr)
Countries= df2019$Country_x %>%  str_to_lower()
Countries=sapply(Countries, simpleCap)
df2019$Country_x=as.vector(Countries)





######################################################################

library(ggplot2)
library(ggrepel)


# Example: create your data 
df <- data.frame(
  Country = df2019$Country_x,
  Group = df2019$Group,
  Energy_per_worker = df2019$Energy_per_emp,
  CO2_per_worker = df2019$CO2_per_emp
)

# Scatter plot
ggplot(df, aes(x = Energy_per_worker, y = CO2_per_worker,
               color = Group, shape = Group, label = Country)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 30, box.padding = 0.3) +
  labs(
    x = "Energy per Worker (TWh per million workers)",
    y = "CO2 per Worker (Kt per million workers)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  ) +
  expand_limits (x=0, y=0)

# Save to EPS
ggsave("./Output/plot-per-worker.pdf", device = "pdf",
       width = 7, height = 7, units = "in")



######################################################################

library(ggplot2)
library(ggrepel)


# Example: create your data 
df <- data.frame(
  Country = df2019$Country_x,
  Group = df2019$Group,
  Energy = log(df2019$`Primary energy consumption (TWh)`),
  CO2 = log(df2019$CO2)
)

# Scatter plot
ggplot(df, aes(x = Energy, y = CO2,
               color = Group, shape = Group, label = Country)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 30, box.padding = 0.3) +
  labs(
    x = "Log of Energy (TWh)",
    y = "Log of CO2 (Kt)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

#+ 
#expand_limits (x=6, y=3)+
#ylim(3, 10)

# Save to EPS
ggsave("./Output/plot.pdf", device = "pdf",
       width = 7, height = 7, units = "in")






