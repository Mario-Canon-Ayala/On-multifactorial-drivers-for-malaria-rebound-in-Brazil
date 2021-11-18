########################################################################
#                                                                      #
#         Bayesian models and vivax incidence prediction (ACRE)        #
#                                                                      #
########################################################################
### Packages
library(foreign)
library(lubridate)
library(rgdal)
library('stringr')
library(maptools)
library(sp)
library(spdep)
library(INLA)
INLA:::inla.dynload.workaround()

setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/INLA data")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Subset Amazonas state
BASE_BRA_AMA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<140000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>129999)
ID.area.ama<-rep(1:62,192)
ID.mun.mes.ama<-1:11904
BASE_BRA_AMA2<-data.frame(BASE_BRA_AMA,ID.area.ama,ID.mun.mes.ama)
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AM")
Ama <- readOGR("AM", "13MUE250GC_SIR")
pol <- poly2nb(Ama)
nb2INLA("AMA.graph", pol) 
#this create a file called "AMA.adj" with the graph for INLA
AMA.adj <- paste(getwd(),"/AMA.graph",sep="")
BASE_BRA_AMAb<-BASE_BRA_AMA[1:9672,3:14]
pred<-data.frame(matrix(NA,2232,12))
names(pred)<-names(BASE_BRA_AMAb)
BASE_BRA_AMA2b<-rbind(BASE_BRA_AMAb,pred)
BASE_BRA_AMA3<-data.frame(BASE_BRA_AMA[,1:2],BASE_BRA_AMA2b,POP=BASE_BRA_AMA$POP,ID.area.ama,ID.mun.mes.ama)

###########################################################################################
#                                                                                         #
#                             Bayesian model                                              #
#                                                                                         #
###########################################################################################
BASE_SIMPLES<-data.frame(BASE_BRA_AMA2,Y2=BASE_BRA_AMA3$Y,year=2002+ceiling(BASE_BRA_AMA2$i/12),month=BASE_BRA_AMA2$i-12*(ceiling(BASE_BRA_AMA2$i/12)-1),ID.area.ama2=BASE_BRA_AMA2$ID.area.ama)
######## Run code to make INLA to 2015 and predict 2016, 2017 and 2018 without covariables
formula.4<- Y2 ~ 1 +  f(month, model = "rw2", constr = T, cyclic = T) + f(year, model = "iid", constr = T) +  f(ID.area.ama, model = "iid")  + f(i,model="rw1")
t <- proc.time()
model.inla.5 <- inla(formula.4,family="poisson",data=BASE_SIMPLES,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))#,waic=TRUE))
proc.time()-t
#save(model.inla.5,file="AM_after_INLA_mes_viv_Ypred2016_2018_simples+spa_month_year_ate_2018.rdata")
load(file="AM_after_INLA_mes_viv_Ypred2016_2018_simples+spa_month_year_ate_2018.rdata")
pred5<-exp(model.inla.5[["summary.linear.predictor"]][["mean"]]) 

############# real incidence vs model prevision
I0<-data.frame(i=BASE_BRA_AMA$i,Y=BASE_BRA_AMA$Y,Yp5=pred5*BASE_BRA_AMA$POP/100000)
Inci<-NULL
for (j in 1:192) {
  base<-subset(I0,I0$i==j)
  Y_<-sum(base$Y)
  Y5<-sum(base$Yp5)
  col<-data.frame(i=j,Y=Y_,Pre=Y5)
  Inci<-rbind(Inci,col)
}
# save results
#save(Inci,file="AM_pre_viv.rdata")
###########################################################################################
#                                                                                         #
#                                 Error map                                               #
#                                                                                         #
###########################################################################################
Error<-data.frame(i=BASE_BRA_AMA$i,ID.area=BASE_BRA_AMA$ID.area,Y=BASE_BRA_AMA$Y,Yp5=pred5*BASE_BRA_AMA$POP/100000,diff=BASE_BRA_AMA$Y-pred5*BASE_BRA_AMA$POP/100000)
matrx_err<-matrix(Error$diff,max(BASE_BRA_AMA2$ID.area.ama),192)
matrx_err2<-matrx_err[,157:192]
BASSER<-NULL
for (i in 1:max(BASE_BRA_AMA2$ID.area.ama)) {
  tab<-matrx_err2[i,]
  Y2016<-mean(tab[1:12])
  Y2017<-mean(tab[13:24])
  Y2018<-mean(tab[25:36])
  record<-data.frame(mun=i,ID.area=BASE_BRA_AMA2$ID.area[i],EP2016=Y2016,EP2017=Y2017,EP2018=Y2018)
  BASSER<-rbind(BASSER,record)
}
# Save results
#save(BASSER,file="MAPER_AM_viv.rdata")
