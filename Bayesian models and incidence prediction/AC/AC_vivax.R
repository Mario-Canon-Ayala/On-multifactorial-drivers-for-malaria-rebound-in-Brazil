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

setwd("~/INLA data")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Subset Acre state
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<130000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>119999)

setwd("~/Bayesian models and incidence prediction/AC")
Acre <- readOGR("AC", "12MUE250GC_SIR")
pol <- poly2nb(Acre)
nb2INLA("AC.graph", pol) 
#this create a file called "PA.adj" with the graph for INLA
AC.adj <- paste(getwd(),"/AC.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(Acre@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1))

###########################################################################################
#                                                                                         #
#                             Bayesian model                                              #
#                                                                                         #
###########################################################################################
BASE_SIMPLES2<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),ID.area.sta2=BASE_SIMPLES$ID.area.sta)
######## Run code to make INLA to 2015 and predict 2016, 2017 and 2018 without covariables
formula.4<- Y2 ~ 1 +  f(month, model = "rw2", constr = T, cyclic = T) + f(year, model = "iid", constr = T) +  f(ID.area.sta, model = "iid") +f(ID.area.sta2,model="bym",graph=AC.adj) + f(i,model="rw1")
t <- proc.time()
model.inla.5 <- inla(formula.4,family="poisson",data=BASE_SIMPLES2,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))#,waic=TRUE))
proc.time()-t

############# model prevision
pred5<-exp(model.inla.5[["summary.linear.predictor"]][["mean"]]) 
############# real incidence vs model prevision
I0<-data.frame(i=BASE_BRA_STA2$i,Y=BASE_BRA_STA2$Y,Yp5=pred5*BASE_BRA_STA2$POP/100000)
Inci<-NULL
for (j in 1:192) {
  base<-subset(I0,I0$i==j)
  Y_<-sum(base$Y)
  Y5<-sum(base$Yp5)
  col<-data.frame(i=j,Y=Y_,Pre=Y5)
  Inci<-rbind(Inci,col)
}
# save results state
#save(Inci,file="AC_pre_viv.rdata")

###########################################################################################
#                                                                                         #
#                                 Error map                                               #
#                                                                                         #
###########################################################################################
Error<-data.frame(i=BASE_BRA_STA2$i,ID.area=BASE_BRA_STA2$ID.area.sta,Y=BASE_BRA_STA2$Y,Yp5=pred5*BASE_BRA_STA2$POP/100000,diff=BASE_BRA_STA2$Y-pred5*BASE_BRA_STA2$POP/100000)
matrx_err<-matrix(Error$diff,max(BASE_BRA_STA2$ID.area.sta),192)
matrx_err2<-matrx_err[,157:192]
BASSER<-NULL
for (i in 1:max(BASE_BRA_STA2$ID.area.sta)) {
  tab<-matrx_err2[i,]
  Y2016<-mean(tab[1:12])
  Y2017<-mean(tab[13:24])
  Y2018<-mean(tab[25:36])
  record<-data.frame(mun=i,ID.area=BASE_BRA_STA2$ID.area[i],EP2016=Y2016,EP2017=Y2017,EP2018=Y2018)
  BASSER<-rbind(BASSER,record)
}
# save error results
#save(BASSER,file="MAPER_AC_viv.rdata")

