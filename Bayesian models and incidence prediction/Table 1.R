# C?digo para hacer la corrida con covariables hasta 2015
### Packages
library(udunits2)
library(units)
library(sf)
library(foreign)
library(lubridate)
library(rgdal)
library('stringr')
library(maptools)
library(sp)
library(spdep)
library(INLA)
INLA:::inla.dynload.workaround() #Para correr y evitar actualizar el linux

#setwd("~/ownCloud2/malaria/Analises Mario/Dados/Brazil_todo")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Acre 1200000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<130000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>119999)

setwd("~/ownCloud2/malaria/Analises Mario/Dados/ACRE")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/ACRE")
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
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0

setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/ACRE")

######## ########

formula.1<- Y ~ 1 +  f(i,model="rw1") + f(ID.mun.ano, model="iid") + f(ID.area.sta,model="bym",graph=AC.adj) +despues*GAR_MIN
  despues*Desma + despues*AGRECU + despues*DOM #+ despues*TUR_V #+ despues*OUT+ despues*GAR_MIN 
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
#beep("mario")
## alpha fixed effects
R_18v<-model.inla.1$summary.fixed
R_18ve<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18v,file="R18v.rdata")
save(R_18ve,file="R18ve.rdata")

#rm(list = ls())

#######################################################################################
#
#                               Falciparum
#
#######################################################################################
#setwd("~/ownCloud2/malaria/Analises Mario/Dados/Brazil_todo")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_fal.rdata", verbose=TRUE) #falciparum data
#### Acre 1200000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<130000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>119999)

setwd("~/ownCloud2/malaria/Analises Mario/Dados/ACRE")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/ACRE")
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
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}

setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/ACRE")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AC.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + despues*DOM
  despues*Desma + despues*AGRECU #+ despues*DOM #+ despues*TUR_V + despues*GAR_MIN + despues*OUT
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
beep("mario")
## alpha fixed effects
R_18f<-model.inla.1$summary.fixed
R_18fe<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18f,file="R18f.rdata")
save(R_18fe,file="R18fe.rdata")
#################################################################################
#
#                          Rondonia
#
#################################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Rondonia 1100000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<120000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>109999)
setwd("~/ownCloud2/malaria/Analises Mario/Dados/RO")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RO")
ron <- readOGR("RO", "11MUE250GC_SIR")
pol <- poly2nb(ron)
nb2INLA("RO.graph", pol) 
#this create a file called "RO.adj" with the graph for INLA
RO.adj <- paste(getwd(),"/RO.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(ron@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RO")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=RO.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + 
  despues*Desma + despues*AGRECU + despues*DOM + despues*TUR_V + despues*GAR_MIN + despues*OUT
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
beep("mario")
## alpha fixed effects
R_18v<-model.inla.1$summary.fixed
R_18ve<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18v,file="R18v.rdata")
save(R_18ve,file="R18ve.rdata")
##########################################################################
#
#                                falciparum
###########################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_fal.rdata", verbose=TRUE) #vivax data
#### Rondonia 1100000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<120000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>109999)
setwd("~/ownCloud2/malaria/Analises Mario/Dados/RO")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RO")
ron <- readOGR("RO", "11MUE250GC_SIR")
pol <- poly2nb(ron)
nb2INLA("RO.graph", pol) 
#this create a file called "RO.adj" with the graph for INLA
RO.adj <- paste(getwd(),"/RO.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(ron@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RO")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=RO.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + 
  despues*Desma + despues*AGRECU #+ despues*DOM #+ despues*TUR_V + despues*GAR_MIN + despues*OUT
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
beep("mario")
## alpha fixed effects
R_18f<-model.inla.1$summary.fixed
R_18fe<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18f,file="R18f.rdata")
save(R_18fe,file="R18fe.rdata")

###################################################################################3
#
#                                 Roraima
#
###################################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Roraima 1400000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<150000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>139999)
setwd("~/ownCloud2/malaria/Analises Mario/Dados/RR")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RR")
RR <- readOGR("RR", "14MUE250GC_SIR")
pol <- poly2nb(RR)
nb2INLA("RR.graph", pol) 
#this create a file called "RR.adj" with the graph for INLA
RR.adj <- paste(getwd(),"/RR.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(RR@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RR")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=RR.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") +
  despues*Desma + despues*AGRECU + despues*DOM  + despues*OUT# + despues*TUR_V + despues*GAR_MIN
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
#beep("mario")
## alpha fixed effects
R_18v<-model.inla.1$summary.fixed
R_18ve<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18v,file="R18v.rdata")
save(R_18ve,file="R18ve.rdata")
#####################################################################3
#
#    Falciparum
#####################################################3
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_fal.rdata", verbose=TRUE) #vivax data
#### Roraima 1400000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<150000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>139999)
setwd("~/ownCloud2/malaria/Analises Mario/Dados/RR")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RR")
RR <- readOGR("RR", "14MUE250GC_SIR")
pol <- poly2nb(RR)
nb2INLA("RR.graph", pol) 
#this create a file called "RR.adj" with the graph for INLA
RR.adj <- paste(getwd(),"/RR.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(RR@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/RR")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=RR.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + 
  despues*Desma #+ despues*AGRECU + despues*DOM + despues*TUR_V + despues*GAR_MIN + despues*OUT
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
#beep("mario")
## alpha fixed effects
R_18f<-model.inla.1$summary.fixed
R_18fe<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18f,file="R18f.rdata")
save(R_18fe,file="R18fe.rdata")
###################################################################################
#
#                               Amapá
#
###################################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Amapa 1600000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<170000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>159999)
setwd("~/ownCloud2/malaria/Analises Mario/Dados/AP")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AP")
Amapa <- readOGR("AP", "16MUE250GC_SIR")
pol <- poly2nb(Amapa)
nb2INLA("AP.graph", pol) 
#this create a file called "PA.adj" with the graph for INLA
AP.adj <- paste(getwd(),"/AP.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(Amapa@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AP")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AP.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") +
  despues*Desma + despues*AGRECU + despues*DOM  + despues*GAR_MIN + despues*OUT # + despues*TUR_V
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
#beep("mario")
## alpha fixed effects
R_18v<-model.inla.1$summary.fixed
R_18ve<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18v,file="R18v.rdata")
save(R_18ve,file="R18ve.rdata")
#####################################################################################
#
#                              Falciparum
###################################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_fal.rdata", verbose=TRUE) #vivax data
#### Amapa 1600000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<170000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>159999)
setwd("~/ownCloud2/malaria/Analises Mario/Dados/AP")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AP")
Amapa <- readOGR("AP", "16MUE250GC_SIR")
pol <- poly2nb(Amapa)
nb2INLA("AP.graph", pol) 
#this create a file called "PA.adj" with the graph for INLA
AP.adj <- paste(getwd(),"/AP.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(Amapa@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AP")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AP.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + 
  despues*Desma + despues*AGRECU + despues*DOM  + despues*GAR_MIN + despues*OUT#+ despues*TUR_V
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
#beep("mario")
## alpha fixed effects
R_18f<-model.inla.1$summary.fixed
R_18fe<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18f,file="R18f.rdata")
save(R_18fe,file="R18fe.rdata")
#########################################################################################
#
#                              Amazonas
#
##########################################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data

BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<140000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>129999)
ID.area.ama<-rep(1:62,192)
ID.mun.mes.ama<-1:11904
BASE_BRA_AMA2<-data.frame(BASE_BRA_AMA,ID.area.ama,ID.mun.mes.ama)
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AM")
setwd("~/ownCloud2/malaria/Analises Mario/Dados/AM")
Ama <- readOGR("AM", "13MUE250GC_SIR")
pol <- poly2nb(Ama)
nb2INLA("AMA.graph", pol) 
#this create a file called "AMA.adj" with the graph for INLA
AMA.adj <- paste(getwd(),"/AMA.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(Ama@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AM")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AMA.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + 
  despues*Desma + despues*AGRECU + despues*DOM + despues*TUR_V + despues*GAR_MIN + despues*OUT
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
beep("mario")
## alpha fixed effects
R_18v<-model.inla.1$summary.fixed
R_18ve<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18v,file="R18v.rdata")
save(R_18ve,file="R18ve.rdata")
################################################################################
#                                   falciparum
################################################################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Brazil_todo")
load("BRASIL_INLA_mes_fal.rdata", verbose=TRUE) #vivax data

BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<140000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>129999)
ID.area.ama<-rep(1:62,192)
ID.mun.mes.ama<-1:11904
BASE_BRA_AMA2<-data.frame(BASE_BRA_AMA,ID.area.ama,ID.mun.mes.ama)
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AM")
setwd("~/ownCloud2/malaria/Analises Mario/Dados/AM")
Ama <- readOGR("AM", "13MUE250GC_SIR")
pol <- poly2nb(Ama)
nb2INLA("AMA.graph", pol) 
#this create a file called "AMA.adj" with the graph for INLA
AMA.adj <- paste(getwd(),"/AMA.graph",sep="")
############# Order data ############
BASE_STATE<-NULL
for (j in 1:192) {
  BASEORD<-subset(BASE_BRA_STA,BASE_BRA_STA$i==j)
  order<-match(Ama@data[["CD_GEOCMU"]],BASEORD$muni_cut.CD_GEOCMU.j.)
  BASEORD<-BASEORD[order,]  
  BASE_STATE<-rbind(BASE_STATE,BASEORD)
}
ID.area.sta<-rep(1:sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0),192)
ID.mun.mes.sta<-1:(sum(table(BASE_BRA_STA$muni_cut.geocod_6.j.)>0)*192)
BASE_BRA_STA2<-data.frame(BASE_STATE,ID.area.sta,ID.mun.mes.sta)
BASE_BRA_STAb<-BASE_BRA_STA2[1:(156*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),3:14]
pred<-data.frame(matrix(NA,(36*sum(table(BASE_BRA_STA2$muni_cut.geocod_6.j.)>0)),12))
names(pred)<-names(BASE_BRA_STAb)
BASE_BRA_STA2b<-rbind(BASE_BRA_STAb,pred)
BASE_BRA_STA3<-data.frame(BASE_BRA_STA2[,1:2],BASE_BRA_STA2b,POP=BASE_BRA_STA2$POP,ID.area.sta,ID.mun.mes.sta)
BASE_PRUPRE<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y)
BASE_SIMPLES<-data.frame(BASE_BRA_STA2,Y2=BASE_BRA_STA3$Y,year=2002+ceiling(BASE_BRA_STA2$i/12),month=BASE_BRA_STA2$i-12*(ceiling(BASE_BRA_STA2$i/12)-1),
                         AGRPECU=BASE_BRA_STA2$Y*BASE_BRA_STA2$AGRI_PECU,DOMT=BASE_BRA_STA2$Y*BASE_BRA_STA2$DOM,TURVIA=BASE_BRA_STA2$Y*BASE_BRA_STA2$TUR_VIAJ,
                         GARIMINE=BASE_BRA_STA2$Y*BASE_BRA_STA2$GARI_MINE,ROADHW=BASE_BRA_STA2$Y*BASE_BRA_STA2$ROAD_HUNT_WOOD)

BASEz<-NULL
for (k in 2003:2018) {
  WDS <- subset(BASE_SIMPLES,BASE_SIMPLES$year==k)
  for (j in 1:max(WDS$ID.area.sta)) {
    QWR <- subset(WDS,WDS$ID.area.sta==j)
    Y <- sum(QWR$Y)
    AGRECU <- sum(QWR$AGRPECU)/Y
    DOM <- sum(QWR$DOMT)/Y
    TUR_V <- sum(QWR$TURVIA)/Y
    GAR_MIN <- sum(QWR$GARIMINE)/Y
    OUT <- sum(QWR$ROADHW)/Y
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT)
    BASEz<-rbind(BASEz,REC)
  }
}
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Manuscrito")
load("dados_desmatamento.rdata")
Base_desmatamento<-BASE
DESMA<-NULL
for (i in 1:nrow(BASEz)) {
  rew<-subset(Base_desmatamento,as.character(Base_desmatamento$CodIbge)==as.character(BASEz$MUN_GEO2[i]) &
                Base_desmatamento$Year==BASEz$Year[i])
  recx<-data.frame(MUN=rew$CodIbge[1],Year=rew$Year[1],Desma=rew$Incremento[1])
  DESMA<-rbind(DESMA,recx)
}
BASE_IN<-data.frame(BASEz,Desma=DESMA$Desma,i=BASEz$Year-2002,antes=as.integer(BASEz$Year<2016),despues=as.integer(BASEz$Year>2015),ID.mun.ano=seq(1,nrow(BASEz)))
BASE_IN[is.na(BASE_IN)] <- 0
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/AM")

######## ########

formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AMA.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + 
  despues*Desma + despues*AGRECU + despues*DOM # + despues*OUT#+ despues*TUR_V + despues*GAR_MIN
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
#beep("mario")
## alpha fixed effects
R_18f<-model.inla.1$summary.fixed
R_18fe<-exp(model.inla.1$summary.fixed)
model.inla.1$dic$dic
save(R_18f,file="R18f.rdata")
save(R_18fe,file="R18fe.rdata")