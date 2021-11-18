# C?digo para hacer la corrida con covariables hasta 2015
### Packages
library(beepr)
library(foreign)
library(lubridate)
library(rgdal)
library('stringr')
library(maptools)
library(sp)
library(spdep)
library(INLA)
INLA:::inla.dynload.workaround() #

##################################################################################################################
#                                                                                                                #
#                                                 Vivax                                                          #
#                                                                                                                #
##################################################################################################################

setwd("~/INLA data")
load("BRASIL_INLA_mes_v.rdata", verbose=TRUE) #vivax data
#### Amapa 1600000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<170000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>159999)

setwd("~/Bayesian models and incidence prediction/AP")
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
    IMPO <- sum(QWR$IMPO)/Y 
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT,IMPO)
    BASEz<-rbind(BASEz,REC)
  }
}

setwd("~/Deforestation data")
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
setwd("~/Bayesian models and incidence prediction/AP")

######## ########
formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AP.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + despues*AGRECU + despues*TUR_V + despues*Desma  + despues*DOM  + despues*GAR_MIN + despues*OUT +despues*IMPO
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
beep("mario")
## alpha fixed effects
model.inla.1$summary.fixed
(exp(model.inla.1$summary.fixed*0.01) - 1)*100 
model.inla.1$dic$dic

## Amapá vivax results
model.inla.1$summary.fixed
model.inla.1$dic$dic
# incidence increase from 1% of variable increase
AP_viv <- model.inla.1$summary.fixed[10:16,] 
#save(AP_viv,file="AP_viv_variavels.Rdata")


##################################333
#
#            Falciparum
#########################3
setwd("~/INLA data")
load("BRASIL_INLA_mes_fal.rdata", verbose=TRUE) #falciparum data
#### Amapa 1600000
BASE_BRA_STA<-subset(BASE_BRA,as.character(BASE_BRA$muni_cut.geocod_6.j.)<170000 & as.character(BASE_BRA$muni_cut.geocod_6.j.)>159999)

setwd("~/Bayesian models and incidence prediction/AP")
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
    IMPO <- sum(QWR$IMPO)/Y 
    REC <- data.frame(MUN_GEO=QWR$muni_cut.geocod_6.j.[1],MUN_GEO2=QWR$muni_cut.CD_GEOCMU.j.[1],
                      ID.area.sta=j,Year=QWR$year[1],Y,POP=QWR$POP[1],AGRECU,DOM,TUR_V,GAR_MIN,OUT,IMPO)
    BASEz<-rbind(BASEz,REC)
  }
}

setwd("~/Deforestation data")
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
setwd("~/Bayesian models and incidence prediction/AP")

######## ########
formula.1<- Y ~ 1 + f(ID.area.sta,model="bym",graph=AP.adj) + f(i,model="rw1") + f(ID.mun.ano, model="iid") + despues*AGRECU + despues*TUR_V + despues*Desma  + despues*DOM  + despues*GAR_MIN + despues*OUT +despues*IMPO
######## Run code to make INLA to 2018

t <- proc.time()
model.inla.1 <- inla(formula.1,family="poisson",data=BASE_IN,E=POP/100000, control.predictor=list(compute=TRUE), control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE))
proc.time()-t
beep("mario")
## alpha fixed effects
model.inla.1$summary.fixed
(exp(model.inla.1$summary.fixed*0.01) - 1)*100 
model.inla.1$dic$dic

## falciparum results
model.inla.1$summary.fixed
model.inla.1$dic$dic
# incidence increase from 1% of variable increase
AP_fal <- model.inla.1$summary.fixed[10:16,] 
#save(AP_fal,file="AP_fal_variavels.Rdata")
