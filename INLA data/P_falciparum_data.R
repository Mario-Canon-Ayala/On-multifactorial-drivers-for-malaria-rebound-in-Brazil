############################################
#               Packages                   #
############################################
library(foreign)
library(lubridate)
library(rgdal)
library(readxl)
library('stringr')
library(readODS)
setwd("~/Incidence data")
# Load depurate data
load("BASE_ABA_042019.rdata", verbose=TRUE)
# Calculate data-base at month
mes<-data.frame(month=12*(year(BASEF3$DT_NOTIF)-2002-1)+month(BASEF3$DT_NOTIF))
BASE<-cbind(BASEF3,mes)
# Falciparum data-base
BASE_fal<-subset(BASE,BASE$RES_EXAM == 2 | BASE$RES_EXAM == 3 | BASE$RES_EXAM == 7) #
# Brazil municipalities
setwd("~/Maps")
# Obtain Amazon basin
muni <- readOGR("br_municipios","BRMUE250GC_SIR")
muni@data$geocod_2 <- str_sub(muni$CD_GEOCMU, 1, 2)
muni@data$geocod_6 <- str_sub(muni$CD_GEOCMU, 1, 6)
muni_cut <- muni[{muni$geocod_2 == 11}|{muni$geocod_2 == 12}|{muni$geocod_2 == 13}|{muni$geocod_2 == 14}|{muni$geocod_2 == 15}|
{muni$geocod_2 == 16}|{muni$geocod_2 == 17}|{muni$geocod_2 == 21}|{muni$geocod_2 == 51},] 
###### Base to INLA #####
plot(muni_cut)
BASEP<-NULL
for (i in 1:192){
  set1<- subset(BASE_fal,BASE_fal$month==i)
  for (j in 1:length(muni_cut$geocod_6)) {
    set2<- subset(set1,set1$MUN_NOTI==muni_cut$geocod_6[j])
    Y<-nrow(set2)##Number of cases per week
    EDAD<-mean(set2$ID_PACIE)
    DELAY<-mean(set2$t_delay)
    PARA<-mean(as.numeric(set2$QTD_CRUZ))
    AGRI_PECU<-(sum(set2$Agri)+sum(set2$Pecu))/Y
    DOM<-sum(set2$Dom)/Y
    TUR_VIAJ<-(sum(set2$Tur)+sum(set2$Viaj))/Y
    GARI_MINE<-(sum(set2$Gari)+sum(set2$Mine))/Y
    ROAD_HUNT_WOOD<-(sum(set2$Veg)+sum(set2$Cac)+sum(set2$Cons))/Y
    IMPO<-sum(set2$IMP_CASE)/Y
    MALE<-sum(set2$MALE)/Y
    BA<-data.frame(muni_cut$geocod_6[j],muni_cut$CD_GEOCMU[j],i,Y,EDAD,DELAY,PARA,AGRI_PECU,DOM,TUR_VIAJ,GARI_MINE,ROAD_HUNT_WOOD,IMPO,MALE)
    BASEP<-rbind(BASEP,BA)
  }}
# save results
#setwd("~/INLA data")
#save(BASEP,file="BRASIL_INLA_fal_mes_2.rdata")
#load("BRASIL_INLA_fal_mes_2.rdata", verbose=TRUE)
#       ID   MUN (area-population)         #
ID.area<-NULL
for (i in 1:192) {
  area<-data.frame(1:808)
  colnames(area)<-c("ID.area")
  ID.area<-rbind(ID.area,area)
}
ID.mun.mes<-data.frame(1:nrow(BASEP))
colnames(ID.mun.mes)<-c("ID.mun.mes")
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Population data")
POP <- read_excel("Populations.xls")
order <- match(muni_cut$CD_GEOCMU,POP$COD)
POP_order2<- NULL
for (i in 1:nrow(muni_cut)){
  faj<-data.frame(POP[order[i],])
  POP_order2<-rbind(POP_order2,faj)
}
POP_order2$COD<-NULL
####By twelve months per year####
Bass<-NULL
for(i in 1:12){
  elf<-data.frame(POP=POP_order2)
  Bass<-rbind(Bass,elf)
}
month_POP<-NULL
for(i in 1:16){
  HHH<-data.frame(POP=Bass[i])
  colnames(HHH)<-c("POP")
  month_POP<-rbind(month_POP,HHH)
}
BASE_BRA<-data.frame(BASEP,ID.area,ID.mun.mes,month_POP)
# save results
#setwd("~/INLA data")
#save(BASE_BRA,file="BRASIL_INLA_mes_fal.rdata")
