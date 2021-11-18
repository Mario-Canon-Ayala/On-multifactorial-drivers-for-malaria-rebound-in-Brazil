############################################
#               Packages                   #
############################################
library(foreign)
library(lubridate)
###########################################
#               Data                      #
###########################################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal")
# Raw data: 5,972,715
load(file="data_2003_2018.Rdata")
# Delete the posterior report after treatment (LVC=1). We obtain new cases only
mala_total<- subset(mala_total_1,mala_total_1$ID_LVC.==2) #4,931,318 obs 
# Extract study variables
mala_total_2<- data.frame(DT_NOTIF=mala_total$DT_NOTIF,MUN_NOTI=mala_total$MUN_NOTI,DT_SINTO=mala_total$DT_SINTO,
                          DT_TRATA=mala_total$DT_TRATA,COD_OCUP=mala_total$COD_OCUP,QTD_CRUZ=mala_total$QTD_CRUZ,
                          ID_PACIE=mala_total$ID_PACIE,MUN_INFE=mala_total$MUN_INFE,PAIS_INF=mala_total$PAIS_INF,
                          RES_EXAM=mala_total$RES_EXAM,GENDER=mala_total$RACA) #4,931,318 obs
# Delete records without complete data
m_total<-na.omit(mala_total_2) # 3,973,439
# Delete records offside year
BASE<-NULL
cdall<-NULL
for (i in 3:18){
  mala_ano_tot<-subset(m_total, mala_total_2$DT_NOTIF> paste0(2000+(i-1),"-12-31") & mala_total_2$DT_NOTIF<paste0(2000+(i+1),"-01-01"))
  mala_ano_cle<-subset(m_total, m_total$DT_NOTIF> paste0(2000+(i-1),"-12-31") & m_total$DT_NOTIF<paste0(2000+(i+1),"-01-01")&
                 m_total$DT_SINTO>paste0(2000+(i-1),"-11-30") & m_total$DT_SINTO<paste0(2000+(i+1),"-01-01")&
                 m_total$DT_TRATA>paste0(2000+(i-1),"-12-31") & m_total$DT_TRATA<paste0(2000+(i+1),"-01-01"))
 BASE<-rbind(BASE,mala_ano_cle)
 cd<-dim(mala_ano_cle)[1]/dim(mala_ano_tot)[1]
 cda<-cbind(paste0(2000+i),cd)
 cdall<-rbind(cdall,cda)
}#3,963,255 valid obs
# Delete invalid data: treatment date before symptomp date
t_delay<-as.numeric(BASE$DT_TRATA-BASE$DT_SINTO)
BAS<-data.frame(BASE,t_delay)
BASE2<-subset(BAS,BAS$t_delay>-1) # 3,963,003

#################################################################
#                                                               #
#           Prepare data-base (prepare variables)               #
#                                                               #
#################################################################
# Obtain occupations code
Agri<-BASE2$COD_OCUP == 1 ### Farmer
Pecu<-BASE2$COD_OCUP == 2 ### Cattle raising
Dom<-BASE2$COD_OCUP == 3 #### Housing
Tur<-BASE2$COD_OCUP == 4 #### tourism
Gari<-BASE2$COD_OCUP == 5 ##  Gold Miner
Veg<-BASE2$COD_OCUP == 6 ###  Vegetable explotation
Cac<-BASE2$COD_OCUP == 7 ###  Hunting-fishing
Cons<- BASE2$COD_OCUP == 8 #  Highway building
Mine<- BASE2$COD_OCUP == 9 #  Miner
Viaj<- BASE2$COD_OCUP == 10#  Traveler
Otr<- BASE2$COD_OCUP == 11  # Others
###### Data ocupation ######
OCUP<-data.frame(Agri,Pecu,Dom,Tur,Gari,Veg,Cac,Cons,Mine,Viaj,Otr)
####### Convert TRUE/FALSE to 1/0 
cols <- sapply(OCUP, is.logical)
OCUP[,cols] <- lapply(OCUP[,cols], as.numeric)
#### Imported cases                  
IMP<-as.data.frame(as.integer(as.character(BASE2$MUN_NOTI) != as.character(BASE2$MUN_INFE)))
colnames(IMP)<-c('IMP_CASE')
#   Gender                              
male<-sapply(BASE2$GENDER == "M",as.numeric)
fema<-sapply(BASE2$GENDER == "F",as.numeric)
GEN<-data.frame(MALE=male,FEM=fema)
#######################################################################
#                                                                     #
#                Prepare data-base (prepare date format)              #
#                                                                     #
#######################################################################
test<-male+fema
BASEF<-cbind(BASE2,OCUP,IMP,GEN,test)
week<-as.numeric(week(ymd(BASEF$DT_NOTIF)))
week1<- data.frame(week)
week2<-as.integer(week1>52) #### data from the week 53 is adding to the week 52
week3<-week1-week2          ####  
colnames(week3)<- c('week')
year<- as.numeric(format(as.Date(BASEF$DT_NOTIF, format="%Y/%m/%d"),"%Y"))
year1<- data.frame(year)
colnames(year1)<- c('year')
cod_week<- week3-(2003-year1)*52
colnames(cod_week)<- c('COD_WEEK')
BASEF2<- cbind(BASEF,week3,year1,cod_week)
BASEF2<-BASEF2[order(BASEF2$DT_NOTIF),] # 3,963,003
BASEF3<-subset(BASEF2, BASEF2$test>0) # 3,963,003
save(BASEF3,file="BASE_ABA_042019.rdata") ###Data for applied bayesian analysis ## 042019
Data_vivax <- subset(BASEF3,BASEF3$RES_EXAM == 4) # 3,232,766 obs
Data_falciparum <- subset(BASEF3,BASEF3$RES_EXAM == 2 | BASEF3$RES_EXAM == 3 | BASEF3$RES_EXAM == 7) # 685,317 obs




