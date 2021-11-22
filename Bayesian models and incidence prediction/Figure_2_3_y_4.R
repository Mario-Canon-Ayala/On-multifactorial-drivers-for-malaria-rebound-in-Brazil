###########################################
#               Packages                   #
############################################
library(foreign)
library(lubridate)
library(rgdal)
library(readxl)
library('stringr')
library(gridExtra)
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Maps")
#############################################
#                BRAZIL MAP                 #
#############################################
muni <- readOGR("br_municipios","BRMUE250GC_SIR")
muni@data$geocod_2 <- str_sub(muni$CD_GEOCMU, 1, 2)
muni@data$geocod_6 <- str_sub(muni$CD_GEOCMU, 1, 6)
muni_cut <- muni[{muni$geocod_2 == 11}|{muni$geocod_2 == 12}|{muni$geocod_2 == 13}|{muni$geocod_2 == 14}|{muni$geocod_2 == 15}|
{muni$geocod_2 == 16}|{muni$geocod_2 == 17}|{muni$geocod_2 == 21}|{muni$geocod_2 == 51},]
municipios<-muni_cut@data
#################################################################################
#############
##  Vivax   #
#############
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AP")
load("MAPER_AP_viv.rdata")
AP_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AM")
load("MAPER_AM_viv.rdata")
AM_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/PA")
load("MAPER_PA_viv.rdata")
PA_viv<-BASSER
names(PA_viv)<-names(AP_viv)
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AC")
load("MAPER_AC_viv.rdata")
AC_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/RR")
load("MAPER_RR_viv.rdata")
RR_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/RO")
load("MAPER_RO_viv.rdata")
RO_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/MT")
load("MAPER_MT_viv.rdata")
MT_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/TO")
load("MAPER_TO_viv.rdata")
TO_viv<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/MA")
load("MAPER_MA_viv.rdata")
MA_viv<-BASSER
BASE_viv<-rbind(AP_viv,AM_viv,PA_viv,AC_viv,RR_viv,RO_viv,MT_viv,TO_viv,MA_viv)
BAS_VO<-NULL
for (i in 1:808) {
  count<-subset(BASE_viv,BASE_viv$ID.area == i)
  BAS_VO<-rbind(BAS_VO,count)
}
################
## falciparum  #
################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AP")
load("MAPER_AP_fal.rdata")
AP_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AM")
load("MAPER_AM_fal.rdata")
AM_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/PA")
load("MAPER_PA_fal.rdata")
PA_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AC")
load("MAPER_AC_fal.rdata")
AC_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/RR")
load("MAPER_RR_fal.rdata")
RR_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/RO")
load("MAPER_RO_fal.rdata")
RO_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/MT")
load("MAPER_MT_fal.rdata")
MT_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/TO")
load("MAPER_TO_fal.rdata")
TO_fal<-BASSER
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/MA")
load("MAPER_MA_fal.rdata")
MA_fal<-BASSER
BASE_fal<-rbind(AP_fal,AM_fal,PA_fal,AC_fal,RR_fal,RO_fal,MT_fal,TO_fal,MA_fal)
BAS_FO<-NULL
for (i in 1:808) {
  count<-subset(BASE_fal,BASE_fal$ID.area == i)
  BAS_FO<-rbind(BAS_FO,count)
}
################################################### Vivax
BAS_VO$mun<-NULL
c1<-BAS_VO$ID.area
zetav.cutoff<- c(-1000,-200,-100,-50,-25,-10,-5,0,5,10,25,50,100,200,1000)
BAS_VO$ID.area=NULL
zeta.v2016=cut(BAS_VO$EP2016,breaks=zetav.cutoff,include.lowest=TRUE)
zeta.v2017=cut(BAS_VO$EP2017,breaks=zetav.cutoff,include.lowest=TRUE)
zeta.v2018=cut(BAS_VO$EP2018,breaks=zetav.cutoff,include.lowest=TRUE)
ERR_v<-data.frame(V2016=zeta.v2016,V2017=zeta.v2017,V2018=zeta.v2018)
ERR_v2<-data.frame(P._vivax_2016=zeta.v2016,P._vivax_2017=zeta.v2017,P._vivax_2018=zeta.v2018)
################################################### falciparum
BAS_FO$mun<-NULL
c2<-BAS_FO$ID.area
zetaf.cutoff<- c(-1000,-50,-20,-10,-5,0,5,10,20,50,1000)
BAS_FO$ID.area=NULL
zeta.f2016=cut(BAS_FO$EP2016,breaks=zetaf.cutoff,include.lowest=TRUE)
zeta.f2017=cut(BAS_FO$EP2017,breaks=zetaf.cutoff,include.lowest=TRUE)
zeta.f2018=cut(BAS_FO$EP2018,breaks=zetaf.cutoff,include.lowest=TRUE)
ERR_f<-data.frame(F2016=zeta.f2016,F2017=zeta.f2017,F2018=zeta.f2018)
ERR_f2<-data.frame(P._falciparum_2016=zeta.f2016,P._falciparum_2017=zeta.f2017,P._falciparum_2018=zeta.f2018)
#names(ERR_f2)<-c("P. falciparum 2016","P. falciparum 2017","P. falciparum 2018")
attr(muni_cut, "data")=data.frame(ERR_v2,ERR_f2)
###########################3
#         MAPS
############################
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction")
#
#      Vivax
#
#########################
#library(lattice)
#trellis.par.set(axis.line=list(col=NA))
labelat = c(1, 2, 3, 4, 5,6,7,8,9,10,11,12,13,14)
zetav.cutoff<- c(-1000,-200,-100,-50,-25,-10,-5,0,5,10,25,50,100,200,1000)
labeltextv = c("<-200", "-200 - -100", "-100 - -50", "-50 - -25", "-25 - -10","-10 - -5","-5 - 0","0 - 5","5 - 10","10 - 25","25 - 50","50 - 100","100 - 200",">200")
pdf(paste("Figure 2a.pdf",sep=""),width=5,height=1.5,paper='special')
mapv<-spplot(obj=muni_cut, zcol= c("P._vivax_2016","P._vivax_2017","P._vivax_2018"),layout=c(3,1), col.regions=colorRampPalette(c('blue', 'gray','red'))(14),par.settings=list(fontsize=list(text=7)),
             col = "transparent", 
             names.attr = c("2016","2017","2018") ,
             colorkey = list(labels=list(at=labelat,labels=labeltextv)),
             main= list( as.expression(bquote(italic("P. vivax"))), x=0.1,just = "left") 
              )
print(mapv)
dev.off()
############################
#
#      falciparum
#
#########################
#library(lattice)
#trellis.par.set(axis.line=list(col=NA))
labelatf = c(1, 2, 3, 4, 5,6,7,8,9,10)
zetaf.cutoff<- c(-1000,-50,-20,-10,-5,0,5,10,20,50,1000)
labeltextf = c("<-50", "-50 - -20", "-20 - -10", "-10 - -5", "-5 - 0","0 - 5","5 - 10","10 - 20","20 - 50",">50")
pdf(paste("Figure 2b.pdf",sep=""),width=5,height=1.5,paper='special')
mapf<-spplot(obj=muni_cut, zcol= c("P._falciparum_2016","P._falciparum_2017","P._falciparum_2018"),layout=c(3,1), col.regions=colorRampPalette(c('blue', 'gray','red'))(10),par.settings=list(fontsize=list(text=7)),
             col = "transparent",
             names.attr = c("2016","2017","2018") ,
             colorkey = list(labels=list(at=labelat,labels=labeltextv)),
             main= list( as.expression(bquote(italic("P. falciparum"))), x=0.1,just = "left"))
print(mapf)
dev.off()
##############################################################
#                                                            #
#                        Figure 2                            #
#                                                            #
##############################################################
library(gridExtra)
pdf(paste("Figure_2.pdf",sep=""),width=5,height=3,paper='special')
grid.arrange(mapv,mapf, nrow=2)
dev.off()

############################################################################
#                                                                          #
#                  Head maps (figure 3 and figure 4)                       #
#                                                                          #
############################################################################
#########
# Vivax #
#########
BASEviv<-data.frame(municipios,ERR_v) #808
################################################
#  (25,50] #
################################################
# 2016
Ev25_50_16<-subset(BASEviv,BASEviv$V2016=="(25,50]")
Ev50_100_16<-subset(BASEviv,BASEviv$V2016=="(50,100]")
Ev100_200_16<-subset(BASEviv,BASEviv$V2016=="(100,200]")
Ev200_mas_16<-subset(BASEviv,BASEviv$V2016=="(200,1e+03]")
BASE_TOTV1<-rbind(Ev25_50_16,Ev50_100_16,Ev100_200_16,Ev200_mas_16) #
# 2017
Ev25_50_17<-subset(BASEviv,BASEviv$V2017=="(25,50]")
Ev50_100_17<-subset(BASEviv,BASEviv$V2017=="(50,100]")
Ev100_200_17<-subset(BASEviv,BASEviv$V2017=="(100,200]")
Ev200_mas_17<-subset(BASEviv,BASEviv$V2017=="(200,1e+03]")
BASE_TOTV2<-rbind(Ev25_50_17,Ev50_100_17,Ev100_200_17,Ev200_mas_17)  #
# 2018 
Ev25_50_18<-subset(BASEviv,BASEviv$V2018=="(25,50]")
Ev50_100_18<-subset(BASEviv,BASEviv$V2018=="(50,100]")
Ev100_200_18<-subset(BASEviv,BASEviv$V2018=="(100,200]")
Ev200_mas_18<-subset(BASEviv,BASEviv$V2018=="(200,1e+03]")
BASE_TOTV3<-rbind(Ev25_50_18,Ev50_100_18,Ev100_200_18,Ev200_mas_18) # 
BASE_TOTV<-rbind(BASE_TOTV1,BASE_TOTV2,BASE_TOTV3) # 
BASE_ERRV<-NULL
for (i in 1:nrow(municipios)) {
  conj<-subset(BASE_TOTV,BASE_TOTV$NM_MUNICIP==municipios$NM_MUNICIP[i])
  red<-conj[1,]
  BASE_ERRV<-rbind(BASE_ERRV,red)
}
BASE_ERRV <- BASE_ERRV[!is.na(BASE_ERRV$NM_MUNICIP),]
BASE_ERRV <- BASE_ERRV[order(BASE_ERRV$geocod_2),]
BASEV_2016<-BASE_ERRV
BASEV_2016[,6:7]<-NULL
BASEV_2016[,2]<-NULL
BASEV_2016[,3]<-NULL
BASEV_2016<-data.frame(year=rep(2016,1,nrow(BASEV_2016)),BASEV_2016)
colnames(BASEV_2016)=c("Year","Municipality","State","Difference")
# 2017
BASEV_2017<-BASE_ERRV
BASEV_2017[,4:5]<-NULL
BASEV_2017[,2]<-NULL
BASEV_2017[,4]<-NULL
BASEV_2017<-data.frame(year=rep(2017,1,nrow(BASEV_2017)),BASEV_2017)
colnames(BASEV_2017)=c("Year","Municipality","State","Difference")
# 2018
BASEV_2018<-BASE_ERRV
BASEV_2018[,4:6]<-NULL
BASEV_2018[,2]<-NULL
BASEV_2018<-data.frame(year=rep(2018,1,nrow(BASEV_2018)),BASEV_2018)
colnames(BASEV_2018)=c("Year","Municipality","State","Difference")
#TOT para head_map
HEAD_V<-rbind(BASEV_2016,BASEV_2017,BASEV_2018)
kkk<-paste0(HEAD_V$State," ",HEAD_V$Municipality)
nam_mun<-c("(RO) ARIQUEMES",                   
           "(RO) CANDEIAS DO JAMARI",          
           "(RO) GUAJARÁ-MIRIM",         
           "(RO) MACHADINHO D'OESTE",          
           "(RO) PORTO VELHO",                 
           "(AC) CRUZEIRO DO SUL",             
           "(AC) MÂNCIO LIMA",                
           "(AM) BARCELOS",                    
           "(AM) CARAUARI",                    
           "(AM) CAREIRO DA VÁRZEA",     
           "(AM) COARI",                       
           "(AM) GUAJARÁ",               
           "(AM) ITACOATIARA",                 
           "(AM) ITAMARATI",                   
           "(AM) JAPURÁ",                
           "(AM) JUTAÍ",                 
           "(AM) MARAÃ",                      
           "(AM) NOVO AIRÃO",                 
           "(AM) PRESIDENTE FIGUEIREDO",       
           "(AM) SANTA ISABEL DO RIO NEGRO",   
           "(AM) SANTO ANTÔNIO DO IÇÁ",
           "(AM) SÃO GABRIEL DA CACHOEIRA",   
           "(AM) TAPAUÁ",                
           "(AM) TEFÉ",                       
           "(RR) ALTO ALEGRE",                 
           "(RR) BOA VISTA",                   
           "(RR) CANTÁ",                 
           "(RR) CARACARAÍ",             
           "(RR) PACARAIMA",                   
           "(RR) RORAINÓPOLIS",               
           "(RR) UIRAMUTÃ",                   
           "(PA) AFUÁ",                  
           "(PA) ANAJÁS",                
           "(PA) BAGRE",                       
           "(PA) BAIÃO",                      
           "(PA) BREVES",                      
           "(PA) CAMETÁ",                
           "(PA) CURRALINHO",                  
           "(PA) ITAITUBA",                    
           "(PA) JACAREACANGA",                
           "(PA) LIMOEIRO DO AJURU",           
           "(PA) MOCAJUBA",                    
           "(PA) NOVO PROGRESSO",              
           "(PA) NOVO REPARTIMENTO",           
           "(PA) OEIRAS DO PARÁ",        
           "(PA) OURILÁNDIA DO NORTE",        
           "(PA) PORTEL",                      
           "(PA) SÃO SEBASTIÃO DA BOA VISTA",
           "(PA) TUCUMÃ",                     
           "(PA) TUCURUÍ",               
           "(AP) MACAPÁ",                
           "(AP) MAZAGÃO",                    
           "(AP) PORTO GRANDE",                
           "(AP) SANTANA",  
           "(RO) ARIQUEMES",                   
           "(RO) CANDEIAS DO JAMARI",          
           "(RO) GUAJARÁ-MIRIM",         
           "(RO) MACHADINHO D'OESTE",          
           "(RO) PORTO VELHO",                 
           "(AC) CRUZEIRO DO SUL",             
           "(AC) MÂNCIO LIMA",                
           "(AM) BARCELOS",                    
           "(AM) CARAUARI",                    
           "(AM) CAREIRO DA VÁRZEA",     
           "(AM) COARI",                       
           "(AM) GUAJARÁ",               
           "(AM) ITACOATIARA",                 
           "(AM) ITAMARATI",                   
           "(AM) JAPURÁ",                
           "(AM) JUTAÍ",                 
           "(AM) MARAÃ",                      
           "(AM) NOVO AIRÃO",                 
           "(AM) PRESIDENTE FIGUEIREDO",       
           "(AM) SANTA ISABEL DO RIO NEGRO",   
           "(AM) SANTO ANTÔNIO DO IÇÁ",
           "(AM) SÃO GABRIEL DA CACHOEIRA",   
           "(AM) TAPAUÁ",                
           "(AM) TEFÉ",                       
           "(RR) ALTO ALEGRE",                 
           "(RR) BOA VISTA",                   
           "(RR) CANTÁ",                 
           "(RR) CARACARAÍ",             
           "(RR) PACARAIMA",                   
           "(RR) RORAINÓPOLIS",               
           "(RR) UIRAMUTÃ",                   
           "(PA) AFUÁ",                  
           "(PA) ANAJÁS",                
           "(PA) BAGRE",                       
           "(PA) BAIÃO",                      
           "(PA) BREVES",                      
           "(PA) CAMETÁ",                
           "(PA) CURRALINHO",                  
           "(PA) ITAITUBA",                    
           "(PA) JACAREACANGA",                
           "(PA) LIMOEIRO DO AJURU",           
           "(PA) MOCAJUBA",                    
           "(PA) NOVO PROGRESSO",              
           "(PA) NOVO REPARTIMENTO",           
           "(PA) OEIRAS DO PARÁ",        
           "(PA) OURILÁNDIA DO NORTE",        
           "(PA) PORTEL",                      
           "(PA) SÃO SEBASTIÃO DA BOA VISTA",
           "(PA) TUCUMÃ",                     
           "(PA) TUCURUÍ",               
           "(AP) MACAPÁ",                
           "(AP) MAZAGÃO",                    
           "(AP) PORTO GRANDE",                
           "(AP) SANTANA",
           "(RO) ARIQUEMES",                   
           "(RO) CANDEIAS DO JAMARI",          
           "(RO) GUAJARÁ-MIRIM",         
           "(RO) MACHADINHO D'OESTE",          
           "(RO) PORTO VELHO",                 
           "(AC) CRUZEIRO DO SUL",             
           "(AC) MÂNCIO LIMA",                
           "(AM) BARCELOS",                    
           "(AM) CARAUARI",                    
           "(AM) CAREIRO DA VÁRZEA",     
           "(AM) COARI",                       
           "(AM) GUAJARÁ",               
           "(AM) ITACOATIARA",                 
           "(AM) ITAMARATI",                   
           "(AM) JAPURÁ",                
           "(AM) JUTAÍ",                 
           "(AM) MARAÃ",                      
           "(AM) NOVO AIRÃO",                 
           "(AM) PRESIDENTE FIGUEIREDO",       
           "(AM) SANTA ISABEL DO RIO NEGRO",   
           "(AM) SANTO ANTÔNIO DO IÇÁ",
           "(AM) SÃO GABRIEL DA CACHOEIRA",   
           "(AM) TAPAUÁ",                
           "(AM) TEFÉ",                       
           "(RR) ALTO ALEGRE",                 
           "(RR) BOA VISTA",                   
           "(RR) CANTÁ",                 
           "(RR) CARACARAÍ",             
           "(RR) PACARAIMA",                   
           "(RR) RORAINÓPOLIS",               
           "(RR) UIRAMUTÃ",                   
           "(PA) AFUÁ",                  
           "(PA) ANAJÁS",                
           "(PA) BAGRE",                       
           "(PA) BAIÃO",                      
           "(PA) BREVES",                      
           "(PA) CAMETÁ",                
           "(PA) CURRALINHO",                  
           "(PA) ITAITUBA",                    
           "(PA) JACAREACANGA",                
           "(PA) LIMOEIRO DO AJURU",           
           "(PA) MOCAJUBA",                    
           "(PA) NOVO PROGRESSO",              
           "(PA) NOVO REPARTIMENTO",           
           "(PA) OEIRAS DO PARÁ",        
           "(PA) OURILÁNDIA DO NORTE",        
           "(PA) PORTEL",                      
           "(PA) SÃO SEBASTIÃO DA BOA VISTA",
           "(PA) TUCUMÃ",                     
           "(PA) TUCURUÍ",               
           "(AP) MACAPÁ",                
           "(AP) MAZAGÃO",                    
           "(AP) PORTO GRANDE",                
           "(AP) SANTANA")                   
HEAD_V<-data.frame(HEAD_V,STA_Municipality=kkk,nam_mun)
##############################################################
#                                                            #
#                         Figure 3                           #
#                                                            #
##############################################################
library(ggplot2) # ggplot() for plotting
pdf(paste("Figure 3.pdf",sep=""),width=8,height=9,paper='special')
 ggplot(HEAD_V,aes(x=Year,y=nam_mun,fill=Difference))+
  geom_tile()+
labs(x="",y="")+
#remove extra space
scale_y_discrete(expand=c(0,0))+
  #set a base size for all fonts
  theme_grey(base_size=12)+
  scale_fill_manual(values=colorRampPalette(c('blue','gray','red'))(12))
dev.off()
#######################################################################
#                          falciparum                                 #
#######################################################################
BASEfal<-data.frame(municipios,ERR_f)
################################################
#     #
################################################
# 2016
Ef50_100_16<-subset(BASEfal,BASEfal$F2016=="(10,20]")
Ef100_200_16<-subset(BASEfal,BASEfal$F2016=="(20,50]")
Ef200_mas_16<-subset(BASEfal,BASEfal$F2016=="(50,1e+03]")
BASE_TOTF1<-rbind(Ef50_100_16,Ef100_200_16,Ef200_mas_16)
# 2017
Ef50_100_17<-subset(BASEfal,BASEfal$F2017=="(10,20]")
Ef100_200_17<-subset(BASEfal,BASEfal$F2017=="(20,50]")
Ef200_mas_17<-subset(BASEfal,BASEfal$F2017=="(50,1e+03]")
BASE_TOTF2<-rbind(Ef50_100_17,Ef100_200_17,Ef200_mas_17)
# 2018
Ef50_100_18<-subset(BASEfal,BASEfal$F2018=="(10,20]")
Ef100_200_18<-subset(BASEfal,BASEfal$F2018=="(20,50]")
Ef200_mas_18<-subset(BASEfal,BASEfal$F2018=="(50,100]")
BASE_TOTF3<-rbind(Ef50_100_18,Ef100_200_18,Ef200_mas_18)
BASE_TOTF<-rbind(BASE_TOTF1,BASE_TOTF2,BASE_TOTF3)
BASE_ERRF<-NULL
for (i in 1:nrow(municipios)) {
  conj<-subset(BASE_TOTF,BASE_TOTF$NM_MUNICIP==municipios$NM_MUNICIP[i])
  red<-conj[1,]
  BASE_ERRF<-rbind(BASE_ERRF,red)
}
BASE_ERRF <- BASE_ERRF[!is.na(BASE_ERRF$NM_MUNICIP),]
BASE_ERRF <- BASE_ERRF[order(BASE_ERRF$geocod_2),]
save(BASE_ERRF,file="BASE_ERRF.rdata")
BASEF_2016<-BASE_ERRF
BASEF_2016[,6:7]<-NULL
BASEF_2016[,2]<-NULL
BASEF_2016[,3]<-NULL
BASEF_2016<-data.frame(year=rep(2016,1,nrow(BASEF_2016)),BASEF_2016)
colnames(BASEF_2016)=c("Year","Municipality","State","Difference")
# 2017
BASEF_2017<-BASE_ERRF
BASEF_2017[,4:5]<-NULL
BASEF_2017[,2]<-NULL
BASEF_2017[,4]<-NULL
BASEF_2017<-data.frame(year=rep(2017,1,nrow(BASEF_2017)),BASEF_2017)
colnames(BASEF_2017)=c("Year","Municipality","State","Difference")
# 2018
BASEF_2018<-BASE_ERRF
BASEF_2018[,4:6]<-NULL
BASEF_2018[,2]<-NULL
BASEF_2018<-data.frame(year=rep(2018,1,nrow(BASEF_2018)),BASEF_2018)
colnames(BASEF_2018)=c("Year","Municipality","State","Difference")
#TOT para head_map
HEAD_F<-rbind(BASEF_2016,BASEF_2017,BASEF_2018)
kkkF<-paste0(HEAD_F$State," ",HEAD_F$Municipality)
name_munf<-c("(RO) CANDEIAS DO JAMARI",       
             "(RO) PORTO VELHO",              
             "(AC) CRUZEIRO DO SUL",          
             "(AC) MÂNCIO LIMA",             
             "(AM) BARCELOS",                 
             "(AM) GUAJARÁ",            
             "(AM) HUMAITÁ",            
             "(AM) JUTAÍ",              
             "(AM) SANTA ISABEL DO RIO NEGRO",
             "(AM) SANTO ANTÔNIO DO IÇÁ",
             "(AM) SÃO GABRIEL DA CACHOEIRA",
             "(AM) TAPAUÁ",             
             "(RR) AMAJARI",                  
             "(RR) RORAINÓPOLIS",            
             "(AP) LARANJAL DO JARI",         
             "(AP) MAZAGÃO",                 
             "(AP) SANTANA", 
             "(RO) CANDEIAS DO JAMARI",       
             "(RO) PORTO VELHO",              
             "(AC) CRUZEIRO DO SUL",          
             "(AC) MÂNCIO LIMA",             
             "(AM) BARCELOS",                 
             "(AM) GUAJARÁ",            
             "(AM) HUMAITÁ",            
             "(AM) JUTAÍ",              
             "(AM) SANTA ISABEL DO RIO NEGRO",
             "(AM) SANTO ANTÔNIO DO IÇÁ",
             "(AM) SÃO GABRIEL DA CACHOEIRA",
             "(AM) TAPAUÁ",             
             "(RR) AMAJARI",                  
             "(RR) RORAINÓPOLIS",            
             "(AP) LARANJAL DO JARI",         
             "(AP) MAZAGÃO",                 
             "(AP) SANTANA",
             "(RO) CANDEIAS DO JAMARI",       
             "(RO) PORTO VELHO",              
             "(AC) CRUZEIRO DO SUL",          
             "(AC) MÂNCIO LIMA",             
             "(AM) BARCELOS",                 
             "(AM) GUAJARÁ",            
             "(AM) HUMAITÁ",            
             "(AM) JUTAÍ",              
             "(AM) SANTA ISABEL DO RIO NEGRO",
             "(AM) SANTO ANTÔNIO DO IÇÁ",
             "(AM) SÃO GABRIEL DA CACHOEIRA",
             "(AM) TAPAUÁ",             
             "(RR) AMAJARI",                  
             "(RR) RORAINÓPOLIS",            
             "(AP) LARANJAL DO JARI",         
             "(AP) MAZAGÃO",                 
             "(AP) SANTANA")                 
HEAD_F<-data.frame(HEAD_F,STA_Municipality=kkkF,name_munf)
##############################################################
#                                                            #
#                         Figure 4                           #
#                                                            #
##############################################################
library(ggplot2) # ggplot() for plotting
pdf(paste("Figure 4.pdf",sep=""),width=8,height=2.85,paper='special')
 ggplot(HEAD_F,aes(x=Year,y=name_munf,fill=Difference))+
  geom_tile()+
  labs(x="",y="")+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+
  #set a base size for all fonts
  theme_grey(base_size=12)+
  scale_fill_manual(values=colorRampPalette(c('blue','gray','red'))(6))
dev.off()
