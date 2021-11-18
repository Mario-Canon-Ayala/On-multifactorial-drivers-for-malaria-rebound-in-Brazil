########################################################################
#                                                                      #
#               Variable analysis with INLA (table 1)                  #
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

# Open results

# Acre
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AC")
load(file="AC_viv_variavels.Rdata")
load(file="AC_fal_variavels.Rdata")
ACRE <-rbind(data.frame(State=rep("AC",7),Variable=row.names(AC_viv),Species=rep("P. vivax",7),AC_viv),
             data.frame(State=rep("AC",7),Variable=row.names(AC_fal),Species=rep("P. falciparum",7),AC_fal))
# Amazonas
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AM")
load(file="AM_viv_variavels.Rdata")
load(file="AM_fal_variavels.Rdata")
AMA <-rbind(data.frame(State=rep("AM",7),Variable=row.names(AM_viv),Species=rep("P. vivax",7),AM_viv),
             data.frame(State=rep("AM",7),Variable=row.names(AM_fal),Species=rep("P. falciparum",7),AM_fal))
# Amapá
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/AP")
load(file="AP_viv_variavels.Rdata")
load(file="AP_fal_variavels.Rdata")
AMAPA <-rbind(data.frame(State=rep("AP",7),Variable=row.names(AP_viv),Species=rep("P. vivax",7),AP_viv),
             data.frame(State=rep("AP",7),Variable=row.names(AP_fal),Species=rep("P. falciparum",7),AP_fal))
# Pará
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/PA")
load(file="PA_viv_variavels.Rdata")
load(file="PA_fal_variavels.Rdata")
PARA <-rbind(data.frame(State=rep("PA",7),Variable=row.names(PA_viv),Species=rep("P. vivax",7),PA_viv),
             data.frame(State=rep("PA",7),Variable=row.names(PA_fal),Species=rep("P. falciparum",7),PA_fal))
# Roraima
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/RR")
load(file="RR_viv_variavels.Rdata")
load(file="RR_fal_variavels.Rdata")
RORAIMA <-rbind(data.frame(State=rep("RR",7),Variable=row.names(RR_viv),Species=rep("P. vivax",7),RR_viv),
             data.frame(State=rep("RR",7),Variable=row.names(RR_fal),Species=rep("P. falciparum",7),RR_fal))
# Rondonia
setwd("C:/Users/Mario Cañon/ownCloud/malaria/Analises Mario/Dados/Artigo malaria Journal/Bayesian models and incidence prediction/RO")
load(file="RO_viv_variavels.Rdata")
load(file="RO_fal_variavels.Rdata")
RONDONIA <-rbind(data.frame(State=rep("RO",7),Variable=row.names(RO_viv),Species=rep("P. vivax",7),RO_viv),
             data.frame(State=rep("RO",7),Variable=row.names(RO_fal),Species=rep("P. falciparum",7),RO_fal))
##########################################################################
#  All data
ALL <- rbind(ACRE,AMA,AMAPA,PARA,RORAIMA,RONDONIA)
rownames(ALL) <- 1:nrow(ALL)
##########################################################################
#  Positive relation
POSITIVE <- subset(ALL, ALL$X0.025quant>0 & ALL$X0.975quant>0)
# Negative relation
NEGATIVE <- subset(ALL, ALL$X0.025quant<0 & ALL$X0.975quant<0)
# With relation
RELT <- rbind(POSITIVE,NEGATIVE)
# DEFORESTATION (1 km2 of increase)
DEF <- subset(RELT,RELT$Variable=="despues:Desma")
DEFT <- cbind(DEF[,1:3],(exp(DEF[,4:10])-1)*100)
# OCCUPATIONS (1% of increase)
OCU <- subset(RELT,RELT$Variable!="despues:Desma")
OCUT <- cbind(OCU[,1:3],(exp(OCU[,4:10]*0.01)-1)*100)
##########################################################################
#
#                         Table 1 results
#
##########################################################################
TABLE <- rbind(DEFT,OCUT)
