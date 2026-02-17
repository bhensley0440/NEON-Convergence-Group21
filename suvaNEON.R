######################################################################################################################## 
#' @title suvaNEON

#' @author Bobby Hensley \email{hensley@battelleecology.org} \cr 

#' @description This script calculates Specific Ultra-Violet Absorbance (SUVA)
#' from NEON surface water or groundwater chemistry data.

#' @return This script produces a .csv file 

#' @references 
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

# changelog and author contributions / copyrights
#   Bobby Hensley (10/21/2025)
#     Original script created
#   Bobby Hensley (02/17/2026)
#     Updated to work for either surface or groundwater data
######################################################################################################################## 
library(neonUtilities)
library(plyr)
library(tidyverse)
library(broom)
library(ggplot2)
######################################################################################################################## 

siteID<-"COMO"

#' Pulls surface water chemistry data from NEON data portal and loads tables into R environment
rawData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteID, startdate="2023-10",enddate="2025-09", 
                                  package="expanded", include.provisional=T, check.size = F)

#' Assigns data table names regardless of swc or gwc prefix
externalLabAbsorbanceScan<-rawData[[grep("externalLabAbsorbanceScan",names(rawData))]]
externalLabDataByAnalyte<-rawData[[grep("externalLabDataByAnalyte",names(rawData))]]

#' Averages absorbance from "A" and "B" spectrum replicates
externalLabAbsorbanceScan$sampleID.wavelength<-paste(externalLabAbsorbanceScan$sampleID, externalLabAbsorbanceScan$wavelength, sep=".")
fullSpecData<-plyr::ddply(externalLabAbsorbanceScan,c("sampleID.wavelength"),summarise,sampleID=unique(sampleID),domainID=unique(domainID),siteID=unique(siteID),
                        collectDate=unique(collectDate),wavelength=unique(wavelength),absorbance=mean(decadicAbsorbance)) 
fullSpecData$sampleID.wavelength<-NULL

#' Combines absorbance and DOC values into wide format using sampleID's
DOC<-externalLabDataByAnalyte[(externalLabDataByAnalyte$analyte=="DOC"),]
DOC<-DOC[,c("sampleID","analyteConcentration")]
colnames(DOC)<-c("sampleID","DOC")
fullSpecData<-merge(fullSpecData, DOC,by.x="sampleID",by.y="sampleID")

#' Calculates specific ultra-violet absorbance (units - L/mg-m)
fullSpecData$SUVA<-fullSpecData$absorbance/fullSpecData$DOC*100 

#' Writes out csv file of full spectrum results
write.csv(fullSpecData,file=paste0(siteID,"_fullSpecData.csv"))

#' Creates secondary table of just SUVA254 and SUVA350
abs254<-fullSpecData[(fullSpecData$wavelength=="254"),]
names(abs254)[names(abs254) == "SUVA"] <- "SUVA254"
abs254[,c("wavelength","absorbance","DOC")]<-NULL
abs350<-fullSpecData[(fullSpecData$wavelength=="350"),]
names(abs350)[names(abs350) == "SUVA"] <- "SUVA350"
abs350<-abs350[,c("sampleID","SUVA350")]
summaryData<-merge(abs254,abs350,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)

#' Calculates spectral slope ratio
fullSpcecClean<-na.omit(fullSpecData) #' Omits samples with missing values to prevent errors
  #' Calculates slope from 275-295 nm
  slopes275 <- fullSpcecClean %>% 
  filter(wavelength >= 275 & wavelength <= 295) %>% 
  group_by(sampleID) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(SUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(sampleID, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::rename("intercept275" = `(Intercept)`,"slope275" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope275 = slope275*-1)
  #' Calculates slope from 350-400 nm
  slopes350 <- fullSpcecClean %>% 
  filter(wavelength >= 350 & wavelength <= 400) %>% 
  group_by(sampleID) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(SUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(sampleID, term, estimate) %>% 
  pivot_wider(names_from = term,values_from = estimate) %>% 
  dplyr::rename("intercept350" = `(Intercept)`,"slope350" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope350 = slope350*-1)
  #' Calculates spectral slope ratio
  SlopeRatios <- merge(x = slopes275, y = slopes350,by = "sampleID") %>%
  mutate(Sr = slope275/slope350)
  #' Merge with SUVA254 and SUVA350 table
  SlopeRatios<-SlopeRatios[,c("sampleID","Sr")]
  summaryData<-merge(summaryData,SlopeRatios,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)

#' Writes out secondary csv file of results
write.csv(summaryData,file=paste0(siteID,"_summaryData.csv"))

##############################################################################################################################################################################
plot<-ggplot(fullSpecData,aes(x=wavelength,y=SUVA, colour=sampleID))+
  geom_line()+scale_x_continuous(limits = c(200, 400))+scale_y_continuous(limits = c(0, 10))
plot
