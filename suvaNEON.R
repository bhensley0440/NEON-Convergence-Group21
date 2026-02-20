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
#   Bobby Hensley (02/19/2026)
#     Updated to include option for Fe correction
######################################################################################################################## 
library(neonUtilities)
library(plyr)
library(tidyverse)
library(broom)
library(ggplot2)
######################################################################################################################## 

#siteID<-c("SUGG","MAYF","TOMB","FLNT")
siteID<-"all"

#' Pulls surface water chemistry data from NEON data portal and loads tables into R environment
rawData<-neonUtilities::loadByProduct(dpID="DP1.20093.001", site=siteID, startdate="2023-10",enddate="2025-09", 
                                  package="expanded", include.provisional=T, check.size = F,
                                  token="eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJoZW5zbGV5QGJhdHRlbGxlZWNvbG9neS5vcmciLCJzY29wZSI6InJhdGU6dW5saW1pdGVkIHJlYWQ6cmVsZWFzZXMgcmVhZDpyZWxlYXNlcy1sYXRlc3QiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTgwMDkyMzAxMiwiaWF0IjoxNjQzMjQzMDEyLCJlbWFpbCI6ImhlbnNsZXlAYmF0dGVsbGVlY29sb2d5Lm9yZyJ9.YmAq_qJdZFiHLZRNrpA0muT3w70pwwukWPOUPJ73_bogfKps-1JSuxnVyp2wTYQuvDAFga7ltDPJj_SUQ9a1Iw")

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
fullSpecClean<-na.omit(fullSpecData) #' Omits samples with missing values to prevent errors
write.csv(fullSpecClean,file="neon_suva_full.csv")

#' Creates table metrics of SUVA254 and SUVA350
abs254<-fullSpecClean[(fullSpecClean$wavelength=="254"),]
names(abs254)[names(abs254) == "SUVA"] <- "SUVA254"
abs254[,c("wavelength","absorbance")]<-NULL
abs350<-fullSpecClean[(fullSpecClean$wavelength=="350"),]
names(abs350)[names(abs350) == "SUVA"] <- "SUVA350"
abs350<-abs350[,c("sampleID","SUVA350")]
metricsData<-merge(abs254,abs350,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)

#' Calculates E2:E3 absorbance ratio 
abs250<-fullSpecClean[(fullSpecClean$wavelength=="250"),]
abs250<-abs250[,c("sampleID","absorbance")]
names(abs250)[names(abs250) == "absorbance"] <- "abs250"
abs365<-fullSpecClean[(fullSpecClean$wavelength=="364"|fullSpecClean$wavelength=="366"),]
abs365<-plyr::ddply(abs365,c("sampleID"),summarise,abs365=mean(absorbance)) 
E2E3<-merge(abs250,abs365,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)
E2E3$E2E3<-E2E3$abs250/E2E3$abs365
E2E3<-E2E3[,c("sampleID","E2E3")]
metricsData<-merge(metricsData,E2E3,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)

#' Calculates spectral slope ratio
  #' Calculates slope from 275-295 nm
  slopes275 <- fullSpecClean %>% 
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
  slopes350 <- fullSpecClean %>% 
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
  mutate(SR = slope275/slope350)
  #' Merge with SUVA254 and SUVA350 table
  SlopeRatios<-SlopeRatios[,c("sampleID","SR")]
  metricsData<-merge(metricsData,SlopeRatios,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)

#' Writes out csv file of SUVA and SR metrics
write.csv(metricsData,file="neon_suva_metrics.csv")

#' Create table of site statistics
statsData<-plyr::ddply(metricsData,c("siteID"),summarise,nSamples=n(),meanDOC=mean(DOC),sdDOC=sd(DOC),meanSUVA254=mean(SUVA254),sdSUVA254=sd(SUVA254),
                       meanSUVA350=mean(SUVA350),sdSUVA350=sd(SUVA350),meanSR=mean(SR),sdSR=sd(SR))
  
#' Writes out csv file of site statistics
write.csv(statsData,file="neon_suva_stats.csv")

####### Option to correct for Fe interference #######
Fe<-externalLabDataByAnalyte[(externalLabDataByAnalyte$analyte=="Fe"),]
Fe<-Fe[,c("sampleID","analyteConcentration")]
colnames(Fe)<-c("sampleID","Fe")
correctedAbs<-merge(fullSpecClean,Fe,by.x="sampleID",by.y="sampleID",all.x=T,all.y=F)
correctedAbs<-na.omit(correctedAbs) #- Removes any samples with missing Fe

#' Keep and correct only wavelengths used in metrics
#' SUVA254
correctedAbs254<-correctedAbs[(correctedAbs$wavelength=="254"),]
correctedAbs254$correctedAbsorbance<-correctedAbs254$absorbance-(0.0653*correctedAbs254$Fe)
correctedAbs254$correctedSUVA254<-correctedAbs254$correctedAbsorbance/correctedAbs254$DOC*100
#' SUVA 350
correctedAbs350<-correctedAbs[(correctedAbs$wavelength=="350"),]
correctedAbs350$correctedAbsorbance<-correctedAbs350$absorbance-(0.0323*correctedAbs350$Fe)
correctedAbs350$correctedSUVA350<-correctedAbs350$correctedAbsorbance/correctedAbs350$DOC*100
#' SR
correctedAbs280<-correctedAbs[(correctedAbs$wavelength>=275),]
correctedAbs280<-correctedAbs280[(correctedAbs280$wavelength<=295),]
correctedAbs280$correctedAbsorbance<-correctedAbs280$absorbance-(0.0570*correctedAbs280$Fe)
correctedAbs400<-correctedAbs[(correctedAbs$wavelength>=350),]
correctedAbs400<-correctedAbs400[(correctedAbs400$wavelength<=400),]
correctedAbs400$correctedAbsorbance<-correctedAbs400$absorbance-(0.0118*correctedAbs400$Fe)
correctedSR<-rbind(correctedAbs280,correctedAbs400)
correctedSR$correctedSUVA<-correctedSR$correctedAbsorbance/correctedSR$DOC*100
correctedSlopes275 <- correctedSR %>% 
  filter(wavelength >= 275 & wavelength <= 295) %>% 
  group_by(sampleID) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(correctedSUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(sampleID, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::rename("intercept275" = `(Intercept)`,"slope275" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope275 = slope275*-1)
#' Calculates slope from 350-400 nm
correctedSlopes350 <- correctedSR %>% 
  filter(wavelength >= 350 & wavelength <= 400) %>% 
  group_by(sampleID) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(correctedSUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(sampleID, term, estimate) %>% 
  pivot_wider(names_from = term,values_from = estimate) %>% 
  dplyr::rename("intercept350" = `(Intercept)`,"slope350" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope350 = slope350*-1)
#' Calculates spectral slope ratio
correctedSlopeRatios <- merge(x = correctedSlopes275, y = correctedSlopes350,by = "sampleID") %>%
  mutate(correctedSR = slope275/slope350)

#' Merges corrected data table
correctedAbs254<-correctedAbs254[,c("sampleID","Fe","correctedSUVA254")]
correctedAbs350<-correctedAbs350[,c("sampleID","correctedSUVA350")]
correctedSlopeRatios<-correctedSlopeRatios[,c("sampleID","correctedSR")]
correctedMetricsData<-merge(correctedAbs254,correctedAbs350,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)
correctedMetricsData<-merge(correctedMetricsData,correctedSlopeRatios,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)

#' Option to merge with uncorrected metrics table
correctedMetricsData<-merge(metricsData,correctedMetricsData,by.x="sampleID",by.y="sampleID",all.x=T,all.y=T)
names(correctedMetricsData)[names(correctedMetricsData) == "SUVA254"] <- "rawSUVA254"
names(correctedMetricsData)[names(correctedMetricsData) == "SUVA350"] <- "rawSUVA350"
names(correctedMetricsData)[names(correctedMetricsData) == "SR"] <- "rawSR"

#' Writes out csv file of corrected SUVA and SR metrics
write.csv(correctedMetricsData,file="neon_suva_metrics_corrected.csv")
