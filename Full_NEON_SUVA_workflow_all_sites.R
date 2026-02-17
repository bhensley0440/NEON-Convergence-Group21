######################################################################################################################## 
#' @title suvaNEON

#' @author Bobby Hensley \email{hensley@battelleecology.org} \cr 

#' @description This script calculates Specific Ultra-Violet Absorbance (SUVA)
#' from NEON surface water chemistry data.

#' @return This script produces a .csv file 

#' @references 
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

# changelog and author contributions / copyrights
#   Bobby Hensley (10/21/2025)
#     Original script created
#   Max Glines
#     Initial spectral slope ratio calculations and time-series generation across sites
#   Alexi Besser
#     Updated spectral slope ratios and time-series across sites
######################################################################################################################## 

# Load required packages

#install.packages("neonUtilities")
library(neonUtilities)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("plyr")
library(plyr)

######################################################################################################################## 

### SURFACE WATER ###

#' Pulls surface water chemistry data from NEON data portal and loads tables into R environment
swc <- neonUtilities::loadByProduct(dpID = "DP1.20093.001",
                                    site = "all",
                                    startdate = "2023-10",
                                    enddate = "2025-12",
                                    package = "expanded",
                                    include.provisional = TRUE,
                                    check.size = FALSE)
list2env(swc,.GlobalEnv)

#' Averages decadicAbsorbance from "A" and "B" spectrum replicates
swc_externalLabAbsorbanceScan$sampleID.wavelength <- paste(swc_externalLabAbsorbanceScan$sampleID,
                                                           swc_externalLabAbsorbanceScan$wavelength,
                                                           sep =".")

swc_fullUV <- ddply(swc_externalLabAbsorbanceScan, c("sampleID.wavelength"), summarize, 
                    sampleID = unique(sampleID), domainID = unique(domainID), siteID = unique(siteID), collectDate = unique(collectDate),
                    wavelength = unique(wavelength), absorbance = mean(decadicAbsorbance)) 

#' Combines DOC and Abs254 values into wide format using sampleIDs
swc_DOC <- swc_externalLabDataByAnalyte[(swc_externalLabDataByAnalyte$analyte == "DOC"), ]
swc_DOC <- swc_DOC[, c("sampleID", "analyteConcentration")]
colnames(swc_DOC) <- c("sampleID", "DOC")
swc_SUVA <- merge(swc_fullUV, swc_DOC, by.x = "sampleID", by.y = "sampleID")

#' Calculates specific ultra-violet absorbance (SUVA) in units of L/mg-m
swc_SUVA$SUVA <- swc_SUVA$absorbance/swc_SUVA$DOC*100 

#' Writes out csv file of results
write.csv(swc_SUVA, file = "Data/NEON_SurfaceWater_SUVA_Data.csv")

#' Plots surface water SUVA for all sites and dates
sw_suva_plot <- ggplot(swc_SUVA, aes(x = wavelength, y = SUVA, color = sampleID)) +
  geom_line() +
  scale_x_continuous(limits = c(200, 400)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme(legend.position = "none")

sw_suva_plot

### GROUND WATER ###

#' Pulls groundwater chemistry data from NEON data portal and loads tables into R environment
gwc <- neonUtilities::loadByProduct(dpID = "DP1.20092.001",
                                    site = "all",
                                    startdate = "2023-10",
                                    enddate = "2025-12",
                                    package = "expanded",
                                    include.provisional = TRUE,
                                    check.size = FALSE)
list2env(gwc,.GlobalEnv)

#' Averages decadicAbsorbance from "A" and "B" spectrum replicates
gwc_externalLabAbsorbanceScan$sampleID.wavelength <- paste(gwc_externalLabAbsorbanceScan$sampleID,
                                                           gwc_externalLabAbsorbanceScan$wavelength,
                                                           sep =".")

gwc_fullUV <- ddply(gwc_externalLabAbsorbanceScan, c("sampleID.wavelength"), summarize, 
                    sampleID = unique(sampleID), domainID = unique(domainID), siteID = unique(siteID), collectDate = unique(collectDate),
                    wavelength = unique(wavelength), absorbance = mean(decadicAbsorbance)) 

#' Combines DOC and Abs254 values into wide format using sampleIDs
gwc_DOC <- gwc_externalLabDataByAnalyte[(gwc_externalLabDataByAnalyte$analyte == "DOC"), ]
gwc_DOC <- gwc_DOC[, c("sampleID", "analyteConcentration")]
colnames(gwc_DOC) <- c("sampleID", "DOC")
gwc_SUVA <- merge(gwc_fullUV, gwc_DOC, by.x = "sampleID", by.y = "sampleID")

#' Calculates specific ultra-violet absorbance (SUVA) in units of L/mg-m
gwc_SUVA$SUVA <- gwc_SUVA$absorbance/gwc_SUVA$DOC*100 

#' Writes out csv file of results
write.csv(gwc_SUVA, file = "Data/NEON_GroundWater_SUVA_Data.csv")

#' Plots groundwater SUVA for all sites and dates
gw_suva_plot <- ggplot(gwc_SUVA, aes(x = wavelength, y = SUVA, color = sampleID)) +
  geom_line() +
  scale_x_continuous(limits = c(200, 400)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme(legend.position = "none")

gw_suva_plot

### SURFACE WATER PLOTS ###

#' Figure 1
#' ### Need to add code in to color streams and lakes differently ###
#' Fixed axes
Fig_1_sw_suva_fixed <- ggplot(swc_SUVA, aes(x = wavelength, y = SUVA, group = collectDate))+
  geom_line() +
  facet_wrap(~siteID, scales = "fixed") +
  theme_linedraw() +
  xlab("Wavelength (nm)") +
  ylab(expression(SUVA~(L/(mg %.% m))))

Fig_1_sw_suva_fixed

ggsave("Figures/Fig_1_SurfaceWater_SUVA_fixed.pdf", plot = last_plot())

#' Free axes
Fig_1_sw_suva_free <- ggplot(swc_SUVA, aes(x = wavelength, y = SUVA, group = collectDate))+
  geom_line() +
  facet_wrap(~siteID, scales = "free") +
  theme_linedraw() +
  xlab("Wavelength (nm)") +
  ylab(expression(SUVA~(L/(mg %.% m))))

Fig_1_sw_suva_free

ggsave("Figures/Fig_1_SurfaceWater_SUVA_free.pdf", plot = last_plot())

### GROUNDWATER PLOTS ###

#' Figure 2
#' ### Need to add code in to color streams and lakes differently ###
#' Fixed axes
Fig_2_gw_suva_fixed <- ggplot(gwc_SUVA, aes(x = wavelength, y = SUVA, group = collectDate))+
  geom_line() +
  facet_wrap(~siteID, scales = "fixed") +
  theme_linedraw() +
  xlab("Wavelength (nm)") +
  ylab(expression(SUVA~(L/(mg %.% m))))

Fig_2_gw_suva_fixed

ggsave("Figures/Fig_2_GroundWater_SUVA_fixed.pdf", plot = last_plot())

#' Free axes
Fig_2_gw_suva_free <- ggplot(gwc_SUVA, aes(x = wavelength, y = SUVA, group = collectDate))+
  geom_line() +
  facet_wrap(~siteID, scales = "free") +
  theme_linedraw() +
  xlab("Wavelength (nm)") +
  ylab(expression(SUVA~(L/(mg %.% m))))

Fig_2_gw_suva_free

ggsave("Figures/Fig_2_GroundWater_SUVA_free.pdf", plot = last_plot())

#spectral slope--------------------
library(broom)

ggplot(suva,aes(wavelength,log(SUVA)))+
  geom_line()+
  facet_wrap(~siteID)

#remove any observations where SUVA is NA at any wavelength
suva_wide <- suva %>% 
  select(-DOC,-absorbance) %>% 
  group_by(siteID,date,wavelength) %>% 
  summarise(SUVA = mean(SUVA,na.rm=TRUE)) %>% 
  pivot_wider(names_from=wavelength,values_from=SUVA) %>% 
  ungroup()

suva_filter <- suva_wide[complete.cases(suva_wide),] %>% 
  pivot_longer(cols=3:ncol(suva_wide),names_to='wavelength',values_to='SUVA') %>% 
  mutate(wavelength=as.numeric(wavelength))

#calculate slope from 275-295 nm
slopes_275 <- suva_filter %>% 
  filter(wavelength >= 275 & wavelength <= 295) %>% 
  group_by(siteID,date) %>% 
  nest() %>% 
  mutate(S = map(data,~lm(log(SUVA)~wavelength,data=.x))) %>% 
  mutate(slope = map(S,~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(siteID,date,term,estimate) %>% 
  pivot_wider(names_from=term,values_from=estimate) %>% 
  rename('intercept350'=`(Intercept)`,
         'slope275'=wavelength) %>% 
  ungroup() %>% 
  mutate(slope275=slope275*-1)

#calculate slope from 350-400 nm
slopes_350 <- suva_filter %>% 
  filter(wavelength >= 350 & wavelength <= 400) %>% 
  group_by(siteID,date) %>% 
  nest() %>% 
  mutate(S = map(data,~lm(log(SUVA)~wavelength,data=.x))) %>% 
  mutate(slope = map(S,~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(siteID,date,term,estimate) %>% 
  pivot_wider(names_from=term,values_from=estimate) %>% 
  rename('intercept350'=`(Intercept)`,
         'slope350'=wavelength) %>% 
  ungroup() %>% 
  mutate(slope350=slope350*-1)

#calculate slope ratio
suva_slopes <- merge(x=slopes_275,y=slopes_350,
                     by=c('siteID','date')) %>% 
  mutate(Sr = slope275/slope350)

#visualize data
ggplot(suva_slopes,aes(Sr,fill=siteID))+
  geom_density(alpha=0.2)

ggplot(suva_slopes,aes(siteID,Sr))+
  geom_boxplot()

ggplot(suva_slopes,aes(date,Sr))+
  geom_point()+
  geom_line()+
  facet_wrap(~siteID,scales='free')

variation <- suva_slopes %>% 
  group_by(siteID) %>% 
  summarise(var = var(Sr),
            mean = mean(Sr))

overall_var = var(suva_slopes$Sr)
overall_mean = mean(suva_slopes$Sr)

var_driver <- merge(x=variation,y=aquatic,
                    by.x='siteID',by.y='site_id')

ggplot(var_driver,aes(mean,var,color=site_subtype))+
  geom_text(aes(label=siteID))+
  geom_vline(xintercept=overall_mean,color='red',linetype='dashed')+
  geom_hline(yintercept=overall_var,color='red',linetype='dashed')

ggplot(var_driver,aes(mean,var))+
  geom_text(aes(label=siteID,color=site_subtype))+
  geom_vline(xintercept=overall_mean,color='red',linetype='dashed')+
  geom_hline(yintercept=overall_var,color='red',linetype='dashed')+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method='lm')

driver_plots <- var_driver %>% 
  select(siteID,site_subtype,var,mean,latitude,longitude,mean_evelation_m,
         mean_annual_precipitation_mm,mean_annual_temperature_C,
         dominant_nlcd_classes,watershed_size_km2,
         lake_depth_mean_m,lake_depth_max_m) %>% 
  mutate(mean_annual_temperature_C = parse_number(mean_annual_temperature_C),
         mean_evelation_m = as.numeric(mean_evelation_m),
         mean_annual_precipitation_mm = as.numeric(mean_annual_precipitation_mm),
         watershed_size_km2 = log(watershed_size_km2)) %>% 
  pivot_longer(cols=c(5:9,11:13),values_to = 'value',names_to = 'variable')

ggplot(driver_plots,aes(value,mean))+
  geom_point(aes(color=site_subtype))+
  geom_smooth(method='lm',color='black')+
  # geom_errorbar(aes(ymin=mean-var,ymax=mean+var),width=0.2)+
  facet_wrap(~variable,scales='free')

ggplot(driver_plots,aes(value,log(var)))+
  geom_point(aes(color=site_subtype))+
  geom_smooth(method='lm',color='black')+
  facet_wrap(~variable,scales='free')
