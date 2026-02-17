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
#   Alexi Besser (02/16/2026)
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


### SURFACE WATER SPECTRAL SLOPE RATIO CALCULATIONS ###

#' Plots log(SUVA) for all sites and dates
ggplot(swc_SUVA, aes(wavelength, log(SUVA))) +
  geom_line() +
  facet_wrap(~siteID)

### I am not sure these next two chunks are doing what is intended ###
#' Removes observations where SUVA is NA at any wavelength
swc_SUVA_wide <- swc_SUVA %>% 
  select(-sampleID.wavelength, -DOC,-absorbance) %>% 
  group_by(siteID, collectDate, wavelength) %>% 
  mutate(SUVA = mean(SUVA, na.rm = TRUE)) %>% 
  pivot_wider(names_from = wavelength, values_from = SUVA) %>% 
  ungroup()

#' Annotate something here...
swc_suva_filter <- swc_SUVA_wide[complete.cases(swc_SUVA_wide), ] 

#' I think this is a simpler way but might not take care of the complete data for each observation issue... 
swc_SUVA_cleaned <- swc_SUVA[!is.na(swc_SUVA$SUVA),]

#' Calculates slope from 275-295 nm
swc_slopes_275 <- swc_SUVA_cleaned %>% 
  filter(wavelength >= 275 & wavelength <= 295) %>% 
  group_by(siteID, collectDate) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(SUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(siteID, collectDate, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::rename("intercept275" = `(Intercept)`,
         "slope275" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope275 = slope275*-1)

#' Calculates slope from 350-400 nm
swc_slopes_350 <- swc_SUVA_cleaned %>% 
  filter(wavelength >= 350 & wavelength <= 400) %>% 
  group_by(siteID, collectDate) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(SUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(siteID, collectDate, term, estimate) %>% 
  pivot_wider(names_from = term,values_from = estimate) %>% 
  dplyr::rename("intercept350" = `(Intercept)`,
         "slope350" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope350 = slope350*-1)

#' Calculates spectral slope ratio
swc_suva_slopes <- merge(x = swc_slopes_275, y = swc_slopes_350,
                         by = c("siteID", "collectDate")) %>%
  mutate(Sr = slope275/slope350)


### GROUNDWATER SPECTRAL SLOPE RATIO CALCULATIONS ###

#' Plots log(SUVA) for all sites and dates
ggplot(gwc_SUVA, aes(wavelength, log(SUVA))) +
  geom_line() +
  facet_wrap(~siteID)

### I am not sure these next two chunks are doing what is intended ###
#' Removes observations where SUVA is NA at any wavelength
gwc_SUVA_wide <- gwc_SUVA %>% 
  select(-sampleID.wavelength, -DOC,-absorbance) %>% 
  group_by(siteID, collectDate, wavelength) %>% 
  mutate(SUVA = mean(SUVA, na.rm = TRUE)) %>% 
  pivot_wider(names_from = wavelength, values_from = SUVA) %>% 
  ungroup()

#' Annotate something here...
gwc_suva_filter <- gwc_SUVA_wide[complete.cases(gwc_SUVA_wide), ] 

#' I think this is a simpler way but might not take care of the complete data for each observation issue... 
gwc_SUVA_cleaned <- gwc_SUVA[!is.na(gwc_SUVA$SUVA),]

#' Calculates slope from 275-295 nm
gwc_slopes_275 <- gwc_SUVA_cleaned %>% 
  filter(wavelength >= 275 & wavelength <= 295) %>% 
  group_by(siteID, collectDate) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(SUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(siteID, collectDate, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  dplyr::rename("intercept275" = `(Intercept)`,
                "slope275" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope275 = slope275*-1)

#' Calculates slope from 350-400 nm
gwc_slopes_350 <- gwc_SUVA_cleaned %>% 
  filter(wavelength >= 350 & wavelength <= 400) %>% 
  group_by(siteID, collectDate) %>% 
  nest() %>% 
  mutate(S = map(data, ~lm(log(SUVA) ~ wavelength, data = .x))) %>% 
  mutate(slope = map(S, ~tidy(.x))) %>% 
  unnest(slope) %>% 
  select(siteID, collectDate, term, estimate) %>% 
  pivot_wider(names_from = term,values_from = estimate) %>% 
  dplyr::rename("intercept350" = `(Intercept)`,
                "slope350" = wavelength) %>% 
  ungroup() %>% 
  mutate(slope350 = slope350*-1)

#' Calculates spectral slope ratio
gwc_suva_slopes <- merge(x = gwc_slopes_275, y = gwc_slopes_350,
                         by = c("siteID", "collectDate")) %>%
  mutate(Sr = slope275/slope350)


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

#' Figure 3
#' Boxplots of spectral slope ratios across sites
Fig_3_sw_sr_boxplots <- ggplot(swc_suva_slopes, aes(siteID, Sr)) +
  geom_boxplot() +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Spectral Slope Ratio")

Fig_3_sw_sr_boxplots

ggsave("Figures/Fig_3_SurfaceWater_Sr_boxplots.pdf", plot = last_plot())

#' Figure 5
#' Time-series of spectral slope ratios across sites
#' Fixed axes
Fig_5_sw_sr_time_fixed <- ggplot(swc_suva_slopes, aes(collectDate, Sr)) +
  geom_point() +
  geom_line() +
  facet_wrap(~siteID, scales = "fixed") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        axis.title.x = element_blank()) +
  ylab("Spectral Slope Ratio")
  
Fig_5_sw_sr_time_fixed

ggsave("Figures/Fig_5_SurfaceWater_Sr_time_fixed.pdf", plot = last_plot())

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

#' Figure 4
#' Boxplots of spectral slope ratios across sites
Fig_4_gw_sr_boxplots <- ggplot(gwc_suva_slopes, aes(siteID, Sr)) +
  geom_boxplot() +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Spectral Slope Ratio")

Fig_4_gw_sr_boxplots

ggsave("Figures/Fig_4_GroundWater_Sr_boxplots.pdf", plot = last_plot())

#' Figure 6
#' Time-series of spectral slope ratios across sites
#' Fixed axes
Fig_6_gw_sr_time_fixed <- ggplot(gwc_suva_slopes, aes(collectDate, Sr)) +
  geom_point() +
  geom_line() +
  facet_wrap(~siteID, scales = "fixed") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        axis.title.x = element_blank()) +
  ylab("Spectral Slope Ratio")

Fig_6_gw_sr_time_fixed

ggsave("Figures/Fig_6_GroundWater_Sr_time_fixed.pdf", plot = last_plot())


#' Groundwater vs. Surface Water Comparisons
#' Boxplots
swc_suva_slopes$type <- "Surface Water"
gwc_suva_slopes$type <- "Groundwater"
suva_slopes <- rbind(swc_suva_slopes, gwc_suva_slopes)

Fig_7_sw_gw_sr_boxplots <- ggplot(suva_slopes, aes(x = siteID, y = Sr,
                                                   color = type)) +
geom_boxplot() +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Spectral Slope Ratio")

Fig_7_sw_gw_sr_boxplots

ggsave("Figures/Fig_7_Ground_Surface_Sr_boxplots.pdf", plot = last_plot())

#' Figure 8
#' Time-series of spectral slope ratios across sites
#' Fixed axes
Fig_8_sw_gw_sr_time_fixed <- ggplot(suva_slopes, aes(x = collectDate, y = Sr,
                                                         color = type)) +
  geom_point() +
  geom_line() +
  facet_wrap(~siteID, scales = "fixed") +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0),
        axis.title.x = element_blank()) +
  ylab("Spectral Slope Ratio")

Fig_8_sw_gw_sr_time_fixed

ggsave("Figures/Fig_8_Ground_Surface_Sr_time_fixed.pdf", plot = last_plot())


#' Other data visualization plots for surface water
ggplot(swc_suva_slopes, aes(Sr, fill = siteID)) +
  geom_density(alpha = 0.2)

ggplot(swc_suva_slopes, aes(collectDate, Sr)) +
  geom_point() +
  geom_line() +
  facet_wrap(~siteID, scales = "free")

variation <- swc_suva_slopes %>% 
  group_by(siteID) %>% 
  summarise(var = var(Sr),
            mean = mean(Sr))

overall_var = var(swc_suva_slopes$Sr)
overall_mean = mean(swc_suva_slopes$Sr)

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
