library(tidyverse)

all_sites <- read.csv('Data/NEON_Field_Site_Metadata_20251125.csv')

aquatic <- all_sites %>% 
  filter(site_type == 'Core Aquatic' | site_type == 'Gradient Aquatic') 


file_list <- list.files('Data/SUVA_data/',full.names = TRUE)

df_list <- list()
for(i in 1:length(file_list)){
  
  df_list[[i]] = read.csv(file_list[i]) %>% 
    select(-X)
  
}

suva_data <- do.call(rbind,df_list)

suva <- suva_data %>% 
  select(-sampleID,-sampleID.wavelength,-domainID) %>% 
  mutate(date = ymd_hms(collectDate)) %>% 
  select(-collectDate)

ggplot(suva,aes(wavelength,SUVA,group=date))+
  geom_line()+
  facet_wrap(~siteID,scales='free')


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
