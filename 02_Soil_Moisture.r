#created: 11/30/2020 last updated: 06/24/2024
#Title: GMDR Site Soil Moisture Analysis
#Authors: Kevin and Ashley

### Set up workspace
library(tidyverse)
library(nlme)
library(car)
library(lubridate)
library(ggthemes)
library(performance)
library(emmeans)
rm(list=ls()) # clean up

# Set working directory and write location
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\soil_moisture\\")
figs_to <- "C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\figures\\"

### Set graphing parameters
source("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\GMDR-NPP-Manuscript\\99_graph_format.R")

### Create Standard Error function 
SE_function<-function(x,na.rm=na.rm){
  SD = sd(x,na.rm=TRUE)
  denom = sqrt(length(x[!is.na(x)]))
  SE=SD/denom
  return(SE)
}

### Read in anova_t3 function
source("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\GMDR-NPP-Manuscript\\99_lme_anova_type3_function.R")

###
### Read in Data and Key, clean up data frames
###
{
Plot_Key <- read.csv("..\\Plot_Treatment_Key_020424.csv") %>%
  rename(PlotID=plot, Site=site)

###
### Thunder BAsin
smoist_tb_raw <- read.csv("soilmoisture_tb_19-23.csv") %>%
  mutate(Date = as.Date(Date, "%m/%d/%Y")) %>%
  mutate(Year=year(Date),
         Month=month(Date),
         Day=day(Date),
         DOY=yday(Date)) %>%
  select(-Drought, -Grazing, -Block, -Paddock, -Plot) %>%
  left_join(Plot_Key, by=c("Site","PlotID"))

## Create DOY2 to merge on to Year and DOY -- this will be used instead of DOY to avoid issues with time points that extended over 2 days -- means and stats that need to have discrete time points
tb_doy2 <- smoist_tb_raw %>%
  group_by(Year, DOY) %>%
  summarize(pmeas=length(Soil.Moisture)/54) %>%
  ungroup() %>%
  group_by(Year) %>%
  mutate(time_period=ceiling(cumsum(pmeas))) %>%
  ungroup() %>%
  group_by(Year, time_period) %>%
  mutate(DOY2 = min(DOY)) %>%
  ungroup() %>%
  dplyr::select(-pmeas, -time_period)

smoist_tb <- smoist_tb_raw %>%
  full_join(tb_doy2, by=c("Year", "DOY"))

### Take a quick look to make sure we have the correct number of values and such
with(smoist_tb_raw, table(Year, DOY))
with(smoist_tb, table(Year, DOY2))
# ggplot(filter(smoist_tb_raw,Date=="2021-07-14"), aes(PlotID, Soil.Moisture, col=factor(rainfall_reduction))) +
#   geom_point()

### Calculate means and standard error for each time point (use DOY2 because it combines sampling points that extended to multiple days)
smoist_tb_means <- smoist_tb %>%
  group_by(Year, DOY2, rainfall_reduction) %>%
  summarize(smoist_mean = mean(Soil.Moisture, na.rm=T),
            smoist_se = SE_function(Soil.Moisture)) %>%
  ungroup()

### Have a quick look
ggplot(smoist_tb_means, aes(DOY2, smoist_mean, ymin=smoist_mean-smoist_se, ymax=smoist_mean+smoist_se, col=factor(rainfall_reduction))) +
  geom_errorbar(width=0.1) +
  geom_point() +
  geom_path() +
  theme_few() +
  facet_grid(~Year)

###
### Fort Keogh
fk_smoist_19_20 <- read.csv("soil moisture fk and tb 2019-2020_2021Jan26.csv") %>% # no soil temperature here as well
  filter(Site=="FK") %>%
  rename(date=Date, block=Block, plotID=Plot, drought=Drought, soil_moist=Soil_Moisture) %>%
  mutate(soil_temp = NA) %>%
  dplyr::select(date, block, plotID, drought, soil_temp, soil_moist) %>%
  mutate(date=as.Date(date, "%Y-%m-%d"),
         drought=as.character(drought))



fk_smoist_21 <- read.csv("soilmoisture_fk_21.csv") %>% ## Note that the soil temp is missing from this -- it's available in the individual sheets but it would just need to be compiled
  dplyr::select(-year) %>%
  pivot_longer(-block:-drought.treatment, names_to="date", values_to="soil_moist") %>%
  mutate(date=substr(date, 2, nchar(date))) %>%
  separate(date, into=c('month','day','year')) %>%
  mutate(date = as.Date(paste(month, day, year, sep="/"), "%m/%d/%Y")) %>%
  rename(drought=drought.treatment) %>%
  mutate(soil_temp=NA) %>%
  dplyr::select(date, block, plotID, drought, soil_temp, soil_moist) %>%
  mutate(drought=as.character(drought))


fk_smoist_22 <- read.csv("soilmoisture_fk_22.csv") %>%
  mutate(date=as.Date(date, "%m/%d/%Y"))

fk_smoist_23 <- read.csv("soilmoisture_fk_23.csv") %>%
  mutate(date=as.Date(date, "%m/%d/%Y"))

fk_smoist_19_23 <- bind_rows(
  fk_smoist_19_20,
  fk_smoist_21,
  fk_smoist_22,
  fk_smoist_23
) %>%
  mutate(Site="FK") %>%
  mutate(Year=year(date),
         Month=month(date),
         Day=day(date),
         DOY=yday(date)) %>%
  dplyr::select(-block) %>%
  rename(PlotID=plotID,
         drought_old=drought,
         Soil.Moisture=soil_moist,
         Date=date) %>%
  mutate(DOY2=DOY) %>%
  left_join(Plot_Key, by=c("Site","PlotID")) %>%
  dplyr::select(Date, Year, Month, Day, DOY, Site, PlotID, block, paddock, rainfall_reduction,
                drought, grazing_category, grazing_treatment, livestock_util_2019, livestock_util_2020,
                livestock_util_2021, Soil.Moisture, DOY2)
  
smoist_all_temp <- smoist_tb %>%
  rename(soil_temp=TemperatureC) %>%
  bind_rows(fk_smoist_19_23)

# Add in time period for model runs
time_period_key <- smoist_all_temp %>%
  dplyr::select(Site, Year, DOY2) %>%
  unique(.) %>%
  group_by(Site, Year) %>%
  mutate(Time_period = 1:length(DOY2)) %>%
  ungroup()

smoist_all <- smoist_all_temp %>%
  full_join(time_period_key, by=c("Site","Year","DOY2")) %>%
  dplyr::select(-drought) %>%
  rename(Plot=PlotID, Grazing=grazing_category, Drought=rainfall_reduction, Block=block, Paddock=paddock)


### Check that we have the correct number of measurements for all dates and years -- looks good
with(smoist_all, table(DOY2, Year, Site))

# ggplot(smoist_all, aes(Date, Soil.Moisture, col=factor(rainfall_reduction))) +
#   geom_point() + facet_grid(Site~.)

### Calculate means and se's for plotting
smoist_means_all <- smoist_all %>%
  group_by(Site, Year, DOY2, Drought) %>%
  summarize(smoist_mean = mean(Soil.Moisture, na.rm=T),
            smoist_se = SE_function(Soil.Moisture)) %>%
  ungroup()

}

###
### Plotting soil moisture and precip (daily)
###
{
# Theme set -- text size, background, and grid lines
theme_set(theme_few())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20),
             legend.position="none")
#Colors
All_Drought <- c('#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603')
Sub_Drought <- c('#feedde','#a63603')
Individuals <- c("#CC79A7","#009E73")
#Drought Symbol
#Drought_Symbol <- c(25,24,23,22,21)
Drought_Symbol <- c(21,22,23,24,25)
Sub_droughtSymbol <- c(21,25)

#designating colors for drought and grazing treatments
droughtColor <- c('#feedde', '#fdbe85', '#fd8d3c', '#e6550d', '#a63603')
grazingColor <- c('#ABDEFF', '#469BEC', '#6D882B') #from HHMMM to MMMMM to MLLMM
droughtSymbol <- c(21,22,23,24,25) # Allows for fill

soil_moisture_figure <- ggplot(smoist_means_all, aes(x=DOY2, y=smoist_mean,
                                                        ymin=smoist_mean-smoist_se,
                                                        ymax=smoist_mean+smoist_se,
                                                        group=Drought, col=as.factor(Drought), 
                                                        fill=as.factor(Drought), shape=as.factor(Drought))) +
  geom_errorbar(col="black",width=1) +
  geom_path(col="black", aes(lty=as.factor(Drought))) +
  geom_point(color="black",size=2 ) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  scale_colour_manual(values=droughtColor, name='Drought\nTreatment') +
  scale_shape_manual(values=droughtSymbol) +
  xlab("Day of year") + ylab("Soil Moisture (%)") +
  theme_update(axis.title.x=element_text(size=14, vjust=-0.35), 
               axis.text.x=element_text(size=14),
               axis.title.y=element_text(size=14, angle=90, vjust=0.5), 
               axis.text.y=element_text(size=14), legend.position= c(5,29)) +
  facet_grid(Site~Year)

pdf(file=paste0(figs_to,"SoilMoisture_BothSitesAllYears_",Sys.Date(),".pdf"), width=14, height=6, useDingbats = F)
print(soil_moisture_figure)
dev.off()

## Add on precip data (plot separately but with same dimensions and will overlay in inkscape)
tb_ppt_daily <- read.csv("..\\precipitation\\TBPrecip_edited.csv") %>% 
  dplyr::select(-X, -Notes) %>% 
  mutate(date=as.Date(date, "%m/%d/%y")) %>%
  mutate(site="TB") %>%
  dplyr::select(site, date, year, month, day, doy, ppt_mm)

fk_ppt_daily <- read.csv("..\\precipitation\\FKPrecip_edited.csv") %>% 
  select (-X, -PRCP_ATTRIBUTES , -Notes) %>%
  mutate(date=as.Date(date, "%m/%d/%Y")) %>%
  mutate(site="FK") %>%
  rename(ppt_mm = PRCP) %>%
  dplyr::select(site, date, year, month, day, doy, ppt_mm)

ppt_daily_all <- tb_ppt_daily %>%
  bind_rows(fk_ppt_daily)

daily_ppt_fig <- ggplot(filter(ppt_daily_all,doy %in% 100:300 & year %in% 2018:2023), aes(x=doy, y=ppt_mm)) +
                    geom_col() +
                    facet_grid(site~year) +
                    theme_few() +
                    ylab("Precipitation (mm)") + xlab("Day of year")

pdf(file=paste0(figs_to,"daily ppt_BothSitesAllYears_",Sys.Date(),".pdf"), width=9.14, height=2, useDingbats = F)
print(daily_ppt_fig)
dev.off()

}

###
### Run models
###
{
smoist_all$Time_period <- as.factor(smoist_all$Time_period)

  ## Analyze sites seperately and by year

# Model with TB alone 2019

  tb_2019_smoist_lme <- lme(Soil.Moisture ~ Time_period*Grazing + Time_period*as.factor(Drought) + Grazing*as.factor(Drought)
                         , data=subset(smoist_all, Site=="TB" & Year==2019)
                         , random = ~1 |Block/Paddock/Plot
                         , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
  #                       , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                         , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)
  
  anova.lme(tb_2019_smoist_lme, type="marginal")
  
  qqnorm(tb_2019_smoist_lme, abline = c(0,1)) ## qqplot
  AIC(tb_2019_smoist_lme)
  plot(tb_2019_smoist_lme)
  hist(tb_2019_smoist_lme$fitted)
  hist(subset(smoist_all, Site=="TB" & Year==2019)$Soil.Moisture)

  emmeans(tb_2019_smoist_lme, pairwise ~ Drought, by="Time_period", adjust="sidak")
  emmeans(abun_grassy_lme, pairwise ~ fxn_group, by=c("topo","year"), adjust="sidak")
  
smoist_model_tb_2019 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="TB" & Year==2019)
)

# Model with TB alone 2020
smoist_model_tb_2020 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="TB" & Year==2020)
)  
# Model with TB alone 2021
smoist_model_tb_2021 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="TB" & Year==2021)
)  
# Model with TB alone 2022
smoist_model_tb_2022 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="TB" & Year==2022)
) 

# Model with TB alone 2023
smoist_model_tb_2023 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="TB" & Year==2023)
)  

# Model with FK alone 2019
smoist_model_fk_2019 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="FK" & Year==2019)
)  

# Model with FK alone 2020
smoist_model_fk_2020 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="FK" & Year==2020)
) 

# Model with FK alone 2021
smoist_model_fk_2021 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="FK" & Year==2021)
) 

# Model with FK alone 2022
smoist_model_fk_2022 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="FK" & Year==2022)
) 

# Model with FK alone 2023
smoist_model_fk_2023 <-   anova_t3(IndVars=c('Time_period','Grazing','Drought'),
                                   DepVar='Soil.Moisture',
                                   RndForm='~1 |Block/Paddock/Plot',
                                   Data=subset(smoist_all, Site=="FK" & Year==2023)
) 

### Comments from site-year models
# No significant effects of grazing main or interactive effects
# Significant effects of drought, time period, and drought:time period at all sites in 2019-2020 (update -- there is also a significatnt drought x time period effect in FK 2021)
# Now to split by time period in these years to see which time periods there are significant differences among drought treatments

### TB 2019 by time period
time_period_vec_tb19 <- levels(factor(subset(smoist_all, Site=="TB" & Year==2019)$Time_period))
smoist_model_tb_2019_byperiod_master <- {}

for(TIME in 1:length(time_period_vec_tb19)){
  data_temp <- smoist_all %>%
    filter(Site=="TB" & Year==2019 & Time_period==time_period_vec_tb19[TIME])
  model_temp <- lme(Soil.Moisture ~ Drought
                    , data=data_temp
                    , random = ~1 |Block/Paddock/Plot
                    , na.action = na.omit)
  model_out_temp <- data.frame(Site="TB", 
                               Year="2019", 
                               Time_period=time_period_vec_tb19[TIME],
                               Model_term=row.names(anova.lme(model_temp, type="marginal")[2,]),
                               anova.lme(model_temp, type="marginal")[2,])
  smoist_model_tb_2019_byperiod_master <- rbind(smoist_model_tb_2019_byperiod_master, model_out_temp)
  rm(data_temp, model_temp, model_out_temp)
}
smoist_model_tb_2019_byperiod_master <- smoist_model_tb_2019_byperiod_master %>%
  mutate(sig_flag = ifelse(p.value<0.05,"*",ifelse(p.value<0.1,".","")))

### TB 2020 by time period
time_period_vec_tb20 <- levels(factor(subset(smoist_all, Site=="TB" & Year==2020)$Time_period))
smoist_model_tb_2020_byperiod_master <- {}

for(TIME in 1:length(time_period_vec_tb20)){
  data_temp <- smoist_all %>%
    filter(Site=="TB" & Year==2020 & Time_period==time_period_vec_tb20[TIME])
  model_temp <- lme(Soil.Moisture ~ Drought
                    , data=data_temp
                    , random = ~1 |Block/Paddock/Plot
                    , na.action = na.omit)
  model_out_temp <- data.frame(Site="TB", 
                               Year="2020", 
                               Time_period=time_period_vec_tb20[TIME],
                               Model_term=row.names(anova.lme(model_temp, type="marginal")[2,]),
                               anova.lme(model_temp, type="marginal")[2,])
  smoist_model_tb_2020_byperiod_master <- rbind(smoist_model_tb_2020_byperiod_master, model_out_temp)
  rm(data_temp, model_temp, model_out_temp)
}
smoist_model_tb_2020_byperiod_master <- smoist_model_tb_2020_byperiod_master %>%
  mutate(sig_flag = ifelse(p.value<0.05,"*",ifelse(p.value<0.1,".","")))

### FK 2019 by time period
time_period_vec_fk19 <- levels(factor(subset(smoist_all, Site=="FK" & Year==2019)$Time_period))
smoist_model_fk_2019_byperiod_master <- {}

for(TIME in 1:length(time_period_vec_fk19)){
  data_temp <- smoist_all %>%
    filter(Site=="FK" & Year==2019 & Time_period==time_period_vec_fk19[TIME])
  model_temp <- lme(Soil.Moisture ~ Drought
                    , data=data_temp
                    , random = ~1 |Block/Paddock/Plot
                    , na.action = na.omit)
  model_out_temp <- data.frame(Site="FK", 
                               Year="2019", 
                               Time_period=time_period_vec_fk19[TIME],
                               Model_term=row.names(anova.lme(model_temp, type="marginal")[2,]),
                               anova.lme(model_temp, type="marginal")[2,])
  smoist_model_fk_2019_byperiod_master <- rbind(smoist_model_fk_2019_byperiod_master, model_out_temp)
  rm(data_temp, model_temp, model_out_temp)
}
smoist_model_fk_2019_byperiod_master <- smoist_model_fk_2019_byperiod_master %>%
  mutate(sig_flag = ifelse(p.value<0.05,"*",ifelse(p.value<0.1,".","")))


### FK 2020 by time period
time_period_vec_fk20 <- levels(factor(subset(smoist_all, Site=="FK" & Year==2020)$Time_period))
smoist_model_fk_2020_byperiod_master <- {}

for(TIME in 1:length(time_period_vec_fk20)){
  data_temp <- smoist_all %>%
    filter(Site=="FK" & Year==2020 & Time_period==time_period_vec_fk20[TIME])
  model_temp <- lme(Soil.Moisture ~ Drought
                    , data=data_temp
                    , random = ~1 |Block/Paddock/Plot
                    , na.action = na.omit)
  model_out_temp <- data.frame(Site="FK", 
                               Year="2020", 
                               Time_period=time_period_vec_fk20[TIME],
                               Model_term=row.names(anova.lme(model_temp, type="marginal")[2,]),
                               anova.lme(model_temp, type="marginal")[2,])
  smoist_model_fk_2020_byperiod_master <- rbind(smoist_model_fk_2020_byperiod_master, model_out_temp)
  rm(data_temp, model_temp, model_out_temp)
}
smoist_model_fk_2020_byperiod_master <- smoist_model_fk_2020_byperiod_master %>%
  mutate(sig_flag = ifelse(p.value<0.05,"*",ifelse(p.value<0.1,".","")))

### FK 2021 by time period
time_period_vec_fk21 <- levels(factor(subset(smoist_all, Site=="FK" & Year==2021)$Time_period))
smoist_model_fk_2021_byperiod_master <- {}

for(TIME in 1:length(time_period_vec_fk21)){
  data_temp <- smoist_all %>%
    filter(Site=="FK" & Year==2021 & Time_period==time_period_vec_fk21[TIME])
  model_temp <- lme(Soil.Moisture ~ Drought
                    , data=data_temp
                    , random = ~1 |Block/Paddock/Plot
                    , na.action = na.omit)
  model_out_temp <- data.frame(Site="FK", 
                               Year="2021", 
                               Time_period=time_period_vec_fk21[TIME],
                               Model_term=row.names(anova.lme(model_temp, type="marginal")[2,]),
                               anova.lme(model_temp, type="marginal")[2,])
  smoist_model_fk_2021_byperiod_master <- rbind(smoist_model_fk_2021_byperiod_master, model_out_temp)
  rm(data_temp, model_temp, model_out_temp)
}
smoist_model_fk_2021_byperiod_master <- smoist_model_fk_2021_byperiod_master %>%
    mutate(sig_flag = ifelse(p.value<0.05,"*",ifelse(p.value<0.1,".","")))

### Combine model results from annual runs 
yearly_model_out <- smoist_model_tb_2019 %>%
  mutate(Site="TB",
         Year=2019,
         model_effect = row.names(smoist_model_tb_2019)) %>%
  bind_rows(smoist_model_tb_2020 %>%
              mutate(Site="TB",
                     Year=2020,
                     model_effect = row.names(smoist_model_tb_2020))
  ) %>%
  bind_rows(smoist_model_tb_2021 %>%
              mutate(Site="TB",
                     Year=2021,
                     model_effect = row.names(smoist_model_tb_2021))
  ) %>%
  bind_rows(smoist_model_tb_2022 %>%
              mutate(Site="TB",
                     Year=2022,
                     model_effect = row.names(smoist_model_tb_2022))
  ) %>%
  bind_rows(smoist_model_tb_2023 %>%
              mutate(Site="TB",
                     Year=2023,
                     model_effect = row.names(smoist_model_tb_2023))
  ) %>%
  bind_rows(smoist_model_fk_2019 %>%
              mutate(Site="FK",
                     Year=2019,
                     model_effect = row.names(smoist_model_fk_2019))
  ) %>%
  bind_rows(smoist_model_fk_2020 %>%
              mutate(Site="FK",
                     Year=2020,
                     model_effect = row.names(smoist_model_fk_2020))
  ) %>%
  bind_rows(smoist_model_fk_2021 %>%
              mutate(Site="FK",
                     Year=2021,
                     model_effect = row.names(smoist_model_fk_2021))
  ) %>%
  bind_rows(smoist_model_fk_2022 %>%
              mutate(Site="FK",
                     Year=2022,
                     model_effect = row.names(smoist_model_fk_2022))
  ) %>%
  bind_rows(smoist_model_fk_2023 %>%
              mutate(Site="FK",
                     Year=2023,
                     model_effect = row.names(smoist_model_fk_2023))
  ) %>%
  rename(p_value="p-value", F_value="F-value") %>%
  dplyr::select(Site, Year, model_effect, numDF, denDF, F_value, p_value) %>%
  mutate(sig_flag = ifelse(p_value<0.05,"*",ifelse(p_value<0.1,".","")))


### Write files to csv format
write.csv(smoist_all, file="soil moisture fk and tb 2019-2023_2024_07_11.csv", row.names=F)

#by site and year models
write.csv(yearly_model_out, file=paste0("model_out\\soil moisture annual model output_both sites all years_",Sys.Date(),".csv"), row.names=F)

# by time period models
write.csv(smoist_model_tb_2019_byperiod_master, file=paste0("model_out\\soil moisture model output_by time period_tb 2019_", Sys.Date(), ".csv"), row.names=F)
write.csv(smoist_model_tb_2020_byperiod_master, file=paste0("model_out\\soil moisture model output_by time period_tb 2020_", Sys.Date(), ".csv"), row.names=F)
write.csv(smoist_model_fk_2019_byperiod_master, file=paste0("model_out\\soil moisture model output_by time period_fk 2019_", Sys.Date(), ".csv"), row.names=F)
write.csv(smoist_model_fk_2020_byperiod_master, file=paste0("model_out\\soil moisture model output_by time period_fk 2020_", Sys.Date(), ".csv"), row.names=F)
write.csv(smoist_model_fk_2021_byperiod_master, file=paste0("model_out\\soil moisture model output_by time period_fk 2021_", Sys.Date(), ".csv"), row.names=F)


}
















smoist_data_2020 <- read.csv("2020_Soil Moisture.csv") %>% 
  dplyr::select(-Notes, -X, -X.1, -X.2, -Block, -Paddock, -Drought, -Grazing) %>%
  full_join(Plot_Key, by=c("Site", "Plot")) %>%
  mutate(Date = as.Date(paste(Month, Day, Year, sep="-"), "%m-%d-%Y")) %>%
  dplyr::select(Site, Year, Month, Day, Date, Time_Period, 
                Site, Block, Paddock, Plot, Drought, Grazing, Soil_Moisture)

smoist_data_2019_fk <- read.csv("2019 Compiled Moisture Data_FortKeogh.csv") %>%
#  mutate(Date=as.Date(Date), format= "%m/%d/%Y")
  separate(Date, c("Month", "Day", "Year"), sep="/") %>%
  mutate(Date = as.Date(paste(Month, Day, Year, sep="-"), "%m-%d-%Y")) %>%
  mutate(Year=as.numeric(Year),
         Month=as.numeric(Month),
         Day=as.numeric(Day)) %>%
  dplyr::select(-Drought) %>%
  left_join(Plot_Key, by=c("Site", "Block", "Plot")) %>%
  rename(Time_Period=Time_period, Soil_Moisture=Soil_moisture) %>%
  dplyr::select(Site, Year, Month, Day, Date, Time_Period, 
                Site, Block, Paddock, Plot, Drought, Grazing, Soil_Moisture)

smoist_data_2019_tb <- read.csv("Soil Moisture data 11jan2021_fromLP.csv") %>%
  filter(Year==2019) %>%
  dplyr::select(-Year) %>%
  separate(Date, c("Month", "Day", "Year"), sep="/") %>%
  mutate(Date = as.Date(paste(Month, Day, Year, sep="-"), "%m-%d-%Y")) %>%
  mutate(Year=as.numeric(Year),
         Month=as.numeric(Month),
         Day=as.numeric(Day)) %>%
  dplyr::select(-Drought, -Grazing) %>%
  left_join(Plot_Key, by=c("Site", "Block", "Paddock", "Plot")) %>%
  rename(Time_Period=Sampling.Period, Soil_Moisture=Soil.Moisture) %>%
  dplyr::select(Site, Year, Month, Day, Date, Time_Period, 
                Site, Block, Paddock, Plot, Drought, Grazing, Soil_Moisture)

### Combine into one data frame
# The drought treatment identifiers are incorrect in the TB 2020 timepoints 11-14; Getting rid of all drought identiers and replace with plot-treatment key
smoist_data_all <- smoist_data_2019_fk %>%
  bind_rows(smoist_data_2019_tb, smoist_data_2020) %>%
  dplyr::select(-Drought, -Grazing) %>%
  full_join(Plot_Key, by=c("Site", "Block", "Paddock", "Plot")) %>%
  mutate(site_plot=paste(Site, Plot, sep="_"))

#Tell R what are the factors
smoist_data_all$Grazing <- factor(smoist_data_all$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))
smoist_data_all$Block <- as.factor(smoist_data_all$Block)
smoist_data_all$Paddock <- as.factor(smoist_data_all$Paddock)
smoist_data_all$Plot <- as.factor(smoist_data_all$Plot)
smoist_data_all$Time_Period <- as.factor(smoist_data_all$Time_Period)

### Plot data
## Raw data by site
ggplot(smoist_data_all, aes(x=Date, y=Soil_Moisture, col=as.factor(Drought), 
                                                          fill=as.factor(Drought), shape=as.factor(Drought))) +
  geom_point(alpha=0.7, color="black") +
  geom_smooth(se=F) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  scale_colour_manual(values=droughtColor, name='Drought\nTreatment') +
  scale_shape_manual(values=droughtSymbol) +
  facet_grid(Site~Year, scales="free_x")
  
## Means and se's
# Create a Date-Time period key with just one date per time period for plotting purposes
date_period_key <- unique(dplyr::select(smoist_data_all, Site, Year, Date, Time_Period)) %>%
  filter(!(Site=="TB"&Year==2020&Date=="2020-06-18" )) # one time period got sampled on two dates, here I remove one of them for data alignment purposes

smoist_drought_means <- smoist_data_all %>%
  group_by(Site, Year, Time_Period, Drought) %>% 
  summarize_at(vars(Soil_Moisture),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(smoist_mean=fn1,smoist_se=fn2) %>%
  full_join(date_period_key, by=c("Site", "Year", "Time_Period"))

smoist_means_plot <- ggplot(smoist_drought_means, aes(x=Date, y=smoist_mean, ymin=smoist_mean-smoist_se, ymax=smoist_mean+smoist_se,
                                 col=as.factor(Drought), fill=as.factor(Drought), shape=as.factor(Drought))) +
  geom_errorbar(col="black",width=0.1) +
  geom_path(col="black", aes(lty=as.factor(Drought))) +
  geom_point(color="black",size=2 ) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  scale_colour_manual(values=droughtColor, name='Drought\nTreatment') +
  scale_shape_manual(values=droughtSymbol) +
  facet_grid(Site~Year, scales="free_x") +
  xlab("Date") + ylab("Soil Moisture (%)")

pdf("..//..//figures//soil_moisture_means_FK_TB_2019-2020.pdf", width=12, height=6, useDingbats = F)
print(smoist_means_plot)
dev.off()




data_temp <- smoist_all %>%
  filter(Site=="TB" & Year==2019 & Time_period==time_period_vec_tb19[TIME])
model_temp <- lme(Soil.Moisture ~ Drought
                  , data=data_temp
                  , random = ~1 |Block/Paddock/Plot
                  , na.action = na.omit)
summary(model_temp)
