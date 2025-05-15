#created: 11/30/2020 last updated: 05/14/2025
#Title: GMDR Site Soil Moisture Analysis
#Authors: Kevin Wilcox (k_wilcox@uncg.edu)

### Set up workspace
library(tidyverse)
library(nlme)
library(car)
library(lubridate)
library(ggthemes)
library(performance)
library(emmeans)
library(lmerTest)
library(ggResidpanel)
rm(list=ls()) # clean up

# Set working directory and write location
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\soil_moisture\\")
figs_to <- "C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\figures\\"
tables_to <- "C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\tables\\"

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

  ## Analyze sites separately and by year
  ## Note, to look at differences among drought treatments within a time period, I am running emtrends and using 
  ## a slope estimate with 95CI not including 0 as a significant time point

  ### TB 2019
  ###
  tb_2019_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                         , data=subset(smoist_all, Site=="TB" & Year==2019)
                         , random = ~1 |Block/Paddock/Plot
                         , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                         , control=lmeControl(returnObject=TRUE)
                         , na.action = na.omit)
  
  # Model diagnostics
  plot(tb_2019_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(tb_2019_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="TB" & Year==2019)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(tb_2019_smoist_lme, type="marginal")
  anova_tb_smoist_2019 <- data.frame(effect=row.names(anova_temp), anova_temp,site="TB",year=2019)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2019_tb <- as.data.frame(emtrends(tb_2019_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="TB", year=2019)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2019_tb), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)

  ### TB 2020 ##
  ###
   tb_2020_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="TB" & Year==2020)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(tb_2020_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(tb_2020_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="TB" & Year==2020)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(tb_2020_smoist_lme, type="marginal")
  anova_tb_smoist_2020 <- data.frame(effect=row.names(anova_temp),anova_temp, site="TB",year=2020)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2020_tb <- as.data.frame(emtrends(tb_2020_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="TB",
           year=2020)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2020_tb), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  
  ### TB 2021 ##
  ###
  tb_2021_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="TB" & Year==2021)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(tb_2021_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(tb_2021_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="TB" & Year==2021)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(tb_2021_smoist_lme, type="marginal")
  anova_tb_smoist_2021 <- data.frame(effect=row.names(anova_temp), anova_temp, site="TB",year=2021)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2021_tb <- as.data.frame(emtrends(tb_2021_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="TB",
           year=2021)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2021_tb), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  ### TB 2022 ##
  ###
  tb_2022_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="TB" & Year==2022)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(tb_2022_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(tb_2022_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="TB" & Year==2022)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(tb_2022_smoist_lme, type="marginal")
  anova_tb_smoist_2022 <- data.frame(effect=row.names(anova_temp), anova_temp, site="TB",year=2022)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2022_tb <- as.data.frame(emtrends(tb_2022_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="TB",
           year=2022)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2022_tb), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  ### TB 2023 ##
  ###
  tb_2023_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="TB" & Year==2023)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(tb_2023_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(tb_2023_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="TB" & Year==2023)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(tb_2023_smoist_lme, type="marginal")
  anova_tb_smoist_2023 <- data.frame(effect=row.names(anova_temp),anova_temp, site="TB",year=2023)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2023_tb <- as.data.frame(emtrends(tb_2023_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="TB", year=2023)
  
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2023_tb), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
    
  ### FK 2019
  ###
  fk_2019_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="FK" & Year==2019)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(fk_2019_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(fk_2019_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="FK" & Year==2019)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(fk_2019_smoist_lme, type="marginal")
  anova_fk_smoist_2019 <- data.frame(effect=row.names(anova_temp), anova_temp,site="FK",year=2019)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2019_fk <- as.data.frame(emtrends(fk_2019_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="FK", year=2019)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2019_fk), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  ### FK 2020 ##
  ###
  fk_2020_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="FK" & Year==2020)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(fk_2020_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(fk_2020_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="FK" & Year==2020)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(fk_2020_smoist_lme, type="marginal")
  anova_fk_smoist_2020 <- data.frame(effect=row.names(anova_temp),anova_temp, site="FK",year=2020)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2020_fk <- as.data.frame(emtrends(fk_2020_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="FK",
           year=2020)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2020_fk), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  
  ### FK 2021 ##
  ###
  fk_2021_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="FK" & Year==2021)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(fk_2021_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(fk_2021_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="FK" & Year==2021)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(fk_2021_smoist_lme, type="marginal")
  anova_fk_smoist_2021 <- data.frame(effect=row.names(anova_temp), anova_temp, site="FK",year=2021)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2021_fk <- as.data.frame(emtrends(fk_2021_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="FK",
           year=2021)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2021_fk), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  ### FK 2022 ##
  ###
  fk_2022_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="FK" & Year==2022)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(fk_2022_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(fk_2022_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="FK" & Year==2022)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(fk_2022_smoist_lme, type="marginal")
  anova_fk_smoist_2022 <- data.frame(effect=row.names(anova_temp), anova_temp, site="FK",year=2022)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2022_fk <- as.data.frame(emtrends(fk_2022_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="FK",
           year=2022)
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2022_fk), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  
  ### FK 2023 ##
  ###
  fk_2023_smoist_lme <- lme(sqrt(Soil.Moisture) ~ Time_period*Grazing + Time_period*Drought + Grazing*Drought
                            , data=subset(smoist_all, Site=="FK" & Year==2023)
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot) #(AR1 AIC 4986.65, CS AIC 4993.83 -- going with AR1)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  
  # Model diagnostics
  plot(fk_2023_smoist_lme, type=c("p","smooth"), col.line=1)
  qqnorm(fk_2023_smoist_lme, abline = c(0,1)) ## qqplot
  hist(sqrt(subset(smoist_all, Site=="FK" & Year==2023)$Soil.Moisture))
  
  # Save ANCOVA results
  anova_temp <- anova.lme(fk_2023_smoist_lme, type="marginal")
  anova_fk_smoist_2023 <- data.frame(effect=row.names(anova_temp),anova_temp, site="FK",year=2023)
  rm(anova_temp)
  
  # Save emtrends
  emtrends_2023_fk <- as.data.frame(emtrends(fk_2023_smoist_lme, "Time_period", var="Drought")) %>%
    mutate(site="FK", year=2023)
  
  # Plot to see significant drought slopes by time period
  ggplot(as.data.frame(emtrends_2023_fk), aes(Time_period, Drought.trend, ymin=lower.CL, ymax=upper.CL)) + geom_errorbar() + geom_point() + geom_hline(yintercept=0)
  

  anova_smoist_all <- anova_fk_smoist_2019 %>%
    bind_rows(anova_fk_smoist_2020, 
              anova_fk_smoist_2021,
              anova_fk_smoist_2022,
              anova_fk_smoist_2023,
              anova_tb_smoist_2020, 
              anova_tb_smoist_2021,
              anova_tb_smoist_2022,
              anova_tb_smoist_2023
              )
  emtrends_smoist_all <- emtrends_2019_fk %>%
    bind_rows(emtrends_2020_fk,
              emtrends_2021_fk,
              emtrends_2022_fk,
              emtrends_2023_fk,
              emtrends_2019_tb,
              emtrends_2020_tb,
              emtrends_2021_tb,
              emtrends_2022_tb,
              emtrends_2023_tb
              )
    
write.table(emtrends_smoist_all, file=paste0(tables_to, "soil moisture emtrends table both sites all years.csv"), sep=",", row.names=F)
write.table(anova_smoist_all, file=paste0(tables_to, "soil moisture ANCOVA table both sites all years.csv"), sep=",", row.names=F)


}