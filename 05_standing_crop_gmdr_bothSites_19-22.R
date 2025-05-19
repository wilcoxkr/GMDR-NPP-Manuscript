### Cleaning plotting and analysis -- Standing crop data
###
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created May 9, 2025; Last updated May 9, 2025

### Set up workspace
library(tidyverse)
library(nlme)
library(car)
library(lubridate)
library(emmeans)
library(performance)
library(ggthemes)
rm(list=ls()) # clean up

# Set working directory
setwd("C:\\Users\\wilco\\OneDrive - UNCG\\Current projects\\GMDR\\data\\standing_crop_roots\\") # wilcox personal laptop
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\standing_crop_roots\\") # wilcox desktop UNCG

### Set graphing parameters
source("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\graph_format.R")

### Read in Anova t3 function
source("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\lme_anova_type3_function.R")

### Create Standard Error function 
SE_function<-function(x,na.rm=na.rm){
  SE=sd(x,na.rm=TRUE)/sqrt(length(x[!is.na(x)]))
  return(SE)
}

### Read in plot key for assigning treatments to plotID
plot_key <- read.csv("..\\Plot_Treatment_Key.csv") %>%
  rename(PlotID=Plot)

###
### Read in data and calculate standing crop
###
{
## Thunder Basin data

tb_stcrop_all <- read.csv("TB_standing crop root_19-22.csv") %>%
  mutate(stcrop_gm2 = live_root_afdm/core_area_m2) %>%
  mutate(dead_root_gm2 = (coarse_dead_afdm+FPM_som_afdm+(FPM_som_afdm/(FPM_som_afdm+FPM_root_afdm))*FPM_remain_afdm)/core_area_m2) %>%
  subset(!(Year==2019 & PlotID==10 & Sub.Sample==2)) %>% ### subsample missing data
  subset(!(Year==2019 & PlotID==21 & Sub.Sample==2)) %>%### subsample missing data
  subset(!(Year==2019 & PlotID==22 & Sub.Sample==2)) %>%### subsample missing data
  subset(!(Year==2019 & PlotID==31 & Sub.Sample==1)) %>%### subsample missing data
  subset(!(Year==2019 & PlotID==36 & Sub.Sample==1)) %>%### subsample missing data
  subset(!(Year==2019 & PlotID==54 & Sub.Sample==3)) %>%### subsample missing data
  group_by(Year, Site, Block, PlotID) %>%
  summarize(stcrop_gm2 = mean(stcrop_gm2, na.rm=T),
            som_gm2 = mean(dead_root_gm2)) %>% # 
  left_join(plot_key, by=c("Site", "Block", "PlotID")) %>%
  ungroup() %>%
  rename(year=Year, site=Site, block=Block, plot=PlotID, live_root_biomass_gm2=stcrop_gm2, som_biomass_gm2=som_gm2) %>%
  dplyr::select(year, block, plot, site, Paddock, Drought, Grazing, live_root_biomass_gm2, som_biomass_gm2)

### Quick look at the data
with(tb_stcrop_all, table(year, block)) ### All looks good!
with(tb_stcrop_all, table(year, block, Drought)) ### All looks good!

hist(tb_stcrop_all$live_root_biomass_gm2)
hist(tb_stcrop_all$som_biomass_gm2)

ggplot(tb_stcrop_all, aes(Drought, live_root_biomass_gm2)) +
  geom_jitter(width=3) +
  facet_grid(.~year) # an odd data point in 2022, see below
ggplot(tb_stcrop_all, aes(Drought, som_biomass_gm2)) +
  geom_jitter(width=3) +
  facet_grid(.~year)

## Fort Keogh data
### Read in data and calculate standing crop root biomass in g/m^2
stcrop_2018 <- read.csv("RootStandingCrop_FK_2018.csv")
stcrop_2019 <- read.csv("RootStandingCrop_FK_2019.csv")
stcrop_2020 <- read.csv("RootStandingCrop_FK_2020.csv")
stcrop_2021 <- read.csv("RootStandingCrop_FK_2021.csv")
stcrop_2022 <- read.csv("RootStandingCrop_FK_2022.csv")

core_diameter_in <- 5/8 
core_radius_cm <- core_diameter_in*2.54/2
num_of_cores <- 3
total_core_area <- pi*(core_radius_cm^2/10000)*num_of_cores

fk_stcrop_all <- stcrop_2018 %>%
  bind_rows(stcrop_2019, stcrop_2020, stcrop_2021, stcrop_2022) %>%
  mutate(site="FK") %>%
  dplyr::select(-Comments) %>%
  left_join(plot_key, by=c("site" = "Site", "block" = "Block", "plot" = "PlotID")) %>%
  dplyr::select(-tin_number:-ash_mass) %>%
  spread(key=sample_type, value=AFDM_g) %>%
  rename(fine_live='Fine live roots', fine_som='Fine SOM', live_roots='Live Roots', dead_roots='Dead Roots') %>%
  mutate(p_live = fine_live/(fine_live+fine_som)) %>%
  mutate(total_live_g = live_roots+fine_live+(Remainder*p_live)) %>%
  mutate(total_som_g = dead_roots+fine_som+(Remainder*(1-p_live))) %>%
  mutate(live_root_biomass_gm2 = total_live_g/total_core_area) %>% 
  mutate(som_biomass_gm2 = total_som_g/total_core_area)

fk_stcrop_plot_means <- fk_stcrop_all %>%
  group_by(year, block, plot, site, Paddock, Drought, Grazing) %>%
  summarize(live_root_biomass_gm2 = mean(live_root_biomass_gm2, na.rm=T),
            som_biomass_gm2 = mean(som_biomass_gm2, na.rm=T)) %>%
  ungroup()

### Quick look at the data
with(fk_stcrop_all, table(year, block)) ### All looks good!
with(fk_stcrop_all, table(year, block, Drought)) ### All looks good!
with(fk_stcrop_plot_means, table(year, block)) ### All looks good!
with(fk_stcrop_plot_means, table(year, block, Drought)) ### All looks good!

ggplot(fk_stcrop_plot_means, aes(Drought, live_root_biomass_gm2)) +
  geom_jitter(width=3) +
  facet_grid(.~year) # an odd data point in 2022, see below
ggplot(fk_stcrop_plot_means, aes(Drought, som_biomass_gm2)) +
  geom_jitter(width=3) +
  facet_grid(.~year)

### strange data point with very high live roots and very low som plot 44 (drought 25, MMMMM grazing) in 2022. Looked back at the raw data and didn't see anything odd... might remove as an outlier anyways
filter(fk_stcrop_plot_means, live_root_biomass_gm2>2000)
filter(fk_stcrop_plot_means, som_biomass_gm2<22)
min(fk_stcrop_plot_means$som_biomass_gm2)
###

### Combine both sites
###

stcrop_all <- fk_stcrop_plot_means %>%
  bind_rows(tb_stcrop_all)

with(stcrop_all, table(site, year, block))

## Take a look
ggplot(stcrop_all, aes(Drought, live_root_biomass_gm2, col=Grazing)) +
  geom_jitter(width=3) +
  geom_smooth(method="lm",se=F) +
  ylim(0,1750) +
  facet_grid(site~year) 

ggplot(stcrop_all, aes(Drought, som_biomass_gm2, col=Grazing)) +
  geom_jitter(width=3) +
  geom_smooth(method="lm",se=F) +
#  ylim(0,1750) +
  facet_grid(site~year) 

}

###
### Calculate response ratios
###
{
stcrop_controls <- stcrop_all %>%
  filter(Drought==0) %>%
  group_by(site, year, block, Paddock) %>%
  summarize(ctrl_live_root=mean(live_root_biomass_gm2, na.rm=T),
            ctrl_som=mean(som_biomass_gm2, na.rm=T)
  )

stcrop_rr <- stcrop_all %>%
  filter(Drought != 0) %>%
  full_join(stcrop_controls, by=c("year","site","block","Paddock")) %>%
  mutate(lnrr_live_root=log((live_root_biomass_gm2+0.01)/(ctrl_live_root+0.01)), pchange_live_root=((live_root_biomass_gm2+0.01)-(ctrl_live_root+0.01))/(ctrl_live_root+0.01)) %>%
  mutate(lnrr_som=log((som_biomass_gm2+0.01)/(ctrl_som+0.01)), pchange_som=((som_biomass_gm2+0.01)-(ctrl_som+0.01))/(ctrl_som+0.01)) %>%
  dplyr::select(year:Grazing, lnrr_live_root, pchange_live_root, lnrr_som, pchange_som) %>%
  filter(!(year==2022 & plot==44 & site=="FK")) ### plot 44 is a clear outlier so I am removing it (probably need to calculate outlier diagnostics on this)
  

stcrop_rr_long <- stcrop_rr %>%
  pivot_longer(cols=lnrr_live_root:pchange_som, values_to="metric_value", names_to="metric_type")

## Take a quick look at response ratios
ggplot(stcrop_rr, aes(Drought, lnrr_live_root, col=Grazing)) +
  geom_hline(yintercept=0,col="grey") +
  geom_jitter(width=3) +
  stat_smooth(method="lm",se=F)+
  facet_grid(site~year)
         
ggplot(stcrop_rr, aes(Drought, lnrr_som, col=Grazing)) +
  geom_hline(yintercept=0,col="grey") +
  geom_jitter(width=3) +
  stat_smooth(method="lm",se=F)+
  facet_grid(site~year)
}         

###
### Run models with response ratios
###
{
  ##
  ## Live root biomass
  
  ## Fort Keogh  
  fk_live_root_rr_lme_full <- lme(lnrr_live_root ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                                    , data=filter(stcrop_rr, year %in% 2019:2022 & site=="FK")
                                    , random = ~1 |block/Paddock/plot
#                                    , correlation=corCompSymm(form = ~1 |block/Paddock/plot) # aic 264.2936
                                    , correlation=corAR1(form = ~1 |block/Paddock/plot) # aic 261.87 -- going with AR1
                                    , control=lmeControl(returnObject=TRUE)
                                    , na.action = na.omit)
  anova.lme(fk_live_root_rr_lme_full, type="marginal")
  summary(fk_live_root_rr_lme_full)
  
  # check model
  qqnorm(fk_live_root_rr_lme_full, abline = c(0,1)) ## qqplot
  check_model(fk_live_root_rr_lme_full) ## residuals and normality of resids

  # because of year*grazing interaction, look at grazing effect by year  
  emmeans(fk_live_root_rr_lme_full, pairwise ~ Grazing, by="year", adjust="sidak")

  ## Thunder Basin
  tb_live_root_rr_lme_full <- lme(lnrr_live_root ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                                    , data=filter(stcrop_rr, year %in% 2019:2022 & site=="TB")
                                    , random = ~1 |block/Paddock/plot
                                    , correlation=corAR1(form = ~1 |block/Paddock/plot)
                                    , control=lmeControl(returnObject=TRUE)
                                    , na.action = na.omit)
  anova.lme(tb_live_root_rr_lme_full, type="marginal")
  summary(tb_live_root_rr_lme_full)
  
  # check model assumptions
  qqnorm(tb_live_root_rr_lme_full, abline = c(0,1)) ## qqplot
  check_model(tb_live_root_rr_lme_full) ## residuals and normality of resids

  # because of year*drought interaction, look at drought regressions by year  

  #2019
  tb_live_root_rr_lme_2019 <- lme(lnrr_live_root ~ Drought*Grazing
                               , data=filter(stcrop_rr, year == 2019 & site=="TB")
                               , random = ~1 |block/Paddock/plot
                               , na.action = na.omit)
  anova.lme(tb_live_root_rr_lme_2019, type="marginal")
  summary(tb_live_root_rr_lme_2019)
  #2020
  tb_live_root_rr_lme_2020 <- lme(lnrr_live_root ~ Drought*Grazing
                               , data=filter(stcrop_rr, year == 2020 & site=="TB")
                               , random = ~1 |block/Paddock/plot
                               , na.action = na.omit)
  anova.lme(tb_live_root_rr_lme_2020, type="marginal")
  summary(tb_live_root_rr_lme_2020)
  #2021
  tb_live_root_rr_lme_2021 <- lme(lnrr_live_root ~ Drought*Grazing
                               , data=filter(stcrop_rr, year == 2021 & site=="TB")
                               , random = ~1 |block/Paddock/plot
                               , na.action = na.omit)
  anova.lme(tb_live_root_rr_lme_2021, type="marginal")
  summary(tb_live_root_rr_lme_2021)
  #2022
  tb_live_root_rr_lme_2022 <- lme(lnrr_live_root ~ Drought*Grazing
                               , data=filter(stcrop_rr, year == 2022 & site=="TB")
                               , random = ~1 |block/Paddock/plot
                               , na.action = na.omit)
  anova.lme(tb_live_root_rr_lme_2022, type="marginal")
  summary(tb_live_root_rr_lme_2022)
  
  ### Nothing really here - weak pattern in 2021
  
  
  ##
  ## Dead roots and SOM
  
  ## Fort Keogh  
  fk_som_rr_lme_full <- lme(lnrr_som ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                               , data=filter(stcrop_rr, year %in% 2019:2022 & site=="FK")
                               , random = ~1 |block/Paddock/plot
                               #                                    , correlation=corCompSymm(form = ~1 |block/Paddock/plot) # aic 264.2936
                               , correlation=corAR1(form = ~1 |block/Paddock/plot) # aic 261.87 -- going with AR1
                               , control=lmeControl(returnObject=TRUE)
                               , na.action = na.omit)
  anova.lme(fk_som_rr_lme_full, type="marginal")
  summary(fk_som_rr_lme_full)
  
  # check model assumptions -- meh...
  qqnorm(fk_som_rr_lme_full, abline = c(0,1)) ## qqplot
  check_model(fk_som_rr_lme_full) ## residuals and normality of resids

  ## nothing significant, stop here
  
  ## Thunder Basin
  tb_som_rr_lme_full <- lme(lnrr_som ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                               , data=filter(stcrop_rr, year %in% 2019:2022 & site=="TB")
                               , random = ~1 |block/Paddock/plot
                               , correlation=corAR1(form = ~1 |block/Paddock/plot)
                               , control=lmeControl(returnObject=TRUE)
                               , na.action = na.omit)
  anova.lme(tb_som_rr_lme_full, type="marginal")
  summary(tb_som_rr_lme_full)
  
  # check model assumptions
  qqnorm(tb_som_rr_lme_full, abline = c(0,1)) ## qqplot
  check_model(tb_som_rr_lme_full) ## residuals and normality of resids
  
  # because of year*grazing interaction, look at grazing effect by year  
  emmeans(tb_som_rr_lme_full, pairwise ~ Grazing, by="year", adjust="sidak")
  
}
  

###
### Run models with non response ratio data
###
{

##
## Live roots
  
## Fort Keogh  
fk_live_root_lme_full <- lme(live_root_biomass_gm2 ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                             , data=filter(stcrop_all, year %in% 2019:2022 & site=="FK")
                             , random = ~1 |block/Paddock/plot
                             #                                    , correlation=corCompSymm(form = ~1 |block/Paddock/plot) # aic 264.2936
                             , correlation=corAR1(form = ~1 |block/Paddock/plot) # aic 261.87 -- going with AR1
                             , control=lmeControl(returnObject=TRUE)
                             , na.action = na.omit)
anova.lme(fk_live_root_lme_full, type="marginal")
summary(fk_live_root_lme_full)
hist(filter(stcrop_all, year %in% 2019:2022 & site=="FK")$live_root_biomass_gm2)
# check model
qqnorm(fk_live_root_lme_full, abline = c(0,1)) ## qqplot
check_model(fk_live_root_lme_full) ## residuals and normality of resids

## Nothing significant except year so we'll stop here

## Thunder Basin 
tb_live_root_lme_full <- lme(live_root_biomass_gm2 ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                             , data=filter(stcrop_all, year %in% 2019:2022 & site=="TB")
                             , random = ~1 |block/Paddock/plot
                             #                                    , correlation=corCompSymm(form = ~1 |block/Paddock/plot) # aic 264.2936
                             , correlation=corAR1(form = ~1 |block/Paddock/plot) # aic 261.87 -- going with AR1
                             , control=lmeControl(returnObject=TRUE)
                             , na.action = na.omit)
anova.lme(tb_live_root_lme_full, type="marginal")
summary(tb_live_root_lme_full)
hist(filter(stcrop_all, year %in% 2019:2022 & site=="TB")$live_root_biomass_gm2)
# check model
qqnorm(tb_live_root_lme_full, abline = c(0,1)) ## qqplot
check_model(tb_live_root_lme_full) ## residuals and normality of resids

## Nothing significant so we'll stop here

##
## Dead roots and SOM

## Fort Keogh  
fk_som_lme_full <- lme(som_biomass_gm2 ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                             , data=filter(stcrop_all, year %in% 2019:2022 & site=="FK"&som_biomass_gm2<1100)
                             , random = ~1 |block/Paddock/plot
                             #                                    , correlation=corCompSymm(form = ~1 |block/Paddock/plot) # aic 264.2936
                             , correlation=corAR1(form = ~1 |block/Paddock/plot) # aic 261.87 -- going with AR1
                             , control=lmeControl(returnObject=TRUE)
                             , na.action = na.omit)
anova.lme(fk_som_lme_full, type="marginal")
summary(fk_som_lme_full)
hist(filter(stcrop_all, year %in% 2019:2022 & site=="FK")$som_biomass_gm2)

# check model
qqnorm(fk_som_lme_full, abline = c(0,1)) ## qqplot
check_model(fk_som_lme_full) ## residuals and normality of resids

emmeans(fk_som_lme_full, pairwise ~ Grazing, by="year", adjust="sidak")
## Nothing significant except year so we'll stop here

## Thunder Basin 
tb_live_root_lme_full <- lme(live_root_biomass_gm2 ~ as.factor(year)*Drought + as.factor(year)*Grazing + Drought*Grazing
                             , data=filter(stcrop_all, year %in% 2019:2022 & site=="TB")
                             , random = ~1 |block/Paddock/plot
                             #                                    , correlation=corCompSymm(form = ~1 |block/Paddock/plot) # aic 264.2936
                             , correlation=corAR1(form = ~1 |block/Paddock/plot) # aic 261.87 -- going with AR1
                             , control=lmeControl(returnObject=TRUE)
                             , na.action = na.omit)
anova.lme(tb_live_root_lme_full, type="marginal")
summary(tb_live_root_lme_full)
hist(filter(stcrop_all, year %in% 2019:2022 & site=="TB")$live_root_biomass_gm2)
# check model
qqnorm(tb_live_root_lme_full, abline = c(0,1)) ## qqplot
check_model(tb_live_root_lme_full) ## residuals and normality of resids

## Nothing significant so we'll stop here


}


###
### Plotting response ratios
###

  stcrop_rr_drought_means <- stcrop_rr_long %>%
    filter(!(year==2022 & plot==44 & site=="FK")) %>% ### plot 44 is a clear outlier
    group_by(site, year, metric_type, Drought) %>%
    summarize_at(.vars=vars(metric_value), .funs=list(mean=mean, se=SE_function), na.rm=T) %>%
    ungroup() %>%
    rename(metric_mean = mean, metric_se = se)
    
  stcrop_fig <- ggplot(filter(stcrop_rr_drought_means,year%in%2019:2022 & metric_type %in% c("lnrr_live_root", "lnrr_som")), aes(x=Drought, y=metric_mean, 
                ymin=metric_mean-metric_se, ymax=metric_mean+metric_se, fill=metric_type, shape=metric_type)) +
    geom_hline(yintercept=0, col="grey") +
    geom_smooth(inherit.aes=F, data=filter(stcrop_rr_long,year%in%2019:2022&metric_type %in% c("lnrr_live_root", "lnrr_som")),
                aes(Drought, metric_value, lty=metric_type), method="lm",se=F, col="black") +
    geom_errorbar(width=3) +
    geom_point(size=2, col="black") +
    scale_shape_manual(values=c(21,24))+
    scale_fill_manual(values=c("white","black")) +
    facet_grid(site~year) +
    theme_few()
  
pdf("..//..//npp_manuscript//figures//stcrop fig_v1.pdf", width=9, height=4.25, useDingbats = F)
print(stcrop_fig)
dev.off()
  
  
  
  
  

stcrop_means <- stcrop_raw %>%
  group_by(Year, Site, Drought) %>%
  summarize(stcrop_mean = mean(stcrop_gm2, na.rm=T),
            stcrop_se = SE_function(stcrop_gm2),
            som_mean = mean(som_gm2, na.rm=T),
            som_se = SE_function(som_gm2))

ggplot(stcrop_means, aes(Drought, stcrop_mean, ymin=stcrop_mean-stcrop_se, ymax=stcrop_mean+stcrop_se)) +
  geom_errorbar(width=0) + geom_point() + geom_smooth(method="lm", se=F) + facet_grid(~Year) +
  ylab("Live roots (g/m2)") 

ggplot(stcrop_means, aes(Drought, som_mean, ymin=som_mean-som_se, ymax=som_mean+som_se)) +
  geom_errorbar(width=0) + geom_point() + geom_smooth(method="lm", se=F) + facet_grid(~Year) +
  ylab("SOM (g/m2)") 

ggsave("..//..//figures//standing_crop//tb_standing crop_19-22.png", width=10, height=4, units="in")

ggplot(stcrop_raw, aes(x=Year, y=live_root_g)) + geom_jitter()
ggplot(stcrop_raw, aes(x=Year, y=live_root_gm2)) + geom_jitter()

stcrop_raw$Grazing <- factor(stcrop_raw$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))
stcrop_raw$Year <- factor(stcrop_raw$Year)
full_model_tb_rm <- lme(stcrop_gm2 ~ Drought*Grazing + Drought*Year + Grazing*Year
                     , data=stcrop_raw
                     , random = as.formula(~1 |Block/Paddock/PlotID)
                     , correlation=corCompSymm(form = as.formula(~1 |Block/Paddock/PlotID))
                     , control=lmeControl(returnObject=TRUE)
                     , na.action = na.omit)
anova.lme(full_model_tb_rm, type="marginal")
tb_slope_means <- emtrends(full_model_tb_rm, "Year", "Drought")
pairs(tb_slope_means)

tb_stcrop_2019_lme <- lme(stcrop_gm2 ~ Drought*Grazing
                        , data=filter(stcrop_raw, Year==2019)
                        , random = as.formula(~1 |Block/Paddock/PlotID)
                        , na.action = na.omit)
anova.lme(tb_stcrop_2019_lme, type="marginal")
tb_stcrop_2020_lme <- lme(stcrop_gm2 ~ Drought*Grazing
                        , data=filter(stcrop_raw, Year==2020)
                        , random = as.formula(~1 |Block/Paddock/PlotID)
                        , na.action = na.omit)
anova.lme(tb_stcrop_2020_lme, type="marginal")
tb_stcrop_2021_lme <- lme(stcrop_gm2 ~ Drought*Grazing
                        , data=filter(stcrop_raw, Year==2021)
                        , random = as.formula(~1 |Block/Paddock/PlotID)
                        , na.action = na.omit)
anova.lme(tb_stcrop_2021_lme, type="marginal")
tb_stcrop_2022_lme <- lme(stcrop_gm2 ~ Drought*Grazing
                        , data=filter(stcrop_raw, Year==2022)
                        , random = as.formula(~1 |Block/Paddock/PlotID)
                        , na.action = na.omit)
anova.lme(tb_stcrop_2022_lme, type="marginal")

### Plot raw data
ggplot(stcrop_raw, aes(x=Drought, y=live_root_gm2)) +
  geom_point() +
  facet_wrap(~Year*Grazing) +
  theme_bw()

### Calculate means per plot
stcrop_plot_means <- stcrop_raw %>%
  mutate(som_gm2 = (CoarseDead_afdm+FPMsom_afdm+FPMremainder_afdm*(1-FPM_live_prop))/total_area_cm2*10000,
         fine_live_gm2 = (FPMroots_afdm+FPMremainder_live_afdm)/total_area_cm2*10000) %>%
  group_by(Site, Year, Block, PlotID, Paddock, Drought, Grazing) %>%
  summarize(live_root_gm2 = mean(live_root_gm2, na.rm=T),
            som_gm2 = mean(som_gm2, na.rm=T),
            fine_live_gm2 = mean(fine_live_gm2, na.rm=T))

#with(bnpp_plot_means, table(PlotID))

### Calculate means by drought alone
stcrop_drt_means <- stcrop_plot_means %>%
  group_by(Site, Year, Drought) %>%
  summarize(live_root_mean = mean(live_root_gm2, na.rm=T), 
            live_root_se=SE_function(live_root_gm2),
            som_mean = mean(som_gm2, na.rm=T), 
            som_se=SE_function(som_gm2),
            fine_live_mean = mean(fine_live_gm2, na.rm=T), 
            fine_live_se=SE_function(fine_live_gm2)
  )

stcrop_drt_noyr_means <- stcrop_plot_means %>%
  group_by(Site, Drought) %>%
  summarize(stcrop_mean = mean(live_root_gm2, na.rm=T), 
            stcrop_se=SE_function(live_root_gm2))

## all live roots
png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop drought only_", Sys.Date(), ".png"), width=8.2, height=7.0, units="in", res=600)
print(
  ggplot(stcrop_drt_means, aes(x=Drought, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se)) +
    geom_errorbar(width=1) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    facet_wrap(~Year) +
    ylab("Standing crop g/m2")
  )
dev.off()

## fine live roots only
ggplot(stcrop_drt_means, aes(x=Drought, y=fine_live_mean, ymin=fine_live_mean-fine_live_se, ymax=fine_live_mean+fine_live_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_wrap(~Year) +
  ylab("Fine live roots g/m2")

## coarse dead roots + fine som
png(file=paste0("..\\..\\figures\\standing_crop\\TB_SOM drought only_", Sys.Date(), ".png"), width=8.2, height=7.0, units="in", res=600)
print(
  ggplot(stcrop_drt_means, aes(x=Drought, y=som_mean, ymin=som_mean-som_se, ymax=som_mean+som_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_wrap(~Year) +
  ylab("SOM g/m2")
)
dev.off()

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop drought only_avg across years_", Sys.Date(), ".png"), width=8.2, height=4.2, units="in", res=600)
print(
  ggplot(stcrop_drt_noyr_means, aes(x=Drought, y=stcrop_mean, ymin=stcrop_mean-stcrop_se, ymax=stcrop_mean+stcrop_se)) +
    geom_errorbar(width=1) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    ylab("Standing crop g/m2")
)
dev.off()

### Calculate means by grazing alone
stcrop_grz_means <- stcrop_plot_means %>%
  group_by(Site, Year, Grazing) %>%
  summarize(live_root_mean = mean(live_root_gm2, na.rm=T), 
            live_root_se=SE_function(live_root_gm2),
            som_mean = mean(som_gm2, na.rm=T), 
            som_se=SE_function(som_gm2),
            fine_live_mean = mean(fine_live_gm2, na.rm=T), 
            fine_live_se=SE_function(fine_live_gm2)
  )

stcrop_grz_means$Grazing <- factor(stcrop_grz_means$Grazing, levels=c("MLLMM", "MMMMM", "HHMMM"))

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop grazing only_", Sys.Date(), ".png"), width=8.2, height=7.0, units="in", res=600)
print(
  ggplot(stcrop_grz_means, aes(x=Grazing, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se)) +
    geom_errorbar(width=0.05) +
    geom_point(size=3) +
    facet_wrap(~Year) +
    ylab("stcrop g/m2") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
)
dev.off()

ggplot(stcrop_grz_means, aes(x=Grazing, y=fine_live_mean, ymin=fine_live_mean-fine_live_se, ymax=fine_live_mean+fine_live_se)) +
  geom_errorbar(width=0.05) +
  geom_point(size=3) +
  facet_wrap(~Year) +
  ylab("fine live roots g/m2") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

png(file=paste0("..\\..\\figures\\standing_crop\\TB_SOM grazing only_", Sys.Date(), ".png"), width=8.2, height=7.0, units="in", res=600)
print(
  ggplot(stcrop_grz_means, aes(x=Grazing, y=som_mean, ymin=som_mean-som_se, ymax=som_mean+som_se)) +
  geom_errorbar(width=0.05) +
  geom_point(size=3) +
  facet_wrap(~Year) +
  ylab("SOM g/m2") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
)
dev.off()

### Calculate means by drought and grazing
stcrop_drt_grz_means <- stcrop_plot_means %>%
  group_by(Site, Year, Drought, Grazing) %>%
  summarize(live_root_mean = mean(live_root_gm2, na.rm=T), 
            live_root_se=SE_function(live_root_gm2))

stcrop_drt_grz_means$Grazing <- factor(stcrop_drt_grz_means$Grazing, levels=c("MLLMM", "MMMMM", "HHMMM"))
q3 <- c("#00AD9A", "#9183E6", "#E16A86")

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop drought and grazing_", Sys.Date(), ".png"), width=9.2, height=8.2, units="in", res=600)
print(
  ggplot(stcrop_drt_grz_means, aes(x=Drought, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se, col=Grazing)) +
    geom_errorbar(width=1) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    facet_wrap(~Year*Grazing, scales="free_y") +
    ylab("stcrop g/m2")
)
dev.off()

###
### Run models
###
{
  # read in anova_t3 function for repeated measures mixed model anova
  source("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Git projects\\Grazing-Management-for-Drought-Resilience\\lme_anova_type3_function.R")

  stcrop_full_model <- anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                DepVar='live_root_gm2',
                                RndForm='~1 |Block/Paddock/PlotID',
                                Data=stcrop_plot_means)
  
  stcrop_drought_model <- anova_t3(IndVars=c('Year','Drought'),
                                   DepVar='live_root_gm2',
                                   RndForm='~1 |Block/Paddock/PlotID',
                                   Data=stcrop_plot_means)
  
  stcrop_grazing_model <- anova_t3(IndVars=c('Year','Grazing'),
                                   DepVar='live_root_gm2',
                                   RndForm='~1 |Block/Paddock/PlotID',
                                   Data=stcrop_plot_means)
  
  grazing_2018_model <- lme(live_root_gm2 ~ Grazing
                            , data=filter(stcrop_plot_means, Year==2018)
                            , random = ~1 |Block/Paddock/PlotID
                            , na.action = na.omit)
  Anova(grazing_2018_model, type=3)
  
  grazing_2019_model <- lme(live_root_gm2 ~ Grazing
                            , data=filter(stcrop_plot_means, Year==2019)
                            , random = ~1 |Block/Paddock/PlotID
                            , na.action = na.omit)
  Anova(grazing_2019_model, type=3)
  
  grazing_2020_model <- lme(live_root_gm2 ~ Grazing
                            , data=filter(stcrop_plot_means, Year==2020)
                            , random = ~1 |Block/Paddock/PlotID
                            , na.action = na.omit)
  Anova(grazing_2020_model, type=3)
  
  grazing_2021_model <- lme(live_root_gm2 ~ Grazing
                            , data=filter(stcrop_plot_means, Year==2021)
                            , random = ~1 |Block/Paddock/PlotID
                            , na.action = na.omit)
  Anova(grazing_2021_model, type=3)
  
  grazing_2022_model <- lme(live_root_gm2 ~ Grazing
                            , data=filter(stcrop_plot_means, Year==2022)
                            , random = ~1 |Block/Paddock/PlotID
                            , na.action = na.omit)
  Anova(grazing_2022_model, type=3)
  

  fine_live_full_model <- anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                DepVar='fine_live_gm2',
                                RndForm='~1 |Block/Paddock/PlotID',
                                Data=stcrop_plot_means)

  som_full_model <- anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                DepVar='som_gm2',
                                RndForm='~1 |Block/Paddock/PlotID',
                                Data=stcrop_plot_means)
  
  }


###
### Looking at just live fine roots
###

### Calculate means per plot
LiveFine_plot_means <- stcrop_raw %>%
  group_by(Site, Year, Block, PlotID, Paddock, Drought, Grazing) %>%
  summarize(live_root_gm2 = mean(live_root_gm2, na.rm=T))

#with(bnpp_plot_means, table(PlotID))

### Calculate means by drought alone
stcrop_drt_means <- stcrop_plot_means %>%
  group_by(Site, Year, Drought) %>%
  summarize(live_root_mean = mean(live_root_gm2, na.rm=T), 
            live_root_se=SE_function(live_root_gm2))

stcrop_drt_noyr_means <- stcrop_plot_means %>%
  group_by(Site, Drought) %>%
  summarize(stcrop_mean = mean(live_root_gm2, na.rm=T), 
            stcrop_se=SE_function(live_root_gm2))

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop drought only_", Sys.Date(), ".png"), width=8.2, height=4.2, units="in", res=600)
print(
  ggplot(stcrop_drt_means, aes(x=Drought, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se)) +
    geom_errorbar(width=1) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    facet_wrap(~Year) +
    ylab("Standing crop g/m2")
)
dev.off()

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop drought only_avg across years_", Sys.Date(), ".png"), width=8.2, height=4.2, units="in", res=600)
print(
  ggplot(stcrop_drt_noyr_means, aes(x=Drought, y=stcrop_mean, ymin=stcrop_mean-stcrop_se, ymax=stcrop_mean+stcrop_se)) +
    geom_errorbar(width=1) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    ylab("Standing crop g/m2")
)
dev.off()

### Calculate means by grazing alone
stcrop_grz_means <- stcrop_plot_means %>%
  group_by(Site, Year, Grazing) %>%
  summarize(live_root_mean = mean(live_root_gm2, na.rm=T), 
            live_root_se=SE_function(live_root_gm2))

stcrop_grz_means$Grazing <- factor(stcrop_grz_means$Grazing, levels=c("MLLMM", "MMMMM", "HHMMM"))

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop grazing only_", Sys.Date(), ".png"), width=8.2, height=4.2, units="in", res=600)
print(
  ggplot(stcrop_grz_means, aes(x=Grazing, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se)) +
    geom_errorbar(width=0.05) +
    geom_point(size=3) +
    facet_wrap(~Year) +
    ylab("stcrop g/m2") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
)
dev.off()

### Calculate means by drought and grazing
stcrop_drt_grz_means <- stcrop_plot_means %>%
  group_by(Site, Year, Drought, Grazing) %>%
  summarize(live_root_mean = mean(live_root_gm2, na.rm=T), 
            live_root_se=SE_function(live_root_gm2))

stcrop_drt_grz_means$Grazing <- factor(stcrop_drt_grz_means$Grazing, levels=c("MLLMM", "MMMMM", "HHMMM"))
q3 <- c("#00AD9A", "#9183E6", "#E16A86")

png(file=paste0("..\\..\\figures\\standing_crop\\TB_stcrop drought and grazing_", Sys.Date(), ".png"), width=9.2, height=8.2, units="in", res=600)
print(
  ggplot(stcrop_drt_grz_means, aes(x=Drought, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se, col=Grazing)) +
    geom_errorbar(width=1) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    facet_wrap(~Year*Grazing, scales="free_y") +
    ylab("stcrop g/m2")
)
dev.off()

drt_model_2019 <- lme(live_root_gm2 ~ Drought*Grazing
                       , data=filter(stcrop_plot_means, Year==2019)
                       , random = ~1 |Block/Paddock/PlotID
                       , na.action = na.omit)


Anova(drt_model_2019, type=3)


drt_model_2020 <- lme(live_root_gm2 ~ Drought*Grazing
                      , data=filter(stcrop_plot_means, Year==2020)
                      , random = ~1 |Block/Paddock/PlotID
                      , na.action = na.omit)


Anova(drt_model_2020, type=3)


drt_model_2021 <- lme(live_root_gm2 ~ Drought*Grazing
                      , data=filter(stcrop_plot_means, Year==2021)
                      , random = ~1 |Block/Paddock/PlotID
                      , na.action = na.omit)


Anova(drt_model_2021, type=3)

stcrop_drt_grz_means <- stcrop_plot_means %>%
  group_by(Year, Drought, Grazing) %>%
  summarize_at(vars(live_root_gm2),.funs = c(mean, SE_function),na.rm=TRUE) %>% 
  rename(live_root_mean=fn1, live_root_se=fn2)

ggplot(stcrop_drt_grz_means, aes(x=Drought, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  ylab("Live root biomass g/m2") +
  facet_wrap(Year~Grazing)




stcrop_grz_means <- stcrop_plot_means %>%
  group_by(Year, Grazing) %>%
  summarize_at(vars(live_root_gm2),.funs = c(mean, SE_function),na.rm=TRUE) %>% 
  rename(live_root_mean=fn1, live_root_se=fn2)

stcrop_grz_means$Grazing <- factor(stcrop_grz_means$Grazing, levels=c("MLLMM", "MMMMM", "HHMMM"))

ggplot(stcrop_grz_means, aes(x=Grazing, y=live_root_mean, ymin=live_root_mean-live_root_se, ymax=live_root_mean+live_root_se)) +
  geom_errorbar(width=0.05) +
  geom_point(size=3) +
  facet_wrap(~Year) +
  ylab("Standing Crop Live Roots g/m2")



###
### Regress standing crop in 2020 with response ratio of ANPP in 2021 -- need to read in anpp_tb_rr from anpp_gmdr_forTBRI2022 presentation.R
### NOTE: nothing really here unfortunately
{
  ### Calculate standing crop response ratio
  stcrop_controls <- stcrop_plot_means %>%
    filter(Drought==0) %>%
    group_by(Site, Year, Block, Paddock) %>%
    summarize(ctrl_stcrop=mean(stcrop_gm2, na.rm=T))
  
  stcrop_rr <- stcrop_plot_means %>%
    filter(Drought != 0) %>%
    full_join(stcrop_controls, by=c("Year","Site","Block","Paddock")) %>%
    mutate(lnrr_stcrop=log(stcrop_gm2+0.01/ctrl_stcrop+0.01), 
           pchange_stcrop=((stcrop_gm2+0.01)-(ctrl_stcrop+0.01))/(ctrl_stcrop+0.01))
           
  stcrop_and_anpp_rr <- anpp_tb_rr %>%
    dplyr::select(Year:PlotID, lnrr_ANPP, pchange_ANPP, pchange_C3P, pchange_C4P) %>%
    filter(Year == 2021) %>%
    rename(anpp_year = Year) %>%
    mutate(PlotID = as.numeric(PlotID)) %>%
    left_join(
      stcrop_rr %>% filter(Year == 2020) %>% rename(stcrop_year = Year), 
      by=c("Site", "Block", "PlotID", "Paddock", "Drought", "Grazing")
      )
  
  ggplot(stcrop_and_anpp_rr, aes(x=pchange_stcrop, y=pchange_ANPP)) +
    geom_point() + geom_smooth(method="lm",se=F)
  ggplot(stcrop_and_anpp_rr, aes(x=pchange_stcrop, y=pchange_C3P)) +
    geom_point() + geom_smooth(method="lm")
  ggplot(subset(stcrop_and_anpp_rr,pchange_C4P<4), aes(x=pchange_stcrop, y=pchange_C4P)) +
    geom_point() + geom_smooth(method="lm",se=F)

  ggplot(stcrop_rr, aes(Drought, pchange_stcrop, col=Grazing)) +
    geom_point() + geom_smooth(method="lm",se=F) + facet_grid(Year~Grazing, scales="free")
  
  c3p_stcrop_regression <- lm(pchange_C3P ~ pchange_stcrop, data=stcrop_and_anpp_rr)
  anova(c3p_stcrop_regression)
  summary(c3p_stcrop_regression)
}




