###
### Precipitation-npp regressions - analyses and plotting
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### Created: June 19, 2024; last updated: May 23, 2025

### Set up workspace
library(tidyverse)
library(nlme)
library(car)
library(lubridate)
library(performance)
library(ggthemes)
#library(broom)
rm(list=ls()) # clean up

# Set working directory 
# setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\GMDR\\data\\anpp\\") # Kevin's laptop
# setwd("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer UNCG
setwd("C:\\Users\\wilco\\OneDrive - UNCG\\Current projects\\GMDR\\data\\anpp\\") # Kevins laptop computer UNCG

# Set writing directory
write_dir <- "C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\"

### Set graphing parameters
# source("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Git projects\\Grazing-Management-for-Drought-Resilience\\graph_format.R")
# source("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Git projects\\Grazing-Management-for-Drought-Resilience\\graph_format.R")
source("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\graph_format.R")
source("C:\\Users\\wilco\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\graph_format.R")

### Create Standard Error function 
SE_function<-function(x,na.rm=na.rm){
  SD = sd(x,na.rm=TRUE)
  denom = sqrt(length(x[!is.na(x)]))
  SE=SD/denom
  return(SE)
}


### Read in anova_t3 function
# source("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Git projects\\Grazing-Management-for-Drought-Resilience\\lme_anova_type3_function.R")
# source("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Git projects\\Grazing-Management-for-Drought-Resilience\\lme_anova_type3_function.R")
source("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\lme_anova_type3_function.R")
source("C:\\Users\\wilco\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\lme_anova_type3_function.R")

###
### Read in npp and precipitation data and merge together
### 
{
  npp_df <- read.csv("ANPP BNPP total NPP_plot level_cleaned and processed from R.csv")
  
  ppt_rainyear <- read.csv("..\\precipitation\\rain year precipitation by treatment_both sites.csv") %>%
    rename(Site=site, Plot=plot, Year=rain_year) %>%
    filter(Year %in% 2018:2023)
  
  npp_ppt_df <- npp_df %>% ### Note: there are some plots in 2018 that don't match up due to the accidental shifting of the fences -- these are creating extra rows when merging. We can remove these extra rows 
    full_join(ppt_rainyear, by=c("Site","Plot", "Year")) %>%
    mutate(trt_on = ifelse(Year%in%2019:2020, "yes","no"))
  
  npp_ppt_means <- npp_ppt_df %>%
    group_by(Year, Site, npp_type, Drought) %>%
    summarize(npp_mean = mean(npp_gm2, na.rm=T),
              npp_se = SE_function(npp_gm2),
              rain_year_ppt_mean = mean(rain_year_ppt, na.rm=T))
## take a quick look at plots
  # ggplot(filter(npp_ppt_df,npp_type=="ANPP"), aes(x=rain_year_ppt, y=npp_gm2, col=trt_on)) +
  #   geom_point() +
  #   facet_wrap(~Site)
  # 
  # ggplot(filter(npp_ppt_means,npp_type=="ANPP"), aes(x=rain_year_ppt_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, col=factor(Year), shape=factor(Drought))) +
  #   geom_errorbar(width=5) +
  #   geom_point() +
  #   geom_smooth(inherit.aes=F, aes(x=rain_year_ppt_mean, y=npp_mean), method="lm", se=F, col="black") +
  #   facet_wrap(~Site)
  # 
  # ggplot(filter(npp_ppt_means,npp_type=="BNPP"), aes(x=rain_year_ppt_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, col=factor(Year), shape=factor(Drought))) +
  #   geom_errorbar(width=5) +
  #   geom_point() +
  #   geom_smooth(inherit.aes=F, aes(x=rain_year_ppt_mean, y=npp_mean), method="lm", se=F, col="black") +
  #   facet_wrap(~Site)
  # 
  # ggplot(filter(npp_ppt_means,npp_type=="totNPP"), aes(x=rain_year_ppt_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, col=factor(Year), shape=factor(Drought))) +
  #   geom_errorbar(width=5) +
  #   geom_point() +
  #   geom_smooth(inherit.aes=F, aes(x=rain_year_ppt_mean, y=npp_mean), method="lm", se=F, col="black") +
  #   facet_wrap(~Site)
}

###
### Read in npp and two year precipitation window data and merge together
### 
{
  npp_df <- read.csv("ANPP BNPP total NPP_plot level_cleaned and processed from R.csv")
  
  ppt_2rainyear <- read.csv("..\\precipitation\\rain year precipitation over two years by treatment_both sites.csv") %>%
    rename(Site=site, Plot=plot, Year=rain_year) %>%
    filter(Year %in% 2019:2023)
  
  npp_ppt_2yr_df <- npp_df %>% ### Note: there are some plots in 2018 that don't match up due to the accidental shifting of the fences -- these are creating extra rows when merging. We can remove these extra rows 
    full_join(ppt_2rainyear, by=c("Site","Plot", "Year")) %>%
    mutate(trt_on = ifelse(Year%in%2019:2020, "yes","no"))
  
  npp_ppt_2yr_means <- npp_ppt_2yr_df %>%
    group_by(Year, Site, npp_type, Drought) %>%
    summarize(npp_mean = mean(npp_gm2, na.rm=T),
              npp_se = SE_function(npp_gm2),
              ppt_2Yr_mean = mean(ppt_2Yr, na.rm=T),
              ppt_prev1 = mean(ppt_prev1)) %>%
    ungroup()
  
}


###
### Testing linear vs non-linear slopes for 1 and 2 year models (STILL NEEDS TO BE CLEANED UP A BIT)
###
{
  ### Current rain year models
  ## Testing transformations
  performance::r2(lm(npp_gm2~rain_year_ppt*Site, data=filter(npp_ppt_df,npp_type=="ANPP")))
  performance::r2(lm(log(npp_gm2)~rain_year_ppt*Site, data=filter(npp_ppt_df,npp_type=="ANPP")))
  performance::r2(lm(sqrt(npp_gm2)~rain_year_ppt*Site, data=filter(npp_ppt_df,npp_type=="ANPP")))
  
  performance::r2(lm(npp_gm2~rain_year_ppt, data=filter(npp_ppt_df,npp_type=="ANPP"&Site=="FK")))
  performance::r2(lm(log(npp_gm2)~rain_year_ppt, data=filter(npp_ppt_df,npp_type=="ANPP"&Site=="FK")))
  performance::r2(lm(sqrt(npp_gm2)~rain_year_ppt, data=filter(npp_ppt_df,npp_type=="ANPP"&Site=="FK")))
  
  current_ppt_anpp_lm <- lm(sqrt(npp_gm2)~rain_year_ppt*Site, data=filter(npp_ppt_df,npp_type=="ANPP"))
  Anova(current_ppt_anpp_lm, type=3)
  current_ppt_anpp_anova <- Anova(current_ppt_anpp_lm, type=3)
  # Significant interaction, split by site
  
  plot(current_ppt_anpp_lm)
  hist(sqrt(filter(npp_ppt_df,npp_type=="ANPP")$npp_gm2))
  
  current_ppt_anova_df <- data.frame(effect=row.names(current_ppt_anpp_anova), current_ppt_anpp_anova, ppt_type="rain_year")
  current_ppt_diff_emtrends <- data.frame(test(emtrends(current_ppt_anpp_lm, "Site", var="rain_year_ppt")), ppt_type="rain_year")
  
  # For trendlines in Fig 3
  fk_current_ppt_lm <- lm(sqrt(npp_gm2)~rain_year_ppt, data=filter(npp_ppt_df, npp_type=="ANPP"& Site=="FK"))
  performance::r2(fk_current_ppt_lm)
  
  fk_current_ppt_xhat <- with(filter(npp_ppt_df, npp_type=="ANPP"& Site=="FK"), seq(min(rain_year_ppt), max(rain_year_ppt), by=1))
  fk_current_ppt_df <- data.frame(rain_year_ppt=fk_current_ppt_xhat)
  fk_current_ppt_yhat <- predict.lm(fk_current_ppt_lm, newdata = fk_current_ppt_df)
  plot(fk_current_ppt_xhat, fk_current_ppt_yhat^2)
  
  fk_current_ppt_trendlines <- data.frame(
    xhat=fk_current_ppt_xhat,
    yhat=fk_current_ppt_yhat,
    Site="FK",
    ppt_type="rain_year"
  )
  
  tb_current_ppt_lm <- lm(sqrt(npp_gm2)~rain_year_ppt, data=filter(npp_ppt_df, npp_type=="ANPP"& Site=="TB"))
  summary(tb_current_ppt_lm)
  performance::r2(tb_current_ppt_lm)
  
  tb_current_ppt_xhat <- with(filter(npp_ppt_df, npp_type=="ANPP"& Site=="TB"), seq(min(rain_year_ppt), max(rain_year_ppt), by=1))
  tb_current_ppt_df <- data.frame(rain_year_ppt=tb_current_ppt_xhat)
  tb_current_ppt_yhat <- predict.lm(tb_current_ppt_lm, newdata = tb_current_ppt_df)
  plot(tb_current_ppt_xhat, tb_current_ppt_yhat^2)
  
  tb_current_ppt_trendlines <- data.frame(
    xhat=tb_current_ppt_xhat,
    yhat=tb_current_ppt_yhat,
    Site="TB",
    ppt_type="rain_year"
  )
  
  current_ppt_trendlines_both <- tb_current_ppt_trendlines %>%
    bind_rows(fk_current_ppt_trendlines)
  
  ### 2 year precip models
  ## Testing transformations - LOG MODELS HAVE HIGHEST R2
  performance::r2(lm(npp_gm2~ppt_2Yr*Site, data=filter(npp_ppt_2yr_df,npp_type=="ANPP")))
  performance::r2(lm(log(npp_gm2)~ppt_2Yr*Site, data=filter(npp_ppt_2yr_df,npp_type=="ANPP")))
  performance::r2(lm(sqrt(npp_gm2)~ppt_2Yr*Site, data=filter(npp_ppt_2yr_df,npp_type=="ANPP")))
  
  performance::r2(lm(npp_gm2~ppt_2Yr, data=filter(npp_ppt_2yr_df,npp_type=="ANPP"&Site=="FK")))
  performance::r2(lm(log(npp_gm2)~ppt_2Yr, data=filter(npp_ppt_2yr_df,npp_type=="ANPP"&Site=="FK")))
  performance::r2(lm(sqrt(npp_gm2)~ppt_2Yr, data=filter(npp_ppt_2yr_df,npp_type=="ANPP"&Site=="FK")))
  
  ppt_2yr_anpp_lm <- lm(log(npp_gm2)~ppt_2Yr*Site, data=filter(npp_ppt_2yr_df,npp_type=="ANPP"))
  Anova(ppt_2yr_anpp_lm, type=3)
  ppt_2yr_anpp_anova <- Anova(ppt_2yr_anpp_lm, type=3)
  # Significant interaction, split by site
  
  plot(ppt_2yr_anpp_lm)
  hist(log(filter(npp_ppt_2yr_df,npp_type=="ANPP")$npp_gm2))
  
  ppt_2yr_anova_df <- data.frame(effect=row.names(ppt_2yr_anpp_anova), ppt_2yr_anpp_anova, ppt_type="2yr")
  ppt_2yr_diff_emtrends <- data.frame(test(emtrends(ppt_2yr_anpp_lm, "Site", var="ppt_2Yr")), ppt_type="2yr")
  
  # For trendlines and r2 in Fig 3
  fk_ppt_2yr_lm <- lm(log(npp_gm2)~ppt_2Yr, data=filter(npp_ppt_2yr_df, npp_type=="ANPP"& Site=="FK"))
  performance::r2(fk_ppt_2yr_lm)
  
  fk_ppt_2yr_xhat <- with(filter(npp_ppt_2yr_df, npp_type=="ANPP"& Site=="FK" & Year %in% 2019:2023), seq(min(ppt_2Yr), max(ppt_2Yr), by=1))
  fk_ppt_2yr_df <- data.frame(ppt_2Yr=fk_ppt_2yr_xhat)
  fk_ppt_2yr_yhat <- predict.lm(fk_ppt_2yr_lm, newdata = fk_ppt_2yr_df)
  plot(fk_ppt_2yr_xhat, exp(fk_ppt_2yr_yhat))
  
  fk_ppt_2yr_trendlines <- data.frame(
    xhat=fk_ppt_2yr_xhat,
    yhat=fk_ppt_2yr_yhat,
    Site="FK",
    ppt_type="2yr"
  )
  
  tb_ppt_2yr_lm <- lm(log(npp_gm2)~ppt_2Yr, data=filter(npp_ppt_2yr_df, npp_type=="ANPP"& Site=="TB"))
  summary(tb_ppt_2yr_lm)
  performance::r2(tb_ppt_2yr_lm)
  
  tb_ppt_2yr_xhat <- with(filter(npp_ppt_2yr_df, npp_type=="ANPP"& Site=="TB"&Year%in%2019:2023), seq(min(ppt_2Yr), max(ppt_2Yr), by=1))
  tb_ppt_2yr_df <- data.frame(ppt_2Yr=tb_ppt_2yr_xhat)
  tb_ppt_2yr_yhat <- predict.lm(tb_ppt_2yr_lm, newdata = tb_ppt_2yr_df)
  plot(tb_ppt_2yr_xhat, exp(tb_ppt_2yr_yhat))
  
  tb_ppt_2yr_trendlines <- data.frame(
    xhat=tb_ppt_2yr_xhat,
    yhat=tb_ppt_2yr_yhat,
    Site="TB",
    ppt_type="2yr"
  )
  
  ppt_2yr_trendlines_both <- tb_ppt_2yr_trendlines %>%
    bind_rows(fk_ppt_2yr_trendlines)
  
}

###
### Plotting 1 year and 2 year vs ANPP and BNPP by site
###
{
  # current year
  ppt1yrPlot <- ggplot(filter(npp_ppt_means,npp_type=="ANPP"), aes(x=rain_year_ppt_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, fill=factor(Drought), shape=factor(Year))) +
    geom_errorbar(width=15, col="black") +
    geom_point(col="black", size=3) +
    geom_line(inherit.aes=F, data=current_ppt_trendlines_both, aes(x=xhat, y=yhat^2), col="black") +
    scale_fill_manual(values=droughtColor) +
    scale_shape_manual(values=c(13,21:25)) +
    facet_wrap(~Site, scales="free_y") +
    theme_few()
  
  pdf(paste0(write_dir, "figures//ppt-anpp regression currenty rainYear",Sys.Date(),".pdf"), width=7, height=3, useDingbats = F)
  print(ppt1yrPlot)
  dev.off()
  
  
  # 2 year window
  ppt2yrPlot <- ggplot(filter(npp_ppt_2yr_means,npp_type=="ANPP"), aes(x=ppt_2Yr_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, fill=factor(Drought), shape=factor(Year))) +
    geom_errorbar(width=15) +
    geom_point(size=3) +
    geom_line(inherit.aes=F, data=ppt_2yr_trendlines_both, aes(x=xhat, y=exp(yhat)), col="black") +
    scale_fill_manual(values=droughtColor) +
    scale_shape_manual(values=c(13,21:25)) +
    facet_wrap(~Site, scales="free_y") +
    theme_few()
  
  pdf(paste0(write_dir, "figures//ppt-anpp regression last 2 rainYears",Sys.Date(),".pdf"), width=7, height=3, useDingbats = F)
  print(ppt2yrPlot)
  dev.off()
  
}


###
### Plotting slopes varying by grazing treatment
###
{
  ## 
  npp_ppt_grz_means <- npp_ppt_df %>%
    group_by(Year, Site, npp_type, Grazing, Drought) %>%
    summarize(npp_mean = mean(npp_gm2, na.rm=T),
              npp_se = SE_function(npp_gm2),
              rain_year_ppt_mean = mean(rain_year_ppt, na.rm=T)) %>%
    ungroup()
  
  npp_ppt_grz_means$Grazing <- factor(npp_ppt_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM")) 
  
  ggplot(filter(npp_ppt_grz_means,npp_type=="ANPP"), aes(x=rain_year_ppt_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, fill=factor(Drought), shape=factor(Year))) +
    geom_errorbar(width=15) +
    geom_point(size=3) +
    geom_smooth(inherit.aes=F, aes(x=rain_year_ppt_mean, y=npp_mean), method="lm", se=F, col="black") +
    scale_fill_manual(values=droughtColor) +
    scale_shape_manual(values=c(13,21:25)) +
    facet_grid(Grazing~Site, scales="free_y") +
    theme_few()
  
  ggsave(paste0(write_dir, "figures//1 year ppt vs npp - split by grazing_",Sys.Date(),".png"), width=6,height=6, units="in")
  
  
  ## 2 year cumulative
  npp_ppt_2yr_grz_means <- npp_ppt_2yr_df %>%
    group_by(Year, Site, npp_type, Drought, Grazing) %>%
    summarize(npp_mean = mean(npp_gm2, na.rm=T),
              npp_se = SE_function(npp_gm2),
              ppt_2Yr_mean = mean(ppt_2Yr, na.rm=T),
              ppt_prev1 = mean(ppt_prev1)) %>%
    ungroup()
  
  npp_ppt_2yr_grz_means$Grazing <- factor(npp_ppt_2yr_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM")) 
  ggplot(filter(npp_ppt_2yr_grz_means,npp_type=="ANPP"), aes(x=ppt_2Yr_mean, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, fill=factor(Drought), shape=factor(Year))) +
    geom_errorbar(width=15) +
    geom_point(size=3) +
    geom_smooth(inherit.aes=F, aes(x=ppt_2Yr_mean, y=npp_mean), method="lm", se=F, col="black") +
    scale_fill_manual(values=droughtColor) +
    scale_shape_manual(values=c(13,21:25)) +
    facet_grid(Grazing~Site, scales="free_y") +
    theme_few()
  
  ggsave(paste0(write_dir, "figures//2yr ppt vs npp - split by grazing_",Sys.Date(),".png"), width=6,height=6, units="in")
  
  

}


###
### Models comparing ppt slopes with ANPP by grazing treatment
###
{
  ## Current year precipitation  
  fk_cuyr_ppt_anpp_grz_model <- lme(sqrt(npp_gm2) ~ rain_year_ppt*Grazing
                                 , data=filter(npp_ppt_2yr_df, Site=="FK" & npp_type=="ANPP")
                                 , random = ~1 |Block/Paddock
                                 , na.action = na.omit)
  anova.lme(fk_cuyr_ppt_anpp_grz_model, type="marginal")  

  tb_cuyr_ppt_anpp_grz_model <- lme(sqrt(npp_gm2) ~ rain_year_ppt*Grazing
                                 , data=filter(npp_ppt_2yr_df, Site=="TB" & npp_type=="ANPP")
                                 , random = ~1 |Block/Paddock
                                 , na.action = na.omit)
  anova.lme(tb_cuyr_ppt_anpp_grz_model, type="marginal")  
  
  ## Current year precipitation  
  fk_2yr_ppt_anpp_grz_model <- lme(log(npp_gm2) ~ ppt_2Yr*Grazing
                                    , data=filter(npp_ppt_2yr_df, Site=="FK" & npp_type=="ANPP")
                                    , random = ~1 |Block/Paddock
                                    , na.action = na.omit)
  anova.lme(fk_2yr_ppt_anpp_grz_model, type="marginal")  
  
  tb_2yr_ppt_anpp_grz_model <- lme(log(npp_gm2) ~ ppt_2Yr*Grazing
                                    , data=filter(npp_ppt_2yr_df, Site=="TB" & npp_type=="ANPP")
                                    , random = ~1 |Block/Paddock
                                    , na.action = na.omit)
  anova.lme(tb_2yr_ppt_anpp_grz_model, type="marginal")  
  
}

### looking at models of grazing vs drought
npp_ppt_df
tb_anpp_grz_model <- lme(npp_gm2 ~ Drought*Grazing
                         , data=filter(npp_ppt_df, Year==2021 & Site=="TB" & npp_type=="ANPP")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
anova(tb_anpp_grz_model, type="marginal")
summary(tb_anpp_grz_model)

### Read in TB Data and Key, clean up data frames -- need to run L69-139 in anpp_gmdr_Apr2024.R script

tb_anpp_2021_grz_model <- lme(ANPP_gm2 ~ Drought*Grazing
                              , data=filter(anpp_tb_plot_means, Year==2021)
                              , random = ~1 |Block/Paddock
                              , na.action = na.omit)
anova(tb_anpp_2021_grz_model, type="marginal")
summary(tb_anpp_2021_grz_model)

tb_C3P_2021_grz_model <- lme(C3P_gm2 ~ Drought*Grazing
                             , data=filter(anpp_tb_plot_means, Year==2021)
                             , random = ~1 |Block/Paddock
                             , na.action = na.omit)
anova(tb_C3P_2021_grz_model, type="marginal")
summary(tb_C3P_2021_grz_model)

tb_C4P_2021_grz_model <- lme(C4P_gm2 ~ Drought*Grazing
                             , data=filter(anpp_tb_plot_means, Year==2021)
                             , random = ~1 |Block/Paddock
                             , na.action = na.omit)
anova(tb_C4P_2021_grz_model, type="marginal")
summary(tb_C4P_2021_grz_model)

tb_Forb_2021_grz_model <- lme(Forb_gm2 ~ Drought*Grazing
                              , data=filter(anpp_tb_plot_means, Year==2021)
                              , random = ~1 |Block/Paddock
                              , na.action = na.omit)
anova(tb_Forb_2021_grz_model, type="marginal")
summary(tb_Forb_2021_grz_model)

tb_AnnualGrass_2021_grz_model <- lme(AnnualGrass_gm2 ~ Drought*Grazing
                                     , data=filter(anpp_tb_plot_means, Year==2021)
                                     , random = ~1 |Block/Paddock
                                     , na.action = na.omit)
anova(tb_AnnualGrass_2021_grz_model, type="marginal")
summary(tb_AnnualGrass_2021_grz_model)


