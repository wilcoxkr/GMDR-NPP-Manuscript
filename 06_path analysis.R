###
### Path analysis 
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### Date created: May 14, 2026 -- last updated: Nov 9, 2026

### Set up workspace
library(tidyverse)
library(piecewiseSEM)
library(GGally)
library(nlme)
rm(list=ls())
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\data_sets\\")

###
### Read in data frames
npp_master <- read.csv("ANPP BNPP total NPP_plot level_cleaned and processed from R2025-07-09.csv")
smoist_master <- read.csv("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\soil_moisture\\smoist_bothSites_allYears.csv")
fxn_npp_master <- read.csv("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\data_sets\\ANPP by fxn group_2018-2023_2025-05-20.csv")

###
### Path analysis using raw data (not response ratios)
###
{
  ### ANPP and BNPP 2nd year of drought, and ANPP during recovery
  {
    
    anpp_drt2020 <- npp_master %>%
      filter(Year==2020 & npp_type=="ANPP") %>%
      mutate(Grazing_num=case_when(Grazing=="MLLMM" ~ 1,
                                   Grazing=="MMMMM" ~ 2,
                                   Grazing=="HHMMM" ~ 3)) %>%
      dplyr::select(Site, Plot, Drought, Grazing_num, npp_gm2) %>%
      rename(ANPPdrt=npp_gm2)
    
    bnpp_drt2020 <- npp_master %>%
      filter(Year==2020 & npp_type=="BNPP") %>%
      dplyr::select(Site, Plot, npp_gm2) %>%
      rename(BNPPdrt=npp_gm2)
    
    anpp_recovery <- npp_master %>%
      filter(Year%in%2021:2023 & npp_type=="ANPP") %>%
      dplyr::select(Site, Year, Plot, npp_gm2) %>%
      rename(ANPPrec=npp_gm2)
    
  }
  
  ###
  ### Soil moisture
  ###
  {
    smoist_april <- smoist_master %>%
      filter(Year %in% 2021:2023) %>%
      group_by(Year, Site, Block, Paddock, Plot, Drought) %>%
      summarize(smoist=mean(Soil.Moisture, na.rm=T)) %>%
      ungroup() %>%
      dplyr::select(Site, Year, Plot, smoist)
    
  }
  
  ###
  ### Comm change (Euclidean distance of plot to control)
  ###
  {
    fxn_npp_rel <- fxn_npp_master %>%
      group_by(Year, Site, Block, Paddock, Drought, Plot) %>%
      mutate(rel_abun = ANPP_gm2/sum(ANPP_gm2)) %>%
      ungroup()
    
    annual_relabun <- fxn_npp_rel %>%
      filter(Year %in% 2021:2023 & fxn_type=="AnnualGrass_gm2") %>%
      dplyr::select(Site, Year, Plot, rel_abun) %>%
      rename(AnnualRelabun=rel_abun)

    annual_relabun2022 <- fxn_npp_rel %>%
      filter(Year == 2022 & fxn_type=="AnnualGrass_gm2") %>%
      dplyr::select(Site, Plot, rel_abun, ANPP_gm2) %>%
      rename(AnnualRelabun2022=rel_abun,
             AnnualBiomass2022=ANPP_gm2)
    C3P_relabun2020 <- fxn_npp_rel %>%
      filter(Year == 2020 & fxn_type=="C3P_gm2") %>%
      dplyr::select(Site, Plot, rel_abun, ANPP_gm2) %>%
      rename(C3PRelabun2020=rel_abun,
             C3PBiomass2020=ANPP_gm2)
    

    
  }
  
  
  ###
  ### Put together full data frame for running path analysis
  ###
  {
    rawData_df <- anpp_drt2020 %>%
      full_join(bnpp_drt2020, by=c("Site","Plot")) %>%
      full_join(anpp_recovery, by=c("Site","Plot")) %>%
      full_join(smoist_april, by=c("Year","Site","Plot")) %>%
      full_join(annual_relabun, by=c("Year","Site","Plot")) %>%
      full_join(annual_relabun2022, by=c("Site","Plot")) %>%
      full_join(C3P_relabun2020, by=c("Site","Plot"))
      
  }
  
  ###
  ### Run path analysis
  ###
  {
    ##
    ## Check for colinearities across variables without arrows
    ggpairs(rawData_df %>% 
              filter(Year==2022) %>% 
              filter(ANPPdrt<400) %>%
              dplyr::select(ANPPrec, ANPPdrt, BNPPdrt, AnnualRelabun, smoist))
    
    ##
    ## Model structure
    
    drt_model <- '
  anppRR_recovery ~ lnrr_smoist + ANPPdrtRR + BNPPdrtRR + fxn_dist
  lnrr_smoist ~ Drought
  ANPPdrtRR ~ Drought
  BNPPdrtRR ~ Drought
  fxn_dist ~ Drought
 '
  }
  
  ### some odd stuff with the path analysis -- specifically strong but variable positive relationship between annual grasses and anpp during recovery.. perhaps this is because in a certain year where drought legacies are present, annual grasses are not fully compensating... ok, what if we tried perennial rel abun
  data4sem <- filter(rawData_df, Site=="FK" & Year==2021) %>% drop_na()
  data4sem_bothSites <- filter(rawData_df, Year==2023) %>% drop_na()
  # just FK
  drt_model_fk_2021_fit <- psem(
    lm(ANPPrec ~ smoist + ANPPdrt + BNPPdrt + AnnualBiomass2022 + Drought, data= data4sem),
    lm(smoist ~ Drought + ANPPdrt, data= data4sem),
    lm(ANPPdrt ~ Drought, data= data4sem),
    lm(BNPPdrt ~ Drought, data= data4sem),
    lm(AnnualBiomass2022 ~ Drought + ANPPdrt, data= data4sem),
    data= data4sem) # Best fit of data to model
  summary(drt_model_fk_2021_fit, standardized=TRUE)
  # Both sites
  drt_model_bothSites_2023_fit <- psem(
    lm(ANPPrec ~ smoist + ANPPdrt + BNPPdrt + AnnualBiomass2022 + Drought, data= data4sem_bothSites),
    lm(smoist ~ Drought + ANPPdrt, data= data4sem_bothSites),
    lm(ANPPdrt ~ Drought, data= data4sem_bothSites),
    lm(BNPPdrt ~ Drought, data= data4sem_bothSites),
    lm(AnnualBiomass2022 ~ Drought + ANPPdrt, data= data4sem_bothSites),
    data= data4sem_bothSites) # Best fit of data to model
  summary(drt_model_bothSites_2023_fit, standardized=TRUE)
  
  
  ggplot(rawData_df, aes(x=AnnualBiomass2022, y=ANPPrec, col=as.factor(Drought))) +
    geom_point(size=2) + facet_grid(Site~Year)
  ggplot(rawData_df, aes(x=C3PBiomass2020, y=ANPPrec, col=as.factor(Drought))) +
    geom_point(size=2) + facet_grid(Site~Year)
  ggplot(rawData_df, aes(x=C3PRelabun2020, y=ANPPrec, col=as.factor(Drought))) +
    geom_point(size=2) + facet_grid(Site~Year)
  
}


###
### Doing everything with response ratios
###
{
  rm(list=ls()) # reset because of name duplicates from above
### ANPP and BNPP 2nd year of drought
{


npp_controls <- npp_master %>%
  filter(Drought==0) %>%
  group_by(Site, Year, Block, Paddock, npp_type) %>%
  summarize(ctrl_npp=mean(npp_gm2, na.rm=T)
  )

npp_rr <- npp_master %>%
  filter(Drought != 0) %>%
  full_join(npp_controls, by=c("Year","Site","Block","Paddock", "npp_type")) %>%
  mutate(lnrr_npp=log((npp_gm2+0.01)/(ctrl_npp+0.01)), pchange_npp=((npp_gm2+0.01)-(ctrl_npp+0.01))/(ctrl_npp+0.01)) %>%
  dplyr::select(Year:npp_type, ctrl_npp, lnrr_npp, pchange_npp)

anpp_drt <- npp_rr %>%
  filter(Year==2020 & npp_type=="ANPP") %>%
  mutate(Grazing_num=case_when(Grazing=="MLLMM" ~ 1,
                        Grazing=="MMMMM" ~ 2,
                        Grazing=="HHMMM" ~ 3)) %>%
  dplyr::select(Site, Plot, Drought, Grazing_num, lnrr_npp) %>%
  rename(ANPPdrtRR=lnrr_npp)
  
bnpp_drt <- npp_rr %>%
  filter(Year==2020 & npp_type=="BNPP") %>%
  mutate(Grazing_num=case_when(Grazing=="MLLMM" ~ 1,
                               Grazing=="MMMMM" ~ 2,
                               Grazing=="HHMMM" ~ 3)) %>%
  dplyr::select(Site, Plot, Drought, Grazing_num, lnrr_npp) %>%
  rename(BNPPdrtRR=lnrr_npp)

anpp_recov <- npp_rr %>%
  filter(Year %in% 2021:2023 & npp_type=="ANPP") %>%
  dplyr::select(Year, Site, Plot, lnrr_npp) %>%
  rename(anppRR_recovery=lnrr_npp)

}

###
### Soil moisture
###
{
  smoist_april <- smoist_master %>%
    filter(Year %in% 2021:2023) %>%
    group_by(Year, Site, Block, Paddock, Plot, Drought) %>%
    summarize(smoist=mean(Soil.Moisture, na.rm=T)) %>%
    ungroup()
  
  smoist_controls <- smoist_april %>%
    filter(Drought==0) %>%
    group_by(Site, Year, Block, Paddock) %>%
    summarize(ctrl_smoist=mean(smoist, na.rm=T)
    )
  
  smoistRR_april <- smoist_april %>%
    filter(Drought != 0) %>%
    full_join(smoist_controls, by=c("Year","Site","Block","Paddock")) %>%
    mutate(lnrr_smoist=log((smoist)/(ctrl_smoist))) %>%
    dplyr::select(Year,Site, Plot, lnrr_smoist)
  
}

###
### Precipitation  
###
  {
  ## I am only looking at ambient precip during the recovery years. There are slight differences due to the window we used for rain year ppt
    ## So we might try to do this by treatment instead in the future... but I think we should start with just ambient precip for now
  ppt_rainyear <- read.csv("C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\GMDR\\data\\precipitation\\rain year precipitation by treatment_both sites.csv") %>%
    rename(Site=site, Plot=plot, Year=rain_year) %>%
    filter(Year %in% 2021:2023, rainfall_reduction==0) %>%
    dplyr::select(Year, Site, rain_year_ppt) %>%
    unique(.)
   
  # ppt_trt_both_sites <- read.csv("C:\\Users\\k_wilcox\\OneDrive - UNCG\\Current projects\\GMDR\\data\\precipitation\\daily ppt cleaned both sites 2026-05-12.csv") %>%
  #   rename(Site=site, Year=year, rainfall_reduction=Drought) %>%
  #   mutate(Rain_year=ifelse(month >= 8, Year+1, Year)) %>% ## Rain year is August 1 - July 31
  #   mutate(three_mo_flag = ifelse(month %in% 5:7,1,0),
  #          six_mo_flag = ifelse(month %in% 2:7,1,0),
  #          nine_mo_flag = ifelse(month %in% c(11:12,1:7),1,0))
  }
  
  
  
  
###
### Comm change (Euclidean distance of plot to control)
###
{
fxn_npp_rel <- fxn_npp_master %>%
  group_by(Year, Site, Block, Paddock, Drought, Plot) %>%
  mutate(rel_abun = ANPP_gm2/sum(ANPP_gm2)) %>%
  ungroup()

fxn_npp_wide <- fxn_npp_rel %>%
  dplyr::select(-slope_simple, -Grazing, -ANPP_gm2) %>%
  pivot_wider(names_from="fxn_type", values_from="rel_abun")

env_df <- fxn_npp_wide %>% dplyr::select(Year:Plot)
fxn_df <- fxn_npp_wide %>% dplyr::select(C4P_gm2:Forb_gm2)

pc_fxn <- prcomp(fxn_df, scale.=T)

fxn_scores <- bind_cols(env_df, as.data.frame(pc_fxn$x))


fxn_controls <- fxn_scores %>%
  filter(Drought == 0) %>%
  select(Year, Site, Block, Paddock, PC1, PC2) %>%
  rename(control_PC1 = PC1,
         control_PC2 = PC2) %>%
  group_by(Year, Site, Block, Paddock) %>%
  summarize(control_PC1=mean(control_PC1),
            control_PC2=mean(control_PC2)) %>%
  ungroup()

fxn_dist <- fxn_scores %>%
  filter(Drought != 0) %>%
  select(Year, Site, Block, Paddock, Plot, Drought, PC1, PC2) %>%
  full_join(fxn_controls, by=c("Year","Site","Block","Paddock")) %>%
  mutate(dist_to_control = sqrt((PC1 - control_PC1)^2 + (PC2 - control_PC2)^2)) %>%
  rename(fxn_dist=dist_to_control) %>%
  select(Year, Site, Plot, fxn_dist)

}

###
### Annual relative abundance 
###
  {
  fxn_npp_rel <- fxn_npp_master %>%
    group_by(Year, Site, Block, Paddock, Drought, Plot) %>%
    mutate(rel_abun = ANPP_gm2/sum(ANPP_gm2)) %>%
    ungroup()
  
  annual_relabun_noCtrl <- fxn_npp_rel %>%
    filter(Year %in% 2021:2023 & fxn_type=="AnnualGrass_gm2" & Drought > 0) %>%
    dplyr::select(Site, Year, Plot, rel_abun) %>%
    rename(AnnualRelabun=rel_abun)
  
  annual_relabun_2020 <- fxn_npp_rel %>%
    filter(Year %in% 2021:2023 & fxn_type=="AnnualGrass_gm2" & Drought > 0) %>%
    dplyr::select(Site, Year, Plot, rel_abun) %>%
    rename(AnnualRelabun=rel_abun)
  
  annual_relabun_mean <- fxn_npp_rel %>%
    filter(Year %in% 2021:2023 & fxn_type=="AnnualGrass_gm2" & Drought > 0) %>%
    dplyr::select(Site, Year, Plot, rel_abun) %>%
    group_by(Site, Plot) %>%
    summarize(AnnualRelabunMean=mean(rel_abun, na.rm=T))

  }
  

###
### Put together full data frame for running path analysis
###
{
full_df_responses <- anpp_drt %>%
  full_join(bnpp_drt, by=c("Site","Plot","Drought","Grazing_num")) %>%
  full_join(anpp_recov, by=c("Site","Plot")) %>%
  full_join(smoistRR_april, by=c("Year","Site","Plot")) %>%
  full_join(filter(fxn_dist,Year%in%2021:2023), by=c("Year","Site","Plot")) %>%
  full_join(annual_relabun_noCtrl, by=c("Year","Site","Plot")) %>%
  full_join(annual_relabun_mean, by=c("Site","Plot")) %>%
  full_join(ppt_rainyear, by=c("Year","Site"))

full_df_responses$Site <- factor(full_df_responses$Site)
full_df_responses$Plot <- factor(paste(full_df_responses$Site, full_df_responses$Plot, sep="_"))

}

###
### Run path analysis
###
{
  ##
  ## Check for colinearities across variables without arrows
  ggpairs(full_df_responses %>% filter(Site=="FK",Year==2023) %>% dplyr::select(anppRR_recovery, ANPPdrtRR, BNPPdrtRR, fxn_dist, lnrr_smoist))
  ggpairs(full_df_responses  %>% filter(Year==2021) %>% dplyr::select(anppRR_recovery, ANPPdrtRR, BNPPdrtRR, fxn_dist, lnrr_smoist, AnnualRelabun))
  ggpairs(full_df_responses  %>% filter(Year==2022) %>% dplyr::select(anppRR_recovery, ANPPdrtRR, BNPPdrtRR, fxn_dist, lnrr_smoist, AnnualRelabun))
  ggpairs(full_df_responses  %>% filter(Year==2023) %>% dplyr::select(anppRR_recovery, ANPPdrtRR, BNPPdrtRR, fxn_dist, lnrr_smoist, AnnualRelabun))
  ggpairs(full_df_responses  %>% dplyr::select(anppRR_recovery, ANPPdrtRR, BNPPdrtRR, lnrr_smoist, AnnualRelabunMean, rain_year_ppt))

  ##
  ## Model structure
  
  drt_model <- '
  anppRR_recovery ~ lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean + rain_year_ppt
  lnrr_smoist ~ Drought
  ANPPdrtRR ~ Drought
  BNPPdrtRR ~ Drought
  AnnualRelabunMean ~ Drought
  rain_year_ppt ~ Drought
 '
}

drt_model_rep <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean + rain_year_ppt,
      random = ~ 1 | Site/Plot, 
      data = full_df_responses),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Site, data=full_df_responses),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Site, data=full_df_responses),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Site, data=full_df_responses),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Site, data=full_df_responses),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  
  data=full_df_responses
)
dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")

summary(drt_model_rep)
dSep(drt_model_rep, direction=dirs) 

drt_model_fk <- drt_model_rep <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean + rain_year_ppt,
      random = ~ 1 | Plot, 
      data = filter(full_df_responses, Site=="FK")),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Plot, data=filter(full_df_responses, Site=="FK")),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Plot, data=filter(full_df_responses, Site=="FK")),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Plot, data=filter(full_df_responses, Site=="FK")),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Plot, data=filter(full_df_responses, Site=="FK")),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  
  data=filter(full_df_responses, Site=="FK")
)
dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")

summary(drt_model_fk)
dSep(drt_model_fk, direction=dirs) 


### splitting by year
drt_model_2021 <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
      random = ~ 1 | Site,
      data = filter(full_df_responses, Year==2021)),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2021)),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2021)),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2021)),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Site, data = filter(full_df_responses, Year==2021)),
  
#  lnrr_smoist %~~% AnnualRelabunMean,
  lnrr_smoist %~~% ANPPdrtRR,
  
  data = filter(full_df_responses, Year==2021)
)
summary(drt_model_2021)

dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")
dSep(drt_model_2021, direction=dirs) 


drt_model_2022 <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
      random = ~ 1 | Site,
      data = filter(full_df_responses, Year==2022)),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  lnrr_smoist %~~% ANPPdrtRR,
  
  data = filter(full_df_responses, Year==2022)
)
summary(drt_model_2022)

dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")
dSep(drt_model_2022, direction=dirs) 

drt_model_2023 <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
      random = ~ 1 | Site,
      data = filter(full_df_responses, Year==2023)),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  lnrr_smoist %~~% ANPPdrtRR,
  
  data = filter(full_df_responses, Year==2023)
)
summary(drt_model_2023)

dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")
dSep(drt_model_2023, direction=dirs) 


drt_model_2022 <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
      random = ~ 1 | Site/Plot, 
      data = filter(full_df_responses, Year==2022)),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Site, data = filter(full_df_responses, Year==2022)),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  
  data = filter(full_df_responses, Year==2022)
)
dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")

summary(drt_model_2022)
dSep(drt_model_2022, direction=dirs) 

drt_model_2023 <- psem(
  lme(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
      random = ~ 1 | Site/Plot, 
      data = filter(full_df_responses, Year==2023)),
  lme(lnrr_smoist ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  lme(ANPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  lme(BNPPdrtRR ~ Drought, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  lme(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, random = ~ 1 | Site, data = filter(full_df_responses, Year==2023)),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  
  data = filter(full_df_responses, Year==2023)
)
dirs <- c("lnrr_smoist <- ANPPdrtRR", "lnrr_smoist <- BNPPdrtRR")

summary(drt_model_2023)
dSep(drt_model_2023, direction=dirs) 


### splitting by year and site
drt_model_2021_FK <- psem(
  lm(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
      data = filter(full_df_responses, Year==2021 & Site=="FK")),
  lm(lnrr_smoist ~ Drought, data = filter(full_df_responses, Year==2021 & Site=="FK")),
  lm(ANPPdrtRR ~ Drought, data = filter(full_df_responses, Year==2021 & Site=="FK")),
  lm(BNPPdrtRR ~ Drought, data = filter(full_df_responses, Year==2021 & Site=="FK")),
  lm(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, data = filter(full_df_responses, Year==2021 & Site=="FK")),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  lnrr_smoist %~~% ANPPdrtRR,
  
  data = filter(full_df_responses, Year==2021 & Site=="FK")
)
summary(drt_model_2021_FK)

drt_model_2022_FK <- psem(
  lm(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
     data = filter(full_df_responses, Year==2022 & Site=="FK")),
  lm(lnrr_smoist ~ Drought, data = filter(full_df_responses, Year==2022 & Site=="FK")),
  lm(ANPPdrtRR ~ Drought, data = filter(full_df_responses, Year==2022 & Site=="FK")),
  lm(BNPPdrtRR ~ Drought, data = filter(full_df_responses, Year==2022 & Site=="FK")),
  lm(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, data = filter(full_df_responses, Year==2022 & Site=="FK")),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  lnrr_smoist %~~% ANPPdrtRR,
  
  data = filter(full_df_responses, Year==2022 & Site=="FK")
)
summary(drt_model_2022_FK)

drt_model_2023_FK <- psem(
  lm(anppRR_recovery ~ Drought + lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean,
     data = filter(full_df_responses, Year==2023 & Site=="FK")),
  lm(lnrr_smoist ~ Drought, data = filter(full_df_responses, Year==2023 & Site=="FK")),
  lm(ANPPdrtRR ~ Drought, data = filter(full_df_responses, Year==2023 & Site=="FK")),
  lm(BNPPdrtRR ~ Drought, data = filter(full_df_responses, Year==2023 & Site=="FK")),
  lm(AnnualRelabunMean ~ Drought + ANPPdrtRR + BNPPdrtRR, data = filter(full_df_responses, Year==2023 & Site=="FK")),
  
  lnrr_smoist %~~% AnnualRelabunMean,
  lnrr_smoist %~~% ANPPdrtRR,
  
  data = filter(full_df_responses, Year==2023 & Site=="FK")
)
summary(drt_model_2023_FK)




drt_model_fit <- psem(
  lm(anppRR_recovery ~ lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean + rain_year_ppt, data=full_df_responses),
  lm(lnrr_smoist ~ Drought, data=full_df_responses),
  lm(ANPPdrtRR ~ Drought, data=full_df_responses),
  lm(BNPPdrtRR ~ Drought, data=full_df_responses),
  lm(AnnualRelabunMean ~ Drought, data=full_df_responses),
  lm(rain_year_ppt ~ Drought, data=full_df_responses),
  data=full_df_responses) # Best fit of data to model
summary(drt_model_fit, standardized=TRUE)

drt_model_fk_2021_fit <- psem(
  lm(anppRR_recovery ~ lnrr_smoist + ANPPdrtRR + BNPPdrtRR + fxn_dist, data=filter(full_df_responses, Site=="FK" & Year==2021)),
  lm(lnrr_smoist ~ Drought, data=filter(full_df_responses, Site=="FK" & Year==2021)),
  lm(ANPPdrtRR ~ Drought, data=filter(full_df_responses, Site=="FK" & Year==2021)),
  lm(BNPPdrtRR ~ Drought, data=filter(full_df_responses, Site=="FK" & Year==2021)),
  data=filter(full_df_responses, Site=="FK" & Year==2021)) # Best fit of data to model
summary(drt_model_fk_2021_fit, standardized=TRUE)

drt_model_2023_fit <- psem(
  lm(anppRR_recovery ~ lnrr_smoist + ANPPdrtRR + BNPPdrtRR + AnnualRelabunMean, data=filter(full_df_responses, Year==2023)),
  lm(lnrr_smoist ~ Drought, data=filter(full_df_responses, Year==2023)),
  lm(ANPPdrtRR ~ Drought, data=filter(full_df_responses, Year==2023)),
  lm(BNPPdrtRR ~ Drought, data=filter(full_df_responses, Year==2023)),
  lm(AnnualRelabunMean ~ Drought, data=filter(full_df_responses, Year==2023)),
  data=filter(full_df_responses, Year==2023)) # Best fit of data to model
summary(drt_model_2023_fit, standardized=TRUE)

}