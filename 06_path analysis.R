###
### Path analysis 
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### Date created: May 14, 2026

### Set up workspace
library(tidyverse)
library(piecewiseSEM)
rm(list=ls())
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\data_sets\\")

###
### Read in data frames
npp_master <- read.csv("ANPP BNPP total NPP_plot level_cleaned and processed from R2025-07-09.csv")
smoist_master <- read.csv("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\soil_moisture\\smoist_bothSites_allYears.csv")
fxn_npp_master <- read.csv("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\data_sets\\ANPP by fxn group_2018-2023_2025-05-20.csv")

### ANPP and BNPP 2nd year of drought
{
### Calculate response ratios

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
  rename(ANPPdrtRR=lnrr_npp)

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

