### Cleaning plotting and analysis -- ANPP data
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### Created January 25, 2022; Last updated May 5, 2025

### Set up workspace
library(tidyverse)
library(nlme)
library(car)
library(lubridate)
library(performance)
library(ggthemes)
library(broom)
library(sjstats)
library(emmeans)
rm(list=ls()) # clean up

# Set working directory 
# setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\GMDR\\data\\anpp\\") # Kevin's laptop
# setwd("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer UNCG
setwd("C:\\Users\\wilco\\OneDrive - UNCG\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer UNCG

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


### Read in Data and Key, clean up data frames
Plot_Key <- read.csv("..\\Plot_Treatment_Key_020424.csv") %>%
  rename(PlotID=plot, Site=site) %>%
  dplyr::select(-drought) %>%
  rename(Drought=rainfall_reduction, Block=block, Paddock=paddock, Grazing=grazing_category)

###
### Data cleaning ANPP
###
{
  ### ANPP
  ###
  
  ###
  ### Thunder Basin
  ###
  
  #BEGIN_ANPP_SOURCE
  
anpp_tb_22 <- read.csv("GMDR_anpp_2018-2022_TB.csv") %>% ## Fixing block, paddock, drought and grazing labels in 2022 data
    filter(Year==2022) %>%
    rename(Clip_strip=Clip.strip, C3P_gm2=C3P.g.m2.,C4P_gm2=C4P.g.m2., Forb_gm2=Forb.g.m2.,
           ANPP_gm2=ANPP.g.m2., StandingDead_gm2=StandingDead.g.m2., C3P_gm2=C3P.g.m2., AnnualGrass_gm2=AnnualGrass.g.m2.,
           Vulpia_gm2=Vulpia.g.m2., Bromus_gm2=Bromus.g.m2.) %>%
    select(Year:Clip_strip, C3P_gm2:Bromus_gm2) %>%
    mutate(Drought=replace(Drought, Drought %in% c(1,2),0)) %>% ## Ambient plots were coded as 1 or 2, changed to 0 here
    select(-Drought, -Grazing, -Block, -Paddock) %>%
    left_join(Plot_Key, by=c("Site","PlotID")) %>%
    dplyr::select(Year, Site, Block, Paddock, Drought, Grazing, PlotID, Clip_strip, slope_simple,
                  C3P_gm2, C4P_gm2, Forb_gm2, ANPP_gm2, StandingDead_gm2, AnnualGrass_gm2, 
                  Vulpia_gm2,Bromus_gm2)
  
anpp_tb_18_21 <- read.csv("GMDR_anpp_2018-2022_TB.csv") %>%
  filter(Year!=2022) %>% ### there is a mistake in the raw data for 2022 where all plots are labeled as 75% drought and MMMMM so here I remove 2022, fix it separately (above) and put them back together.
  rename(Clip_strip=Clip.strip, C3P_gm2=C3P.g.m2.,C4P_gm2=C4P.g.m2., Forb_gm2=Forb.g.m2.,
         ANPP_gm2=ANPP.g.m2., StandingDead_gm2=StandingDead.g.m2., C3P_gm2=C3P.g.m2., AnnualGrass_gm2=AnnualGrass.g.m2.,
         Vulpia_gm2=Vulpia.g.m2., Bromus_gm2=Bromus.g.m2.) %>%
  select(Year:Clip_strip, C3P_gm2:Bromus_gm2) %>%
#  mutate(Drought=replace(Drought, Drought %in% c(1,2),0)) %>% ## Ambient plots were coded as 1 or 2, changed to 0 here
  select(-Drought, -Grazing, -Block, -Paddock) %>%
  left_join(Plot_Key, by=c("Site","PlotID")) %>%
  dplyr::select(Year, Site, Block, Paddock, Drought, Grazing, PlotID, Clip_strip, slope_simple,
                C3P_gm2, C4P_gm2, Forb_gm2, ANPP_gm2, StandingDead_gm2, AnnualGrass_gm2, 
                Vulpia_gm2,Bromus_gm2)

anpp_tb_18_22 <- anpp_tb_18_21 %>%
  bind_rows(anpp_tb_22)

### There are some mistakes in the raw data with plot identifiers, so I drop all plot identifiers from this dataset except plot ID and then merge together with the PlotID key to fix this
### However, plot numbers in 2018 are different so I need to clean 2018 separate then bind back with the rest of the data
anpp_tb_18 <- read.csv("GMDR_anpp_2018-2022_TB.csv") %>%
  filter(Year==2018) %>% 
  rename(Clip_strip=Clip.strip, C3P_gm2=C3P.g.m2.,C4P_gm2=C4P.g.m2., Forb_gm2=Forb.g.m2.,
         ANPP_gm2=ANPP.g.m2., StandingDead_gm2=StandingDead.g.m2., C3P_gm2=C3P.g.m2., AnnualGrass_gm2=AnnualGrass.g.m2.,
         Vulpia_gm2=Vulpia.g.m2., Bromus_gm2=Bromus.g.m2.) %>%
  select(Year:Clip_strip, C3P_gm2:Bromus_gm2) %>%
#  mutate(Drought=replace(Drought, Drought %in% c(1,2),0)) ## Ambient plots were coded as 1 or 2, changed to 0 here
  select(-Drought, -Grazing, -Block, -Paddock) %>%
  left_join(Plot_Key, by=c("Site","PlotID")) %>%
  dplyr::select(Year, Site, Block, Paddock, Drought, Grazing, PlotID, Clip_strip, slope_simple,
                C3P_gm2, C4P_gm2, Forb_gm2, ANPP_gm2, StandingDead_gm2, AnnualGrass_gm2, 
                Vulpia_gm2,Bromus_gm2)

anpp_tb_19_23 <- read.csv("GMDR_anpp_2018-2023_TB.csv") %>%
  filter(Year != 2018) %>%
  rename(Clip_strip=Clip.strip, C3P_gm2=C3P.g.m2.,C4P_gm2=C4P.g.m2., Forb_gm2=Forb.g.m2.,
         ANPP_gm2=ANPP.g.m2., StandingDead_gm2=StandingDead.g.m2., C3P_gm2=C3P.g.m2., AnnualGrass_gm2=AnnualGrass.g.m2.,
         Vulpia_gm2=Vulpia.g.m2., Bromus_gm2=Bromus.g.m2.) %>%
  select(Year, Site, PlotID, Clip_strip, C3P_gm2:Bromus_gm2) %>%
  left_join(Plot_Key, by=c("Site", "PlotID")) %>%
  dplyr::select(Year, Site, Block, Paddock, Drought, Grazing, PlotID, Clip_strip, slope_simple,
                C3P_gm2, C4P_gm2, Forb_gm2, ANPP_gm2, StandingDead_gm2,
                AnnualGrass_gm2, Vulpia_gm2, Bromus_gm2)

anpp_tb_18_23 <- anpp_tb_18 %>%
  bind_rows(anpp_tb_19_23) %>%
  filter(!is.na(slope_simple)) # remove pretreatment data for plots that weren't included in fence

## check data structure and ANPP values
with(anpp_tb_18_23, table(Year, Drought, Grazing))
with(anpp_tb_18_23, hist(ANPP_gm2))
filter(anpp_tb_18_23, ANPP_gm2>400) ### plot had a large fabaceae in it which is why the forb biomass is so high.. could remove this as an outlier

### REmove outlier clip strip
anpp_tb_18_23 <- anpp_tb_18_23[!(anpp_tb_18_23$Year==2019 & 
                                   anpp_tb_18_23$PlotID==28 & 
                                   anpp_tb_18_23$Clip_strip==2),]

### Average over clip strip to get plot level means
anpp_tb_plot_means <- anpp_tb_18_23 %>%
  group_by(Year, Site, Block, Paddock, Drought, Grazing,  PlotID, slope_simple) %>%
  summarize_at(vars(C3P_gm2:Bromus_gm2),.funs = c(mean),na.rm=TRUE) %>%
  rename(Plot=PlotID) %>%
  ungroup()
with(anpp_tb_plot_means, table(Year, Drought, Grazing))
with(anpp_tb_plot_means, hist(ANPP_gm2))

  ###
  ### Fort Keogh - Data cleaning
  ###

  anpp_fk_18_23 <- read.csv("DxG_ANPP_FK_2018.csv") %>%
    bind_rows(read.csv("DxG_ANPP_FK_2019.csv")) %>%
    bind_rows(read.csv("DxG_ANPP_FK_2020.csv")) %>%
    bind_rows(read.csv("DxG_ANPP_FK_2021_v2.csv")) %>%
    bind_rows(read.csv("DxG_ANPP_FK_2022.csv")) %>%
    bind_rows(read.csv("DxG_ANPP_FK_2023.csv")) %>%
    replace(is.na(.), 0) %>%
    mutate(C4_perennial=C4_perennial*10,
           C3_perennial=C3_perennial*10,
           annual_grass=annual_grass*10,
           forb=forb*10,
           subshrub_shrub=subshrub_shrub*10
    ) %>%
    mutate(anpp_gm2=C4_perennial+C3_perennial+annual_grass+forb+subshrub_shrub) %>%
    rename(Clip_strip=clip_subplot, C3P_gm2=C3_perennial,C4P_gm2=C4_perennial, Forb_gm2=forb,
           ANPP_gm2=anpp_gm2, AnnualGrass_gm2=annual_grass, Year=date, Site=site, Block=block, PlotID=plot) %>%
    dplyr::select(Site, Year, PlotID, Clip_strip, C4P_gm2:comments) %>%
    left_join(dplyr::select(Plot_Key, Site, PlotID, Block, Paddock, Drought, Grazing, slope_simple), 
              by=c("Site", "PlotID"))
  
  
  
  ### Check for missing data and that data values are reasonable -- LOOKS GOOD
  with(anpp_fk_18_23, table(Year, Drought, Grazing))
  with(anpp_fk_18_23,hist(ANPP_gm2))
  
  ### removing subplot outlier
  filter(anpp_fk_18_23,ANPP_gm2>600) 
  anpp_fk_18_23 <- anpp_fk_18_23 %>%
    filter(!(Site=="FK" & Year==2020 & PlotID==44 & Clip_strip==2))
  
  ### Average over clip strip to get plot level means
  anpp_fk_plot_means <- anpp_fk_18_23 %>%
    group_by(Year, Site, Block, Paddock, Drought, Grazing, PlotID, slope_simple) %>%
    summarize_at(vars(c(C4P_gm2:litter, ANPP_gm2)),.funs = c(mean),na.rm=TRUE) %>%
    rename(Plot=PlotID) %>%
    ungroup()
  
  #END_ANPP_SOURCE
  
  ### Check for missing data and that data values are reasonable -- LOOKS GOOD
  with(anpp_fk_plot_means, table(Year, Drought, Grazing))
  with(anpp_fk_plot_means, hist(ANPP_gm2))
  }

###
### Data cleaning BNPP
###
{
  ### Thunder basin
 
  ### Read in data and calculate BNPP
  bnpp_tb_19_22_raw <- read.csv("..//bnpp//TB_BNPP_2019-2022.csv") %>%
    rename(Subsample = Sub.Sample,
           crs_roots_mass = Coarse.Roots..Dry.wt..g.,
           crs_roots_tin_plus_ash = Coarse.Roots.Tin...Ash.wt..g.,
           crs_roots_tin_only = Coarse.Roots.Tin.wt..g.,
           fpm_roots_mass = FPMroots..Dry.wt..g.,
           fpm_roots_tin_plus_ash = FPMroots.Tin...Ash.wt..g.,
           fpm_roots_tin_only = FPMroots.Tin.wt..g.,
           fpm_som_mass = FPMsom..Dry.wt..g.,
           fpm_som_tin_plus_ash = FPMsom.Tin...Ash.wt..g.,
           fpm_som_tin_only = FPMsom.Tin.wt..g.,
           fpm_remainder_mass = FPMremainder..Dry.wt..g.,
           fpm_remainder_tin_plus_ash = FPMremainder.Tin...Ash.wt..g.,
           fpm_remainder_tin_only = FPMremainder.Tin.wt..g.
           ) %>%
    mutate(
      crs_roots_afdm = crs_roots_mass - (crs_roots_tin_plus_ash - crs_roots_tin_only),
      fpm_roots_afdm = fpm_roots_mass - (fpm_roots_tin_plus_ash - fpm_roots_tin_only),
      fpm_som_afdm = fpm_som_mass - (fpm_som_tin_plus_ash - fpm_som_tin_only),
      fpm_remainder_afdm = fpm_remainder_mass - (fpm_remainder_tin_plus_ash - fpm_remainder_tin_only)
           )
  
### There are a couple of missing values for ash mass so below I create regressions relating dry mass to AFDM to fill in those AFDM values from their dry mass
###  -- this was an issue for fpm SOM in plot 33 subplot 1 in 2022 and for fpm roots in plot 47 subplot 1 in 2022
ggplot(bnpp_tb_19_22_raw, aes(fpm_roots_mass, fpm_roots_afdm)) +geom_point() + geom_smooth(method="lm",se=F)
fpm_roots_afdm_summary <- summary(lm(fpm_roots_afdm ~ fpm_roots_mass, data=bnpp_tb_19_22_raw))
fpm_root_intercept <- fpm_roots_afdm_summary$coefficients[1,1]
fpm_root_slope <- fpm_roots_afdm_summary$coefficients[2,1]

### Doing 2022 alone since its a much tigher fit than all together (and the data we're filling is from 2022)... not sure what was going on in 2019 with ash mass
ggplot(filter(bnpp_tb_19_22_raw,Year==2022), aes(fpm_som_mass, fpm_som_afdm, col=factor(Year), label=PlotID)) +geom_text() + geom_smooth(method="lm",se=F)
fpm_som_afdm_summary <- summary(lm(fpm_som_afdm ~ fpm_som_mass, data=filter(bnpp_tb_19_22_raw, Year==2022)))
fpm_som_intercept <- fpm_som_afdm_summary$coefficients[1,1]
fpm_som_slope <- fpm_som_afdm_summary$coefficients[2,1]

plot33_1_som <- bnpp_tb_19_22_raw %>% filter(Year == 2022 & PlotID== 33 & Subsample == 1) %>% pull(fpm_som_mass)
plot47_1_root <- bnpp_tb_19_22_raw %>% filter(Year == 2022 & PlotID== 47 & Subsample == 1) %>% pull(fpm_roots_mass)

bnpp_tb_19_22 <- bnpp_tb_19_22_raw %>%
  mutate(fpm_som_afdm = replace(fpm_som_afdm, Year == 2022 & PlotID== 33 & Subsample == 1, plot33_1_som*fpm_som_slope + fpm_som_intercept)) %>%
  mutate(fpm_roots_afdm = replace(fpm_roots_afdm, Year == 2022 & PlotID== 47 & Subsample == 1, plot47_1_root*fpm_root_slope + fpm_root_intercept))
  
bnpp_tb_19_22 %>% filter(Year == 2022 & PlotID== 33 & Subsample == 1)
bnpp_tb_19_22 %>% filter(Year == 2022 & PlotID== 47 & Subsample == 1)

## just double check that nothing weird happened
hist(bnpp_tb_19_22$fpm_roots_afdm)
hist(bnpp_tb_19_22_raw$fpm_roots_afdm)
hist(bnpp_tb_19_22$fpm_som_afdm)
hist(bnpp_tb_19_22_raw$fpm_som_afdm)


## Back to formatting TB bnpp data
## Also, remove samples where 50% or more of the sample was lost during extraction -- plot 39 (sub 1 and 2), 37 (sub 1), 5 (sub 1) all in 2020
bnpp_tb_19_22 <- bnpp_tb_19_22 %>%
  filter(!(Year==2020 & PlotID==39 & Subsample==1),
         !(Year==2020 & PlotID==39 & Subsample==2),
         !(Year==2020 & PlotID==37 & Subsample==2),
         !(Year==2020 & PlotID==5 & Subsample==1)
         ) %>%
  dplyr::select(Year:Subsample, crs_roots_afdm, fpm_roots_afdm, fpm_som_afdm, fpm_remainder_afdm) %>%
  mutate(fpm_root_perc = fpm_roots_afdm/(fpm_roots_afdm + fpm_som_afdm)) %>%
  mutate(fpm_roots_from_remainder = fpm_remainder_afdm*fpm_root_perc) %>%
  mutate(total_root_afdm = crs_roots_afdm + fpm_roots_afdm + fpm_roots_from_remainder) %>%
  mutate(bnpp_gm2 = total_root_afdm/(pi*2.5^2/10000)) %>%
  left_join(Plot_Key, by=c("Site", "Block", "PlotID"))

bnpp_tb_plot_means <- bnpp_tb_19_22 %>%
  group_by(Year, Site, Block, Paddock, PlotID, Drought, Grazing, slope, slope_simple) %>%
  summarize(bnpp_gm2=mean(bnpp_gm2, na.rm=T)) %>%
  ungroup()
  

hist(bnpp_tb_19_22$bnpp_gm2)
ggplot(bnpp_tb_19_22, aes(x=Drought, y=bnpp_gm2)) + geom_jitter(width=2) + geom_smooth(method="lm",se=F) + facet_wrap(~Year)
#ggplot(bnpp_tb_rr, aes(x=Grazing, y=pchange_bnpp)) + geom_jitter(width=0.1) + facet_wrap(~Year)



###
### Fort Keogh BNPP
bnpp_fk_19_22 <- read.csv("..//bnpp//FK bnpp_19-22.csv") %>%
  group_by(Year, Block, PlotID, Subsample, root_cat) %>%
  summarize(AFDM=sum(AFDM)) %>% # This step is because they split the remainder measurements into two parts in 2022
  ungroup() %>%
  filter(!(Year==2022 & PlotID==51 & Subsample==1)) %>% # The fine live root entry here is incorrect but I couldn't figure out why, so we need to remove this sample altogether I think (KW May 2024)
  pivot_wider(names_from = root_cat, values_from = AFDM) %>%
  rename(
    crs_roots_afdm='Live Roots',
    fpm_roots_afdm='Fine live roots',
    fpm_som_afdm='Fine SOM',
    fpm_remainder_afdm='Remainder'
         ) %>%
  mutate(fpm_root_perc = fpm_roots_afdm/(fpm_roots_afdm + fpm_som_afdm)) %>%
  mutate(fpm_roots_from_remainder = fpm_remainder_afdm*fpm_root_perc) %>%
  mutate(total_root_afdm = crs_roots_afdm + fpm_roots_afdm + fpm_roots_from_remainder) %>%
  mutate(bnpp_gm2 = total_root_afdm/(pi*2.5^2/10000)) %>%
  mutate(Site="FK") %>%
  left_join(Plot_Key, by=c("Site", "Block", "PlotID"))

bnpp_fk_plot_means <- bnpp_fk_19_22 %>%
  group_by(Year, Site, Block, Paddock, PlotID, Drought, Grazing, slope, slope_simple) %>%
  summarize(bnpp_gm2=mean(bnpp_gm2, na.rm=T)) %>%
  ungroup()


hist(bnpp_fk_19_22$bnpp_gm2)

ggplot(bnpp_fk_19_22, aes(x=Drought, y=bnpp_gm2)) + geom_jitter(width=.1) + geom_smooth(method="lm",se=F) + facet_wrap(~Year)
}

###
### Calculate Total NPP and combine with ANPP and BNPP and both sites
###
{
totnpp_tb_plot_means <- anpp_tb_plot_means %>%
  dplyr::select(Year, Site, Block, Paddock, Drought, Grazing, Plot, slope_simple, ANPP_gm2) %>%
  full_join(
    bnpp_tb_plot_means %>% dplyr::select(Year, Site, Block, Paddock, PlotID, Drought, Grazing, slope_simple, bnpp_gm2) %>%
      rename(Plot=PlotID),
    by=c("Year","Site","Block","Paddock","Drought","Grazing","Plot","slope_simple")
  ) %>%
  mutate(totnpp=ANPP_gm2 + bnpp_gm2)

  ggplot(totnpp_tb_plot_means, aes(x=Drought, y=totnpp)) + geom_jitter(width=1) + geom_smooth(method="lm",se=F) + facet_wrap(~Year)
  
totnpp_fk_plot_means <- anpp_fk_plot_means %>%
  dplyr::select(Year, Site, Block, Paddock, Drought, Grazing, Plot, slope_simple, ANPP_gm2) %>%
  full_join(
    bnpp_fk_plot_means %>% dplyr::select(Year, Site, Block, Paddock, PlotID, Drought, Grazing, slope_simple, bnpp_gm2) %>%
      rename(Plot=PlotID),
    by=c("Year","Site","Block","Paddock","Drought","Grazing","Plot","slope_simple")
  ) %>%
  mutate(totnpp=ANPP_gm2 + bnpp_gm2)

ggplot(totnpp_fk_plot_means, aes(x=Drought, y=totnpp)) + geom_jitter(width=1) + geom_smooth(method="lm",se=F) + facet_wrap(~Year)


npp_master <- totnpp_tb_plot_means %>%
  bind_rows(totnpp_fk_plot_means) %>%
  rename(ANPP=ANPP_gm2, BNPP=bnpp_gm2, totNPP=totnpp) %>%
  pivot_longer(cols=ANPP:totNPP, names_to="npp_type", values_to="npp_gm2")

write.csv(npp_master, file=paste0(write_dir,"data_sets\\ANPP BNPP total NPP_plot level_cleaned and processed from R",Sys.Date(),".csv"), row.names=F)


}

###
### Calculate response ratios
###
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
  
}

###
### Run models and testing for normality - Response ratios
###
{
  
  ### Fort Keogh
  ###
  
  ### ANPP model all years
  fk_anpp_model_full <- lme(lnrr_npp ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                                   , data=filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="ANPP")
                                   , random = ~1 |Block/Paddock/Plot
                                   , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit)
  fk_anpp_anova <- anova.lme(fk_anpp_model_full, type="marginal")
  plot(fk_anpp_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(fk_anpp_model_full, abline = c(0,1)) ## qqplot
  hist(filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="ANPP")$lnrr_npp)
  
  emtrends(fk_anpp_model_full, "Year", var="Drought")
  test(emtrends(fk_anpp_model_full, "Year", var="Drought"))
  emmeans(fk_anpp_model_full, "Grazing", by="Year")
  pairs(emmeans(fk_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison

  # Save to writable tables
  fk_anpp_anova_df <- data.frame(effect=row.names(fk_anpp_anova), fk_anpp_anova, site="FK")
  fk_anpp_emtrends <- data.frame(test(emtrends(fk_anpp_model_full, "Year", var="Drought")),
                                 site="FK")
  
  ### Split by year to get R2 values for significant regressions (I DON'T THINK R2 IS WORKING, ITS GIVING ME 0.99 FOR R2, AND I'M PRETTY SURE THAT IS INCORRECT)
    # 2019
  fk_anpp_2019_lme <- lme(lnrr_npp ~ Drought
                          , data=filter(npp_rr, Year==2019 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_2019_lme, type="marginal")
  performance::r2(fk_anpp_2019_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
    # 2020
  fk_anpp_2020_lme <- lme(lnrr_npp ~ Drought
                          , data=filter(npp_rr, Year==2020 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_2020_lme, type="marginal")
  performance::r2(fk_anpp_2020_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2021
  fk_anpp_2021_lme <- lme(lnrr_npp ~ Drought
                          , data=filter(npp_rr, Year==2021 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_2021_lme, type="marginal")
  performance::r2(fk_anpp_2021_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2023
  fk_anpp_2023_lme <- lme(lnrr_npp ~ Drought
                          , data=filter(npp_rr, Year==2023 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_2023_lme, type="marginal")
  performance::r2(fk_anpp_2023_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  

  ### Thunder Basin
  ###
  
  ### ANPP model all years
  tb_anpp_model_full <- lme(lnrr_npp ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="ANPP")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  tb_anpp_anova <- anova.lme(tb_anpp_model_full, type="marginal")
  plot(tb_anpp_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(tb_anpp_model_full, abline = c(0,1)) ## qqplot
  hist(filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="ANPP")$lnrr_npp)
  
  emtrends(tb_anpp_model_full, "Year", var="Drought")
  test(emtrends(tb_anpp_model_full, "Year", var="Drought"))
  
  # Save to writable tables
  tb_anpp_anova_df <- data.frame(effect=row.names(tb_anpp_anova), tb_anpp_anova, site="TB")
  tb_anpp_emtrends <- data.frame(test(emtrends(tb_anpp_model_full, "Year", var="Drought")),
                                 site="TB")
  
  ### Split by year to get R2 values for significant regressions
  # 2020
  tb_anpp_2020_lme <- lme(lnrr_npp ~ Drought
                          , data=filter(npp_rr, Year==2020 & Site=="TB" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(tb_anpp_2020_lme, type="marginal")
  performance::r2(tb_anpp_2020_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2021
  tb_anpp_2021_lme <- lme(lnrr_npp ~ Drought
                          , data=filter(npp_rr, Year==2021 & Site=="TB" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock
                          , na.action = na.omit)
  anova.lme(tb_anpp_2021_lme, type="marginal")
  performance::r2(tb_anpp_2021_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  
  ### WRite tables of model output for both sites
  anova_anpp_df <- tb_anpp_anova_df %>% bind_rows(fk_anpp_anova_df)
  emtrends_anpp_df <- tb_anpp_emtrends %>% bind_rows(fk_anpp_emtrends)
  write.csv(anova_anpp_df, file=paste0(write_dir,"tables\\anpp lme ANCOVA output_both sites",Sys.Date(),".csv"), row.names=F)
  write.csv(emtrends_anpp_df, file=paste0(write_dir,"tables\\anpp emtrends_both sites",Sys.Date(),".csv"), row.names=F)
  
 ###
 ### BNPP models
 ###
 
  ### Fort Keogh
  ###
  
  ### BNPP model all years
  fk_bnpp_model_full <- lme(lnrr_npp ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="BNPP")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  fk_bnpp_anova <- anova.lme(fk_bnpp_model_full, type="marginal")
  plot(fk_bnpp_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(fk_bnpp_model_full, abline = c(0,1)) ## qqplot
  hist(filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="BNPP")$lnrr_npp)
  
  # Only a significant year by grazing interaction, so split by grazing
  emmeans(fk_bnpp_model_full, "Grazing", by="Year")
  pairs(emmeans(fk_bnpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison except at p=.09 contrast between H and L grazing in 2019
  
  # Save to writable tables
  fk_bnpp_anova_df <- data.frame(effect=row.names(fk_bnpp_anova), fk_bnpp_anova, site="FK")

  
  ### Thunder Basin
  ###
  
  ### BNPP model all years
  tb_bnpp_model_full <- lme(lnrr_npp ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="BNPP")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  tb_bnpp_anova <- anova.lme(tb_bnpp_model_full, type="marginal")
  plot(tb_bnpp_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(tb_bnpp_model_full, abline = c(0,1)) ## qqplot
  hist(filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="BNPP")$lnrr_npp)
  
  ### NO significant interactions so we stop here

  
  
  ### WRite tables of model output for both sites
  anova_bnpp_df <- tb_bnpp_anova_df %>% bind_rows(fk_bnpp_anova_df)
  write.csv(anova_bnpp_df, file=paste0(write_dir,"tables\\bnpp lme ANCOVA output_both sites",Sys.Date(),".csv"), row.names=F)

   
}
 
###
### Run models and testing for normality - Raw data
###
{
  ### ANPP data
  ###

    ## Fort Keogh
    ##
  fk_anpp_raw_model_full <- lme(log(npp_gm2) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(npp_master, Year %in% 2019:2023 & Site=="FK" & npp_type=="ANPP")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  fk_anpp_raw_anova <- anova.lme(fk_anpp_raw_model_full, type="marginal")
  plot(fk_anpp_raw_model_full, type=c("p","smooth"), col.line=1) # NEEDS TO BE LOG TRANSFORMED -- LOG LOOKS MUCH BETTER
  qqnorm(fk_anpp_raw_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(npp_master, Year %in% 2019:2023 & Site=="FK" & npp_type=="ANPP")$npp_gm2))
  
  emtrends(fk_anpp_raw_model_full, "Year", var="Drought")
  test(emtrends(fk_anpp_raw_model_full, "Year", var="Drought"))

  # Save to writable tables
  fk_anpp_raw_anova_df <- data.frame(effect=row.names(fk_anpp_raw_anova), fk_anpp_raw_anova, site="FK")
  fk_anpp_raw_emtrends <- data.frame(test(emtrends(fk_anpp_raw_model_full, "Year", var="Drought")),
                                 site="FK")
  
  ### Split by year to get R2 values for significant regressions (I DON'T THINK R2 IS WORKING, ITS GIVING ME 0.99 FOR R2, AND I'M PRETTY SURE THAT IS INCORRECT)
  # 2019
  fk_anpp_raw_2019_lme <- lme(log(npp_gm2) ~ Drought
                          , data=filter(npp_master, Year==2019 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_raw_2019_lme, type="marginal")
  performance::r2(fk_anpp_raw_2019_lme) #  Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2020
  fk_anpp_raw_2020_lme <- lme(log(npp_gm2) ~ Drought
                          , data=filter(npp_master, Year==2020 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_raw_2020_lme, type="marginal")
  performance::r2(fk_anpp_raw_2020_lme) # R2 OF 1... NOT RIGHT.. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2021
  fk_anpp_raw_2021_lme <- lme(log(npp_gm2) ~ Drought
                          , data=filter(npp_master, Year==2021 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_raw_2021_lme, type="marginal")
  performance::r2(fk_anpp_raw_2021_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2023
  fk_anpp_raw_2023_lme <- lme(log(npp_gm2) ~ Drought
                          , data=filter(npp_master, Year==2023 & Site=="FK" & npp_type=="ANPP")
                          , random = ~1 |Block/Paddock/Plot
                          , na.action = na.omit)
  anova.lme(fk_anpp_raw_2023_lme, type="marginal")
  performance::r2(fk_anpp_raw_2023_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  
    ##
    ## Thunder Basin
  
  tb_anpp_raw_model_full <- lme(log(npp_gm2) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                                , data=filter(npp_master, Year %in% 2019:2023 & Site=="TB" & npp_type=="ANPP")
                                , random = ~1 |Block/Paddock/Plot
                                , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
  tb_anpp_raw_anova <- anova.lme(tb_anpp_raw_model_full, type="marginal")
  plot(tb_anpp_raw_model_full, type=c("p","smooth"), col.line=1) # NEEDS TO BE LOG TRANSFORMED -- LOG LOOKS MUCH BETTER
  qqnorm(tb_anpp_raw_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(npp_master, Year %in% 2019:2023 & Site=="TB" & npp_type=="ANPP")$npp_gm2))
  
  emtrends(tb_anpp_raw_model_full, "Year", var="Drought")
  test(emtrends(tb_anpp_raw_model_full, "Year", var="Drought"))
  
  # Save to writable tables
  tb_anpp_raw_anova_df <- data.frame(effect=row.names(tb_anpp_raw_anova), tb_anpp_raw_anova, site="TB")
  tb_anpp_raw_emtrends <- data.frame(test(emtrends(tb_anpp_raw_model_full, "Year", var="Drought")),
                                     site="TB")
  
  ### Split by year to get R2 values for significant regressions 

    # 2020
  tb_anpp_raw_2020_lme <- lme(log(npp_gm2) ~ Drought
                              , data=filter(npp_master, Year==2020 & Site=="TB" & npp_type=="ANPP")
                              , random = ~1 |Block/Paddock/Plot
                              , na.action = na.omit)
  anova.lme(tb_anpp_raw_2020_lme, type="marginal")
  performance::r2(tb_anpp_raw_2020_lme) # R2 OF 1... NOT RIGHT.. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2021
  tb_anpp_raw_2021_lme <- lme(log(npp_gm2) ~ Drought
                              , data=filter(npp_master, Year==2021 & Site=="TB" & npp_type=="ANPP")
                              , random = ~1 |Block/Paddock/Plot
                              , na.action = na.omit)
  anova.lme(tb_anpp_raw_2021_lme, type="marginal")
  performance::r2(tb_anpp_raw_2021_lme) # Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  
  ### Write tables of model output for both sites
  anova_anpp_raw_df <- tb_anpp_raw_anova_df %>% bind_rows(fk_anpp_raw_anova_df)
  emtrends_anpp_raw_df <- tb_anpp_raw_emtrends %>% bind_rows(fk_anpp_raw_emtrends)
  write.csv(anova_anpp_raw_df, file=paste0(write_dir,"tables\\anpp raw lme ANCOVA output_both sites",Sys.Date(),".csv"), row.names=F)
  write.csv(emtrends_anpp_raw_df, file=paste0(write_dir,"tables\\anpp raw emtrends_both sites",Sys.Date(),".csv"), row.names=F)
  
  ### BNPP data
  ###
  
  ## Fort Keogh
  ##
  fk_bnpp_raw_model_full <- lme(log(npp_gm2) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                                , data=filter(npp_master, Year %in% 2019:2023 & Site=="FK" & npp_type=="BNPP")
                                , random = ~1 |Block/Paddock/Plot
                                , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
  fk_bnpp_raw_anova <- anova.lme(fk_bnpp_raw_model_full, type="marginal")
  plot(fk_bnpp_raw_model_full, type=c("p","smooth"), col.line=1) # NEEDS TO BE LOG TRANSFORMED -- LOG LOOKS MUCH BETTER
  qqnorm(fk_bnpp_raw_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(npp_master, Year %in% 2019:2023 & Site=="FK" & npp_type=="BNPP")$npp_gm2))
  
  ## NO SIGNIFICANT INTERACTIONS SO STOP HERE
  
  # Save to writable tables
  fk_bnpp_raw_anova_df <- data.frame(effect=row.names(fk_bnpp_raw_anova), fk_bnpp_raw_anova, site="FK")
  
  
  ##
  ## Thunder Basin
  
  tb_bnpp_raw_model_full <- lme(log(npp_gm2) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                                , data=filter(npp_master, Year %in% 2019:2023 & Site=="TB" & npp_type=="BNPP")
                                , random = ~1 |Block/Paddock/Plot
                                , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
  tb_bnpp_raw_anova <- anova.lme(tb_bnpp_raw_model_full, type="marginal")
  plot(tb_bnpp_raw_model_full, type=c("p","smooth"), col.line=1) # NEEDS TO BE LOG TRANSFORMED -- LOG LOOKS MUCH BETTER
  qqnorm(tb_bnpp_raw_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(npp_master, Year %in% 2019:2023 & Site=="TB" & npp_type=="BNPP")$npp_gm2))
  
  emmeans(tb_bnpp_raw_model_full, "Grazing", by="Year")
  pairs(emmeans(tb_bnpp_raw_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  emtrends(tb_bnpp_raw_model_full, "Year", var="Drought") ## No significant interaction, but I'm curious -- did come out with a sig negative drought trend in 2020... but I think we stop here since no significant interaction
  test(emtrends(tb_bnpp_raw_model_full, "Year", var="Drought"))
  
  # Save to writable tables
  tb_bnpp_raw_anova_df <- data.frame(effect=row.names(tb_bnpp_raw_anova), tb_bnpp_raw_anova, site="TB")
                                     
  
  
  ### Write tables of model output for both sites
  anova_bnpp_raw_df <- tb_bnpp_raw_anova_df %>% bind_rows(fk_bnpp_raw_anova_df)
  write.csv(anova_bnpp_raw_df, file=paste0(write_dir,"tables\\bnpp raw lme ANCOVA output_both sites",Sys.Date(),".csv"), row.names=F)
}

###
### Plotting ANPP and BNPP ln(RR) by drought
###
{
  npp_rr_drought_means <- npp_rr %>%
    filter(Year %in% 2019:2023) %>% ### plot 44 is a clear outlier
    group_by(Site, Year, npp_type, Drought) %>%
    summarize_at(.vars=vars(lnrr_npp), .funs=list(mean=mean, se=SE_function), na.rm=T) %>%
    ungroup() %>%
    rename(npp_mean = mean, npp_se = se)
  
  npp_fig <- ggplot(filter(npp_rr_drought_means,Year%in%2019:2023 & npp_type %in% c("ANPP","BNPP")), 
                    aes(x=Drought, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se, fill=npp_type, shape=npp_type))+
    scale_shape_manual(values=c(21,24))+
    scale_fill_manual(values=c("springgreen3","saddlebrown")) +
    scale_color_manual(values=c("springgreen3","saddlebrown")) +
    geom_hline(yintercept=0, col="grey") +
    geom_smooth(inherit.aes=F, data=filter(npp_rr,Year%in%2019:2023&npp_type %in% c("ANPP", "BNPP")),
                aes(Drought, lnrr_npp, col=npp_type), method="lm",se=F) +
    geom_errorbar(width=3) +
    geom_point(size=2, col="black") +
    facet_grid(Site~Year) +
    xlab("Rainfall reduction (%)") + ylab("Net primary productivity (g/m2)") +
    theme_few()
  
  pdf(paste0(write_dir,"figures//anpp and bnpp response ratios_v1_",Sys.Date(),".pdf"), width=10.5, height=4.5, useDingbats = F)
  print(npp_fig)
  dev.off()

}


###
### Plotting raw ANPP and BNPP by drought
###
{

### Calculate means by drought only
npp_drt_means <- npp_master %>%
  group_by(Year, Site, npp_type, Drought) %>%
  summarize_at(vars(c(npp_gm2)),.funs = c(mean=mean, se=SE_function),na.rm=TRUE) %>% 
  ungroup() %>%
  rename(npp_mean=mean,
         npp_se=se)

anpp_raw_fig <- ggplot(filter(npp_drt_means,Year %in%2019:2023 & npp_type == "ANPP"),
                      aes(x=Drought, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se)) +
                geom_errorbar(width=1) +
                geom_point() +
                geom_smooth(method="lm",se=F) +
                facet_grid(Site~Year) +
                ylab("ANPP g/m2") + xlab("Rainfall reduction (%)") +
                theme_few()

pdf(paste0(write_dir, "figures//raw anpp figure", Sys.Date(), ".pdf"), width=9, height=4.5, useDingbats = F)
print(anpp_raw_fig)
dev.off()

ggsave(filename=paste0(write_dir, "figures//raw anpp figure", Sys.Date(), ".jpeg"), plot=anpp_raw_fig, width=9, height=4.5, units="in", dpi=300)

bnpp_raw_fig <- ggplot(filter(npp_drt_means,Year %in%2019:2023 & npp_type == "BNPP"),
                       aes(x=Drought, y=npp_mean, ymin=npp_mean-npp_se, ymax=npp_mean+npp_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(Site~Year, scales="free_y") +
  ylab("BNPP g/m2") + xlab("Rainfall reduction (%)") +
  theme_few()

pdf(paste0(write_dir, "figures//raw bnpp figure", Sys.Date(), ".pdf"), width=9, height=4.5, useDingbats = F)
print(bnpp_raw_fig)
dev.off()

ggsave(filename=paste0(write_dir, "figures//raw bnpp figure", Sys.Date(), ".jpeg"), plot=bnpp_raw_fig, width=9, height=4.5, units="in", dpi=300)
}

