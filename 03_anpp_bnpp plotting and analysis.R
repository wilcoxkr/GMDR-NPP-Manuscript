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
# ggplot(bnpp_fk_19_22_raw, aes(x=Year, y=AFDM)) + geom_jitter(width=.1) + facet_wrap(~root_cat, scales="free_y")
# filter(bnpp_fk_19_22_raw, root_cat=="Fine SOM" & AFDM==0)
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
### Run models and testing for normality
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
  performance::r2(fk_anpp_2019_lme) # THIS DOESN'T SEEM RIGHT - Marginal R2 considers only the variance of the fixed effects, which is what I want
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
  write.csv(tb_anpp_emtrends, file=paste0(write_dir,"tables\\anpp emtrends_both sites",Sys.Date(),".csv"), row.names=F)
  
 ###
 ### BNPP models
 ###
 
 ###
 ### Thunder Basin
 ###
 
 tb_bnpp_model_full <-   anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='lnrr_npp',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="BNPP")
 ) 
 
 
 
 
 ## just checking for normality and testing how similar model results are to anovat3 function -- good news is that it is identical
 tb_bnpp_model_noanovaT3fxn <- lme(lnrr_npp ~ Year*Drought*Grazing
                                   , data=filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="BNPP")
                                   , random = ~1 |Block/Paddock/Plot
                                   , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit)
 anova(tb_bnpp_model_noanovaT3fxn, type="marginal")
 summary(tb_bnpp_model_noanovaT3fxn)
 check_model(tb_bnpp_model_noanovaT3fxn)
 
 ### Split by year
 # 2019
 tb_bnpp_2019_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2019 & Site=="TB" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_bnpp_2019_lme, type="marginal")
 
 # 2020
 tb_bnpp_2020_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2020 & Site=="TB" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_bnpp_2020_lme, type="marginal")
 
 # 2021
 tb_bnpp_2021_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2021 & Site=="TB" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_bnpp_2021_lme, type="marginal")
 
 # 2022
 tb_bnpp_2022_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2022 & Site=="TB" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_bnpp_2022_lme, type="marginal")
 

 ### Create model table
 tb_bnpp_model_out <- 
   as.data.frame(anova(tb_bnpp_2019_lme, type="marginal")) %>% mutate(Year=2019, Data="BNPP", Site="TB", Term=rownames(anova(tb_bnpp_2019_lme))) %>%
   bind_rows(
     as.data.frame(anova(tb_bnpp_2020_lme, type="marginal")) %>% mutate(Year=2020, Data="BNPP", Site="TB", Term=rownames(anova(tb_bnpp_2019_lme))),
     as.data.frame(anova(tb_bnpp_2021_lme, type="marginal")) %>% mutate(Year=2021, Data="BNPP", Site="TB", Term=rownames(anova(tb_bnpp_2019_lme))),
     as.data.frame(anova(tb_bnpp_2022_lme, type="marginal")) %>% mutate(Year=2022, Data="BNPP", Site="TB", Term=rownames(anova(tb_bnpp_2019_lme))),
   ) %>%
   dplyr::select(Data, Site, Year, Term, numDF, denDF, "F-value", "p-value")
 
 write.csv(tb_bnpp_model_out, file=paste0(write_dir,"tables\\","tb bnpp lme model output by year_",Sys.Date(),".csv"), row.names=F)
 
 ### Create drought slope table - not using this atm (June 17, 2024)
 tb_bnpp_drought_slopes <- data.frame(
   Year=2019:2022, 
   Data="BNPP", 
   Site="TB",
   pchange_slope = c(
     tb_bnpp_2019_lme$coefficients$fixed[2],
     tb_bnpp_2020_lme$coefficients$fixed[2],
     tb_bnpp_2021_lme$coefficients$fixed[2],
     tb_bnpp_2022_lme$coefficients$fixed[2]
   ),
   pchange_slope_se = c(
     summary(tb_bnpp_2019_lme)$tTable[2,2],
     summary(tb_bnpp_2020_lme)$tTable[2,2],
     summary(tb_bnpp_2021_lme)$tTable[2,2],
     summary(tb_bnpp_2022_lme)$tTable[2,2]
   )
 )
 
 write.csv(tb_bnpp_drought_slopes, file=paste0(write_dir,"tables\\","tb bnpp drought slope and se from lme_",Sys.Date(),".csv"), row.names=F)
 
 
 ###
 ### Fort Keogh
 ###
 
 fk_bnpp_model_full <-   anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='lnrr_npp',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="BNPP")
 ) 
 
 
 
 
 ## just checking for normality and testing how similar model results are to anovat3 function -- good news is that it is identical
 fk_bnpp_model_noanovaT3fxn <- lme(lnrr_npp ~ Year*Drought*Grazing
                                   , data=filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="BNPP")
                                   , random = ~1 |Block/Paddock/Plot
                                   , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit)
 anova(fk_bnpp_model_noanovaT3fxn, type="marginal")
 summary(fk_bnpp_model_noanovaT3fxn)
 check_model(fk_bnpp_model_noanovaT3fxn)
 
 ### Split by year
 # 2019
 fk_bnpp_2019_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2019 & Site=="FK" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_bnpp_2019_lme, type="marginal")
 
 # 2020
 fk_bnpp_2020_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2020 & Site=="FK" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_bnpp_2020_lme, type="marginal")
 
 # 2021
 fk_bnpp_2021_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2021 & Site=="FK" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_bnpp_2021_lme, type="marginal")
 
 # 2022
 fk_bnpp_2022_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2022 & Site=="FK" & npp_type=="BNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_bnpp_2022_lme, type="marginal")
 
 
 ### Create model table
 fk_bnpp_model_out <- 
   as.data.frame(anova(fk_bnpp_2019_lme, type="marginal")) %>% mutate(Year=2019, Data="BNPP", Site="FK", Term=rownames(anova(fk_bnpp_2019_lme))) %>%
   bind_rows(
     as.data.frame(anova(fk_bnpp_2020_lme, type="marginal")) %>% mutate(Year=2020, Data="BNPP", Site="FK", Term=rownames(anova(fk_bnpp_2019_lme))),
     as.data.frame(anova(fk_bnpp_2021_lme, type="marginal")) %>% mutate(Year=2021, Data="BNPP", Site="FK", Term=rownames(anova(fk_bnpp_2019_lme))),
     as.data.frame(anova(fk_bnpp_2022_lme, type="marginal")) %>% mutate(Year=2022, Data="BNPP", Site="FK", Term=rownames(anova(fk_bnpp_2019_lme))),
   ) %>%
   dplyr::select(Data, Site, Year, Term, numDF, denDF, "F-value", "p-value")
 
 write.csv(fk_bnpp_model_out, file=paste0(write_dir,"tables\\","fk bnpp lme model output by year_",Sys.Date(),".csv"), row.names=F)
 
 
 ### Split by grazing in 2022 du to significant interaction
 # light grazign
 fk_bnpp_2022_lg_lme <- lme(lnrr_npp ~ Drought
                         , data=filter(npp_rr, Year==2022 & Site=="FK" & npp_type=="BNPP" & Grazing == "MLLMM")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_bnpp_2022_lg_lme, type="marginal")
 
 # moderate grazing
 fk_bnpp_2022_mg_lme <- lme(lnrr_npp ~ Drought
                            , data=filter(npp_rr, Year==2022 & Site=="FK" & npp_type=="BNPP" & Grazing == "MMMMM")
                            , random = ~1 |Block/Paddock/Plot
                            , na.action = na.omit)
 anova(fk_bnpp_2022_mg_lme, type="marginal")
 
 # heavy grazing
 fk_bnpp_2022_hg_lme <- lme(lnrr_npp ~ Drought
                            , data=filter(npp_rr, Year==2022 & Site=="FK" & npp_type=="BNPP" & Grazing == "HHMMM")
                            , random = ~1 |Block/Paddock/Plot
                            , na.action = na.omit)
 anova(fk_bnpp_2022_hg_lme, type="marginal")
 
 ### Create drought slope table
 fk_bnpp_drought_slopes <- data.frame(
   Year=2019:2022, 
   Data="BNPP", 
   Site="fk",
   pchange_slope = c(
     fk_bnpp_2019_lme$coefficients$fixed[2],
     fk_bnpp_2020_lme$coefficients$fixed[2],
     fk_bnpp_2021_lme$coefficients$fixed[2],
     fk_bnpp_2022_lme$coefficients$fixed[2]
   ),
   pchange_slope_se = c(
     summary(fk_bnpp_2019_lme)$tTable[2,2],
     summary(fk_bnpp_2020_lme)$tTable[2,2],
     summary(fk_bnpp_2021_lme)$tTable[2,2],
     summary(fk_bnpp_2022_lme)$tTable[2,2]
   )
 )
 
 write.csv(fk_bnpp_drought_slopes, file=paste0(write_dir,"tables\\","fk bnpp drought slope and se from lme_",Sys.Date(),".csv"), row.names=F)
 
   

 
 
 ##### tot npp
 
 ###
 ### Thunder Basin
 ###
 
 tb_totnpp_model_full <-   anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='lnrr_npp',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="totNPP")
 ) 
 
 
 
 
 ## just checking for normality and testing how similar model results are to anovat3 function -- good news is that it is identical
 tb_totnpp_model_noanovaT3fxn <- lme(lnrr_npp ~ Year*Drought*Grazing
                                   , data=filter(npp_rr, Year %in% 2019:2023 & Site=="TB" & npp_type=="totNPP")
                                   , random = ~1 |Block/Paddock/Plot
                                   , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit)
 anova(tb_totnpp_model_noanovaT3fxn, type="marginal")
 summary(tb_totnpp_model_noanovaT3fxn)
 check_model(tb_totnpp_model_noanovaT3fxn)
 
 ### Split by year
 # 2019
 tb_totnpp_2019_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2019 & Site=="TB" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_totnpp_2019_lme, type="marginal")
 
 # 2020
 tb_totnpp_2020_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2020 & Site=="TB" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_totnpp_2020_lme, type="marginal")
 
 # 2021
 tb_totnpp_2021_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2021 & Site=="TB" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_totnpp_2021_lme, type="marginal")
 
 # 2022
 tb_totnpp_2022_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2022 & Site=="TB" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(tb_totnpp_2022_lme, type="marginal")
 
 
 ### Create model table
 tb_totnpp_model_out <- 
   as.data.frame(anova(tb_totnpp_2019_lme, type="marginal")) %>% mutate(Year=2019, Data="totNPP", Site="TB", Term=rownames(anova(tb_totnpp_2019_lme))) %>%
   bind_rows(
     as.data.frame(anova(tb_totnpp_2020_lme, type="marginal")) %>% mutate(Year=2020, Data="totNPP", Site="TB", Term=rownames(anova(tb_totnpp_2019_lme))),
     as.data.frame(anova(tb_totnpp_2021_lme, type="marginal")) %>% mutate(Year=2021, Data="totNPP", Site="TB", Term=rownames(anova(tb_totnpp_2019_lme))),
     as.data.frame(anova(tb_totnpp_2022_lme, type="marginal")) %>% mutate(Year=2022, Data="totNPP", Site="TB", Term=rownames(anova(tb_totnpp_2019_lme))),
   ) %>%
   dplyr::select(Data, Site, Year, Term, numDF, denDF, "F-value", "p-value")
 
 write.csv(tb_totnpp_model_out, file=paste0(write_dir,"tables\\","tb totnpp lme model output by year_",Sys.Date(),".csv"), row.names=F)
 
 ### Create drought slope table - not using this atm (June 17, 2024)
 tb_totnpp_drought_slopes <- data.frame(
   Year=2019:2022, 
   Data="totNPP", 
   Site="TB",
   pchange_slope = c(
     tb_totnpp_2019_lme$coefficients$fixed[2],
     tb_totnpp_2020_lme$coefficients$fixed[2],
     tb_totnpp_2021_lme$coefficients$fixed[2],
     tb_totnpp_2022_lme$coefficients$fixed[2]
   ),
   pchange_slope_se = c(
     summary(tb_totnpp_2019_lme)$tTable[2,2],
     summary(tb_totnpp_2020_lme)$tTable[2,2],
     summary(tb_totnpp_2021_lme)$tTable[2,2],
     summary(tb_totnpp_2022_lme)$tTable[2,2]
   )
 )
 
 write.csv(tb_totnpp_drought_slopes, file=paste0(write_dir,"tables\\","tb totnpp drought slope and se from lme_",Sys.Date(),".csv"), row.names=F)
 
 
 ###
 ### Fort Keogh
 ###
 
 fk_totnpp_model_full <-   anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='lnrr_npp',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="totNPP")
 ) 
 
 
 
 
 ## just checking for normality and testing how similar model results are to anovat3 function -- good news is that it is identical
 fk_totnpp_model_noanovaT3fxn <- lme(lnrr_npp ~ Year*Drought*Grazing
                                   , data=filter(npp_rr, Year %in% 2019:2023 & Site=="FK" & npp_type=="totNPP")
                                   , random = ~1 |Block/Paddock/Plot
                                   , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit)
 anova(fk_totnpp_model_noanovaT3fxn, type="marginal")
 summary(fk_totnpp_model_noanovaT3fxn)
 check_model(fk_totnpp_model_noanovaT3fxn)
 
 ### Split by year
 # 2019
 fk_totnpp_2019_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2019 & Site=="FK" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_totnpp_2019_lme, type="marginal")
 
 # 2020
 fk_totnpp_2020_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2020 & Site=="FK" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_totnpp_2020_lme, type="marginal")
 
 # 2021
 fk_totnpp_2021_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2021 & Site=="FK" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_totnpp_2021_lme, type="marginal")
 
 # 2022
 fk_totnpp_2022_lme <- lme(lnrr_npp ~ Drought*Grazing
                         , data=filter(npp_rr, Year==2022 & Site=="FK" & npp_type=="totNPP")
                         , random = ~1 |Block/Paddock/Plot
                         , na.action = na.omit)
 anova(fk_totnpp_2022_lme, type="marginal")
 
 
 ### Create model table
 fk_totnpp_model_out <- 
   as.data.frame(anova(fk_totnpp_2019_lme, type="marginal")) %>% mutate(Year=2019, Data="totNPP", Site="FK", Term=rownames(anova(fk_totnpp_2019_lme))) %>%
   bind_rows(
     as.data.frame(anova(fk_totnpp_2020_lme, type="marginal")) %>% mutate(Year=2020, Data="totNPP", Site="FK", Term=rownames(anova(fk_totnpp_2019_lme))),
     as.data.frame(anova(fk_totnpp_2021_lme, type="marginal")) %>% mutate(Year=2021, Data="totNPP", Site="FK", Term=rownames(anova(fk_totnpp_2019_lme))),
     as.data.frame(anova(fk_totnpp_2022_lme, type="marginal")) %>% mutate(Year=2022, Data="totNPP", Site="FK", Term=rownames(anova(fk_totnpp_2019_lme))),
   ) %>%
   dplyr::select(Data, Site, Year, Term, numDF, denDF, "F-value", "p-value")
 
 write.csv(fk_totnpp_model_out, file=paste0(write_dir,"tables\\","fk totnpp lme model output by year_",Sys.Date(),".csv"), row.names=F)
}


###
### Plotting raw NPP data all together for each site
###

ggplot(filter(npp_rr, npp_type=="ANPP"), aes(x=Drought, y=lnrr_npp)) + geom_point() + geom_smooth(method="lm") + facet_grid(Site~Year)


 
 ###
 ### BNPP:ANPP ratios - calculate, figures, and models
 ###
{
rs_ratios <- npp_master %>%
  filter(Year %in% 2019:2022) %>%
  filter(npp_type != "totNPP") %>%
  pivot_wider(names_from="npp_type", values_from="npp_gm2") %>%
  mutate(root_shoot = BNPP/ANPP) %>%
  mutate(ln_rs = log(root_shoot))


  rsratio_means <- rs_ratios %>%
    group_by(Year, Site, Drought) %>%
    summarize(
      rs_mean=mean(root_shoot, na.rm=T),
      rs_se = SE_function(root_shoot)
    )
  
# ggplot(filter(rsratio_means, Site=="TB"), aes(x=Drought, y=rs_mean, ymin=rs_mean-rs_se, ymax=rs_mean+rs_se)) +
#           geom_errorbar(width=2) +
#           geom_point(size=3) +
#           stat_smooth(method="lm",se=F) +
#           facet_grid(.~Year)
# 
# ggsave(file=paste0(write_dir,"figures\\tb root shoot ratios by year_",Sys.Date(),".png"), width=11, height=3.5)
# 
# ggplot(filter(rsratio_means, Site=="FK"), aes(x=Drought, y=rs_mean, ymin=rs_mean-rs_se, ymax=rs_mean+rs_se)) +
#   geom_errorbar(width=2) +
#   geom_point(size=3) +
#   stat_smooth(method="lm",se=F) +
#   facet_grid(.~Year)
# 
# ggsave(file=paste0(write_dir,"figures\\fk root shoot ratios by year_",Sys.Date(),".png"), width=11, height=3.5)
# 
# ggplot(rsratio_means, aes(x=Drought, y=rs_mean, ymin=rs_mean-rs_se, ymax=rs_mean+rs_se)) +
#   geom_errorbar(width=2) +
#   geom_point(size=3) +
#   stat_smooth(method="lm",se=F) +
#   facet_grid(Site~Year)
# ggsave(file=paste0(write_dir,"figures\\both sites root shoot ratios by year_",Sys.Date(),".png"), width=11, height=7)

rsratio_fig_short_y <- ggplot(rsratio_means, aes(x=Drought, y=rs_mean, ymin=rs_mean-rs_se, ymax=rs_mean+rs_se)) +
      geom_errorbar(width=2) +
      geom_point(size=3) +
      stat_smooth(method="lm",se=F) +
      facet_wrap(Site~Year, ncol=4) +
      ylim(0,10)
pdf(file=paste0(write_dir,"figures\\both sites root shoot ratios by year small ylim_",Sys.Date(),".pdf"), width=11, height=6.5)
print(rsratio_fig_short_y)
dev.off()

rsratio_fig_long_y <- ggplot(rsratio_means, aes(x=Drought, y=rs_mean, ymin=rs_mean-rs_se, ymax=rs_mean+rs_se)) +
  geom_errorbar(width=2) +
  geom_point(size=3) +
  stat_smooth(method="lm",se=F) +
  facet_wrap(Site~Year, ncol=4) +
  ylim(0,36)
pdf(file=paste0(write_dir,"figures\\both sites root shoot ratios by year large ylim_",Sys.Date(),".pdf"), width=11, height=6.5)
print(rsratio_fig_long_y)
dev.off()


### Running models
tb_rs_model_full <-   anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='ln_rs',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=filter(rs_ratios, Year %in% 2019:2022 & Site=="TB")
)




## just checking for normality and testing how similar model results are to anovat3 function 

### Thunder basin

tb_rs_model_noanovaT3fxn <- lme(log(root_shoot) ~ Year*Drought*Grazing
                                  , data=filter(rs_ratios, Year %in% 2019:2022 & Site=="FK")
                                  , random = ~1 |Block/Paddock/Plot
                                  , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                                  , control=lmeControl(returnObject=TRUE)
                                  , na.action = na.omit)
anova(tb_rs_model_noanovaT3fxn)
summary(tb_rs_model_noanovaT3fxn)
check_model(tb_rs_model_noanovaT3fxn)

### Split by year
# 2019
tb_rs_2019_lme <- lme(ln_rs ~ Drought*Grazing
                        , data=filter(rs_ratios, Year==2019 & Site=="FK")
                        , random = ~1 |Block/Paddock/Plot
                        , na.action = na.omit)

anova(tb_rs_2019_lme)
anova(tb_rs_2019_lme, type="marginal")
Anova(tb_rs_2019_lme, type=3)

# 2020
tb_rs_2020_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2020 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(tb_rs_2020_lme)
anova(tb_rs_2020_lme, type="marginal")
Anova(tb_rs_2020_lme, type=3)

# 2021
tb_rs_2021_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2021 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(tb_rs_2021_lme)
anova(tb_rs_2021_lme, type="marginal")
Anova(tb_rs_2021_lme, type=3)

# 2022
tb_rs_2022_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2022 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(tb_rs_2022_lme)
anova(tb_rs_2022_lme, type="marginal")
Anova(tb_rs_2022_lme, type=3)


### Create model table
tb_rs_model_out <- 
  as.data.frame(anova(tb_rs_2019_lme, type="marginal")) %>% mutate(Year=2019, Data="ANPP", Site="FK", Term=rownames(anova(tb_rs_2019_lme))) %>%
  bind_rows(
    as.data.frame(anova(tb_rs_2020_lme, type="marginal")) %>% mutate(Year=2020, Data="ANPP", Site="FK", Term=rownames(anova(tb_rs_2019_lme))),
    as.data.frame(anova(tb_rs_2021_lme, type="marginal")) %>% mutate(Year=2021, Data="ANPP", Site="FK", Term=rownames(anova(tb_rs_2019_lme))),
    as.data.frame(anova(tb_rs_2022_lme, type="marginal")) %>% mutate(Year=2022, Data="ANPP", Site="FK", Term=rownames(anova(tb_rs_2019_lme))),
  ) %>%
  dplyr::select(Data, Site, Year, Term, numDF, denDF, "F-value", "p-value")

write.csv(tb_rs_model_out, file=paste0(write_dir,"tables\\","tb rs lme model output by year_",Sys.Date(),".csv"), row.names=F)


### Fort Keogh

fk_rs_model_noanovaT3fxn <- lme(log(root_shoot) ~ Year*Drought*Grazing
                                , data=filter(rs_ratios, Year %in% 2019:2022 & Site=="FK")
                                , random = ~1 |Block/Paddock/Plot
                                , correlation=corCompSymm(form = ~1 |Block/Paddock/Plot)
                                , control=lmeControl(returnObject=TRUE)
                                , na.action = na.omit)
anova(fk_rs_model_noanovaT3fxn)
summary(fk_rs_model_noanovaT3fxn)
check_model(fk_rs_model_noanovaT3fxn)

### Split by year
# 2019
fk_rs_2019_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2019 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(fk_rs_2019_lme)
anova(fk_rs_2019_lme, type="marginal")
Anova(fk_rs_2019_lme, type=3)

# 2020
fk_rs_2020_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2020 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(fk_rs_2020_lme)
anova(fk_rs_2020_lme, type="marginal")
Anova(fk_rs_2020_lme, type=3)

# 2021
fk_rs_2021_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2021 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(fk_rs_2021_lme)
anova(fk_rs_2021_lme, type="marginal")
Anova(fk_rs_2021_lme, type=3)

# 2022
fk_rs_2022_lme <- lme(ln_rs ~ Drought*Grazing
                      , data=filter(rs_ratios, Year==2022 & Site=="FK")
                      , random = ~1 |Block/Paddock/Plot
                      , na.action = na.omit)

anova(fk_rs_2022_lme)
anova(fk_rs_2022_lme, type="marginal")
Anova(fk_rs_2022_lme, type=3)


### Create model table
fk_rs_model_out <- 
  as.data.frame(anova(fk_rs_2019_lme, type="marginal")) %>% mutate(Year=2019, Data="ANPP", Site="FK", Term=rownames(anova(fk_rs_2019_lme))) %>%
  bind_rows(
    as.data.frame(anova(fk_rs_2020_lme, type="marginal")) %>% mutate(Year=2020, Data="ANPP", Site="FK", Term=rownames(anova(fk_rs_2019_lme))),
    as.data.frame(anova(fk_rs_2021_lme, type="marginal")) %>% mutate(Year=2021, Data="ANPP", Site="FK", Term=rownames(anova(fk_rs_2019_lme))),
    as.data.frame(anova(fk_rs_2022_lme, type="marginal")) %>% mutate(Year=2022, Data="ANPP", Site="FK", Term=rownames(anova(fk_rs_2019_lme))),
  ) %>%
  dplyr::select(Data, Site, Year, Term, numDF, denDF, "F-value", "p-value")

write.csv(fk_rs_model_out, file=paste0(write_dir,"tables\\","fk rs lme model output by year_",Sys.Date(),".csv"), row.names=F)





}

###
### ANPP by functional groups - models and plotting
###
{
  anpp_fxn_grps <- anpp_fk_plot_means %>%
    dplyr::select(-dead, -litter, -subshrub_shrub, -ANPP_gm2) %>%
    pivot_longer(cols=C4P_gm2:Forb_gm2,names_to="fxn_type") %>%
  bind_rows(
    anpp_tb_plot_means %>%
      dplyr::select(-StandingDead_gm2, -Vulpia_gm2, -Bromus_gm2, -ANPP_gm2) %>%
      pivot_longer(cols=C3P_gm2:AnnualGrass_gm2,names_to="fxn_type")
  ) %>%
    rename(ANPP_gm2=value)
  
  ## Running models for different functional groups (NEED TO PUT THESE IN A TABLE)
  # Thunder basin overall models
  
  fxn_vec <- unique(anpp_fxn_grps$fxn_type)
  
  for(FXN in 1:length(fxn_vec)){
    df_temp
    
  }
  fk_C3P_full_model <- lme(ANPP_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                           , data=filter(anpp_fxn_grps, Year != 2018 & fxn_type=="C3P_gm2" & Site=="FK")
                           , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(fk_C3P_full_model, type="marginal")

  tb_C3P_full_model <- lme(C3P_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                               , data=filter(anpp_tb_plot_means, Year != 2018)
                               , random = ~1 |Block/Paddock
                               , na.action = na.omit)
  anova(tb_C3P_full_model, type="marginal")
  
  ### Testing anova_t3 function again -- similar output but main effects are weaker here -- potentially because den df are reduced by including the 3 way interactions
  anpp_tb_plot_means$Year <- as.factor(anpp_tb_plot_means$Year)
  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
           DepVar='C3P_gm2',
           RndForm='~1 |Block/Paddock',
           Data=anpp_tb_plot_means
  )
  
  tb_C4P_full_model <- lme(C4P_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                           , data=filter(anpp_tb_plot_means, Year != 2018)
                           , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(tb_C4P_full_model, type="marginal")
  
  tb_Forb_full_model <- lme(Forb_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                            , data=filter(anpp_tb_plot_means, Year != 2018)
                            , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(tb_Forb_full_model, type="marginal")
  
  tb_AnnualGrass_full_model <- lme(AnnualGrass_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                            , data=filter(anpp_tb_plot_means, Year != 2018)
                            , random = ~1 |Block/Paddock
                            , na.action = na.omit)
  anova(tb_AnnualGrass_full_model, type="marginal")
  
  
  ggplot(anpp_tb_plot_means, aes(x=Drought, y=AnnualGrass_gm2, col=Grazing)) +
    geom_point() +
    geom_smooth(method="lm",se=F) +
    facet_wrap(~Year)

  # Fort Keogh overall models
  fk_C3P_full_model <- lme(ANPP_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                           , data=filter(anpp_fxn_grps, Year != 2018 & fxn_type=="C3P_gm2" & Site=="FK")
                           , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(fk_C3P_full_model, type="marginal")
  
  fk_C4P_full_model <- lme(ANPP_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                           , data=filter(anpp_fxn_grps, Year != 2018 & fxn_type=="C4P_gm2" & Site=="FK")
                           , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(fk_C4P_full_model, type="marginal")

  fk_AnnualGrass_full_model <- lme(ANPP_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                           , data=filter(anpp_fxn_grps, Year != 2018 & fxn_type=="AnnualGrass_gm2" & Site=="FK")
                           , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(fk_AnnualGrass_full_model, type="marginal")
  
  fk_Forb_full_model <- lme(ANPP_gm2 ~ Drought*Grazing + Drought*as.factor(Year) + Grazing*as.factor(Year)
                           , data=filter(anpp_fxn_grps, Year != 2018 & fxn_type=="Forb_gm2" & Site=="FK")
                           , random = ~1 |Block/Paddock
                           , na.action = na.omit)
  anova(fk_Forb_full_model, type="marginal")
  
  
  
  
  ## Plotting functional group rel cov
  anpp_relanpp_fxn <- anpp_fxn_grps %>%
    group_by(Year, Site, Block, Paddock, Drought, Grazing, Plot) %>%
    mutate(rel_ANPP = ANPP_gm2/sum(ANPP_gm2))

  anpp_relanpp_fxn_drt_means <- anpp_fxn_grps %>%
    group_by(Year, Site, Drought, fxn_type) %>%
    summarize(ANPP_gm2 = mean(ANPP_gm2, na.rm=T)) %>%
    group_by(Year, Site, Drought) %>%
    mutate(rel_ANPP = ANPP_gm2/sum(ANPP_gm2))

  fxn_grp_plot <- ggplot(anpp_relanpp_fxn_drt_means, aes(x=factor(Drought), y=rel_ANPP, fill=fxn_type)) +
    geom_col() +
    theme_few() +
    facet_grid(Site ~ Year)

  pdf(file=paste0(write_dir, "figures\\fxn group rel abundance plot", Sys.Date(),".pdf"), width=10, height=3.5)
  print(fxn_grp_plot)
  dev.off()

  
  annual_plot_means <- anpp_fxn_grps %>%
    group_by(Year, Site, Drought, fxn_type) %>%
    summarize_at(vars(ANPP_gm2), .funs=list(mean = mean, se = SE_function)) %>%
    ungroup()
  
  
  ggplot(filter(annual_plot_means, fxn_type=="AnnualGrass_gm2"&Site=="FK"), aes(x=Drought, y=mean, ymin=mean-se, ymax=mean+se)) +
    geom_errorbar(width=1) +
    geom_point(size=2) +
    theme_few() +
    facet_grid(.~Year)

  ggplot(filter(annual_plot_means, fxn_type=="C3P_gm2"&Site=="FK"), aes(x=Drought, y=mean, ymin=mean-se, ymax=mean+se)) +
    geom_errorbar(width=1) +
    geom_point(size=2) +
    theme_few() +
    facet_grid(.~Year)

  ggplot(filter(annual_plot_means, fxn_type=="C4P_gm2"&Site=="FK"), aes(x=Drought, y=mean, ymin=mean-se, ymax=mean+se)) +
    geom_errorbar(width=1) +
    geom_point(size=2) +
    theme_few() +
    facet_grid(.~Year)

  ggplot(filter(annual_plot_means, fxn_type=="Forb_gm2"&Site=="FK"), aes(x=Drought, y=mean, ymin=mean-se, ymax=mean+se)) +
    geom_errorbar(width=1) +
    geom_point(size=2) +
    theme_few() +
    facet_grid(.~Year)

  ggplot(filter(annual_plot_means, fxn_type=="AnnualGrass_gm2"&Site=="FK"), aes(x=Drought, y=mean)) +
    geom_jitter(width=2) +
    theme_few() +
    facet_grid(.~Year)
  
  ggplot(filter(annual_plot_means, Site=="FK" & Year==2023), aes(x=Drought, y=mean, ymin=mean-se, ymax=mean+se)) +
    geom_errorbar(width=1) +
    geom_point(size=2) +
    theme_few() 

  ggplot(filter(annual_plot_means, Site=="FK" & Year==2023), aes(x=Drought, y=mean, ymin=mean-se, ymax=mean+se)) +
    geom_errorbar(width=1) +
    geom_point(size=2) +
    theme_few() 
  
    with(filter(anpp_fxn_grps, fxn_type=="AnnualGrass_gm2"&Site=="FK"), table(Year, Drought))
    
    
### Calculating Functional group response ratios
    anpp_fxn_controls <- anpp_fxn_grps %>%
      filter(Drought==0) %>%
      group_by(Site, Year, Block, Paddock, fxn_type) %>%
      summarize(ctrl_anpp=mean(ANPP_gm2, na.rm=T)
      )
    
    anpp_fxn_rr <- anpp_fxn_grps %>%
      dplyr::select(-slope_simple) %>%
      filter(Drought != 0) %>%
      full_join(anpp_fxn_controls, by=c("Year","Site","Block","Paddock", "fxn_type")) %>%
      mutate(lnrr_anpp=log((ANPP_gm2+0.01)/(ctrl_anpp+0.01)), 
             pchange_anpp=((ANPP_gm2+0.01)-(ctrl_anpp+0.01))/(ctrl_anpp+0.01),
             diff_anpp = ANPP_gm2-ctrl_anpp) 

  ## quick plots  
    ggplot(filter(anpp_fxn_rr,fxn_type=="C4P_gm2"), aes(x=Drought, y=lnrr_anpp)) +
      geom_hline(yintercept=0, col="grey") +
      geom_point() +
      geom_smooth(method="lm",se=F) +
      facet_grid(Site~Year)
    ggplot(filter(anpp_fxn_rr,fxn_type=="C3P_gm2"), aes(x=Drought, y=diff_anpp)) +
      geom_hline(yintercept=0, col="grey") +
      geom_point() +
      geom_smooth(method="lm",se=F) +
      facet_grid(Site~Year)
    ggplot(filter(anpp_fxn_rr,fxn_type=="Forb_gm2"), aes(x=Drought, y=lnrr_anpp)) +
      geom_hline(yintercept=0, col="grey") +
      geom_point() +
      geom_smooth(method="lm",se=F) +
      facet_grid(Site~Year)
    ggplot(filter(anpp_fxn_rr,fxn_type=="AnnualGrass_gm2"), aes(x=Drought, y=diff_anpp)) +
      geom_hline(yintercept=0, col="grey") +
      geom_point() +
      geom_smooth(method="lm",se=F) +
      facet_grid(Site~Year)
    
    # Calculate means and se for response ratios and difference measures
    anpp_rr_fxn_means <- anpp_fxn_rr %>%
      group_by(Year, Site, Drought, fxn_type) %>%
      summarize_at(vars(lnrr_anpp, pchange_anpp, diff_anpp), .funs=list(mean = mean, se = SE_function)) %>%
      ungroup()
    
    ## Plot means and ses
    ggplot(filter(anpp_rr_fxn_means, fxn_type=="AnnualGrass_gm2"&Year!=2018), aes(x=Drought, y=diff_anpp_mean, ymin=diff_anpp_mean-diff_anpp_se, ymax=diff_anpp_mean+diff_anpp_se)) +
      geom_hline(yintercept=0, col="grey") +
      geom_errorbar(width=1) +
      geom_point(size=2) +
      geom_smooth(method="lm",se=F) +
      theme_few() +
      facet_grid(Site~Year)
    
  fxn_diff_plot <-  ggplot(filter(anpp_rr_fxn_means,Year!=2018), aes(x=Drought, y=diff_anpp_mean, ymin=diff_anpp_mean-diff_anpp_se, ymax=diff_anpp_mean+diff_anpp_se, col=fxn_type, fill=fxn_type, pch=fxn_type)) +
            geom_hline(yintercept=0, col="grey") +
            geom_errorbar(width=1) +
            geom_point(size=2, col="black") +
            geom_smooth(method="lm",se=F) +
            scale_shape_manual(values=21:24) +
            theme_few() +
            facet_grid(Site~Year)
  
  pdf(paste0(write_dir,"figures//anpp diff by functional group_",Sys.Date(),".pdf"), useDingbats = F, width=11.5, height=4.25)
  print(fxn_diff_plot)
  dev.off()
}

# response ratio
ggplot(filter(anpp_rr_fxn_means,Year!=2018), aes(x=Drought, y=lnrr_anpp_mean, ymin=lnrr_anpp_mean-lnrr_anpp_se, ymax=lnrr_anpp_mean+lnrr_anpp_se, col=fxn_type, fill=fxn_type, pch=fxn_type)) +
  geom_hline(yintercept=0, col="grey") +
  geom_errorbar(width=1) +
  geom_point(size=2, col="black") +
  geom_smooth(method="lm",se=F) +
  scale_shape_manual(values=21:24) +
  theme_few() +
  facet_grid(Site~Year)

# percent change
ggplot(filter(anpp_rr_fxn_means,Year!=2018), aes(x=Drought, y=pchange_anpp_mean, ymin=pchange_anpp_mean-pchange_anpp_se, ymax=pchange_anpp_mean+pchange_anpp_se, col=fxn_type, fill=fxn_type, pch=fxn_type)) +
  geom_hline(yintercept=0, col="grey") +
  geom_errorbar(width=1) +
  geom_point(size=2, col="black") +
  geom_smooth(method="lm",se=F) +
  scale_shape_manual(values=21:24) +
  theme_few() +
  ylim(-200,200) +
  facet_grid(Site~Year)
### Calculate and plot compensation metric (annuals compensating for C3 reductions)
anpp_fxn_grps <- anpp_fxn_grps %>%
  mutate(Block_Paddock=paste(Block, Paddock, sep="_"))

install.packages("groupedstats")
library(groupedstats)

fxn_grp_drt_slopes <- anpp_fxn_grps %>% 
  group_by(Site, Year, Block, Paddock, fxn_type) %>%
  do(tidy(lm(ANPP_gm2 ~ Drought, .))) %>%
  filter(term=="Drought") %>%
  ungroup()

ggplot(filter(fxn_grp_drt_slopes,fxn_type=="AnnualGrass_gm2"), aes(x=estimate))+
  geom_histogram() +
  facet_wrap(Site~Year) +
  xlim(-1,1) +
  geom_vline(xintercept=0)

ggplot(filter(fxn_grp_drt_slopes,fxn_type=="C3P_gm2"), aes(x=estimate))+
  geom_histogram() +
  facet_wrap(Site~Year) +
  xlim(-1,1) +
  geom_vline(xintercept=0)

compensation_df <- fxn_grp_drt_slopes %>%
  dplyr::select(-std.error, -statistic, -p.value) %>%
  pivot_wider(names_from=fxn_type, values_from=estimate) %>%
  mutate(ann_c3p_comp = (AnnualGrass_gm2/C3P_gm2)*(-1)) %>%
  mutate(high_low_comp = (pmax(AnnualGrass_gm2,C3P_gm2,C4P_gm2,Forb_gm2)/
           pmin(AnnualGrass_gm2,C3P_gm2,C4P_gm2,Forb_gm2))*(-1))

compensation_means <- compensation_df %>%
  group_by(Site, Year) %>%
  summarize_at(.vars=vars("ann_c3p_comp","high_low_comp"), .funs=list(mean=mean, se=SE_function)) %>%
  ungroup()

ann_c3p_comp_fig <- ggplot(filter(compensation_means,Year!=2018), 
                   aes(Year, ann_c3p_comp_mean, ymin=ann_c3p_comp_mean-ann_c3p_comp_se, 
                       ymax=ann_c3p_comp_mean+ann_c3p_comp_se)) +
  geom_hline(yintercept=1, col="grey") +
  geom_path() +
  geom_errorbar(width=0.1) +
  geom_point(size=2) +
  facet_wrap(~Site) +
  theme_few() +
  ylab("Fxn Group Compensation -(mAnn/mC3P)")

high_low_comp_fig <- ggplot(filter(compensation_means,Year %in% 2021:2023), 
                           aes(Year, high_low_comp_mean, ymin=high_low_comp_mean-high_low_comp_se, 
                               ymax=high_low_comp_mean+high_low_comp_se)) +
  geom_hline(yintercept=1, col="grey") +
  geom_path() +
  geom_errorbar(width=0.1) +
  geom_point(size=2) +
  facet_wrap(~Site) +
  theme_few() +
  ylab("Fxn Group Compensation -(mhigh/mlow)")

pdf(file=paste0(write_dir,"figures\\compensation fig_",Sys.Date(),".pdf"), width=6, height=3.5, useDingbats = F)
print(comp_fig)
dev.off()


ggplot(compensation_df, aes(Year, ann_c3p_comp)) +
  geom_point() +
  facet_wrap(~Site) +
#  ylim(-10,10) +
  geom_hline(yintercept=0)

ggplot(filter(anpp_fxn_grps, fxn_type=="AnnualGrass_gm2" & Site=="FK"), aes(x=Drought, y=ANPP_gm2, col=Block_Paddock)) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(.~Year) +
  ylim(0,150)

ggplot(filter(fxn_anpp_compens, Site=="FK"), aes(x=Drought, y=c3_ann_comp)) +
  geom_jitter(width=2) +
  theme_few() +
  facet_grid(.~Year)

hist(log(fxn_anpp_compens$c3_ann_comp))




### Simulation for assessing how compensation metric varies across a gradient of slopes
max(abs(compensation_df$AnnualGrass_gm2))

sim_master <- {}
slope_vec <- seq(-2.2,2.2,by=.1)

for(SLOPE in 1:length(slope_vec)){
sim_temp <- data.frame(
  slope_x = slope_vec[SLOPE],
  slope_y = slope_vec
  )
sim_master <- rbind(sim_master, sim_temp)
rm(sim_temp)
}

sim_master <- sim_master %>%
  mutate(compensation_metric = -(slope_y/slope_x))

ggplot(sim_master, aes(x=slope_x, y=slope_y)) +
  geom_tile(aes(fill=compensation_metric)) +
#  scale_fill_gradient(low = "white", high = "orange", name="Compensation")
  scale_fill_gradient(low = "blue", high = "yellow", limits=c(-5,5), oob=scales::squish, name="Compensation")



### Calculate means by drought only
anpp_tb_drt_means <- anpp_tb_plot_means %>%
  group_by(Year, Site, Drought) %>%
  summarize_at(vars(c(C3P_gm2:ANPP_gm2, AnnualGrass_gm2)),.funs = c(mean, SE_function),na.rm=TRUE) %>% 
  rename(C3P_mean=C3P_gm2_fn1, C3P_se=C3P_gm2_fn2,
         C4P_mean=C4P_gm2_fn1, C4P_se=C4P_gm2_fn2,
         Forb_mean=Forb_gm2_fn1, Forb_se=Forb_gm2_fn2,
         ANPP_mean=ANPP_gm2_fn1, ANPP_se=ANPP_gm2_fn2,
         AnnualGrass_mean=AnnualGrass_gm2_fn1, AnnualGrass_se=AnnualGrass_gm2_fn2)

ggplot(anpp_tb_drt_means, aes(x=Drought, y=ANPP_mean, ymin=ANPP_mean-ANPP_se, ymax=ANPP_mean+ANPP_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free_y") +
  xlim(0,100) +
  ylab("ANPP g/m2") + xlab("Drought magnitude (% ppt reduction)")
ggsave("..\\..\\figures\\ANPP\\tb_TotANPP_18-23.png", width=9, height=3, units="in")

### Plots of different functional groups
ggplot(anpp_tb_drt_means, aes(x=Drought, y=C3P_mean, ymin=C3P_mean-C3P_se, ymax=C3P_mean+C3P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("C3 Grass (g/m2)")
ggsave("..\\..\\figures\\ANPP\\tb_C3Grass_18-22.png", width=9, height=3, units="in")

ggplot(anpp_tb_drt_means, aes(x=Drought, y=C4P_mean, ymin=C4P_mean-C4P_se, ymax=C4P_mean+C4P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("C4 Grass (g/m2)")
ggsave("..\\..\\figures\\ANPP\\tb_c4Grass_18-22.png", width=9, height=3, units="in")

ggplot(anpp_tb_drt_means, aes(x=Drought, y=Forb_mean, ymin=Forb_mean-Forb_se, ymax=Forb_mean+Forb_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("Forbs (g/m2)")
ggsave("..\\..\\figures\\ANPP\\tb_Forb_18-22.png", width=9, height=3, units="in")

ggplot(anpp_tb_drt_means, aes(x=Drought, y=AnnualGrass_mean, ymin=AnnualGrass_mean-AnnualGrass_se, ymax=AnnualGrass_mean+AnnualGrass_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("Annual Grass (g/m2)")
ggsave("..\\..\\figures\\ANPP\\tb_annualGrass_18-22.png", width=9, height=3, units="in")

###
### Calculate and plot means by drought and grazing treatments
###
anpp_tb_drt_grz_means <- anpp_tb_plot_means %>%
  group_by(Year, Site, Drought, Grazing) %>%
  summarize_at(vars(c(C3P_gm2:ANPP_gm2, AnnualGrass_gm2)),.funs = c(mean, SE_function),na.rm=TRUE) %>% 
  rename(C3P_mean=C3P_gm2_fn1, C3P_se=C3P_gm2_fn2,
         C4P_mean=C4P_gm2_fn1, C4P_se=C4P_gm2_fn2,
         Forb_mean=Forb_gm2_fn1, Forb_se=Forb_gm2_fn2,
         ANPP_mean=ANPP_gm2_fn1, ANPP_se=ANPP_gm2_fn2,
         AnnualGrass_mean=AnnualGrass_gm2_fn1, AnnualGrass_se=AnnualGrass_gm2_fn2)

anpp_tb_drt_grz_means$Grazing <- factor(anpp_tb_drt_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))


ggplot(anpp_tb_drt_grz_means, aes(x=Drought, y=ANPP_mean, ymin=ANPP_mean-ANPP_se, ymax=ANPP_mean+ANPP_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("ANPP g/m2")

ggsave("..\\..\\figures\\ANPP\\tb_drtbygrz_ANPP_18-23.png", width=9, height=7, units="in")

ggplot(anpp_tb_drt_grz_means, aes(x=Drought, y=C3P_mean, ymin=C3P_mean-C3P_se, ymax=C3P_mean+C3P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("C3 ANPP g/m2")
ggsave("..\\..\\figures\\ANPP\\tb_drtbygrz_c3_18-23.png", width=9, height=7, units="in")

ggplot(anpp_tb_drt_grz_means, aes(x=Drought, y=C4P_mean, ymin=C4P_mean-C4P_se, ymax=C4P_mean+C4P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("C4 ANPP g/m2")

ggsave("..\\..\\figures\\ANPP\\tb_drtbygrz_c4_18-23.png", width=9, height=7, units="in")


ggplot(anpp_tb_drt_grz_means, aes(x=Drought, y=Forb_mean, ymin=Forb_mean-Forb_se, ymax=Forb_mean+Forb_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("Forb ANPP g/m2")

ggsave("..\\..\\figures\\ANPP\\tb_drtbygrz_Forbs_18-23.png", width=9, height=7, units="in")

ggplot(anpp_tb_drt_grz_means, aes(x=Drought, y=AnnualGrass_mean, ymin=AnnualGrass_mean-AnnualGrass_se, ymax=AnnualGrass_mean+AnnualGrass_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("AnGrass ANPP g/m2")
ggsave("..\\..\\figures\\ANPP\\tb_drtbygrz_AnnGrass_18-23.png", width=9, height=7, units="in")
}

###
### FK - Data cleaning and plots
###
{
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
  left_join(dplyr::select(Plot_Key, Site, PlotID, Block, Paddock, Drought, Grazing), 
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
  group_by(Year, Site, Block, Paddock, Drought, Grazing, PlotID) %>%
  summarize_at(vars(c(C4P_gm2:litter, ANPP_gm2)),.funs = c(mean),na.rm=TRUE)

### Check for missing data and that data values are reasonable -- LOOKS GOOD
with(anpp_fk_plot_means, table(Year, Drought, Grazing))
with(anpp_fk_plot_means, hist(ANPP_gm2))


### Plot to take a quick look
ggplot(anpp_fk_18_23, aes(x=Drought, y=ANPP_gm2)) +
  geom_jitter(width=1, alpha=0.5) +
  facet_grid(rows=vars(Year), cols=vars(Grazing), scales="free")

ggplot(anpp_fk_18_23, aes(x=Drought, y=Forb_gm2)) +
  geom_jitter(width=1, alpha=0.5) +
  facet_grid(rows=vars(Year), cols=vars(Grazing), scales="free")

ggplot(anpp_fk_18_23, aes(x=Drought, y=AnnualGrass_gm2)) +
  geom_jitter(width=1, alpha=0.5) +
  facet_grid(rows=vars(Year), cols=vars(Grazing), scales="free")

### Calculate means by drought only
anpp_fk_drt_means <- anpp_fk_plot_means %>%
  group_by(Year, Site, Drought) %>%
  summarize_at(vars(c(C4P_gm2:litter, ANPP_gm2)),.funs = c(mean, SE_function),na.rm=TRUE) %>% 
  rename(C3P_mean=C3P_gm2_fn1, C3P_se=C3P_gm2_fn2,
         C4P_mean=C4P_gm2_fn1, C4P_se=C4P_gm2_fn2,
         Forb_mean=Forb_gm2_fn1, Forb_se=Forb_gm2_fn2,
         ANPP_mean=ANPP_gm2_fn1, ANPP_se=ANPP_gm2_fn2,
         AnnualGrass_mean=AnnualGrass_gm2_fn1, AnnualGrass_se=AnnualGrass_gm2_fn2)

### Plot and save to file
png(file="..\\..\\figures\\TBRI2022\\anpp drought only.png", width=12.6, height=4.2, units="in", res=600)
print(
ggplot(anpp_fk_drt_means, aes(x=Drought, y=ANPP_mean, ymin=ANPP_mean-ANPP_se, ymax=ANPP_mean+ANPP_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free_y") +
  xlim(0,100) +
  ylab("ANPP g/m2") + xlab("Drought magnitude (% ppt reduction)")
)
dev.off()

ggplot(anpp_fk_drt_means, aes(x=Drought, y=ANPP_mean, ymin=ANPP_mean-ANPP_se, ymax=ANPP_mean+ANPP_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free_y") +
  xlim(0,100) +
  ylab("ANPP g/m2") + xlab("Drought magnitude (% ppt reduction)")
ggsave("..\\..\\figures\\ANPP\\fk_TotANPP_18-23.png", width=9, height=3, units="in")

### Plots of different functional groups
ggplot(anpp_fk_drt_means, aes(x=Drought, y=C3P_mean, ymin=C3P_mean-C3P_se, ymax=C3P_mean+C3P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("C3 Grass (g/m2)")
ggsave("..\\..\\figures\\ANPP\\fk_C3Grass_18-23.png", width=9, height=3, units="in")

ggplot(anpp_fk_drt_means, aes(x=Drought, y=C4P_mean, ymin=C4P_mean-C4P_se, ymax=C4P_mean+C4P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("C4 Grass (g/m2)")
ggsave("..\\..\\figures\\ANPP\\fk_c4Grass_18-23.png", width=9, height=3, units="in")

ggplot(anpp_fk_drt_means, aes(x=Drought, y=Forb_mean, ymin=Forb_mean-Forb_se, ymax=Forb_mean+Forb_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("Forbs (g/m2)")
ggsave("..\\..\\figures\\ANPP\\fk_Forb_18-23.png", width=9, height=3, units="in")

ggplot(anpp_fk_drt_means, aes(x=Drought, y=AnnualGrass_mean, ymin=AnnualGrass_mean-AnnualGrass_se, ymax=AnnualGrass_mean+AnnualGrass_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(cols=vars(Year), scales="free") +
  xlim(0,100) + ylab("Annual Grass (g/m2)")
ggsave("..\\..\\figures\\ANPP\\fk_annualGrass_18-23.png", width=9, height=3, units="in")

###
### Calculate and plot means by drought and grazing treatments
###
anpp_fk_drt_grz_means <- anpp_fk_plot_means %>%
  group_by(Year, Site, Drought, Grazing) %>%
  summarize_at(vars(c(C4P_gm2:litter, ANPP_gm2)),.funs = c(mean, SE_function),na.rm=TRUE) %>% 
  rename(C3P_mean=C3P_gm2_fn1, C3P_se=C3P_gm2_fn2,
         C4P_mean=C4P_gm2_fn1, C4P_se=C4P_gm2_fn2,
         Forb_mean=Forb_gm2_fn1, Forb_se=Forb_gm2_fn2,
         ANPP_mean=ANPP_gm2_fn1, ANPP_se=ANPP_gm2_fn2,
         AnnualGrass_mean=AnnualGrass_gm2_fn1, AnnualGrass_se=AnnualGrass_gm2_fn2)

anpp_fk_drt_grz_means$Grazing <- factor(anpp_fk_drt_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))



ggplot(anpp_fk_drt_grz_means, aes(x=Drought, y=ANPP_mean, ymin=ANPP_mean-ANPP_se, ymax=ANPP_mean+ANPP_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("ANPP g/m2")

ggsave("..\\..\\figures\\ANPP\\fk_drtbygrz_ANPP_18-23.png", width=9, height=7, units="in")

ggplot(anpp_fk_drt_grz_means, aes(x=Drought, y=C3P_mean, ymin=C3P_mean-C3P_se, ymax=C3P_mean+C3P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("C3 ANPP g/m2")
ggsave("..\\..\\figures\\ANPP\\fk_drtbygrz_c3_18-23.png", width=9, height=7, units="in")

ggplot(anpp_fk_drt_grz_means, aes(x=Drought, y=C4P_mean, ymin=C4P_mean-C4P_se, ymax=C4P_mean+C4P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("C4 ANPP g/m2")

ggsave("..\\..\\figures\\ANPP\\fk_drtbygrz_c4_18-23.png", width=9, height=7, units="in")


ggplot(anpp_fk_drt_grz_means, aes(x=Drought, y=Forb_mean, ymin=Forb_mean-Forb_se, ymax=Forb_mean+Forb_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("Forb ANPP g/m2")

ggsave("..\\..\\figures\\ANPP\\fk_drtbygrz_Forbs_18-23.png", width=9, height=7, units="in")

ggplot(anpp_fk_drt_grz_means, aes(x=Drought, y=AnnualGrass_mean, ymin=AnnualGrass_mean-AnnualGrass_se, ymax=AnnualGrass_mean+AnnualGrass_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(method="lm",se=F) +
  facet_grid(rows=vars(Grazing), cols=vars(Year), scales="free") +
  xlim(0,100) +
  ylab("AnGrass ANPP g/m2")
ggsave("..\\..\\figures\\ANPP\\fk_drtbygrz_AnnGrass_18-23.png", width=9, height=7, units="in")


### Splines for c3p in 2022
ggplot(subset(anpp_fk_drt_grz_means,Year==2022), aes(x=Drought, y=C3P_mean, col=Grazing)) +
  geom_smooth(span=1, se=F) +
  xlim(0,80) +
  scale_color_manual(values=c("forestgreen", "gold","orangered2")) +
  ylab("C3 ANPP g/m2") + xlab("Drought intensity (% reduction)")
ggsave("..\\..\\figures\\ANPP\\fk_drtbygrz_c3_18-22.png", width=9, height=7, units="in")

ggplot(subset(anpp_fk_drt_grz_means,Year==2022), aes(x=Drought, y=C3P_mean, ymin=C3P_mean-C3P_se, ymax=C3P_mean+C3P_se)) +
  geom_errorbar(width=1) +
  geom_point() +
  geom_smooth(span=1, se=F) +
  facet_wrap(~Grazing) +
  xlim(0,80) +
  ylab("C3 ANPP g/m2")


}

###
### Response Ratios
###
{
###
### Thunder Basin
###
anpp_tb_controls <- anpp_tb_plot_means %>%
  filter(Drought==0) %>%
  group_by(Site, Year, Block, Paddock) %>%
  summarize(ctrl_ANPP=mean(ANPP_gm2, na.rm=T),
            ctrl_C3P=mean(C3P_gm2, na.rm=T),
            ctrl_C4P=mean(C4P_gm2, na.rm=T),
            ctrl_Forb=mean(Forb_gm2, na.rm=T),
            ctrl_StandingDead=mean(StandingDead_gm2, na.rm=T),
            ctrl_AnnualGrass=mean(AnnualGrass_gm2, na.rm=T),
            ctrl_Vulpia=mean(Vulpia_gm2, na.rm=T),
            ctrl_Bromus=mean(Bromus_gm2, na.rm=T)
  )
  
  ## calculating anpp means for controls for site description
  anpp_tb_ctrl_means <- anpp_tb_controls %>%
    group_by(Site) %>%
    summarize_at(.vars=vars(ctrl_ANPP:ctrl_Bromus), .funs=list(mean=mean)) %>%
    ungroup() %>%
    dplyr::select(-Site)
  
  with(anpp_tb_ctrl_means, ctrl_C3P_mean/ctrl_ANPP_mean)
  with(anpp_tb_ctrl_means, ctrl_C4P_mean/ctrl_ANPP_mean)
  with(anpp_tb_ctrl_means, ctrl_Forb_mean/ctrl_ANPP_mean)
  with(anpp_tb_ctrl_means, ctrl_AnnualGrass_mean/ctrl_ANPP_mean)

### IMPORTANT!!! lnrr is not working ... need to add parentheses around the +.01
### I believe this got fixed at some point (KW comment Feb 2024)
  anpp_tb_rr <- anpp_tb_plot_means %>%
  filter(Drought != 0) %>%
  full_join(anpp_tb_controls, by=c("Year","Site","Block","Paddock")) %>%
  mutate(lnrr_ANPP=log((ANPP_gm2+0.01)/(ctrl_ANPP+0.01)), pchange_ANPP=((ANPP_gm2+0.01)-(ctrl_ANPP+0.01))/(ctrl_ANPP+0.01),
         lnrr_C3P=log((C3P_gm2+0.01)/(ctrl_C3P+0.01)), pchange_C3P=((C3P_gm2+0.01)-(ctrl_C3P+0.01))/(ctrl_C3P+0.01),
         lnrr_C4P=log((C4P_gm2+0.01)/(ctrl_C4P+0.01)), pchange_C4P=((C4P_gm2+0.01)-(ctrl_C4P+0.01))/(ctrl_C4P+0.01),
         lnrr_Forb=log((Forb_gm2+0.01)/(ctrl_Forb+0.01)), pchange_Forb=((Forb_gm2+0.01)-(ctrl_Forb+0.01))/(ctrl_Forb+0.01),
         lnrr_StandingDead=log((StandingDead_gm2+0.01)/(ctrl_StandingDead+0.01)), pchange_StandingDead=((StandingDead_gm2+0.01)-(ctrl_StandingDead+0.01))/(ctrl_StandingDead+0.01),
         lnrr_AnnualGrass=log((AnnualGrass_gm2+0.01)/(ctrl_AnnualGrass+0.01)), pchange_AnnualGrass=((AnnualGrass_gm2+0.01)-(ctrl_AnnualGrass+0.01))/(ctrl_AnnualGrass+0.01),
         lnrr_Vulpia=log((Vulpia_gm2+0.01)/(ctrl_Vulpia+0.01)), pchange_Vulpia=((Vulpia_gm2+0.01)-(ctrl_Vulpia+0.01))/(ctrl_Vulpia+0.01),
         lnrr_Bromus=log((Bromus_gm2+0.01)/(ctrl_Bromus+0.01)), pchange_Bromus=((Bromus_gm2+0.01)-(ctrl_Bromus+0.01))/(ctrl_Bromus+0.01)
        )

###
### Generating predictions of percent change of anpp at 99% drought from models

### Thunder Basin
### total ANPP
# full model just to check for whether transformations are necessary
tb_anpp_pchange_full_model <- lme(pchange_ANPP ~ Year*Drought + Drought*Grazing + Year*Grazing
                              , data = filter(anpp_tb_rr, Year!=2018)
                              , random = ~1|Block/Paddock
                              , correlation=corCompSymm(form = ~1|Block/Paddock)
                              , control=lmeControl(returnObject=TRUE)
                              , na.action = na.omit
)

anova(tb_anpp_pchange_full_model)      
anova_t3(IndVars=c('Year','Drought','Grazing'),
         DepVar='pchange_ANPP',
         RndForm='~1 |Block/Paddock',
         Data=filter(anpp_tb_rr, Year!=2018)
)  

check_model(tb_anpp_pchange_full_model)
##

### get predictions of npp at 99% drought for each year
# 2018
tb_anpp_pchange_2018_model <-   lme(pchange_ANPP ~ Drought*Grazing
                            , data = filter(anpp_tb_rr, Year==2018)
                            , random = ~1|Block/Paddock/Plot
)

Anova(tb_anpp_pchange_2018_model, type=3)    
check_model(tb_anpp_pchange_2018_model)    
tb_2018_drought_int <- tb_anpp_pchange_2018_model$coefficients$fixed[1]
tb_2018_drought_slope <- tb_anpp_pchange_2018_model$coefficients$fixed[2]
tb_2018_yhat_99 <- tb_2018_drought_slope*99+tb_2018_drought_int

ggplot(filter(anpp_tb_rr, Year==2018), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2019
tb_anpp_pchange_2019_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2019)
                                    , random = ~1|Block/Paddock/Plot
)

anova(tb_anpp_pchange_2019_model)    
check_model(tb_anpp_pchange_2019_model)    
tb_2019_drought_int <- tb_anpp_pchange_2019_model$coefficients$fixed[1]
tb_2019_drought_slope <- tb_anpp_pchange_2019_model$coefficients$fixed[2]
tb_2019_yhat_99 <- tb_2019_drought_slope*99+tb_2019_drought_int

anova_t3(IndVars=c('Drought','Grazing'),
         DepVar='pchange_ANPP',
         RndForm='~1 |Block/Paddock',
         Data=filter(anpp_tb_rr, Year==2019)
)  

ggplot(filter(anpp_tb_rr, Year==2019), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2020
tb_anpp_pchange_2020_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2020)
                                    , random = ~1|Block/Paddock
)

anova(tb_anpp_pchange_2020_model)    
check_model(tb_anpp_pchange_2020_model)    
tb_2020_drought_int <- tb_anpp_pchange_2020_model$coefficients$fixed[1]
tb_2020_drought_slope <- tb_anpp_pchange_2020_model$coefficients$fixed[2]
tb_2020_yhat_99 <- tb_2020_drought_slope*99+tb_2020_drought_int
anova_t3(IndVars=c('Drought','Grazing'),
         DepVar='pchange_ANPP',
         RndForm='~1 |Block/Paddock',
         Data=filter(anpp_tb_rr, Year==2020)
)  

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2021
tb_anpp_pchange_2021_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2021)
                                    , random = ~1|Block/Paddock
)

anova(tb_anpp_pchange_2021_model)    
check_model(tb_anpp_pchange_2021_model)    
tb_2021_drought_int <- tb_anpp_pchange_2021_model$coefficients$fixed[1]
tb_2021_drought_slope <- tb_anpp_pchange_2021_model$coefficients$fixed[2]
tb_2021_yhat_99 <- tb_2021_drought_slope*99+tb_2021_drought_int
anova_t3(IndVars=c('Drought','Grazing'),
         DepVar='pchange_ANPP',
         RndForm='~1 |Block/Paddock',
         Data=filter(anpp_tb_rr, Year==2021)
)  

ggplot(filter(anpp_tb_rr, Year==2021), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2022
tb_anpp_pchange_2022_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2022)
                                    , random = ~1|Block/Paddock/Plot
)

anova(tb_anpp_pchange_2022_model)    
check_model(tb_anpp_pchange_2022_model)    
tb_2022_drought_int <- tb_anpp_pchange_2022_model$coefficients$fixed[1]
tb_2022_drought_slope <- tb_anpp_pchange_2022_model$coefficients$fixed[2]
tb_2022_yhat_99 <- tb_2022_drought_slope*99+tb_2022_drought_int
anova_t3(IndVars=c('Drought','Grazing'),
         DepVar='pchange_ANPP',
         RndForm='~1 |Block/Paddock',
         Data=filter(anpp_tb_rr, Year==2022)
)  

ggplot(filter(anpp_tb_rr, Year==2022), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2023
tb_anpp_pchange_2023_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2023)
                                    , random = ~1|Block/Paddock/Plot
)

anova(tb_anpp_pchange_2023_model)    
check_model(tb_anpp_pchange_2023_model)    
tb_2023_drought_int <- tb_anpp_pchange_2023_model$coefficients$fixed[1]
tb_2023_drought_slope <- tb_anpp_pchange_2023_model$coefficients$fixed[2]
tb_2023_yhat_99 <- tb_2023_drought_slope*99+tb_2023_drought_int
anova_t3(IndVars=c('Drought','Grazing'),
         DepVar='pchange_ANPP',
         RndForm='~1 |Block/Paddock',
         Data=filter(anpp_tb_rr, Year==2023)
)  

ggplot(filter(anpp_tb_rr, Year==2023), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
tb_anpp_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "Total ANPP",
  drt_magnitude = 99,
  pchange = c(tb_2018_yhat_99,
              tb_2019_yhat_99,
              tb_2020_yhat_99,
              tb_2021_yhat_99,
              tb_2022_yhat_99)
  )

ggplot(filter(tb_anpp_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### C3P
# full model just to check for whether transformations are necessary
tb_C3P_pchange_full_model <- lme(pchange_C3P ~ Year*Drought*Grazing
                                  , data = anpp_tb_rr
                                  , random = ~1|Block/Paddock/Plot
                                  , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                  , control=lmeControl(returnObject=TRUE)
                                  , na.action = na.omit
)

Anova(tb_C3P_pchange_full_model, type=3)      
check_model(tb_C3P_pchange_full_model)

##

### get predictions of npp at 99% drought for each year
# 2018
tb_C3P_pchange_2018_model <-   lme(pchange_C3P ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2018)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(tb_C3P_pchange_2018_model, type=3)    
check_model(tb_C3P_pchange_2018_model)    
tb_2018_C3P_drought_int <- tb_C3P_pchange_2018_model$coefficients$fixed[1]
tb_2018_C3P_drought_slope <- tb_C3P_pchange_2018_model$coefficients$fixed[2]
tb_2018_C3P_yhat_99 <- tb_2018_C3P_drought_slope*99+tb_2018_C3P_drought_int

ggplot(filter(anpp_tb_rr, Year==2018), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2019
tb_C3P_pchange_2019_model <-   lme(pchange_C3P ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2019 & pchange_C3P < 6)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(tb_C3P_pchange_2019_model, type=3)    
check_model(tb_C3P_pchange_2019_model)    
tb_2019_C3P_drought_int <- tb_C3P_pchange_2019_model$coefficients$fixed[1]
tb_2019_C3P_drought_slope <- tb_C3P_pchange_2019_model$coefficients$fixed[2]
tb_2019_C3P_yhat_99 <- tb_2019_C3P_drought_slope*99+tb_2019_C3P_drought_int

ggplot(filter(anpp_tb_rr, Year==2019 & pchange_C3P<6), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2020
tb_C3P_pchange_2020_model <-   lme(pchange_C3P ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2020)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(tb_C3P_pchange_2020_model, type=3)    
check_model(tb_C3P_pchange_2020_model)    
tb_2020_C3P_drought_int <- tb_C3P_pchange_2020_model$coefficients$fixed[1]
tb_2020_C3P_drought_slope <- tb_C3P_pchange_2020_model$coefficients$fixed[2]
tb_2020_C3P_yhat_99 <- tb_2020_C3P_drought_slope*99+tb_2020_C3P_drought_int

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2021
tb_C3P_pchange_2021_model <-   lme(pchange_C3P ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2021)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(tb_C3P_pchange_2021_model, type=3)    
check_model(tb_C3P_pchange_2021_model)    
tb_2021_C3P_drought_int <- tb_C3P_pchange_2021_model$coefficients$fixed[1]
tb_2021_C3P_drought_slope <- tb_C3P_pchange_2021_model$coefficients$fixed[2]
tb_2021_C3P_yhat_99 <- tb_2021_C3P_drought_slope*99+tb_2021_C3P_drought_int

ggplot(filter(anpp_tb_rr, Year==2021), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2022
tb_C3P_pchange_2022_model <-   lme(pchange_C3P ~ Drought*Grazing
                                    , data = filter(anpp_tb_rr, Year==2022)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(tb_C3P_pchange_2022_model, type=3)    
check_model(tb_C3P_pchange_2022_model)    
tb_2022_C3P_drought_int <- tb_C3P_pchange_2022_model$coefficients$fixed[1]
tb_2022_C3P_drought_slope <- tb_C3P_pchange_2022_model$coefficients$fixed[2]
tb_2022_C3P_yhat_99 <- tb_2022_C3P_drought_slope*99+tb_2022_C3P_drought_int

ggplot(filter(anpp_tb_rr, Year==2022), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
tb_C3P_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "C3P",
  drt_magnitude = 99,
  pchange = c(tb_2018_C3P_yhat_99,
              tb_2019_C3P_yhat_99,
              tb_2020_C3P_yhat_99,
              tb_2021_C3P_yhat_99,
              tb_2022_C3P_yhat_99)
)

ggplot(filter(tb_C3P_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### C4P
# full model just to check for whether transformations are necessary
anpp_tb_rr$pchange_C4P_log <- log(anpp_tb_rr$pchange_C4P - min(anpp_tb_rr$pchange_C4P) +.01)

tb_C4P_pchange_full_model <- lme(pchange_C4P_log ~ Year*Drought*Grazing
                                 , data = filter(anpp_tb_rr, pchange_C4P<6)
                                 , random = ~1|Block/Paddock/Plot
                                 , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit
)
#filter(anpp_tb_rr, pchange_C4P >6) %>% dplyr::select(Year:Plot, C4P_gm2, ctrl_C4P, pchange_C4P)
Anova(tb_C4P_pchange_full_model, type=3)      
check_model(tb_C4P_pchange_full_model)
check_outliers(tb_C4P_pchange_full_model)

##

log(.47)
exp(-.7550226)

### get predictions of npp at 99% drought for each year
# 2018
tb_C4P_pchange_2018_model <-   lme(pchange_C4P_log ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2018 & pchange_C4P < 6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_C4P_pchange_2018_model, type=3)    
check_model(tb_C4P_pchange_2018_model)    
tb_2018_C4P_drought_int <- tb_C4P_pchange_2018_model$coefficients$fixed[1]
tb_2018_C4P_drought_slope <- tb_C4P_pchange_2018_model$coefficients$fixed[2]
tb_2018_C4P_yhat_99_notback <- tb_2018_C4P_drought_slope*99+tb_2018_C4P_drought_int
tb_2018_C4P_yhat_99 <- exp(tb_2018_C4P_yhat_99_notback) +min(anpp_tb_rr$pchange_C4P) +.01



#2019
tb_C4P_pchange_2019_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2019 & pchange_C4P < 6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_C4P_pchange_2019_model, type=3)    
check_model(tb_C4P_pchange_2019_model)    
tb_2019_C4P_drought_int <- tb_C4P_pchange_2019_model$coefficients$fixed[1]
tb_2019_C4P_drought_slope <- tb_C4P_pchange_2019_model$coefficients$fixed[2]
tb_2019_C4P_yhat_99 <- tb_2019_C4P_drought_slope*99+tb_2019_C4P_drought_int

ggplot(filter(anpp_tb_rr, Year==2019 & pchange_C4P<6), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

#2020
tb_C4P_pchange_2020_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2020)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_C4P_pchange_2020_model, type=3)    
check_model(tb_C4P_pchange_2020_model)    
tb_2020_C4P_drought_int <- tb_C4P_pchange_2020_model$coefficients$fixed[1]
tb_2020_C4P_drought_slope <- tb_C4P_pchange_2020_model$coefficients$fixed[2]
tb_2020_C4P_yhat_99 <- tb_2020_C4P_drought_slope*99+tb_2020_C4P_drought_int

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

#2021
tb_C4P_pchange_2021_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2021 & pchange_C4P<6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_C4P_pchange_2021_model, type=3)    
check_model(tb_C4P_pchange_2021_model)    
tb_2021_C4P_drought_int <- tb_C4P_pchange_2021_model$coefficients$fixed[1]
tb_2021_C4P_drought_slope <- tb_C4P_pchange_2021_model$coefficients$fixed[2]
tb_2021_C4P_yhat_99 <- tb_2021_C4P_drought_slope*99+tb_2021_C4P_drought_int

ggplot(filter(anpp_tb_rr, Year==2021), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

#2022
tb_C4P_pchange_2022_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2022)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_C4P_pchange_2022_model, type=3)    
check_model(tb_C4P_pchange_2022_model)
tb_2022_C4P_drought_int <- tb_C4P_pchange_2022_model$coefficients$fixed[1]
tb_2022_C4P_drought_slope <- tb_C4P_pchange_2022_model$coefficients$fixed[2]
tb_2022_C4P_yhat_99 <- tb_2022_C4P_drought_slope*99+tb_2022_C4P_drought_int

ggplot(filter(anpp_tb_rr, Year==2022), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
tb_C4P_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "C4P",
  drt_magnitude = 99,
  pchange = c(tb_2018_C4P_yhat_99,
              tb_2019_C4P_yhat_99,
              tb_2020_C4P_yhat_99,
              tb_2021_C4P_yhat_99,
              tb_2022_C4P_yhat_99)
)

ggplot(filter(tb_C4P_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### Forb
# full model just to check for whether transformations are necessary

tb_Forb_pchange_full_model <- lme(pchange_Forb ~ Year*Drought*Grazing
                                 , data = filter(anpp_tb_rr, pchange_Forb < 6)
                                 , random = ~1|Block/Paddock/Plot
                                 , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit
)

Anova(tb_Forb_pchange_full_model, type=3)      
check_model(tb_Forb_pchange_full_model)
hist(anpp_tb_rr$pchange_Forb)
##

### get predictions of npp at 99% drought for each year
# 2018
tb_Forb_pchange_2018_model <-   lme(pchange_Forb ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2018 & pchange_Forb<6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_Forb_pchange_2018_model, type=3)    
check_model(tb_Forb_pchange_2018_model)    
tb_2018_Forb_drought_int <- tb_Forb_pchange_2018_model$coefficients$fixed[1]
tb_2018_Forb_drought_slope <- tb_Forb_pchange_2018_model$coefficients$fixed[2]
tb_2018_Forb_yhat_99 <- tb_2018_Forb_drought_slope*99+tb_2018_Forb_drought_int

ggplot(filter(anpp_tb_rr, Year==2018 & pchange_Forb<6), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2019
tb_Forb_pchange_2019_model <-   lme(pchange_Forb ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2019 & pchange_Forb < 6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_Forb_pchange_2019_model, type=3)    
check_model(tb_Forb_pchange_2019_model)    
tb_2019_Forb_drought_int <- tb_Forb_pchange_2019_model$coefficients$fixed[1]
tb_2019_Forb_drought_slope <- tb_Forb_pchange_2019_model$coefficients$fixed[2]
tb_2019_Forb_yhat_99 <- tb_2019_Forb_drought_slope*99+tb_2019_Forb_drought_int

ggplot(filter(anpp_tb_rr, Year==2019 & pchange_Forb<6), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2020
tb_Forb_pchange_2020_model <-   lme(pchange_Forb ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2020)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_Forb_pchange_2020_model, type=3)    
check_model(tb_Forb_pchange_2020_model)    
tb_2020_Forb_drought_int <- tb_Forb_pchange_2020_model$coefficients$fixed[1]
tb_2020_Forb_drought_slope <- tb_Forb_pchange_2020_model$coefficients$fixed[2]
tb_2020_Forb_yhat_99 <- tb_2020_Forb_drought_slope*99+tb_2020_Forb_drought_int

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2021
tb_Forb_pchange_2021_model <-   lme(pchange_Forb ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2021)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_Forb_pchange_2021_model, type=3)    
check_model(tb_Forb_pchange_2021_model)    
tb_2021_Forb_drought_int <- tb_Forb_pchange_2021_model$coefficients$fixed[1]
tb_2021_Forb_drought_slope <- tb_Forb_pchange_2021_model$coefficients$fixed[2]
tb_2021_Forb_yhat_99 <- tb_2021_Forb_drought_slope*99+tb_2021_Forb_drought_int

ggplot(filter(anpp_tb_rr, Year==2021), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2022
tb_Forb_pchange_2022_model <-   lme(pchange_Forb ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2022)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_Forb_pchange_2022_model, type=3)    
check_model(tb_Forb_pchange_2022_model)    
tb_2022_Forb_drought_int <- tb_Forb_pchange_2022_model$coefficients$fixed[1]
tb_2022_Forb_drought_slope <- tb_Forb_pchange_2022_model$coefficients$fixed[2]
tb_2022_Forb_yhat_99 <- tb_2022_Forb_drought_slope*99+tb_2022_Forb_drought_int

ggplot(filter(anpp_tb_rr, Year==2022), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
tb_Forb_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "Forb",
  drt_magnitude = 99,
  pchange = c(tb_2018_Forb_yhat_99,
              tb_2019_Forb_yhat_99,
              tb_2020_Forb_yhat_99,
              tb_2021_Forb_yhat_99,
              tb_2022_Forb_yhat_99)
)

ggplot(filter(tb_Forb_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### AnnualGrass
# full model just to check for whether transformations are necessary
tb_AnnualGrass_pchange_full_model <- lme(pchange_AnnualGrass ~ Year*Drought*Grazing
                                 , data = filter(anpp_tb_rr, pchange_AnnualGrass<6)
                                 , random = ~1|Block/Paddock/Plot
                                 , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit
)

Anova(tb_AnnualGrass_pchange_full_model, type=3)      
check_model(tb_AnnualGrass_pchange_full_model)
##

### get predictions of npp at 99% drought for each year
# 2018
tb_AnnualGrass_pchange_2018_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2018 & pchange_AnnualGrass<6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_AnnualGrass_pchange_2018_model, type=3)    
check_model(tb_AnnualGrass_pchange_2018_model)    
tb_2018_AnnualGrass_drought_int <- tb_AnnualGrass_pchange_2018_model$coefficients$fixed[1]
tb_2018_AnnualGrass_drought_slope <- tb_AnnualGrass_pchange_2018_model$coefficients$fixed[2]
tb_2018_AnnualGrass_yhat_99 <- tb_2018_AnnualGrass_drought_slope*99+tb_2018_AnnualGrass_drought_int

ggplot(filter(anpp_tb_rr, Year==2018& pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2019
tb_AnnualGrass_pchange_2019_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2019 & pchange_AnnualGrass < 6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_AnnualGrass_pchange_2019_model, type=3)    
check_model(tb_AnnualGrass_pchange_2019_model)    
tb_2019_AnnualGrass_drought_int <- tb_AnnualGrass_pchange_2019_model$coefficients$fixed[1]
tb_2019_AnnualGrass_drought_slope <- tb_AnnualGrass_pchange_2019_model$coefficients$fixed[2]
tb_2019_AnnualGrass_yhat_99 <- tb_2019_AnnualGrass_drought_slope*99+tb_2019_AnnualGrass_drought_int

ggplot(filter(anpp_tb_rr, Year==2019 & pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2020
tb_AnnualGrass_pchange_2020_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2020& pchange_AnnualGrass<6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_AnnualGrass_pchange_2020_model, type=3)    
check_model(tb_AnnualGrass_pchange_2020_model)    
tb_2020_AnnualGrass_drought_int <- tb_AnnualGrass_pchange_2020_model$coefficients$fixed[1]
tb_2020_AnnualGrass_drought_slope <- tb_AnnualGrass_pchange_2020_model$coefficients$fixed[2]
tb_2020_AnnualGrass_yhat_99 <- tb_2020_AnnualGrass_drought_slope*99+tb_2020_AnnualGrass_drought_int

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2021
tb_AnnualGrass_pchange_2021_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2021& pchange_AnnualGrass<6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_AnnualGrass_pchange_2021_model, type=3)    
check_model(tb_AnnualGrass_pchange_2021_model)    
tb_2021_AnnualGrass_drought_int <- tb_AnnualGrass_pchange_2021_model$coefficients$fixed[1]
tb_2021_AnnualGrass_drought_slope <- tb_AnnualGrass_pchange_2021_model$coefficients$fixed[2]
tb_2021_AnnualGrass_yhat_99 <- tb_2021_AnnualGrass_drought_slope*99+tb_2021_AnnualGrass_drought_int

ggplot(filter(anpp_tb_rr, Year==2021& pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2022
tb_AnnualGrass_pchange_2022_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                   , data = filter(anpp_tb_rr, Year==2022& pchange_AnnualGrass<6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(tb_AnnualGrass_pchange_2022_model, type=3)    
check_model(tb_AnnualGrass_pchange_2022_model)    
tb_2022_AnnualGrass_drought_int <- tb_AnnualGrass_pchange_2022_model$coefficients$fixed[1]
tb_2022_AnnualGrass_drought_slope <- tb_AnnualGrass_pchange_2022_model$coefficients$fixed[2]
tb_2022_AnnualGrass_yhat_99 <- tb_2022_AnnualGrass_drought_slope*99+tb_2022_AnnualGrass_drought_int

ggplot(filter(anpp_tb_rr, Year==2022& pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
tb_AnnualGrass_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "AnnualGrass",
  drt_magnitude = 99,
  pchange = c(tb_2018_AnnualGrass_yhat_99,
              tb_2019_AnnualGrass_yhat_99,
              tb_2020_AnnualGrass_yhat_99,
              tb_2021_AnnualGrass_yhat_99,
              tb_2022_AnnualGrass_yhat_99)
)

ggplot(filter(tb_AnnualGrass_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

year_key <- data.frame(Year=2018:2022, trt_group=c("no_drt","drought","drought","no_drt","no_drt"))

tb_pchange_99_full <- bind_rows(
  tb_anpp_pchange_99,
  tb_C3P_pchange_99,
  tb_C4P_pchange_99,
  tb_AnnualGrass_pchange_99,
) %>%
  full_join(year_key, by="Year")

tb_pchange_99_full$variable <- factor(tb_pchange_99_full$variable, levels=c("Total ANPP", "C3P","C4P","AnnualGrass"))

ggplot(filter(tb_pchange_99_full,Year!=2018), aes(x=Year, y=pchange*100, fill=trt_group)) +
  geom_col(col='black') +
  scale_fill_manual(values=c("firebrick2","dodgerblue")) +
  ylab("Drought response (% change)") +
  facet_wrap(~variable)


### Just Drought
### Note: We need to deal with the infinity values calculated in here for lnRR (0s causing this), but I'm just using pchange currently so I'm ignoring
anpp_tb_rr_drt_means <- anpp_tb_rr %>% 
  group_by(Year, Drought) %>%
  summarize_at(vars(lnrr_ANPP:pchange_Bromus),.funs = c(mean,SE_function, median),na.rm=TRUE) %>% 
  rename(
         lnrr_ANPP_mean=lnrr_ANPP_fn1, lnrr_ANPP_se=lnrr_ANPP_fn2,lnrr_ANPP_median=lnrr_ANPP_fn3,
         pchange_ANPP_mean=pchange_ANPP_fn1, pchange_ANPP_se=pchange_ANPP_fn2,pchange_ANPP_median=pchange_ANPP_fn3,
         lnrr_C3P_mean=lnrr_C3P_fn1, lnrr_C3P_se=lnrr_C3P_fn2,lnrr_C3P_median=lnrr_C3P_fn3,
         pchange_C3P_mean=pchange_C3P_fn1, pchange_C3P_se=pchange_C3P_fn2,pchange_C3P_median=pchange_C3P_fn3,
         lnrr_C4P_mean=lnrr_C4P_fn1, lnrr_C4P_se=lnrr_C4P_fn2,lnrr_C4P_median=lnrr_C4P_fn3,
         pchange_C4P_mean=pchange_C4P_fn1, pchange_C4P_se=pchange_C4P_fn2,pchange_C4P_median=pchange_C4P_fn3,
         lnrr_Forb_mean=lnrr_Forb_fn1, lnrr_Forb_se=lnrr_Forb_fn2,lnrr_Forb_median=lnrr_Forb_fn3,
         pchange_Forb_mean=pchange_Forb_fn1, pchange_Forb_se=pchange_Forb_fn2,pchange_Forb_median=pchange_Forb_fn3,
         lnrr_StandingDead_mean=lnrr_StandingDead_fn1, lnrr_StandingDead_se=lnrr_StandingDead_fn2,lnrr_StandingDead_median=lnrr_StandingDead_fn3,
         pchange_StandingDead_mean=pchange_StandingDead_fn1, pchange_StandingDead_se=pchange_StandingDead_fn2,pchange_StandingDead_median=pchange_StandingDead_fn3,
         lnrr_AnnualGrass_mean=lnrr_AnnualGrass_fn1, lnrr_AnnualGrass_se=lnrr_AnnualGrass_fn2,lnrr_AnnualGrass_median=lnrr_AnnualGrass_fn3,
         pchange_AnnualGrass_mean=pchange_AnnualGrass_fn1, pchange_AnnualGrass_se=pchange_AnnualGrass_fn2,pchange_AnnualGrass_median=pchange_AnnualGrass_fn3,
         lnrr_Vulpia_mean=lnrr_Vulpia_fn1, lnrr_Vulpia_se=lnrr_Vulpia_fn2,lnrr_Vulpia_median=lnrr_Vulpia_fn3,
         pchange_Vulpia_mean=pchange_Vulpia_fn1, pchange_Vulpia_se=pchange_Vulpia_fn2,pchange_Vulpia_median=pchange_Vulpia_fn3,
         lnrr_Bromus_mean=lnrr_Bromus_fn1, lnrr_Bromus_se=lnrr_Bromus_fn2,lnrr_Bromus_median=lnrr_Bromus_fn3,
         pchange_Bromus_mean=pchange_Bromus_fn1, pchange_Bromus_se=pchange_Bromus_fn2, pchange_Bromus_median=pchange_Bromus_fn3
         )

png(file="..\\..\\figures\\ANPP\\anpp pchange drought only19-23.png", width=14.0, height=4.2, units="in", res=600)
print(
ggplot(filter(anpp_tb_rr_drt_means,Year!=2018), aes(x= Drought, y=pchange_ANPP_mean*100,
                             ymin=(pchange_ANPP_mean-pchange_ANPP_se)*100, ymax=(pchange_ANPP_mean+pchange_ANPP_se)*100 ))+
  geom_point(size=3) +
  geom_errorbar(width=0) +
  #  geom_smooth(method="lm",se=F, col="red") +
  geom_smooth(method="lm", col="red", se=F) +
  geom_hline(yintercept=0) +
  facet_grid(cols=vars(Year)) +
  ylab("Drought effect (% ANPP change)") +
  xlab("Drought magnitude (% ppt reduction)") +
  xlim(25,100)
)
dev.off()

## Plot only 99% drought treatments
ggplot(filter(anpp_tb_rr_drt_means, Drought==99&Year%in%2019:2023), aes(x=Year, y=pchange_ANPP_median*100, 
                                                      ymin=(pchange_ANPP_median-pchange_ANPP_se)*100, 
                                                      ymax=(pchange_ANPP_median+pchange_ANPP_se)*100)) +
  geom_col(col='black', fill='darkseagreen') +
  geom_errorbar(width=.05) +
  ylab("Drought response (% change)")
  
ggsave("..//..//figures//ANPP//tb_anpp_droughtResponse_99only_median-19-23.png", width=6, height=4.5, units="in")

ggplot(filter(anpp_tb_rr_drt_means, Drought==99&Year%in%2019:2023), aes(x=Year, y=pchange_C3P_median*100, 
                                                                        ymin=(pchange_C3P_median-pchange_C3P_se)*100, 
                                                                        ymax=(pchange_C3P_median+pchange_C3P_se)*100)) +
  geom_col(col='black', fill='darkseagreen') +
  geom_errorbar(width=.05) +
  ylab("Drought C3P response (% change)")

ggsave("..//..//figures//ANPP//tb_C3P_droughtResponse_99only_median.png", width=6, height=4.5, units="in")

ggplot(filter(anpp_tb_rr_drt_means, Drought==99&Year%in%2019:2023), aes(x=Year, y=pchange_C4P_median*100, 
                                                                        ymin=(pchange_C4P_median-pchange_C4P_se)*100, 
                                                                        ymax=(pchange_C4P_median+pchange_C4P_se)*100)) +
  geom_col(col='black', fill='darkseagreen') +
  geom_errorbar(width=.05) +
  ylab("Drought C4P response (% change)") +
  ylim(-100,100)

ggsave("..//..//figures//ANPP//tb_C4P_droughtResponse_99only_median.png", width=6, height=4.5, units="in")
ggplot(filter(anpp_tb_rr_drt_means, Drought==99&Year%in%2019:2023), aes(x=Year, y=pchange_Forb_median*100, 
                                                                        ymin=(pchange_Forb_median-pchange_Forb_se)*100, 
                                                                        ymax=(pchange_Forb_median+pchange_Forb_se)*100)) +
  geom_col(col='black', fill='darkseagreen') +
  geom_errorbar(width=.05) +
  ylab("Drought Forb response (% change)") + ylim(-100,100)

ggsave("..//..//figures//ANPP//tb_Forb_droughtResponse_99only_median.png", width=6, height=4.5, units="in")

ggplot(filter(anpp_tb_rr_drt_means, Drought==99&Year%in%2019:2023), aes(x=Year, y=pchange_AnnualGrass_median*100, 
                                                                        ymin=(pchange_AnnualGrass_median-pchange_AnnualGrass_se)*100, 
                                                                        ymax=(pchange_AnnualGrass_median+pchange_AnnualGrass_se)*100)) +
  geom_col(col='black', fill='darkseagreen') +
  geom_errorbar(width=.05) +
  ylab("Drought AnnualGrass response (% change)") + ylim(-200,200)

ggsave("..//..//figures//ANPP//tb_AnnualGrass_droughtResponse_99only_median.png", width=6, height=4.5, units="in")

### Grazing and Drought treatments
### Note: We need to deal with the infinity values calculated in here for lnRR (0s causing this), but I'm just using pchange currently so I'm ignoring
anpp_tb_rr_drt_grz_means <- anpp_tb_rr %>% 
  group_by(Year, Drought, Grazing) %>%
  summarize_at(vars(lnrr_ANPP:pchange_Bromus),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(
    lnrr_ANPP_mean=lnrr_ANPP_fn1, lnrr_ANPP_se=lnrr_ANPP_fn2,
    pchange_ANPP_mean=pchange_ANPP_fn1, pchange_ANPP_se=pchange_ANPP_fn2,
    lnrr_C3P_mean=lnrr_C3P_fn1, lnrr_C3P_se=lnrr_C3P_fn2,
    pchange_C3P_mean=pchange_C3P_fn1, pchange_C3P_se=pchange_C3P_fn2,
    lnrr_C4P_mean=lnrr_C4P_fn1, lnrr_C4P_se=lnrr_C4P_fn2,
    pchange_C4P_mean=pchange_C4P_fn1, pchange_C4P_se=pchange_C4P_fn2,
    lnrr_Forb_mean=lnrr_Forb_fn1, lnrr_Forb_se=lnrr_Forb_fn2,
    pchange_Forb_mean=pchange_Forb_fn1, pchange_Forb_se=pchange_Forb_fn2,
    lnrr_StandingDead_mean=lnrr_StandingDead_fn1, lnrr_StandingDead_se=lnrr_StandingDead_fn2,
    pchange_StandingDead_mean=pchange_StandingDead_fn1, pchange_StandingDead_se=pchange_StandingDead_fn2,
    lnrr_AnnualGrass_mean=lnrr_AnnualGrass_fn1, lnrr_AnnualGrass_se=lnrr_AnnualGrass_fn2,
    pchange_AnnualGrass_mean=pchange_AnnualGrass_fn1, pchange_AnnualGrass_se=pchange_AnnualGrass_fn2,
    lnrr_Vulpia_mean=lnrr_Vulpia_fn1, lnrr_Vulpia_se=lnrr_Vulpia_fn2,
    pchange_Vulpia_mean=pchange_Vulpia_fn1, pchange_Vulpia_se=pchange_Vulpia_fn2,
    lnrr_Bromus_mean=lnrr_Bromus_fn1, lnrr_Bromus_se=lnrr_Bromus_fn2,
    pchange_Bromus_mean=pchange_Bromus_fn1, pchange_Bromus_se=pchange_Bromus_fn2
  )

anpp_tb_rr_drt_grz_means$Grazing <- factor(anpp_tb_rr_drt_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))

q3 <- c("#00AD9A", "#9183E6", "#E16A86")

png(file="..\\..\\figures\\TBRI2022\\anpp pchange drought and grazing.png", width=13.6, height=4.2, units="in", res=600)
print(
ggplot(anpp_tb_rr_drt_grz_means, aes(x= Drought, y=pchange_ANPP_mean*100, col=Grazing,
                                     ymin=(pchange_ANPP_mean-pchange_ANPP_se)*100, ymax=(pchange_ANPP_mean+pchange_ANPP_se)*100 ))+
   geom_point(size=3) +
   geom_errorbar(width=1) +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  facet_grid(cols=vars(Year), scales="free") +
  scale_colour_manual(values=q3, name='Grazing\nTreatment') +
  ylab("Drought effect (% ANPP change)") +
  xlab("Drought magnitude (% ppt reduction)") +
  xlim(25,100) +ylim(-75, 80)
)
dev.off()

png(file="..\\..\\figures\\TBRI2022\\c3p pchange drought and grazing.png", width=13.6, height=4.2, units="in", res=600)
print(
  ggplot(anpp_tb_rr_drt_grz_means, aes(x= Drought, y=pchange_C3P_mean*100, col=Grazing,
                                       ymin=(pchange_C3P_mean-pchange_C3P_se)*100, ymax=(pchange_C3P_mean+pchange_C3P_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
    facet_grid(cols=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (%)") +
    xlab("Precipitation reduction (%)") +
    xlim(100,25)
)
dev.off()

png(file="..\\..\\figures\\TBRI2022\\c3p pchange drought and grazing 2021 only.png", width=6.2, height=4.2, units="in", res=600)
print(
  ggplot(filter(anpp_tb_rr_drt_grz_means,Year==2021), aes(x= Drought, y=pchange_C3P_mean*100, col=Grazing,
                                       ymin=(pchange_C3P_mean-pchange_C3P_se)*100, ymax=(pchange_C3P_mean+pchange_C3P_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
  #  facet_grid(cols=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (%)") +
    xlab("Precipitation reduction (%)") +
    xlim(100,24)
)
dev.off()

ggplot(anpp_tb_rr_drt_grz_means, aes(x= Drought, y=pchange_ANPP_mean*100, col=Grazing,
                                     ymin=(pchange_ANPP_mean-pchange_ANPP_se)*100, ymax=(pchange_ANPP_mean+pchange_ANPP_se)*100 ))+
  geom_point(size=3) +
  geom_errorbar(width=1) +
  geom_smooth(method="lm", se=F) +
  geom_hline(yintercept=0) +
  facet_grid(cols=vars(Year), scales="free") +
  scale_colour_manual(values=q3, name='Grazing\nTreatment') +
  ylab("Drought effect (% ANPP change)") +
  xlab("Drought magnitude (% ppt reduction)") +
  xlim(25,100) +ylim(-75, 80)
)

anpp_tb_rr_drt_grz_means$Grazing <- factor(anpp_tb_rr_drt_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))
ggplot(subset(anpp_tb_rr_drt_grz_means,Drought==99&Year==2021), aes(x=factor(Grazing), y=pchange_ANPP_mean*100, fill=factor(Grazing)))+
  geom_col(col="black") +
#  geom_hline(yintercept=0) +
  scale_fill_manual(values=c("forestgreen","gold","darkorange2")) +
  ylab("Drt Response % change") +
  xlab("Grazing scenario") +
  ylim(-50, 50)
)

## C4P and AnnualGrass not all that interesting
pdf(file="..\\..\\figures\\TBRI2022\\c4p pchange drought and grazing.pdf", width=6.5, height=12.6)
print(
  ggplot(anpp_tb_rr_drt_grz_means, aes(x= Drought, y=pchange_AnnualGrass_mean*100, col=Grazing,
                                       ymin=(pchange_AnnualGrass_mean-pchange_AnnualGrass_se)*100, ymax=(pchange_AnnualGrass_mean+pchange_AnnualGrass_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
    facet_grid(rows=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (%)") +
    xlab("Precipitation reduction (%)") +
    xlim(100,25) + ylim(-200,200)
)
dev.off()

###
### Fort Keogh response ratios
###
anpp_fk_controls <- anpp_fk_plot_means %>%
  filter(Drought==0) %>%
  group_by(Site, Year, Block, Grazing) %>%
  summarize(ctrl_ANPP=mean(ANPP_gm2, na.rm=T),
            ctrl_C3P=mean(C3P_gm2, na.rm=T),
            ctrl_C4P=mean(C4P_gm2, na.rm=T),
            ctrl_Forb=mean(Forb_gm2, na.rm=T),
            ctrl_dead=mean(dead, na.rm=T),
            ctrl_AnnualGrass=mean(AnnualGrass_gm2, na.rm=T)
  )
## calculating anpp means for controls for site description
anpp_fk_ctrl_means <- anpp_fk_controls %>%
  group_by(Site) %>%
  summarize_at(.vars=vars(ctrl_ANPP:ctrl_AnnualGrass), .funs=list(mean=mean)) %>%
  ungroup() %>%
  dplyr::select(-Site)

with(anpp_fk_ctrl_means, ctrl_C3P_mean/ctrl_ANPP_mean)
with(anpp_fk_ctrl_means, ctrl_C4P_mean/ctrl_ANPP_mean)
with(anpp_fk_ctrl_means, ctrl_Forb_mean/ctrl_ANPP_mean)
with(anpp_fk_ctrl_means, ctrl_AnnualGrass_mean/ctrl_ANPP_mean)

### IMPORTANT -- RR is not currently working, need toput parentheses around the ANPP_gm2+.01 and all other similar values because of order of math operations
anpp_fk_rr <- anpp_fk_plot_means %>%
  filter(Drought != 0) %>%
  full_join(anpp_fk_controls, by=c("Year","Site","Block","Grazing")) %>%
  mutate(lnrr_ANPP=log(ANPP_gm2+0.01/ctrl_ANPP+0.01), pchange_ANPP=((ANPP_gm2+0.01)-(ctrl_ANPP+0.01))/(ctrl_ANPP+0.01),
         lnrr_C3P=log(C3P_gm2+0.01/ctrl_C3P+0.01), pchange_C3P=((C3P_gm2+0.01)-(ctrl_C3P+0.01))/(ctrl_C3P+0.01),
         lnrr_C4P=log(C4P_gm2+0.01/ctrl_C4P+0.01), pchange_C4P=((C4P_gm2+0.01)-(ctrl_C4P+0.01))/(ctrl_C4P+0.01),
         lnrr_Forb=log(Forb_gm2+0.01/ctrl_Forb+0.01), pchange_Forb=((Forb_gm2+0.01)-(ctrl_Forb+0.01))/(ctrl_Forb+0.01),
         lnrr_dead=log(dead+0.01/ctrl_dead+0.01), pchange_dead=((dead+0.01)-(ctrl_dead+0.01))/(ctrl_dead+0.01),
         lnrr_AnnualGrass=log(AnnualGrass_gm2+0.01/ctrl_AnnualGrass+0.01), pchange_AnnualGrass=((AnnualGrass_gm2+0.01)-(ctrl_AnnualGrass+0.01))/(ctrl_AnnualGrass+0.01)
  )

###
### Generating predictions of percent change of anpp at 99% drought from models

### Fort Keogh
### total ANPP
# full model just to check for whether transformations are necessary
fk_anpp_pchange_full_model <- lme(pchange_ANPP ~ Year*Drought*Grazing
                                  , data = anpp_fk_rr
                                  , random = ~1|Block/Paddock/Plot
                                  , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                  , control=lmeControl(returnObject=TRUE)
                                  , na.action = na.omit
)

Anova(fk_anpp_pchange_full_model, type=3)      
check_model(fk_anpp_pchange_full_model)
##

### get predictions of npp at 99% drought for each year
# 2018
fk_anpp_pchange_2018_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2018)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_anpp_pchange_2018_model, type=3)    
check_model(fk_anpp_pchange_2018_model)    
fk_2018_drought_int <- fk_anpp_pchange_2018_model$coefficients$fixed[1]
fk_2018_drought_slope <- fk_anpp_pchange_2018_model$coefficients$fixed[2]
fk_2018_yhat_99 <- fk_2018_drought_slope*99+fk_2018_drought_int

ggplot(filter(anpp_fk_rr, Year==2018), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2019
fk_anpp_pchange_2019_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2019)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_anpp_pchange_2019_model, type=3)    
check_model(fk_anpp_pchange_2019_model)    
fk_2019_drought_int <- fk_anpp_pchange_2019_model$coefficients$fixed[1]
fk_2019_drought_slope <- fk_anpp_pchange_2019_model$coefficients$fixed[2]
fk_2019_yhat_99 <- fk_2019_drought_slope*99+fk_2019_drought_int

ggplot(filter(anpp_fk_rr, Year==2019), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2020
fk_anpp_pchange_2020_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2020)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_anpp_pchange_2020_model, type=3)    
check_model(fk_anpp_pchange_2020_model)    
fk_2020_drought_int <- fk_anpp_pchange_2020_model$coefficients$fixed[1]
fk_2020_drought_slope <- fk_anpp_pchange_2020_model$coefficients$fixed[2]
fk_2020_yhat_99 <- fk_2020_drought_slope*99+fk_2020_drought_int

ggplot(filter(anpp_fk_rr, Year==2020), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2021
fk_anpp_pchange_2021_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2021)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_anpp_pchange_2021_model, type=3)    
check_model(fk_anpp_pchange_2021_model)    
fk_2021_drought_int <- fk_anpp_pchange_2021_model$coefficients$fixed[1]
fk_2021_drought_slope <- fk_anpp_pchange_2021_model$coefficients$fixed[2]
fk_2021_yhat_99 <- fk_2021_drought_slope*99+fk_2021_drought_int

ggplot(filter(anpp_fk_rr, Year==2021), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

#2022
fk_anpp_pchange_2022_model <-   lme(pchange_ANPP ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2022)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_anpp_pchange_2022_model, type=3)    
check_model(fk_anpp_pchange_2022_model)    
fk_2022_drought_int <- fk_anpp_pchange_2022_model$coefficients$fixed[1]
fk_2022_drought_slope <- fk_anpp_pchange_2022_model$coefficients$fixed[2]
fk_2022_yhat_99 <- fk_2022_drought_slope*99+fk_2022_drought_int

ggplot(filter(anpp_fk_rr, Year==2022), aes(x=Drought, y=pchange_ANPP)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
fk_anpp_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "Total ANPP",
  drt_magnitude = 99,
  pchange = c(fk_2018_yhat_99,
              fk_2019_yhat_99,
              fk_2020_yhat_99,
              fk_2021_yhat_99,
              fk_2022_yhat_99)
)

ggplot(filter(fk_anpp_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### C3P
# full model just to check for whether transformations are necessary
fk_C3P_pchange_full_model <- lme(pchange_C3P ~ Year*Drought*Grazing
                                 , data = anpp_fk_rr
                                 , random = ~1|Block/Paddock/Plot
                                 , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit
)

Anova(fk_C3P_pchange_full_model, type=3)      
check_model(fk_C3P_pchange_full_model)
##

### get predictions of npp at 99% drought for each year
# 2018
fk_C3P_pchange_2018_model <-   lme(pchange_C3P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2018)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C3P_pchange_2018_model, type=3)    
check_model(fk_C3P_pchange_2018_model)    
fk_2018_C3P_drought_int <- fk_C3P_pchange_2018_model$coefficients$fixed[1]
fk_2018_C3P_drought_slope <- fk_C3P_pchange_2018_model$coefficients$fixed[2]
fk_2018_C3P_yhat_99 <- fk_2018_C3P_drought_slope*99+fk_2018_C3P_drought_int

ggplot(filter(anpp_fk_rr, Year==2018), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2019
fk_C3P_pchange_2019_model <-   lme(pchange_C3P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2019 & pchange_C3P < 6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C3P_pchange_2019_model, type=3)    
check_model(fk_C3P_pchange_2019_model)    
fk_2019_C3P_drought_int <- fk_C3P_pchange_2019_model$coefficients$fixed[1]
fk_2019_C3P_drought_slope <- fk_C3P_pchange_2019_model$coefficients$fixed[2]
fk_2019_C3P_yhat_99 <- fk_2019_C3P_drought_slope*99+fk_2019_C3P_drought_int

ggplot(filter(anpp_fk_rr, Year==2019 & pchange_C3P<6), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2020
fk_C3P_pchange_2020_model <-   lme(pchange_C3P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2020)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C3P_pchange_2020_model, type=3)    
check_model(fk_C3P_pchange_2020_model)    
fk_2020_C3P_drought_int <- fk_C3P_pchange_2020_model$coefficients$fixed[1]
fk_2020_C3P_drought_slope <- fk_C3P_pchange_2020_model$coefficients$fixed[2]
fk_2020_C3P_yhat_99 <- fk_2020_C3P_drought_slope*99+fk_2020_C3P_drought_int

ggplot(filter(anpp_fk_rr, Year==2020), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2021
fk_C3P_pchange_2021_model <-   lme(pchange_C3P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2021)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C3P_pchange_2021_model, type=3)    
check_model(fk_C3P_pchange_2021_model)    
fk_2021_C3P_drought_int <- fk_C3P_pchange_2021_model$coefficients$fixed[1]
fk_2021_C3P_drought_slope <- fk_C3P_pchange_2021_model$coefficients$fixed[2]
fk_2021_C3P_yhat_99 <- fk_2021_C3P_drought_slope*99+fk_2021_C3P_drought_int

ggplot(filter(anpp_fk_rr, Year==2021), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

#2022
fk_C3P_pchange_2022_model <-   lme(pchange_C3P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2022)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C3P_pchange_2022_model, type=3)    
check_model(fk_C3P_pchange_2022_model)    
fk_2022_C3P_drought_int <- fk_C3P_pchange_2022_model$coefficients$fixed[1]
fk_2022_C3P_drought_slope <- fk_C3P_pchange_2022_model$coefficients$fixed[2]
fk_2022_C3P_yhat_99 <- fk_2022_C3P_drought_slope*99+fk_2022_C3P_drought_int

ggplot(filter(anpp_fk_rr, Year==2022), aes(x=Drought, y=pchange_C3P)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
fk_C3P_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "C3P",
  drt_magnitude = 99,
  pchange = c(fk_2018_C3P_yhat_99,
              fk_2019_C3P_yhat_99,
              fk_2020_C3P_yhat_99,
              fk_2021_C3P_yhat_99,
              fk_2022_C3P_yhat_99)
)

ggplot(filter(fk_C3P_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### C4P -- THIS SECTION NEEDS A LOT OF TRANSFORMATIONS OF THE DATA -- I AM NOT DOING ANY OF THAT FOR THE HANDOUT
# full model just to check for whether transformations are necessary
anpp_fk_rr$pchange_C4P_log <- log(anpp_fk_rr$pchange_C4P - min(anpp_fk_rr$pchange_C4P) +.01)

fk_C4P_pchange_full_model <- lme(pchange_C4P_log ~ Year*Drought*Grazing
                                 , data = filter(anpp_fk_rr, pchange_C4P<6)
                                 , random = ~1|Block/Paddock/Plot
                                 , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                 , control=lmeControl(returnObject=TRUE)
                                 , na.action = na.omit
)
#filter(anpp_fk_rr, pchange_C4P >6) %>% dplyr::select(Year:Plot, C4P_gm2, ctrl_C4P, pchange_C4P)
Anova(fk_C4P_pchange_full_model, type=3)      
check_model(fk_C4P_pchange_full_model)
check_outliers(fk_C4P_pchange_full_model)

##

log(.47)
exp(-.7550226)

### get predictions of npp at 99% drought for each year
# 2018
fk_C4P_pchange_2018_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2018 & pchange_C4P < 2)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C4P_pchange_2018_model, type=3)    
check_model(fk_C4P_pchange_2018_model)    
fk_2018_C4P_drought_int <- fk_C4P_pchange_2018_model$coefficients$fixed[1]
fk_2018_C4P_drought_slope <- fk_C4P_pchange_2018_model$coefficients$fixed[2]
fk_2018_C4P_yhat_99 <- fk_2018_C4P_drought_slope*99+fk_2018_C4P_drought_int

ggplot(filter(anpp_fk_rr, Year==2018 & pchange_C4P<2), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")


#2019
fk_C4P_pchange_2019_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2019 & pchange_C4P < 6)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C4P_pchange_2019_model, type=3)    
check_model(fk_C4P_pchange_2019_model)    
fk_2019_C4P_drought_int <- fk_C4P_pchange_2019_model$coefficients$fixed[1]
fk_2019_C4P_drought_slope <- fk_C4P_pchange_2019_model$coefficients$fixed[2]
fk_2019_C4P_yhat_99 <- fk_2019_C4P_drought_slope*99+fk_2019_C4P_drought_int

ggplot(filter(anpp_fk_rr, Year==2019 & pchange_C4P<6), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

#2020
fk_C4P_pchange_2020_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2020 & pchange_C4P < 4)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C4P_pchange_2020_model, type=3)    
check_model(fk_C4P_pchange_2020_model)    
fk_2020_C4P_drought_int <- fk_C4P_pchange_2020_model$coefficients$fixed[1]
fk_2020_C4P_drought_slope <- fk_C4P_pchange_2020_model$coefficients$fixed[2]
fk_2020_C4P_yhat_99 <- fk_2020_C4P_drought_slope*99+fk_2020_C4P_drought_int

ggplot(filter(anpp_fk_rr, Year==2020&pchange_C4P<4), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

#2021
fk_C4P_pchange_2021_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2021 & pchange_C4P<4)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C4P_pchange_2021_model, type=3)    
check_model(fk_C4P_pchange_2021_model)    
fk_2021_C4P_drought_int <- fk_C4P_pchange_2021_model$coefficients$fixed[1]
fk_2021_C4P_drought_slope <- fk_C4P_pchange_2021_model$coefficients$fixed[2]
fk_2021_C4P_yhat_99 <- fk_2021_C4P_drought_slope*99+fk_2021_C4P_drought_int

ggplot(filter(anpp_fk_rr, Year==2021& pchange_C4P<6), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

#2022
fk_C4P_pchange_2022_model <-   lme(pchange_C4P ~ Drought*Grazing
                                   , data = filter(anpp_fk_rr, Year==2022& pchange_C4P<4)
                                   , random = ~1|Block/Paddock/Plot
)

Anova(fk_C4P_pchange_2022_model, type=3)    
check_model(fk_C4P_pchange_2022_model)
fk_2022_C4P_drought_int <- fk_C4P_pchange_2022_model$coefficients$fixed[1]
fk_2022_C4P_drought_slope <- fk_C4P_pchange_2022_model$coefficients$fixed[2]
fk_2022_C4P_yhat_99 <- fk_2022_C4P_drought_slope*99+fk_2022_C4P_drought_int

ggplot(filter(anpp_fk_rr, Year==2022& pchange_C4P<4), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
fk_C4P_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "C4P",
  drt_magnitude = 99,
  pchange = c(fk_2018_C4P_yhat_99,
              fk_2019_C4P_yhat_99,
              fk_2020_C4P_yhat_99,
              fk_2021_C4P_yhat_99,
              fk_2022_C4P_yhat_99)
)

ggplot(filter(fk_C4P_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### Forb
# full model just to check for whether transformations are necessary

fk_Forb_pchange_full_model <- lme(pchange_Forb ~ Year*Drought*Grazing
                                  , data = filter(anpp_fk_rr, pchange_Forb < 6)
                                  , random = ~1|Block/Paddock/Plot
                                  , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                  , control=lmeControl(returnObject=TRUE)
                                  , na.action = na.omit
)

Anova(fk_Forb_pchange_full_model, type=3)      
check_model(fk_Forb_pchange_full_model)
hist(anpp_fk_rr$pchange_Forb)
##

### get predictions of npp at 99% drought for each year
# 2018
fk_Forb_pchange_2018_model <-   lme(pchange_Forb ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2018 & pchange_Forb<6)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_Forb_pchange_2018_model, type=3)    
check_model(fk_Forb_pchange_2018_model)    
fk_2018_Forb_drought_int <- fk_Forb_pchange_2018_model$coefficients$fixed[1]
fk_2018_Forb_drought_slope <- fk_Forb_pchange_2018_model$coefficients$fixed[2]
fk_2018_Forb_yhat_99 <- fk_2018_Forb_drought_slope*99+fk_2018_Forb_drought_int

ggplot(filter(anpp_fk_rr, Year==2018 & pchange_Forb<6), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2019
fk_Forb_pchange_2019_model <-   lme(pchange_Forb ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2019 & pchange_Forb < 6)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_Forb_pchange_2019_model, type=3)    
check_model(fk_Forb_pchange_2019_model)    
fk_2019_Forb_drought_int <- fk_Forb_pchange_2019_model$coefficients$fixed[1]
fk_2019_Forb_drought_slope <- fk_Forb_pchange_2019_model$coefficients$fixed[2]
fk_2019_Forb_yhat_99 <- fk_2019_Forb_drought_slope*99+fk_2019_Forb_drought_int

ggplot(filter(anpp_fk_rr, Year==2019 & pchange_Forb<6), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2020
fk_Forb_pchange_2020_model <-   lme(pchange_Forb ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2020)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_Forb_pchange_2020_model, type=3)    
check_model(fk_Forb_pchange_2020_model)    
fk_2020_Forb_drought_int <- fk_Forb_pchange_2020_model$coefficients$fixed[1]
fk_2020_Forb_drought_slope <- fk_Forb_pchange_2020_model$coefficients$fixed[2]
fk_2020_Forb_yhat_99 <- fk_2020_Forb_drought_slope*99+fk_2020_Forb_drought_int

ggplot(filter(anpp_fk_rr, Year==2020), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2021
fk_Forb_pchange_2021_model <-   lme(pchange_Forb ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2021)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_Forb_pchange_2021_model, type=3)    
check_model(fk_Forb_pchange_2021_model)    
fk_2021_Forb_drought_int <- fk_Forb_pchange_2021_model$coefficients$fixed[1]
fk_2021_Forb_drought_slope <- fk_Forb_pchange_2021_model$coefficients$fixed[2]
fk_2021_Forb_yhat_99 <- fk_2021_Forb_drought_slope*99+fk_2021_Forb_drought_int

ggplot(filter(anpp_fk_rr, Year==2021), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

#2022
fk_Forb_pchange_2022_model <-   lme(pchange_Forb ~ Drought*Grazing
                                    , data = filter(anpp_fk_rr, Year==2022)
                                    , random = ~1|Block/Paddock/Plot
)

Anova(fk_Forb_pchange_2022_model, type=3)    
check_model(fk_Forb_pchange_2022_model)    
fk_2022_Forb_drought_int <- fk_Forb_pchange_2022_model$coefficients$fixed[1]
fk_2022_Forb_drought_slope <- fk_Forb_pchange_2022_model$coefficients$fixed[2]
fk_2022_Forb_yhat_99 <- fk_2022_Forb_drought_slope*99+fk_2022_Forb_drought_int

ggplot(filter(anpp_fk_rr, Year==2022), aes(x=Drought, y=pchange_Forb)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
fk_Forb_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "Forb",
  drt_magnitude = 99,
  pchange = c(fk_2018_Forb_yhat_99,
              fk_2019_Forb_yhat_99,
              fk_2020_Forb_yhat_99,
              fk_2021_Forb_yhat_99,
              fk_2022_Forb_yhat_99)
)

ggplot(filter(fk_Forb_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

### AnnualGrass
# full model just to check for whether transformations are necessary
fk_AnnualGrass_pchange_full_model <- lme(pchange_AnnualGrass ~ Year*Drought*Grazing
                                         , data = filter(anpp_fk_rr, pchange_AnnualGrass<6)
                                         , random = ~1|Block/Paddock/Plot
                                         , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                         , control=lmeControl(returnObject=TRUE)
                                         , na.action = na.omit
)

Anova(fk_AnnualGrass_pchange_full_model, type=3)      
check_model(fk_AnnualGrass_pchange_full_model)
##

### get predictions of npp at 99% drought for each year
# 2018
fk_AnnualGrass_pchange_2018_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                           , data = filter(anpp_fk_rr, Year==2018 & pchange_AnnualGrass<6)
                                           , random = ~1|Block/Paddock/Plot
)

Anova(fk_AnnualGrass_pchange_2018_model, type=3)    
check_model(fk_AnnualGrass_pchange_2018_model)    
fk_2018_AnnualGrass_drought_int <- fk_AnnualGrass_pchange_2018_model$coefficients$fixed[1]
fk_2018_AnnualGrass_drought_slope <- fk_AnnualGrass_pchange_2018_model$coefficients$fixed[2]
fk_2018_AnnualGrass_yhat_99 <- fk_2018_AnnualGrass_drought_slope*99+fk_2018_AnnualGrass_drought_int

ggplot(filter(anpp_fk_rr, Year==2018& pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2019
fk_AnnualGrass_pchange_2019_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                           , data = filter(anpp_fk_rr, Year==2019 & pchange_AnnualGrass < 6)
                                           , random = ~1|Block/Paddock/Plot
)

Anova(fk_AnnualGrass_pchange_2019_model, type=3)    
check_model(fk_AnnualGrass_pchange_2019_model)    
fk_2019_AnnualGrass_drought_int <- fk_AnnualGrass_pchange_2019_model$coefficients$fixed[1]
fk_2019_AnnualGrass_drought_slope <- fk_AnnualGrass_pchange_2019_model$coefficients$fixed[2]
fk_2019_AnnualGrass_yhat_99 <- fk_2019_AnnualGrass_drought_slope*99+fk_2019_AnnualGrass_drought_int

ggplot(filter(anpp_fk_rr, Year==2019 & pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2020
fk_AnnualGrass_pchange_2020_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                           , data = filter(anpp_fk_rr, Year==2020& pchange_AnnualGrass<6)
                                           , random = ~1|Block/Paddock/Plot
)

Anova(fk_AnnualGrass_pchange_2020_model, type=3)    
check_model(fk_AnnualGrass_pchange_2020_model)    
fk_2020_AnnualGrass_drought_int <- fk_AnnualGrass_pchange_2020_model$coefficients$fixed[1]
fk_2020_AnnualGrass_drought_slope <- fk_AnnualGrass_pchange_2020_model$coefficients$fixed[2]
fk_2020_AnnualGrass_yhat_99 <- fk_2020_AnnualGrass_drought_slope*99+fk_2020_AnnualGrass_drought_int

ggplot(filter(anpp_fk_rr, Year==2020), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2021
fk_AnnualGrass_pchange_2021_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                           , data = filter(anpp_fk_rr, Year==2021& pchange_AnnualGrass<4)
                                           , random = ~1|Block/Paddock/Plot
)

Anova(fk_AnnualGrass_pchange_2021_model, type=3)    
check_model(fk_AnnualGrass_pchange_2021_model)    
fk_2021_AnnualGrass_drought_int <- fk_AnnualGrass_pchange_2021_model$coefficients$fixed[1]
fk_2021_AnnualGrass_drought_slope <- fk_AnnualGrass_pchange_2021_model$coefficients$fixed[2]
fk_2021_AnnualGrass_yhat_99 <- fk_2021_AnnualGrass_drought_slope*99+fk_2021_AnnualGrass_drought_int

ggplot(filter(anpp_fk_rr, Year==2021& pchange_AnnualGrass<4), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

#2022
fk_AnnualGrass_pchange_2022_model <-   lme(pchange_AnnualGrass ~ Drought*Grazing
                                           , data = filter(anpp_fk_rr, Year==2022& pchange_AnnualGrass<6)
                                           , random = ~1|Block/Paddock/Plot
)

Anova(fk_AnnualGrass_pchange_2022_model, type=3)    
check_model(fk_AnnualGrass_pchange_2022_model)    
fk_2022_AnnualGrass_drought_int <- fk_AnnualGrass_pchange_2022_model$coefficients$fixed[1]
fk_2022_AnnualGrass_drought_slope <- fk_AnnualGrass_pchange_2022_model$coefficients$fixed[2]
fk_2022_AnnualGrass_yhat_99 <- fk_2022_AnnualGrass_drought_slope*99+fk_2022_AnnualGrass_drought_int

ggplot(filter(anpp_fk_rr, Year==2022& pchange_AnnualGrass<6), aes(x=Drought, y=pchange_AnnualGrass)) +geom_point() + geom_smooth(method="lm")

# Put all yhats together
fk_AnnualGrass_pchange_99 <- data.frame(
  Year = 2018:2022,
  variable = "AnnualGrass",
  drt_magnitude = 99,
  pchange = c(fk_2018_AnnualGrass_yhat_99,
              fk_2019_AnnualGrass_yhat_99,
              fk_2020_AnnualGrass_yhat_99,
              fk_2021_AnnualGrass_yhat_99,
              fk_2022_AnnualGrass_yhat_99)
)

ggplot(filter(fk_AnnualGrass_pchange_99,Year!=2018), aes(x=Year, y=pchange)) +
  geom_col(col='black', fill=c('firebrick1', 'firebrick1', 'dodgerblue','dodgerblue')) +
  ylab("Drought response (% change)")

year_key <- data.frame(Year=2018:2022, trt_group=c("no_drt","drought","drought","no_drt","no_drt"))

fk_pchange_99_full <- bind_rows(
  fk_anpp_pchange_99,
  fk_C3P_pchange_99,
  fk_C4P_pchange_99,
  fk_AnnualGrass_pchange_99,
) %>%
  full_join(year_key, by="Year") %>%
  mutate(pchange=replace(pchange, pchange<=-1,-1))

fk_pchange_99_full$variable <- factor(fk_pchange_99_full$variable, levels=c("Total ANPP", "C3P","C4P","AnnualGrass"))

ggplot(filter(fk_pchange_99_full,Year!=2018), aes(x=Year, y=pchange*100, fill=trt_group)) +
  geom_col(col='black') +
  scale_fill_manual(values=c("firebrick2","dodgerblue")) +
  ylab("Drought response (% change)") +
  facet_wrap(~variable)



### Just Drought
### Note: We need to deal with the infinity values calculated in here for lnRR (0s causing this), but I'm just using pchange currently so I'm ignoring
anpp_fk_rr_drt_means <- anpp_fk_rr %>% 
  group_by(Year, Drought) %>%
  summarize_at(vars(lnrr_ANPP:pchange_AnnualGrass),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(
    lnrr_ANPP_mean=lnrr_ANPP_fn1, lnrr_ANPP_se=lnrr_ANPP_fn2,
    pchange_ANPP_mean=pchange_ANPP_fn1, pchange_ANPP_se=pchange_ANPP_fn2,
    lnrr_C3P_mean=lnrr_C3P_fn1, lnrr_C3P_se=lnrr_C3P_fn2,
    pchange_C3P_mean=pchange_C3P_fn1, pchange_C3P_se=pchange_C3P_fn2,
    lnrr_C4P_mean=lnrr_C4P_fn1, lnrr_C4P_se=lnrr_C4P_fn2,
    pchange_C4P_mean=pchange_C4P_fn1, pchange_C4P_se=pchange_C4P_fn2,
    lnrr_Forb_mean=lnrr_Forb_fn1, lnrr_Forb_se=lnrr_Forb_fn2,
    pchange_Forb_mean=pchange_Forb_fn1, pchange_Forb_se=pchange_Forb_fn2,
    lnrr_dead_mean=lnrr_dead_fn1, lnrr_dead_se=lnrr_dead_fn2,
    pchange_dead_mean=pchange_dead_fn1, pchange_dead_se=pchange_dead_fn2,
    lnrr_AnnualGrass_mean=lnrr_AnnualGrass_fn1, lnrr_AnnualGrass_se=lnrr_AnnualGrass_fn2,
    pchange_AnnualGrass_mean=pchange_AnnualGrass_fn1, pchange_AnnualGrass_se=pchange_AnnualGrass_fn2
  )

png(file="..\\..\\figures\\TBRI2022\\anpp pchange drought only_18-22_fk.png", width=12.6, height=4.2, units="in", res=600)
print(
  ggplot(anpp_fk_rr_drt_means, aes(x= Drought, y=pchange_ANPP_mean*100,
                                   ymin=(pchange_ANPP_mean-pchange_ANPP_se)*100, ymax=(pchange_ANPP_mean+pchange_ANPP_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=0) +
    #  geom_smooth(method="lm",se=F, col="red") +
    geom_smooth(method="lm", col="red", se=F) +
    geom_hline(yintercept=0) +
    facet_grid(cols=vars(Year)) +
    ylab("Drought effect (% ANPP change)") +
    xlab("Drought magnitude (% ppt reduction)") +
    xlim(25,100)
)
dev.off()

### Grazing and Drought treatments
### Note: We need to deal with the infinity values calculated in here for lnRR (0s causing this), but I'm just using pchange currently so I'm ignoring
anpp_fk_rr_drt_grz_means <- anpp_fk_rr %>% 
  group_by(Year, Drought, Grazing) %>%
  summarize_at(vars(lnrr_ANPP:pchange_AnnualGrass),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(
    lnrr_ANPP_mean=lnrr_ANPP_fn1, lnrr_ANPP_se=lnrr_ANPP_fn2,
    pchange_ANPP_mean=pchange_ANPP_fn1, pchange_ANPP_se=pchange_ANPP_fn2,
    lnrr_C3P_mean=lnrr_C3P_fn1, lnrr_C3P_se=lnrr_C3P_fn2,
    pchange_C3P_mean=pchange_C3P_fn1, pchange_C3P_se=pchange_C3P_fn2,
    lnrr_C4P_mean=lnrr_C4P_fn1, lnrr_C4P_se=lnrr_C4P_fn2,
    pchange_C4P_mean=pchange_C4P_fn1, pchange_C4P_se=pchange_C4P_fn2,
    lnrr_Forb_mean=lnrr_Forb_fn1, lnrr_Forb_se=lnrr_Forb_fn2,
    pchange_Forb_mean=pchange_Forb_fn1, pchange_Forb_se=pchange_Forb_fn2,
    lnrr_dead_mean=lnrr_dead_fn1, lnrr_dead_se=lnrr_dead_fn2,
    pchange_dead_mean=pchange_dead_fn1, pchange_dead_se=pchange_dead_fn2,
    lnrr_AnnualGrass_mean=lnrr_AnnualGrass_fn1, lnrr_AnnualGrass_se=lnrr_AnnualGrass_fn2,
    pchange_AnnualGrass_mean=pchange_AnnualGrass_fn1, pchange_AnnualGrass_se=pchange_AnnualGrass_fn2
  )

anpp_fk_rr_drt_grz_means$Grazing <- factor(anpp_fk_rr_drt_grz_means$Grazing, levels=c("MLLMM","MMMMM","HHMMM"))

q3 <- c("#00AD9A", "#9183E6", "#E16A86")

png(file="..\\..\\figures\\TBRI2022\\anpp pchange drought and grazing.png", width=13.6, height=4.2, units="in", res=600)
print(
  ggplot(anpp_fk_rr_drt_grz_means, aes(x= Drought, y=pchange_ANPP_mean*100, col=Grazing,
                                       ymin=(pchange_ANPP_mean-pchange_ANPP_se)*100, ymax=(pchange_ANPP_mean+pchange_ANPP_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
    facet_grid(cols=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (% ANPP change)") +
    xlab("Drought magnitude (% ppt reduction)") +
    xlim(25,100) +ylim(-75, 80)
)
dev.off()

png(file="..\\..\\figures\\TBRI2022\\c3p pchange drought and grazing.png", width=13.6, height=4.2, units="in", res=600)
print(
  ggplot(anpp_fk_rr_drt_grz_means, aes(x= Drought, y=pchange_C3P_mean*100, col=Grazing,
                                       ymin=(pchange_C3P_mean-pchange_C3P_se)*100, ymax=(pchange_C3P_mean+pchange_C3P_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
    facet_grid(cols=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (%)") +
    xlab("Precipitation reduction (%)") +
    xlim(25,100)
)
dev.off()

png(file="..\\..\\figures\\TBRI2022\\c3p pchange drought and grazing 2021 only.png", width=6.2, height=4.2, units="in", res=600)
print(
  ggplot(filter(anpp_fk_rr_drt_grz_means,Year==2021), aes(x= Drought, y=pchange_C3P_mean*100, col=Grazing,
                                                          ymin=(pchange_C3P_mean-pchange_C3P_se)*100, ymax=(pchange_C3P_mean+pchange_C3P_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
    #  facet_grid(cols=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (%)") +
    xlab("Precipitation reduction (%)") +
    xlim(24,100)
)
dev.off()

## C4P and AnnualGrass -- annual grass is higher in heavy grazed plots under any drought magnitude... could be because it was already a dry year and there was a threshold there that was crossed even in the 25% plots... I'm a little skeptical though
pdf(file="..\\..\\figures\\TBRI2022\\c4p pchange drought and grazing.pdf", width=6.5, height=12.6)
print(
  ggplot(filter(anpp_fk_rr_drt_grz_means, Year==2021), aes(x= Drought, y=pchange_AnnualGrass_mean*100, col=Grazing,
                                       ymin=(pchange_AnnualGrass_mean-pchange_AnnualGrass_se)*100, ymax=(pchange_AnnualGrass_mean+pchange_AnnualGrass_se)*100 ))+
    geom_point(size=3) +
    geom_errorbar(width=1) +
    geom_smooth(method="lm", se=F) +
    geom_hline(yintercept=0) +
    facet_grid(rows=vars(Year), scales="free") +
    scale_colour_manual(values=q3, name='Grazing\nTreatment') +
    ylab("Drought effect (%)") +
    xlab("Precipitation reduction (%)") +
    xlim(25,100) 
)
dev.off()

}

###
### Models - Fort Keogh
###
{

  ### Total ANPP
anpp_fk_plot_means$logANPP <- log(anpp_fk_plot_means$ANPP_gm2)  
 fk_anpp_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='logANPP',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=anpp_fk_plot_means
                        )  

  
  fk_anpp_full_model <- lme(log(ANPP_gm2) ~ Year*Drought*Grazing
                        , data = anpp_fk_plot_means
                        , random = ~1|Block/Paddock/Plot
                        , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                        , control=lmeControl(returnObject=TRUE)
                        , na.action = na.omit
        )

    Anova(fk_anpp_full_model, type=3)      
    check_model(fk_anpp_full_model)

  fk_anpp_2018_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2018)
                              , random = ~1|Block/Paddock/Plot
  )

  Anova(fk_anpp_2018_model, type=3)    
  
  fk_anpp_2019_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2019)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_anpp_2019_model, type=3)    

  fk_anpp_2020_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2020)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_anpp_2020_model, type=3)    
  
  fk_anpp_2021_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2021)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_anpp_2021_model, type=3)    
  
  fk_anpp_2022_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2022)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_anpp_2022_model, type=3)    
  
  ### C3Perennial
  anpp_fk_plot_means$logC3P <- log(anpp_fk_plot_means$C3P_gm2+.1)  
  fk_C3P_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='logC3P',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=anpp_fk_plot_means
  )  
  
  
  fk_C3P_full_model <- lme(log(C3P_gm2+.1) ~ Year*Drought*Grazing
                            , data = anpp_fk_plot_means
                            , random = ~1|Block/Paddock/Plot
                            , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit
  )
  
  Anova(fk_C3P_full_model, type=3)      
  check_model(fk_C3P_full_model)
  
  fk_C3P_2018_model <-   lme(log(C3P_gm2+.1) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2018)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C3P_2018_model, type=3)    
  
  fk_C3P_2019_model <-   lme(log(C3P_gm2+.1) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2019)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C3P_2019_model, type=3)    
  
  fk_C3P_2020_model <-   lme(log(C3P_gm2+.1) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2020)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C3P_2020_model, type=3)    
  
  fk_C3P_2021_model <-   lme(log(C3P_gm2+.1) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2021)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C3P_2021_model, type=3)    
  
  fk_C3P_2022_model <-   lme(log(C3P_gm2+.1) ~ Drought*Grazing
                              , data = filter(anpp_fk_plot_means, Year==2022)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C3P_2022_model, type=3) 
  
  ### C4Perennial
  anpp_fk_plot_means$logC4P <- log(anpp_fk_plot_means$C4P_gm2+.1)  
  fk_C4P_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='logC4P',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=anpp_fk_plot_means
  )  
  
  
  fk_C4P_full_model <- lme(log(C4P_gm2+.1) ~ Year*Drought*Grazing
                           , data = anpp_fk_plot_means
                           , random = ~1|Block/Paddock/Plot
                           , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit
  )
  
  Anova(fk_C4P_full_model, type=3)      
  check_model(fk_C4P_full_model)
  hist(anpp_fk_plot_means$C4P_gm2)

  # Split by year
  fk_C4P_2018_model <-   lme(log(C4P_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2018)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C4P_2018_model, type=3)    
  
  fk_C4P_2019_model <-   lme(log(C4P_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2019)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C4P_2019_model, type=3)    
  
  fk_C4P_2020_model <-   lme(log(C4P_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2020)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C4P_2020_model, type=3)    
  
  fk_C4P_2021_model <-   lme(log(C4P_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2021)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C4P_2021_model, type=3)    
  
  fk_C4P_2022_model <-   lme(log(C4P_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2022)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_C4P_2022_model, type=3) 
  
  ### Forberennial
  anpp_fk_plot_means$sqrtForb <- sqrt(anpp_fk_plot_means$Forb_gm2)  
  fk_Forb_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='sqrtForb',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=anpp_fk_plot_means
  )  
  
  
  fk_Forb_full_model <- lme(sqrt(Forb_gm2) ~ Year*Drought*Grazing
                           , data = anpp_fk_plot_means
                           , random = ~1|Block/Paddock/Plot
                           , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit
  )
  
  check_model(fk_Forb_full_model)
  Anova(fk_Forb_full_model, type=3)      
  hist(anpp_fk_plot_means$Forb_gm2)
  
  # Split by year
  fk_Forb_2018_model <-   lme(sqrt(Forb_gm2) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2018)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_Forb_2018_model, type=3)    
  
  fk_Forb_2019_model <-   lme(sqrt(Forb_gm2) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2019)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_Forb_2019_model, type=3)    
  
  fk_Forb_2020_model <-   lme(sqrt(Forb_gm2) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2020)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_Forb_2020_model, type=3)    
  
  fk_Forb_2021_model <-   lme(sqrt(Forb_gm2) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2021)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_Forb_2021_model, type=3)    
  
  fk_Forb_2022_model <-   lme(sqrt(Forb_gm2) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2022)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_Forb_2022_model, type=3) 
  
  ### AnnualGrasserennial
  anpp_fk_plot_means$logAnnualGrass <- log(anpp_fk_plot_means$AnnualGrass_gm2+.1)  
  fk_AnnualGrass_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='logAnnualGrass',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=anpp_fk_plot_means
  )  
  
  
  fk_AnnualGrass_full_model <- lme(log(AnnualGrass_gm2+.1) ~ Year*Drought*Grazing
                           , data = anpp_fk_plot_means
                           , random = ~1|Block/Paddock/Plot
                           , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit
  )
  
  Anova(fk_AnnualGrass_full_model, type=3)      
  check_model(fk_AnnualGrass_full_model)
  hist(anpp_fk_plot_means$AnnualGrass_gm2)
  
  # Split by year
  fk_AnnualGrass_2018_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2018)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_AnnualGrass_2018_model, type=3)    
  
  fk_AnnualGrass_2019_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2019)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_AnnualGrass_2019_model, type=3)    
  
  fk_AnnualGrass_2020_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2020)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_AnnualGrass_2020_model, type=3)    
  
  fk_AnnualGrass_2021_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2021)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_AnnualGrass_2021_model, type=3)    
  
  fk_AnnualGrass_2022_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                             , data = filter(anpp_fk_plot_means, Year==2022)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(fk_AnnualGrass_2022_model, type=3) 
}
  
###
### Models - Thunder Basin
###
{
  
  ### Total ANPP
  anpp_tb_plot_means$logANPP <- log(anpp_tb_plot_means$ANPP_gm2+.1)
  tb_anpp_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='logANPP',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=anpp_tb_plot_means
  )  
  
  tb_anpp_full_model <- lme(log(ANPP_gm2+.1) ~ Year*Drought*Grazing
                            , data = anpp_tb_plot_means
                            , random = ~1|Block/Paddock/Plot
                            , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit
  )
  
  Anova(tb_anpp_full_model, type=3)      
  check_model(tb_anpp_full_model)
  
  tb_anpp_2018_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2018)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_anpp_2018_model, type=3)    
  
  tb_anpp_2019_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2019)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_anpp_2019_model, type=3)    
  
  tb_anpp_2020_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2020)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_anpp_2020_model, type=3)    
  predict(tb_anpp_2020_model)    
  predict(tb_anpp_2020_model, newdata = data.frame(Drought=c(0,99)))    
fitted(tb_anpp_2022_model, level=c(0,1))
  tb_anpp_2021_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2021)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_anpp_2021_model, type=3)    
  
  tb_anpp_2022_model <-   lme(log(ANPP_gm2) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2022)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_anpp_2022_model, type=3)    
  
  ### C3Perennial
  anpp_tb_plot_means$sqrtC3P <- sqrt(anpp_tb_plot_means$C3P_gm2)  
  tb_C3P_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='sqrtC3P',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=anpp_tb_plot_means
  )  
  
  
  tb_C3P_full_model <- lme(sqrt(C3P_gm2) ~ Year*Drought*Grazing
                           , data = anpp_tb_plot_means
                           , random = ~1|Block/Paddock/Plot
                           , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit
  )
  
  Anova(tb_C3P_full_model, type=3)      
  check_model(tb_C3P_full_model)
  
  tb_C3P_2018_model <-   lme(sqrt(C3P_gm2) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2018)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C3P_2018_model, type=3)    
  
  tb_C3P_2019_model <-   lme(sqrt(C3P_gm2) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2019)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C3P_2019_model, type=3)    
  
  tb_C3P_2020_model <-   lme(sqrt(C3P_gm2) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2020)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C3P_2020_model, type=3)    
  
  tb_C3P_2021_model <-   lme(sqrt(C3P_gm2) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2021)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C3P_2021_model, type=3)    
  
  tb_C3P_2022_model <-   lme(sqrt(C3P_gm2) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2022)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C3P_2022_model, type=3) 
  
  ### C4Perennial
  anpp_tb_plot_means$cbrtC4P <- anpp_tb_plot_means$C4P_gm2^(1/3)  
  tb_C4P_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                 DepVar='cbrtC4P',
                                 RndForm='~1 |Block/Paddock/Plot',
                                 Data=anpp_tb_plot_means
  )  
  
  
  tb_C4P_full_model <- lme(C4P_gm2^(1/3) ~ Year*Drought*Grazing
                           , data = anpp_tb_plot_means
                           , random = ~1|Block/Paddock/Plot
                           , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit
  )
  
  Anova(tb_C4P_full_model, type=3)      
  check_model(tb_C4P_full_model)
  hist(anpp_tb_plot_means$C4P_gm2)
  
  # Split by year
  tb_C4P_2018_model <-   lme(C4P_gm2^(1/3) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2018)
                             , random = ~1|Block/Paddock/Plot
  )
  
  check_model(tb_C4P_2018_model)    
  Anova(tb_C4P_2018_model, type=3)    
  
  tb_C4P_2019_model <-   lme(C4P_gm2^(1/3) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2019)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C4P_2019_model, type=3)    
  
  tb_C4P_2020_model <-   lme(C4P_gm2^(1/3) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2020)
                             , random = ~1|Block/Paddock/Plot
  )
  
  check_model(tb_C4P_2020_model, type=3)    
  Anova(tb_C4P_2020_model, type=3)    
  
  tb_C4P_2021_model <-   lme(C4P_gm2^(1/3) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2021)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C4P_2021_model, type=3)    
  
  tb_C4P_2022_model <-   lme(C4P_gm2^(1/3) ~ Drought*Grazing
                             , data = filter(anpp_tb_plot_means, Year==2022)
                             , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_C4P_2022_model, type=3) 
  
  ### Forberennial
  anpp_tb_plot_means$cbrtForb <- anpp_tb_plot_means$Forb_gm2^(1/3)  
  tb_Forb_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                  DepVar='cbrtForb',
                                  RndForm='~1 |Block/Paddock/Plot',
                                  Data=anpp_tb_plot_means
  )  
  
  
  tb_Forb_full_model <- lme(Forb_gm2^(1/3) ~ Year*Drought*Grazing
                            , data = anpp_tb_plot_means
                            , random = ~1|Block/Paddock/Plot
                            , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit
  )
  
  check_model(tb_Forb_full_model)
  Anova(tb_Forb_full_model, type=3)      
  hist(anpp_tb_plot_means$Forb_gm2^(1/3))
  
  # Split by year
  tb_Forb_2018_model <-   lme(Forb_gm2^(1/3) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2018)
                              , random = ~1|Block/Paddock/Plot
  )
  
  check_model(tb_Forb_2018_model, type=3)    
  Anova(tb_Forb_2018_model, type=3)    
  
  tb_Forb_2019_model <-   lme(Forb_gm2^(1/3) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2019)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_Forb_2019_model, type=3)    
  
  tb_Forb_2020_model <-   lme(Forb_gm2^(1/3) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2020)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_Forb_2020_model, type=3)    
  
  tb_Forb_2021_model <-   lme(Forb_gm2^(1/3) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2021)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_Forb_2021_model, type=3)    
  
  tb_Forb_2022_model <-   lme(Forb_gm2^(1/3) ~ Drought*Grazing
                              , data = filter(anpp_tb_plot_means, Year==2022)
                              , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_Forb_2022_model, type=3) 
  
  ### AnnualGrasserennial
  anpp_tb_plot_means$logAnnualGrass <- log(anpp_tb_plot_means$AnnualGrass_gm2+.1)  
  tb_AnnualGrass_full_model <-  anova_t3(IndVars=c('Year','Drought', 'Grazing'),
                                         DepVar='logAnnualGrass',
                                         RndForm='~1 |Block/Paddock/Plot',
                                         Data=anpp_tb_plot_means
  )  
  
  
  tb_AnnualGrass_full_model <- lme(log(AnnualGrass_gm2+.1) ~ Year*Drought*Grazing
                                   , data = anpp_tb_plot_means
                                   , random = ~1|Block/Paddock/Plot
                                   , correlation=corCompSymm(form = ~1|Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit
  )
  
  Anova(tb_AnnualGrass_full_model, type=3)      
  check_model(tb_AnnualGrass_full_model)
  hist(anpp_tb_plot_means$AnnualGrass_gm2)
  
  # Split by year
  tb_AnnualGrass_2018_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                                     , data = filter(anpp_tb_plot_means, Year==2018)
                                     , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_AnnualGrass_2018_model, type=3)    
  
  tb_AnnualGrass_2019_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                                     , data = filter(anpp_tb_plot_means, Year==2019)
                                     , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_AnnualGrass_2019_model, type=3)    
  
  tb_AnnualGrass_2020_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                                     , data = filter(anpp_tb_plot_means, Year==2020)
                                     , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_AnnualGrass_2020_model, type=3)    
  
  tb_AnnualGrass_2021_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                                     , data = filter(anpp_tb_plot_means, Year==2021)
                                     , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_AnnualGrass_2021_model, type=3)    
  
  tb_AnnualGrass_2022_model <-   lme(log(AnnualGrass_gm2+.1) ~ Drought*Grazing
                                     , data = filter(anpp_tb_plot_means, Year==2022)
                                     , random = ~1|Block/Paddock/Plot
  )
  
  Anova(tb_AnnualGrass_2022_model, type=3) 
  
  
}



########### Extra code below ####################
{

### means by grazing treatment
anpp_tb_grz_means <- anpp_tb_plotmeans %>%
  group_by(site, date, Grazing) %>%
  summarize_at(vars(C3_perennial:ARFR,anpp),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(C3_perennial_mean=C3_perennial_fn1,C3_perennial_se=C3_perennial_fn2,
         C4_perennial_mean=C4_perennial_fn1,C4_perennial_se=C4_perennial_fn2,
         forb_mean=forb_fn1,forb_se=forb_fn2,
         annual_grass_mean=annual_grass_fn1,annual_grass_se=annual_grass_fn2,
         vulpia_mean=vulpia_fn1,vulpia_se=vulpia_fn2,
         bromus_mean=bromus_fn1,bromus_se=bromus_fn2,
         GUSA_mean=GUSA_fn1,GUSA_se=GUSA_fn2,
         ARFR_mean=ARFR_fn1,ARFR_se=ARFR_fn2,
         anpp_mean=anpp_fn1,anpp_se=anpp_fn2) %>%
  mutate(Grazing=factor(Grazing, levels=c("MLLMM","MMMMM","HHMMM")))




### Plot drought means
drt_plot_tb <- ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=anpp_mean*10, fill=factor(Drought),
                                             ymin=(anpp_mean-anpp_se)*10, ymax=(anpp_mean+anpp_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)") +
  ylim(0,225)

pdf("..//..//figures//anpp_drought_means_TB_2018-2020.pdf", width=12, height=6, useDingbats = F)
print(drt_plot_tb)
dev.off()

### Plot grazing means
ggplot(anpp_tb_grz_means, aes(x= factor(Grazing), y=anpp_mean*10, fill=factor(Grazing),
                              ymin=(anpp_mean-anpp_se)*10, ymax=(anpp_mean+anpp_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=grazingColor, name='Grazing\nTreatment') +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")

### Run models
# anpp_tb_plotmeans <- anpp_tb_plotmeans %>%
#   left_join(Plot_Key, by=c("site", "block", "Paddock", "plot"))


anpp_model_tb_2019_20 <-   anova_t3(IndVars=c('date','Grazing','Drought'),
                                    DepVar='anpp',
                                    RndForm='~1 |block/Paddock/plot',
                                    Data=subset(anpp_tb_plotmeans, date %in% c(2019,2020))
)  


anpp_model_tb_2019 <-   lme(anpp ~ Drought*Grazing
                            , data=subset(anpp_tb_plotmeans, date == 2019)
                            , random = ~1 |block/Paddock/plot
                            , na.action = na.omit)
anova.lme(anpp_model_tb_2019, type="marginal") ## This provides different results than anova (even when reordering variables)



anpp_model_tb_2020 <-   lme(anpp ~ Drought*Grazing
                            , data=subset(anpp_tb_plotmeans, date == 2020)
                            , random = ~1 |block/Paddock/plot
                            , na.action = na.omit)

anpp_regression_tb_2020 <-   lm(anpp ~ Drought
                                , data=subset(anpp_tb_plotmeans, date == 2020)
                                , na.action = na.omit)

anova.lme(anpp_model_tb_2020, type="marginal")
anova.lme(anpp_model_tb_2020)
anova(anpp_model_tb_2020)
Anova(anpp_model_tb_2020, type=3)
Anova(anpp_model_tb_2020, type=2, test.statistic="F")
?Anova

# C3 perennials
ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=C3_perennial_mean*10, fill=factor(Drought),
                              ymin=(C3_perennial_mean-C3_perennial_se)*10, ymax=(C3_perennial_mean+C3_perennial_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")

# C4 perennials
ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=C4_perennial_mean*10, fill=factor(Drought),
                              ymin=(C4_perennial_mean-C4_perennial_se)*10, ymax=(C4_perennial_mean+C4_perennial_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")
# Annual grass
ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=annual_grass_mean*10, fill=factor(Drought),
                              ymin=(annual_grass_mean-annual_grass_se)*10, ymax=(annual_grass_mean+annual_grass_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")
# vulpia
ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=vulpia_mean*10, fill=factor(Drought),
                              ymin=(vulpia_mean-vulpia_se)*10, ymax=(vulpia_mean+vulpia_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")

ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=bromus_mean*10, fill=factor(Drought),
                              ymin=(bromus_mean-bromus_se)*10, ymax=(bromus_mean+bromus_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")
# GUSA
ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=GUSA_mean*10, fill=factor(Drought),
                              ymin=(GUSA_mean-GUSA_se)*10, ymax=(GUSA_mean+GUSA_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")
#ARFR
ggplot(anpp_tb_drt_means, aes(x= factor(Drought), y=ARFR_mean*10, fill=factor(Drought),
                              ymin=(ARFR_mean-ARFR_se)*10, ymax=(ARFR_mean+ARFR_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")

### Fort Keogh
### Calculate means and SE for ANPP
anpp_fk_plotmeans <- anpp_fk %>%
  group_by(site, date, block, plot) %>%
  summarize_at(vars(C4_perennial:litter,anpp),.funs="mean", na.rm=TRUE) %>%
  ungroup() %>%
  left_join(Plot_Key, by=c("site", "block", "plot"))


### means by drought treatment
anpp_fk_drt_means <- anpp_fk_plotmeans %>%
  group_by(site, date, Drought) %>%
  summarize_at(vars(C4_perennial:anpp),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(C3_perennial_mean=C3_perennial_fn1,C3_perennial_se=C3_perennial_fn2,
         C4_perennial_mean=C4_perennial_fn1,C4_perennial_se=C4_perennial_fn2,
         forb_mean=forb_fn1,forb_se=forb_fn2,
         annual_grass_mean=annual_grass_fn1,annual_grass_se=annual_grass_fn2,
         anpp_mean=anpp_fn1,anpp_se=anpp_fn2)

### means by grazing treatment
anpp_fk_grz_means <- anpp_fk_plotmeans %>%
  group_by(site, date, Grazing) %>%
  summarize_at(vars(C4_perennial:anpp),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(C3_perennial_mean=C3_perennial_fn1,C3_perennial_se=C3_perennial_fn2,
         C4_perennial_mean=C4_perennial_fn1,C4_perennial_se=C4_perennial_fn2,
         forb_mean=forb_fn1,forb_se=forb_fn2,
         annual_grass_mean=annual_grass_fn1,annual_grass_se=annual_grass_fn2,
         anpp_mean=anpp_fn1,anpp_se=anpp_fn2) %>%
  mutate(Grazing=factor(Grazing, levels=c("MLLMM","MMMMM","HHMMM")))

### Means by grazing and drought treatments
anpp_fk_drt_grz_means <- anpp_fk_plotmeans %>%
  group_by(site, date, Drought, Grazing) %>%
  summarize_at(vars(C4_perennial:anpp),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(C3_perennial_mean=C3_perennial_fn1,C3_perennial_se=C3_perennial_fn2,
         C4_perennial_mean=C4_perennial_fn1,C4_perennial_se=C4_perennial_fn2,
         forb_mean=forb_fn1,forb_se=forb_fn2,
         annual_grass_mean=annual_grass_fn1,annual_grass_se=annual_grass_fn2,
         anpp_mean=anpp_fn1,anpp_se=anpp_fn2) %>%
  mutate(Grazing=factor(Grazing, levels=c("MLLMM","MMMMM","HHMMM")))


### Plot drought means
drt_plot_fk <- ggplot(anpp_fk_drt_means, aes(x= factor(Drought), y=anpp_mean*10, fill=factor(Drought),
                                             ymin=(anpp_mean-anpp_se)*10, ymax=(anpp_mean+anpp_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)") +
  ylim(0,225)

pdf("..//..//figures//anpp_drought_means_FK_2018-2020.pdf", width=12, height=6, useDingbats = F)
print(drt_plot_fk)
dev.off()

### Plot grazing means
grz_plot_fk <- ggplot(anpp_fk_grz_means, aes(x= factor(Grazing), y=anpp_mean*10, fill=factor(Grazing),
                                             ymin=(anpp_mean-anpp_se)*10, ymax=(anpp_mean+anpp_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=grazingColor, name='Grazing\nTreatment') +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")

pdf("..//..//figures//anpp_grazing_means_FK_2018-2020.pdf", width=12, height=6, useDingbats = F)
print(grz_plot_fk)
dev.off()

### Plot drought x grazing means
grz_drt_plot_fk <- ggplot(anpp_fk_drt_grz_means, aes(x= factor(Drought), y=anpp_mean*10, fill=factor(Drought),
                                                     ymin=(anpp_mean-anpp_se)*10, ymax=(anpp_mean+anpp_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  facet_grid(Grazing~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)")

pdf("..//..//figures//anpp_grazing x drought_means_FK_2018-2020.pdf", width=12, height=12, useDingbats = F)
print(grz_drt_plot_fk)
dev.off()

### Run models
# Full models
anpp_model_fk_2019_20 <-   anova_t3(IndVars=c('date','Grazing','Drought'),
                                    DepVar='anpp',
                                    RndForm='~1 |block/Paddock/plot',
                                    Data=subset(anpp_fk_plotmeans, date %in% c(2019,2020))
)

anpp_model_fk_2018_20 <-   anova_t3(IndVars=c('date','Grazing','Drought'), ### Drought x Date (year) is significant so split by year
                                    DepVar='anpp',
                                    RndForm='~1 |block/Paddock/plot',
                                    Data=anpp_fk_plotmeans
)  

# Models by year
anpp_model_fk_2018 <-   lme(anpp ~ Grazing*Drought
                            , data=subset(anpp_fk_plotmeans, date == 2018)
                            , random = ~1 |block/Paddock/plot
                            , na.action = na.omit)
anova.lme(anpp_model_fk_2018, type="marginal")
anova(anpp_model_fk_2018)#????

anpp_model_fk_2019 <-   lme(anpp ~ Grazing*Drought
                            , data=subset(anpp_fk_plotmeans, date == 2019)
                            , random = ~1 |block/Paddock/plot
                            , na.action = na.omit)
anova.lme(anpp_model_fk_2019, type="marginal")
anova(anpp_model_fk_2019)#????

anpp_model_fk_2020 <-   lme(anpp ~ Drought*Grazing
                            , data=subset(anpp_fk_plotmeans, date == 2020)
                            , random = ~1 |block/Paddock/plot
                            , na.action = na.omit)

anova.lme(anpp_model_fk_2020, type="marginal")
anova(anpp_model_fk_2020)#????
Anova(anpp_model_fk_2020, type=3)



###
### Response Ratios
###

### Fort Keogh
anpp_fk_controls <- anpp_fk_plotmeans %>%
  filter(Drought==0) %>%
  group_by(site, date, block, Paddock) %>%
  summarize(control_anpp=mean(anpp, na.rm=T))

anpp_fk_rr <- anpp_fk_plotmeans %>%
  filter(Drought != 0) %>%
  dplyr::select(site:plot, anpp:Grazing) %>%
  full_join(anpp_fk_controls, by=c("site","date","block","Paddock")) %>%
  mutate(ln_rr=log(anpp/control_anpp))

anpp_fk_rr_means <- anpp_fk_rr %>%
  group_by(site, date, Drought) %>%
  summarize_at(vars(ln_rr),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(ln_rr_mean=fn1,ln_rr_se=fn2)

ggplot(anpp_fk_rr_means, aes(x= Drought, y=ln_rr_mean,
                             ymin=ln_rr_mean-ln_rr_se, ymax=ln_rr_mean+ln_rr_se))+
  geom_point(size=3) +
  geom_errorbar(width=0) +
  #  geom_smooth(method="lm",se=F, col="red") +
  geom_smooth(method="lm") +
  facet_grid(~date) +
  ylab("ln(ANPPm/ANPPc)") +
  xlab("Precipitation reduction (%)") +
  geom_hline(yintercept=0, lty=2)

rr_fk_plot <- ggplot(subset(anpp_fk_rr_means,date%in%c(2019,2020)), aes(x= Drought, y=ln_rr_mean, fill=factor(date), col=factor(date),lty=factor(date),
                                                                        ymin=ln_rr_mean-ln_rr_se, ymax=ln_rr_mean+ln_rr_se))+
  geom_errorbar(col="black",width=1, lty=1) +
  geom_point(col="black", pch=21, size=3) +
  scale_fill_manual(values=c("black","white")) +
  #  geom_smooth(method="lm",se=F, col="red") +
  geom_smooth(method="lm", se=F, col="black") +
  ylab("ln(ANPPm/ANPPc)") +
  xlab("Precipitation reduction (%)") +
  geom_hline(yintercept=0, col="grey")+
  ylim(-0.9, 0.35)



### Thunder Basin
anpp_tb_controls <- anpp_tb_plotmeans %>%
  filter(Drought==0) %>%
  group_by(site, date, block, Paddock) %>%
  summarize(control_anpp=mean(anpp, na.rm=T))

anpp_tb_rr <- anpp_tb_plotmeans %>%
  filter(Drought != 0) %>%
  dplyr::select(site:plot, anpp:Grazing) %>%
  full_join(anpp_tb_controls, by=c("site","date","block","Paddock")) %>%
  mutate(ln_rr=log(anpp/control_anpp))

anpp_tb_rr_means <- anpp_tb_rr %>%
  group_by(site, date, Drought) %>%
  summarize_at(vars(ln_rr),.funs = c(mean,SE_function),na.rm=TRUE) %>% 
  rename(ln_rr_mean=fn1,ln_rr_se=fn2)

ggplot(anpp_tb_rr_means, aes(x= Drought, y=ln_rr_mean,
                             ymin=ln_rr_mean-ln_rr_se, ymax=ln_rr_mean+ln_rr_se))+
  geom_point(size=3) +
  geom_errorbar(width=0) +
  #  geom_smooth(method="lm",se=F, col="red") +
  geom_smooth(method="lm", col="red") +
  facet_grid(~date) +
  ylab("ln(ANPPm/ANPPc)") +
  xlab("Precipitation reduction (%)")



rr_tb_plot <- ggplot(subset(anpp_tb_rr_means,date%in%c(2019,2020)), aes(x= Drought, y=ln_rr_mean, fill=factor(date), col=factor(date),lty=factor(date),
                                                                        ymin=ln_rr_mean-ln_rr_se, ymax=ln_rr_mean+ln_rr_se))+
  geom_errorbar(col="black",width=1, lty=1) +
  geom_point(col="black", pch=21, size=3) +
  scale_fill_manual(values=c("black","white")) +
  #  geom_smooth(method="lm",se=F, col="red") +
  geom_smooth(method="lm", se=F, col="black") +
  ylab("ln(ANPPm/ANPPc)") +
  xlab("Precipitation reduction (%)") +
  geom_hline(yintercept=0, col="grey")+
  ylim(-0.9, 0.35)


pdf("..//..//figures//ANPP//anpp_drought_response ratios_TB_2019-2020.pdf", width=8, height=6, useDingbats = F)
print(rr_tb_plot)
dev.off()



###
### Calculating relative functional group biomass
###
tb_prop_fxn_anpp <- anpp_tb_plot_means %>%
  group_by(Year, Grazing) %>%
  summarize(ANPP_gm2=mean(ANPP_gm2, na.rm=T),
            C3P_gm2=mean(C3P_gm2, na.rm=T),
            C4P_gm2=mean(C4P_gm2, na.rm=T),
            Forb_gm2=mean(Forb_gm2, na.rm=T),
            AnnualGrass_gm2=mean(AnnualGrass_gm2, na.rm=T)) %>%
  ungroup() %>%
  mutate(p_C3P = C3P_gm2/(C3P_gm2+C4P_gm2+Forb_gm2+AnnualGrass_gm2),
         p_C4P = C4P_gm2/(C3P_gm2+C4P_gm2+Forb_gm2+AnnualGrass_gm2),
         p_Forb = Forb_gm2/(C3P_gm2+C4P_gm2+Forb_gm2+AnnualGrass_gm2),
         p_AnnualGrass = AnnualGrass_gm2/(C3P_gm2+C4P_gm2+Forb_gm2+AnnualGrass_gm2)) %>%
  dplyr::select(-ANPP_gm2:-AnnualGrass_gm2) %>%
  pivot_longer(cols=p_C3P:p_AnnualGrass, names_to="fxn_type", values_to = "rel_anpp")
    
tb_prop_fxn_anpp$Grazing <- factor(tb_prop_fxn_anpp$Grazing, levels=c("MLLMM", "MMMMM", "HHMMM"))

ggplot(tb_prop_fxn_anpp, aes(x=Year, y=rel_anpp, fill=fxn_type)) + geom_bar(stat="identity") + facet_grid(~Grazing)

fk_prop_fxn_anpp <- anpp_fk_plot_means %>%
  group_by(Year, Drought) %>%
  summarize(ANPP_gm2=mean(ANPP_gm2, na.rm=T),
            C3P_gm2=mean(C3P_gm2, na.rm=T),
            C4P_gm2=mean(C4P_gm2, na.rm=T),
            Forb_gm2=mean(Forb_gm2, na.rm=T),
            AnnualGrass_gm2=mean(AnnualGrass_gm2, na.rm=T)) %>%
  ungroup() %>%
  mutate(p_C3P = C3P_gm2/ANPP_gm2,
         p_C4P = C4P_gm2/ANPP_gm2,
         p_Forb = Forb_gm2/ANPP_gm2,
         p_AnnualGrass = AnnualGrass_gm2/ANPP_gm2) %>%
  dplyr::select(-ANPP_gm2:-AnnualGrass_gm2) %>%
  pivot_longer(cols=p_C3P:p_AnnualGrass, names_to="fxn_type", values_to = "rel_anpp")



ggplot(fk_prop_fxn_anpp, aes(x=factor(Drought), y=rel_anpp, fill=fxn_type)) + geom_bar(stat="identity") + facet_grid(~Year) +
  scale_fill_brewer(palette="Set1")



?geom_bar


### extra

ggplot(anpp_fk_rr_means, aes(x= Drought, y=ln_rr_mean, fill=factor(Drought),
                             ymin=(anpp_mean-anpp_se)*10, ymax=(anpp_mean+anpp_se)*10))+
  geom_col(col="black") +
  geom_errorbar(width=0) +
  scale_fill_manual(values=droughtColor, name='Drought\nTreatment') +
  #  geom_smooth(method="lm",se=F) +
  facet_grid(~date) +
  ylab("Aboveground Net Primary Productivity (g/m2)") +
  xlab("Precipitation reduction (%)") +
  ylim(0,225)





exp(tb_2018_C4P_yhat_99)



xhat <- 1:99
yhat <- tb_2018_C4P_drought_slope*xhat+tb_2018_C4P_drought_int
yhat_back <- exp(yhat)+min(anpp_tb_rr$pchange_C4P) +.01
hats <- bind_cols(xhat=xhat, yhat=yhat, yhat_back=yhat_back)
plot(xhat, exp(yhat)+min(anpp_tb_rr$pchange_C4P) +.01)

x_vec <- 1:99
y_vec <- log(x_vec)*(-.1)+0.2+runif(99, -.01,.01)

dummy_df <- data.frame(x=x_vec, y=y_vec)

#dummy_model <- lm(log(y+min(y)+.01)~x, data=dummy_df)
dummy_model <- lm(log(y_vec-min(y_vec)+.01)~x, data=dummy_df)
summary(dummy_model)

int_hat <- dummy_model$coefficients[1]
slope_hat <- dummy_model$coefficients[2]

exp(int_hat+slope_hat*40)
xhat <- 1:99
yhat <- slope_hat*xhat+int_hat
yhat_back <- exp(yhat)+min(y_vec)+.01
hats <- bind_cols(xhat=xhat, yhat=yhat, yhat_back=yhat_back)

test1 <- int_hat + slope_hat*99
exp(test1) + min(y_vec)+.01  
ggplot(dummy_df, aes(x_vec, y_vec)) + geom_point() +
  geom_line(inherit.aes = F, data=hats, aes(x=xhat, y=yhat_back))


ggplot(filter(anpp_tb_rr, Year==2018 & pchange_C4P < 6), aes(x=Drought, y=pchange_C4P)) +geom_point() + geom_smooth(method="lm") + ylim(-2,2) +
  geom_line(inherit.aes = F, data=hats, aes(x=xhat, y=yhat_back))
}






tb_anpp_2019_lme$coefficients$fixed[2]
z <- summary(tb_anpp_2019_lme)
summary(tb_anpp_2019_lme)$coef[,"Std. Error"]
z$coefficients
z[20]

str(z)
str(tb_anpp_2019_lme)

idy(tb_anpp_2019_lme)

anova(tb_anpp_2023_lme)
test <- as.data.frame(anova(tb_anpp_2023_lme))
}

#### Extra code from messing round May 13 2024
ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=exp(lnrr_ANPP))) + geom_point()


#ln trial
lnrr_mod <- lm(lnrr_ANPP ~ Drought, data=filter(anpp_tb_rr, Year==2020))
ln_trendline <- data.frame(xhat=filter(anpp_tb_rr, Year==2020)$Drought, yhat=lnrr_mod$fitted.values)
summary(lnrr_mod)

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=lnrr_ANPP)) + geom_point() +
  geom_path(data=ln_trendline, inherit.aes=F, aes(x=xhat, y=yhat), col="red")

#non ln trial
rr_mod <- lm(exp(lnrr_ANPP) ~ Drought, data=filter(anpp_tb_rr, Year==2020))
rr_trendline <- data.frame(xhat=filter(anpp_tb_rr, Year==2020)$Drought, yhat=rr_mod$fitted.values)
summary(rr_mod)

ggplot(filter(anpp_tb_rr, Year==2020), aes(x=Drought, y=exp(lnrr_ANPP))) + geom_point() +
  geom_path(data=rr_trendline, inherit.aes=F, aes(x=xhat, y=yhat), col="red")


exp(summary(lnrr_mod)$coef[2,1])
