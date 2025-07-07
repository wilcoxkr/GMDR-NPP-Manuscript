### Generating site-level data for methods
###
### Kevin Wilcox (k_wilcox@uncg.edu)
### Created July 3, 2025

###
library(tidyverse)

setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\")

### Read in plot keys from Sally's script "cleaning fk and tb ppt data_v2.R"
PlotKey_FK<-read.csv("..//Plot_Treatment_Key_020424.csv") %>% 
  filter(site=="FK") %>% 
  mutate(RainAmount=100-rainfall_reduction)
PlotKey_TB<-read.csv("..//Plot_Treatment_Key_020424.csv") %>% 
  filter(site=="TB") %>% 
  mutate(RainAmount=100-rainfall_reduction)

###
### Soil texture
###

### Thunder Basin
tb_soil_texture_df <- read.csv("soil_texture\\tb_soilTexture_2019.csv")

tb_soil_texture_block_means <- tb_soil_texture_df %>%
  group_by(site_code, block) %>%
  summarize(sand_perc = mean(sand_perc, na.rm=T),
            silt_perc = mean(silt_perc, na.rm=T),
            clay_perc = mean(clay_perc, na.rm=T),
            .groups='drop'
  )

### Block 1: Sandy loam
### Block 2: Sandy clay loam
### Block 3: Clay

### Fort Keogh
fk_soil_texture_df <- read.csv("soil_texture\\fk_soilTexture.csv") %>%
  rename(plot=plot_num) %>%
  left_join(PlotKey_FK, by=c("site", "plot"))

fk_soil_texture_block_means <- fk_soil_texture_df %>%
  group_by(site, block) %>%
  summarize(sand_perc = mean(sand_perc, na.rm=T),
            silt_perc = mean(silt_perc, na.rm=T),
            clay_perc = mean(clay_perc, na.rm=T),
            .groups='drop'
  )

### Block 1: Sandy loam
### Block 2:Sandy loam
### Block 3: Sandy loam

###
### Precipitation and temperature stats for both sites
###


  tb_daymet <- read.csv("TB_Daymet_1980-2023.csv", skip=6) %>%
    dplyr::select(year, yday, prcp..mm.day.,tmax..deg.c.,tmin..deg.c.) %>%
    rename(YEAR=year, DOY=yday, PRCP=prcp..mm.day., TMAX=tmax..deg.c., TMIN=tmin..deg.c.) %>%
    mutate(TMED=(TMIN+TMAX)/2) # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  tb_monthly <- tb_daymet %>%
    mutate(DATE= as.Date(DOY-1, origin=paste0(YEAR, "-01-01")),
           MONTH= as.numeric(strftime(DATE, "%m")), 
           DAY=as.numeric(strftime(DATE,"%d"))) %>%
    group_by(YEAR, MONTH) %>%
    summarize(TMAX=mean(TMAX, na.rm=T),
              TMIN=mean(TMIN, na.rm=T),
              TMED=mean(TMED, na.rm=T),
              PRCP=sum(PRCP, na.rm=T)
    )%>%
    ungroup()
  
  tb_annual <- tb_monthly %>%
    group_by(MONTH, .drop=T) %>%
    summarize(TMAX = mean(TMAX, na.rm=T),
              TMIN = mean(TMIN, na.rm=T),
              PRCP = mean(PRCP, na.rm=T)
    )
  
  ggplot(tb_annual, aes(as.factor(MONTH), TMAX))+ geom_col()
  max(tb_annual$TMAX)
  ggplot(tb_annual, aes(as.factor(MONTH), TMIN))+ geom_col()
  min(tb_annual$TMIN)
  ggplot(tb_annual, aes(as.factor(MONTH), PRCP))+ geom_col()
  min(tb_annual$TMIN)
  
  fk_daymet <- read.csv("FK_Daymet_1980-2023.csv", skip=6) %>%
    dplyr::select(year, yday, prcp..mm.day.,tmax..deg.c.,tmin..deg.c.) %>%
    rename(YEAR=year, DOY=yday, PRCP=prcp..mm.day., TMAX=tmax..deg.c., TMIN=tmin..deg.c.) %>%
    mutate(TMED=(TMIN+TMAX)/2) # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  fk_monthly <- fk_daymet %>%
    mutate(DATE= as.Date(DOY-1, origin=paste0(YEAR, "-01-01")),
           MONTH= as.numeric(strftime(DATE, "%m")), 
           DAY=as.numeric(strftime(DATE,"%d"))) %>%
    group_by(YEAR, MONTH) %>%
    summarize(TMAX=mean(TMAX, na.rm=T),
              TMIN=mean(TMIN, na.rm=T),
              TMED=mean(TMED, na.rm=T),
              PRCP=sum(PRCP, na.rm=T)
    )%>%
    ungroup()
  
  fk_annual <- fk_monthly %>%
    group_by(MONTH, .drop=T) %>%
    summarize(TMAX = mean(TMAX, na.rm=T),
              TMIN = mean(TMIN, na.rm=T),
              PRCP = mean(PRCP, na.rm=T)
    )
  
  ggplot(fk_annual, aes(as.factor(MONTH), TMAX))+ geom_col()
  max(fk_annual$TMAX)
  ggplot(fk_annual, aes(as.factor(MONTH), TMIN))+ geom_col()
  min(fk_annual$TMIN)
  ggplot(fk_annual, aes(as.factor(MONTH), PRCP))+ geom_col()
  min(fk_annual$TMIN)

