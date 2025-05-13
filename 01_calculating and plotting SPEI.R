### 
### Calculating SPEI for two sites and drought treatements
###
### Author: Kevin Wilcox (k_wilcox@uncg.edu)
### Created May 12 2025; last updated May 12 2025

### Set up workspace
library(tidyverse)
library(lubridate)
library(ggthemes)
library(meteoland)
library(SPEI)
library(ggnewscale)
library(ggfortify)
library(ggridges)
rm(list=ls()) # clean up

# Set working directory
setwd("C:\\Users\\wilco\\OneDrive - UNCG\\Current projects\\GMDR\\data\\precipitation\\") # wilcox personal laptop
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\precipitation\\") # wilcox desktop UNCG

### REad in plot keys from Sally's script "cleaning fk and tb ppt data_v2.R"
PlotKey_FK<-read.csv("..//Plot_Treatment_Key_020424.csv") %>% 
  filter(site=="FK") %>% 
  mutate(RainAmount=100-rainfall_reduction)
PlotKey_TB<-read.csv("..//Plot_Treatment_Key_020424.csv") %>% 
  filter(site=="TB") %>% 
  mutate(RainAmount=100-rainfall_reduction)

###
### Calculate SPEI for each site with Daymet data
###
{
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

### Calculate PET and balance between PET and PPT
tb_spei <- tb_monthly %>%
  mutate(PET=thornthwaite(TMED, 43.30)) %>%
  mutate(BAL=PRCP-PET) #%>%
  
fk_spei<- fk_monthly %>%
  mutate(PET=thornthwaite(TMED, 43.30)) %>%
  mutate(BAL=PRCP-PET)

### Calculate SPEI
tb_spei_12 <- spei(tb_spei$BAL, 12)
fk_spei_12 <- spei(fk_spei$BAL, 12)

### Turn SPEI object into data frame for plotting
tb_spei_df <- data.frame(
  tb_monthly[,c("YEAR","MONTH")],
  SPEI = tb_spei_12$fitted
  ) %>%
  mutate(DATE=as.Date(paste(YEAR, MONTH, "01", sep="-"))) %>%
  mutate(direction = ifelse(SPEI>0,"above","below"))

fk_spei_df <- data.frame(
  fk_monthly[,c("YEAR","MONTH")],
  SPEI = fk_spei_12$fitted
  ) %>%
  mutate(DATE=as.Date(paste(YEAR, MONTH, "01", sep="-"))) %>%
  mutate(direction = ifelse(SPEI>0,"above","below"))

### Plot SPEI from 1980-2025 (just for fun at this point)
par(mfrow=c(2,1))
ggplot(tb_spei_df, aes(x=DATE, y=SPEI, fill=direction)) +
  geom_col(size=2, width=31) +
  scale_fill_manual(values=c("blue","red"))
ggplot(fk_spei_df, aes(x=DATE, y=SPEI, fill=direction)) +
  geom_col(size=2, width=31) +
  scale_fill_manual(values=c("blue","red"))

  }

###
### Calculate SPEI for our drought treatments
###
{
### Thunder basin - Calculating SPEI
tb_trt_ppt_raw <- read.csv("TBPrecip_edited.csv") %>% 
  select (-X, -Notes) %>% 
  pivot_longer(Block_1:Block_3, names_to="Block", values_to = "RainoutOn") %>% 
  separate(Block, into=c("Trash", "block"), sep="_") %>% 
  select(-Trash) %>% 
  mutate(block=as.integer(block)) %>% 
  left_join(PlotKey_TB) %>% 
  mutate(ppt_trt=ifelse(RainoutOn=="On", ppt_mm*(RainAmount/100), ppt_mm)) 

## reorganize to get one ppt value for each treatment and each day
tb_trt_ppt <- tb_trt_ppt_raw %>%
  mutate(DATE=as.Date(date, "%m/%d/%y")) %>%
  dplyr::select(site, DATE, year, month, day, doy, ppt_mm, block, paddock, plot, rainfall_reduction, ppt_trt) %>%
  rename(ppt_amb = ppt_mm, Drought=rainfall_reduction) %>%
  group_by(site, DATE, year, month, day, doy, Drought) %>%
  summarize(ppt_amb = mean(ppt_amb, na.rm=T),
            ppt_trt = mean(ppt_trt, na.rm=T)) %>%
  ungroup()

ggplot(filter(tb_trt_ppt,year==2020), aes(DATE, ppt_trt, col=as.factor(Drought))) +
  geom_path()

tb_trt_monthly_ppt <- tb_trt_ppt %>%
  group_by(year, month, Drought) %>%
  summarize(ppt_trt = sum(ppt_trt)) %>%
  ungroup() %>%
  mutate(DATE=as.Date(paste(year, month, "01", sep="-"))) 
  
ggplot(tb_trt_monthly_ppt, aes(DATE, ppt_trt, fill=as.factor(Drought))) +
  geom_col(position=position_dodge())

tb_trt_monthly_ppt_temp <- tb_monthly %>%
  mutate(DATE=as.Date(paste(YEAR, MONTH, "01", sep="-"))) %>%
#  filter(DATE >= "2017-01-01" & DATE <= "2023-12-31") %>%
  full_join(tb_trt_monthly_ppt, by="DATE") %>%
  mutate(PRCP_2 = ifelse(is.na(ppt_trt), PRCP, ppt_trt)) %>%
  mutate(Drought_2 = ifelse(is.na(Drought), "met_station", Drought)) %>%
  dplyr::select(YEAR, MONTH, DATE, Drought_2, TMAX, TMIN, TMED, PRCP_2) %>%
  mutate(PET=thornthwaite(TMED, 43.30)) %>%
  mutate(BAL=PRCP_2-PET)

### Calculate SPEI for each Drought trt -- LOOP
tb_trt_vec <- as.character(unique(tb_trt_monthly_ppt$Drought))
tb_spei_master <- {}

for(TRT in 1:length(tb_trt_vec)){
  met_temp <- filter(tb_trt_monthly_ppt_temp, Drought_2 %in% c("met_station",tb_trt_vec[TRT]))
  
  spei_12_temp <- spei(met_temp$BAL, 12)
  
  # Turn SPEI object into data frame for plotting
  spei_df_temp <- data.frame(
    site="TB",
    met_temp[,c("YEAR","MONTH","DATE")],
    Drought=as.numeric(tb_trt_vec[TRT]),
    SPEI = spei_12_temp$fitted
    ) %>%
    mutate(direction = ifelse(SPEI>0,"above","below"))
  
  tb_spei_master <- rbind(tb_spei_master, spei_df_temp)
  rm(met_temp, spei_12_temp, spei_df_temp)
  }
  
### End LOOP


###
### Fort Keogh - calculating SPEI
fk_trt_ppt_raw <- read.csv("FKPrecip_edited.csv") %>% 
  select (-X, -PRCP_ATTRIBUTES , -Notes) %>% 
  pivot_longer(Block_1:Block_3, names_to="Block", values_to = "RainoutOn") %>% 
  separate(Block, into=c("Trash", "block"), sep="_") %>% 
  select(-Trash) %>% 
  mutate(block=as.integer(block)) %>% 
  left_join(PlotKey_FK) %>% 
  mutate(ppt_trt=ifelse(RainoutOn=="On", PRCP*(RainAmount/100), PRCP)) 

## reorganize to get one ppt value for each treatment and each day
fk_trt_ppt <- fk_trt_ppt_raw %>%
  mutate(DATE=as.Date(date, "%m/%d/%y")) %>%
  rename(ppt_mm=PRCP) %>%
  dplyr::select(site, DATE, year, month, day, doy, ppt_mm, block, paddock, plot, rainfall_reduction, ppt_trt) %>%
  rename(ppt_amb = ppt_mm, Drought=rainfall_reduction) %>%
  group_by(site, DATE, year, month, day, doy, Drought) %>%
  summarize(ppt_amb = mean(ppt_amb, na.rm=T),
            ppt_trt = mean(ppt_trt, na.rm=T)) %>%
  ungroup()

ggplot(filter(fk_trt_ppt,year==2020), aes(DATE, ppt_trt, col=as.factor(Drought))) +
  geom_path()

fk_trt_monthly_ppt <- fk_trt_ppt %>%
  group_by(year, month, Drought) %>%
  summarize(ppt_trt = sum(ppt_trt)) %>%
  ungroup() %>%
  mutate(DATE=as.Date(paste(year, month, "01", sep="-"))) 

ggplot(fk_trt_monthly_ppt, aes(DATE, ppt_trt, fill=as.factor(Drought))) +
  geom_col(position=position_dodge())

fk_trt_monthly_ppt_temp <- fk_monthly %>%
  mutate(DATE=as.Date(paste(YEAR, MONTH, "01", sep="-"))) %>%
  #  filter(DATE >= "2017-01-01" & DATE <= "2023-12-31") %>%
  full_join(fk_trt_monthly_ppt, by="DATE") %>%
  mutate(PRCP_2 = ifelse(is.na(ppt_trt), PRCP, ppt_trt)) %>%
  mutate(Drought_2 = ifelse(is.na(Drought), "met_station", Drought)) %>%
  dplyr::select(YEAR, MONTH, DATE, Drought_2, TMAX, TMIN, TMED, PRCP_2) %>%
  mutate(PET=thornthwaite(TMED, 43.30)) %>%
  mutate(BAL=PRCP_2-PET)

### Calculate SPEI for each Drought trt -- LOOP
fk_trt_vec <- as.character(unique(fk_trt_monthly_ppt$Drought))
fk_spei_master <- {}

for(TRT in 1:length(fk_trt_vec)){
  met_temp <- filter(fk_trt_monthly_ppt_temp, Drought_2 %in% c("met_station",fk_trt_vec[TRT]))
  
  spei_12_temp <- spei(met_temp$BAL, scale=12)

  # Turn SPEI object into data frame for plotting
  spei_df_temp <- data.frame(
    site="FK",
    met_temp[,c("YEAR","MONTH","DATE")],
    Drought=as.numeric(fk_trt_vec[TRT]),
    SPEI = as.vector(spei_12_temp$fitted)
  ) %>%
    mutate(direction = ifelse(SPEI>0,"above","below"))
  
  fk_spei_master <- rbind(fk_spei_master, spei_df_temp)
  rm(met_temp, spei_12_temp, spei_df_temp)
}

### End LOOP

}

###
### Plot long-term SPEI with control ppt and then treatment SPEI on top
###
{
  
### Combine spei datasets for both sites
spei_df_all <- fk_spei_master %>%
  bind_rows(tb_spei_master) %>%
  year_month  

droughtColor <- c('#feedde', '#fdbe85', '#fd8d3c', '#e6550d', '#a63603')

spei_fig <- ggplot(filter(spei_df_all, Drought==0 & YEAR %in% 2018:2023), aes(x=DATE, y=SPEI, fill=direction)) +
  geom_col(size=2, width=30) +
  scale_fill_manual(values=c("dodgerblue","indianred")) +
  new_scale_fill() +
  geom_path(inherit.aes=F, data=filter(spei_df_all, DATE>="2019-05-01" & DATE<="2021-12-28"), 
             aes(x=DATE, y=SPEI, col=as.factor(Drought), linetype=as.factor(Drought))) +
  geom_point(inherit.aes=F, data=filter(spei_df_all, DATE>="2019-05-01" & DATE<="2021-12-28"), 
             aes(x=DATE, y=SPEI, fill=as.factor(Drought), shape=as.factor(Drought)), col="black") +
  scale_shape_manual(values=21:25) +
  scale_fill_manual(values=droughtColor[1:5]) +
  scale_color_manual(values=droughtColor[1:5]) +
  facet_grid(site~.) +
  scale_x_date(date_breaks = "month", date_labels="%b") +
  theme_few() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=8))
  
pdf("..//..//npp_manuscript//figures//SPEI figure_v1.pdf", width=9, height=5.4, useDingbats = F)
print(spei_fig)
dev.off()

}

###
### Calculating and plotting treatment precipitation within long-term distribution
###
{
  
  ### Fort Keogh
  ###
  
  ### Calculate long-term annual ppt distribution
  fk_rainyr_ppt <- fk_daymet %>%
    mutate(RAIN_YEAR=ifelse(DOY >= 213, YEAR+1, YEAR)) %>% ## Rain year is August 1 - July 31 (DOY 213 is Aug 1)
    group_by(RAIN_YEAR) %>%
    summarize(PRCP_RAINYR=sum(PRCP, na.rm=T)) %>%
    ungroup()
  
  ### Calculate rain year amount in 5th and 95th percentiles
  fk_p5 <- qnorm(0.05, mean=mean(fk_rainyr_ppt$PRCP_RAINYR), sd=sd(fk_rainyr_ppt$PRCP_RAINYR))
  fk_p95 <- qnorm(0.95, mean=mean(fk_rainyr_ppt$PRCP_RAINYR), sd=sd(fk_rainyr_ppt$PRCP_RAINYR))
  
  ### Create column for plotting geom_rug colors
  fk_rainyr_ppt <- fk_rainyr_ppt %>%
    mutate(yr_type = ifelse(PRCP_RAINYR <= fk_p5, "extreme dry",
                            ifelse(PRCP_RAINYR >= fk_p95, "extreme wet",
                                   "non-extreme")))
  fk_rainyr_ppt$yr_type <- factor(fk_rainyr_ppt$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))
  
  ### calculate probability of each drought treatment
  prob_key_for_plotting <- data.frame(
    rain_year = 2018:2023,
    prob_4_plot = seq(.2, .02, length.out=6)
  )
  
  fk_trt_rainyr <- fk_trt_ppt %>%
    mutate(rain_year=ifelse(month %in% 8:12, year+1, year)) %>% ## Rain year is August 1 - July 31
    group_by(rain_year, Drought) %>%
    summarize(ppt_trt_rainyr=sum(ppt_trt, na.rm=T)) %>%
    ungroup() %>%
    left_join(prob_key_for_plotting, by="rain_year") %>%
    filter(rain_year %in% 2018:2023)

  ### Plot up probabilities

fk_dnorm <- data.frame(prcp_hat = seq(0,750,1),
                      probability = dnorm(seq(0,750,1), mean=mean(fk_rainyr_ppt$PRCP_RAINYR), sd=sd(fk_rainyr_ppt$PRCP_RAINYR))*100
                      ) %>%
                      mutate(yr_type = ifelse(prcp_hat <= fk_p5, "extreme dry",
                                              ifelse(prcp_hat >= fk_p95, "extreme wet",
                                              "non-extreme")))
fk_dnorm$yr_type <- factor(fk_dnorm$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

fk_prob_fig <-  ggplot(fk_dnorm, aes(prcp_hat, probability, fill=yr_type)) +
                    geom_col(width = 1) +
                  #  geom_line() +
                    theme_few(base_size=16) +
                    scale_fill_manual(values=c("indianred","grey75","dodgerblue"))+
                    scale_color_manual(values=c("indianred2","grey60","dodgerblue3")) +
                    geom_rug(data=fk_rainyr_ppt, sides="b", inherit.aes = F, aes(x=PRCP_RAINYR, col=yr_type)) +
                    new_scale_fill() +
                    geom_point(data=fk_trt_rainyr, inherit.aes=F, aes(ppt_trt_rainyr, prob_4_plot, shape=as.factor(Drought), fill=as.factor(Drought)), col="black", size=3) +
                    scale_fill_manual(values=droughtColor) +
                    scale_shape_manual(values=21:25) +
                    xlim(0,675)
  
pdf("..//..//npp_manuscript//figures//fk rainyear distribution fig_v1.pdf", width=6, height=3.3, useDingbats = F)
print(fk_prob_fig)
dev.off()

### Thunder Basin
###

### Calculate long-term annual ppt distribution
tb_rainyr_ppt <- tb_daymet %>%
  mutate(RAIN_YEAR=ifelse(DOY >= 213, YEAR+1, YEAR)) %>% ## Rain year is August 1 - July 31 (DOY 213 is Aug 1)
  group_by(RAIN_YEAR) %>%
  summarize(PRCP_RAINYR=sum(PRCP, na.rm=T)) %>%
  ungroup()

### Calculate rain year amount in 5th and 95th percentiles
tb_p5 <- qnorm(0.05, mean=mean(tb_rainyr_ppt$PRCP_RAINYR), sd=sd(tb_rainyr_ppt$PRCP_RAINYR))
tb_p95 <- qnorm(0.95, mean=mean(tb_rainyr_ppt$PRCP_RAINYR), sd=sd(tb_rainyr_ppt$PRCP_RAINYR))

### Create column for plotting geom_rug colors
tb_rainyr_ppt <- tb_rainyr_ppt %>%
  mutate(yr_type = ifelse(PRCP_RAINYR <= tb_p5, "extreme dry",
                          ifelse(PRCP_RAINYR >= tb_p95, "extreme wet",
                                 "non-extreme")))
tb_rainyr_ppt$yr_type <- factor(tb_rainyr_ppt$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

### calculate probability of each drought treatment
prob_key_for_plotting <- data.frame(
  rain_year = 2018:2023,
  prob_4_plot = seq(.2, .02, length.out=6)
)

tb_trt_rainyr <- tb_trt_ppt %>%
  mutate(rain_year=ifelse(month %in% 8:12, year+1, year)) %>% ## Rain year is August 1 - July 31
  group_by(rain_year, Drought) %>%
  summarize(ppt_trt_rainyr=sum(ppt_trt, na.rm=T)) %>%
  ungroup() %>%
  left_join(prob_key_for_plotting, by="rain_year") %>%
  filter(rain_year %in% 2018:2023)

### Plot up probabilities

tb_dnorm <- data.frame(prcp_hat = seq(0,750,1),
                       probability = dnorm(seq(0,750,1), mean=mean(tb_rainyr_ppt$PRCP_RAINYR), sd=sd(tb_rainyr_ppt$PRCP_RAINYR))*100
) %>%
  mutate(yr_type = ifelse(prcp_hat <= tb_p5, "extreme dry",
                          ifelse(prcp_hat >= tb_p95, "extreme wet",
                                 "non-extreme")))
tb_dnorm$yr_type <- factor(tb_dnorm$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

tb_prob_fig <-  ggplot(tb_dnorm, aes(prcp_hat, probability, fill=yr_type)) +
  geom_col(width = 1) +
  #  geom_line() +
  theme_few(base_size=16) +
  scale_fill_manual(values=c("indianred","grey75","dodgerblue"))+
  scale_color_manual(values=c("indianred2","grey60","dodgerblue3")) +
  geom_rug(data=tb_rainyr_ppt, sides="b", inherit.aes = F, aes(x=PRCP_RAINYR, col=yr_type)) +
  new_scale_fill() +
  geom_point(data=tb_trt_rainyr, inherit.aes=F, aes(ppt_trt_rainyr, prob_4_plot, shape=as.factor(Drought), fill=as.factor(Drought)), col="black", size=3) +
  scale_fill_manual(values=droughtColor) +
  scale_shape_manual(values=21:25) +
  xlim(0,700)

pdf("..//..//npp_manuscript//figures//tb rainyear distribution fig_v1.pdf", width=6, height=3.3, useDingbats = F)
print(tb_prob_fig)
dev.off()


###
### Calculate two-year rain ppt

### Fort Keogh
###

fk_2yr_ppt <- fk_rainyr_ppt %>%
  mutate(ppt_prev1 = c(NA,PRCP_RAINYR[1:(length(PRCP_RAINYR)-1)])) %>%
  mutate(ppt_2yr = PRCP_RAINYR+ppt_prev1) %>%
  dplyr::select(-yr_type)
  
### Calculate 2-year rain year amount in 5th and 95th percentiles
fk_2yr_p5 <- qnorm(0.05, mean=mean(fk_2yr_ppt$ppt_2yr,na.rm=T), sd=sd(fk_2yr_ppt$ppt_2yr,na.rm=T))
fk_2yr_p95 <- qnorm(0.95, mean=mean(fk_2yr_ppt$ppt_2yr,na.rm=T), sd=sd(fk_2yr_ppt$ppt_2yr,na.rm=T))

### Create column for plotting geom_rug colors
fk_2yr_ppt <- fk_2yr_ppt %>%
  mutate(yr_type = ifelse(ppt_2yr <= fk_2yr_p5, "extreme dry",
                          ifelse(ppt_2yr >= fk_2yr_p95, "extreme wet",
                                 "non-extreme")))
fk_2yr_ppt$yr_type <- factor(fk_2yr_ppt$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

fk_trt_2yr <- fk_trt_ppt %>%
  mutate(rain_year=ifelse(month %in% 8:12, year+1, year)) %>% ## Rain year is August 1 - July 31
  group_by(rain_year, Drought) %>%
  summarize(ppt_trt_rainyr=sum(ppt_trt, na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(names_from=Drought, values_from=ppt_trt_rainyr, names_prefix="ppt_") %>%
  mutate(prev1_0 = c(NA,ppt_0[1:(length(ppt_0)-1)]),
         prev1_25 = c(NA,ppt_25[1:(length(ppt_25)-1)]),
         prev1_50 = c(NA,ppt_50[1:(length(ppt_50)-1)]),
         prev1_75 = c(NA,ppt_75[1:(length(ppt_75)-1)]),
         prev1_99 = c(NA,ppt_99[1:(length(ppt_99)-1)])
         ) %>%
  mutate(ppt_2yr_0 = ppt_0+prev1_0,
         ppt_2yr_25 = ppt_25+prev1_25,
         ppt_2yr_50 = ppt_50+prev1_50,
         ppt_2yr_75 = ppt_75+prev1_75,
         ppt_2yr_99 = ppt_99+prev1_99
         ) %>%
  dplyr::select(rain_year,ppt_2yr_0:ppt_2yr_99) %>%
  pivot_longer(cols=ppt_2yr_0:ppt_2yr_99, names_to="Drought", values_to="ppt_2yr", names_prefix="ppt_2yr_")

### Plot up probabilities
fk_2yr_dnorm <- data.frame(prcp_hat = seq(200,1200,1),
                       probability = dnorm(seq(200,1200,1), mean=mean(fk_2yr_ppt$ppt_2yr, na.rm=T), sd=sd(fk_2yr_ppt$ppt_2yr, na.rm=T))*100
) %>%
  mutate(yr_type = ifelse(prcp_hat <= fk_2yr_p5, "extreme dry",
                          ifelse(prcp_hat >= fk_2yr_p95, "extreme wet",
                                 "non-extreme"))) %>%
  mutate(site="FK")
fk_2yr_dnorm$yr_type <- factor(fk_2yr_dnorm$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

### Just looking at the distribution
ggplot(fk_2yr_dnorm, aes(prcp_hat, probability, fill=yr_type)) +
  geom_col(width = 1) +
  #  geom_line() +
  theme_few(base_size=16) +
  scale_fill_manual(values=c("indianred","grey75","dodgerblue"))+
  scale_color_manual(values=c("indianred2","grey60","dodgerblue3")) +
  geom_rug(data=fk_2yr_ppt, sides="b", inherit.aes = F, aes(x=ppt_2yr, col=yr_type)) +
  new_scale_fill() +
  geom_point(data=filter(fk_trt_2yr,rain_year==2020), inherit.aes=F, aes(ppt_2yr, 0.02, shape=as.factor(Drought), fill=as.factor(Drought)), col="black", size=3) +
  scale_fill_manual(values=droughtColor) +
  scale_shape_manual(values=21:25) #+
#  xlim(0,675)

### Thunder Basin
###
tb_2yr_ppt <- tb_rainyr_ppt %>%
  mutate(ppt_prev1 = c(NA,PRCP_RAINYR[1:(length(PRCP_RAINYR)-1)])) %>%
  mutate(ppt_2yr = PRCP_RAINYR+ppt_prev1) %>%
  dplyr::select(-yr_type)

### Calculate 2-year rain year amount in 5th and 95th percentiles
tb_2yr_p5 <- qnorm(0.05, mean=mean(tb_2yr_ppt$ppt_2yr,na.rm=T), sd=sd(tb_2yr_ppt$ppt_2yr,na.rm=T))
tb_2yr_p95 <- qnorm(0.95, mean=mean(tb_2yr_ppt$ppt_2yr,na.rm=T), sd=sd(tb_2yr_ppt$ppt_2yr,na.rm=T))

### Create column for plotting geom_rug colors
tb_2yr_ppt <- tb_2yr_ppt %>%
  mutate(yr_type = ifelse(ppt_2yr <= tb_2yr_p5, "extreme dry",
                          ifelse(ppt_2yr >= tb_2yr_p95, "extreme wet",
                                 "non-extreme")))
tb_2yr_ppt$yr_type <- factor(tb_2yr_ppt$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

tb_trt_2yr <- tb_trt_ppt %>%
  mutate(rain_year=ifelse(month %in% 8:12, year+1, year)) %>% ## Rain year is August 1 - July 31
  group_by(rain_year, Drought) %>%
  summarize(ppt_trt_rainyr=sum(ppt_trt, na.rm=T)) %>%
  ungroup() %>%
  pivot_wider(names_from=Drought, values_from=ppt_trt_rainyr, names_prefix="ppt_") %>%
  mutate(prev1_0 = c(NA,ppt_0[1:(length(ppt_0)-1)]),
         prev1_25 = c(NA,ppt_25[1:(length(ppt_25)-1)]),
         prev1_50 = c(NA,ppt_50[1:(length(ppt_50)-1)]),
         prev1_75 = c(NA,ppt_75[1:(length(ppt_75)-1)]),
         prev1_99 = c(NA,ppt_99[1:(length(ppt_99)-1)])
  ) %>%
  mutate(ppt_2yr_0 = ppt_0+prev1_0,
         ppt_2yr_25 = ppt_25+prev1_25,
         ppt_2yr_50 = ppt_50+prev1_50,
         ppt_2yr_75 = ppt_75+prev1_75,
         ppt_2yr_99 = ppt_99+prev1_99
  ) %>%
  dplyr::select(rain_year,ppt_2yr_0:ppt_2yr_99) %>%
  pivot_longer(cols=ppt_2yr_0:ppt_2yr_99, names_to="Drought", values_to="ppt_2yr", names_prefix="ppt_2yr_")

### Plot up probabilities
tb_2yr_dnorm <- data.frame(prcp_hat = seq(200,1200,1),
                           probability = dnorm(seq(200,1200,1), mean=mean(tb_2yr_ppt$ppt_2yr, na.rm=T), sd=sd(tb_2yr_ppt$ppt_2yr, na.rm=T))*100
) %>%
  mutate(yr_type = ifelse(prcp_hat <= tb_2yr_p5, "extreme dry",
                          ifelse(prcp_hat >= tb_2yr_p95, "extreme wet",
                                 "non-extreme"))) %>%
  mutate(site="TB")

tb_2yr_dnorm$yr_type <- factor(tb_2yr_dnorm$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

### Just looking at the distribution
ggplot(tb_2yr_dnorm, aes(prcp_hat, probability, fill=yr_type)) +
  geom_col(width = 1) +
  #  geom_line() +
  theme_few(base_size=16) +
  scale_fill_manual(values=c("indianred","grey75","dodgerblue"))+
  scale_color_manual(values=c("indianred2","grey60","dodgerblue3")) +
  geom_rug(data=tb_2yr_ppt, sides="b", inherit.aes = F, aes(x=ppt_2yr, col=yr_type)) +
  new_scale_fill() +
  geom_point(data=filter(tb_trt_2yr,rain_year==2020), inherit.aes=F, aes(ppt_2yr, 0.02, shape=as.factor(Drought), fill=as.factor(Drought)), col="black", size=3) +
  scale_fill_manual(values=droughtColor) +
  scale_shape_manual(values=21:25) #+

dnorm_2yr_both <- tb_2yr_dnorm %>%
  bind_rows(fk_2yr_dnorm) %>%
  mutate(site_num = ifelse(site=="FK",0.6, 0.3))
dnorm_2yr_both$yr_type <- factor(dnorm_2yr_both$yr_type, levels=c("extreme dry","non-extreme","extreme wet"))

trt_2yr_both <- fk_trt_2yr %>%
  mutate(site="FK") %>%
  bind_rows(
    tb_trt_2yr %>% mutate(site="TB")
  ) %>%
  mutate(site_num = ifelse(site=="FK",0.625, 0.325))

prob_2yr_fig <- ggplot(dnorm_2yr_both, aes(x=prcp_hat, y=site_num, height=probability, fill=yr_type, group=site_num)) +
                  geom_ridgeline_gradient() +
                  theme_few(base_size=16) +
                  scale_fill_manual(values=c("indianred","grey75","dodgerblue")) +
                  new_scale_fill() +
                  geom_point(data=filter(trt_2yr_both,rain_year==2020), inherit.aes=F, aes(ppt_2yr, site_num, shape=as.factor(Drought), fill=as.factor(Drought)), col="black", size=3) +
                  scale_fill_manual(values=droughtColor) +
                  scale_shape_manual(values=21:25) #+

pdf("..//..//npp_manuscript//figures//two year distribution fig_both sites_v1.pdf", width=6, height=3.3, useDingbats = F)
print(prob_2yr_fig)
dev.off()

}


