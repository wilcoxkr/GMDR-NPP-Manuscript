### Analyzing and plotting functional group ANPP 
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
library(sjPlot)
library(emmeans)
rm(list=ls()) # clean up

# Set working directory 
# setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\GMDR\\data\\anpp\\") # Kevin's laptop
# setwd("C:\\Users\\kwilcox4\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer
setwd("C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer UNCG
setwd("C:\\Users\\wilco\\OneDrive - UNCG\\Current projects\\GMDR\\data\\anpp\\") # Kevins office computer UNCG

# Set writing directory
figure_dir <- "C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\figures\\"
table_dir <- "C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\tables\\"
df_dir <- "C:\\Users\\K_WILCOX\\OneDrive - UNCG\\Current projects\\GMDR\\npp_manuscript\\data_sets\\"

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

### Read in ANPP cleaning script(has to have appropriate tags surrounding section to be sourced)
sourcePartial <- function(fn,startTag='#from here',endTag='#to here') {
  lines <- scan(fn, what=character(), sep="\n", quiet=TRUE)
  st<-grep(startTag,lines)
  en<-grep(endTag,lines)
  tc <- textConnection(lines[(st+1):(en-1)])
  source(tc)
  close(tc)
}

sourcePartial("C:\\Users\\k_wilcox\\OneDrive - UNCG\\Git projects\\Grazing-Management-for-Drought-Resilience\\GMDR-NPP-Manuscript\\03_anpp_bnpp plotting and analysis.R"
              , startTag='#BEGIN_ANPP_SOURCE', endTag='#END_ANPP_SOURCE')

rm(anpp_fk_18_23, anpp_tb_18, anpp_tb_18_21, anpp_tb_18_22, anpp_tb_18_23, anpp_tb_19_23, anpp_tb_22)

###
### Combine data into one df and write cleaned data to file
###

anpp_fxn_grps <- anpp_fk_plot_means %>%
  dplyr::select(-dead, -litter, -subshrub_shrub, -ANPP_gm2) %>%
  pivot_longer(cols=C4P_gm2:Forb_gm2,names_to="fxn_type") %>%
  bind_rows(
    anpp_tb_plot_means %>%
      dplyr::select(-StandingDead_gm2, -Vulpia_gm2, -Bromus_gm2, -ANPP_gm2) %>%
      pivot_longer(cols=C3P_gm2:AnnualGrass_gm2,names_to="fxn_type")
  ) %>%
  rename(ANPP_gm2=value)

with(anpp_fxn_grps, table(Year, Site, Drought))
with(anpp_fxn_grps, table(fxn_type, Year, Site, Block, Paddock)) ## Looks good

## take a look at the data
ggplot(filter(anpp_fxn_grps, fxn_type=="C4P_gm2"), aes(x=Drought, y=ANPP_gm2)) + geom_jitter(width=2) + geom_smooth(method="lm",se=F) + facet_grid(Site~Year)
ggplot(filter(anpp_fxn_grps, fxn_type=="C3P_gm2"), aes(x=Drought, y=ANPP_gm2)) + geom_jitter(width=2) + geom_smooth(method="lm",se=F) + facet_grid(Site~Year)
ggplot(filter(anpp_fxn_grps, fxn_type=="Forb_gm2"), aes(x=Drought, y=ANPP_gm2)) + geom_jitter(width=2) + geom_smooth(method="lm",se=F) + facet_grid(Site~Year)
ggplot(filter(anpp_fxn_grps, fxn_type=="AnnualGrass_gm2"), aes(x=Drought, y=ANPP_gm2)) + geom_jitter(width=2) + geom_smooth(method="lm",se=F) + facet_grid(Site~Year)

write.csv(anpp_fxn_grps, file=paste0(df_dir, "ANPP by fxn group_2018-2023_", Sys.Date(), ".csv"), row.names=F)  


###
### Running models -- raw ANPP by functional groups
###
{

  ### Fort Keogh
  ###
  
  ## C3p
  ##
  fk_c3p_model_full <- lme(log(ANPP_gm2+0.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="C3P_gm2")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  fk_c3p_anova <- anova.lme(fk_c3p_model_full, type="marginal")
  plot(fk_c3p_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(fk_c3p_model_full, abline = c(0,1)) ## qqplot
  hist(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="C3P_gm2")$ANPP_gm2)
  
  emtrends(fk_c3p_model_full, "Year", var="Drought")
  test(emtrends(fk_c3p_model_full, "Year", var="Drought"))
  # emmeans(fk_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(fk_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  
  # Save to writable tables
  fk_c3p_anova_df <- data.frame(effect=row.names(fk_c3p_anova), fk_c3p_anova, fxn_type="C3P", site="FK")
  fk_c3p_emtrends <- data.frame(test(emtrends(fk_c3p_model_full, "Year", var="Drought")), fxn_type="C3P", site="FK")
  
  ### Split by year to get R2 values for significant regressions 
  # 2019
  fk_c3p_2019_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                          , data=filter(anpp_fxn_grps, Year==2019 & Site=="FK" & fxn_type=="C3P_gm2")
                          , random = ~1 |Block/Paddock
                          , na.action = na.omit)
  anova.lme(fk_c3p_2019_lme, type="marginal")
  performance::r2(fk_c3p_2019_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT - Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2020
  fk_c3p_2020_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2020 & Site=="FK" & fxn_type=="C3P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c3p_2020_lme, type="marginal")
  performance::r2(fk_c3p_2020_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT.. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2021
  fk_c3p_2021_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2021 & Site=="FK" & fxn_type=="C3P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c3p_2021_lme, type="marginal")
  performance::r2(fk_c3p_2021_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2022
  fk_c3p_2022_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2022 & Site=="FK" & fxn_type=="C3P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c3p_2022_lme, type="marginal")
  performance::r2(fk_c3p_2022_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT .. Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2023
  fk_c3p_2023_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2023 & Site=="FK" & fxn_type=="C3P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c3p_2023_lme, type="marginal")
  performance::r2(fk_c3p_2023_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT .. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  ## C4p
  ##
  fk_c4p_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                           , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="C4P_gm2")
                           , random = ~1 |Block/Paddock/Plot
                           , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit)
  fk_c4p_anova <- anova.lme(fk_c4p_model_full, type="marginal")
  plot(fk_c4p_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(fk_c4p_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="C4P_gm2")$ANPP_gm2+.1))
  
  emtrends(fk_c4p_model_full, "Year", var="Drought")
  test(emtrends(fk_c4p_model_full, "Year", var="Drought"))
  # emmeans(fk_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(fk_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  
  # Save to writable tables
  fk_c4p_anova_df <- data.frame(effect=row.names(fk_c4p_anova), fk_c4p_anova, fxn_type="C4P", site="FK")
  fk_c4p_emtrends <- data.frame(test(emtrends(fk_c4p_model_full, "Year", var="Drought")), fxn_type="C4P", site="FK")
  
  # 2019
  fk_c4p_2019_lme <- lme(log(ANPP_gm2+.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2019 & Site=="FK" & fxn_type=="C4P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c4p_2019_lme, type="marginal")
  performance::r2(fk_c4p_2019_lme) # CAN'T HAVE PLOT IN RANDOM EFFECT .. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2020 (not significant)

  # 2021
  fk_c4p_2021_lme <- lme(log(ANPP_gm2+.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2021 & Site=="FK" & fxn_type=="C4P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c4p_2021_lme, type="marginal")
  performance::r2(fk_c4p_2021_lme) # 0.94 R2 DOESN'T SEEM CORRECT .. Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2022
  fk_c4p_2022_lme <- lme(log(ANPP_gm2+.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2022 & Site=="FK" & fxn_type=="C4P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c4p_2022_lme, type="marginal")
  performance::r2(fk_c4p_2022_lme) # 0.58 R2 DOESN'T SEEM CORRECT .. Marginal R2 considers only the variance of the fixed effects, which is what I want
  # 2023
  fk_c4p_2023_lme <- lme(log(ANPP_gm2+.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2023 & Site=="FK" & fxn_type=="C4P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_c4p_2023_lme, type="marginal")
  performance::r2(fk_c4p_2023_lme) # 0.75 R2 DOESN'T SEEM CORRECT .. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  ## Forb
  ##
  fk_forb_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                           , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="Forb_gm2")
                           , random = ~1 |Block/Paddock/Plot
                           , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit)
  fk_forb_anova <- anova.lme(fk_forb_model_full, type="marginal")
  plot(fk_forb_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(fk_forb_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="Forb_gm2")$ANPP_gm2+.1))
  
  emtrends(fk_forb_model_full, "Year", var="Drought")
  test(emtrends(fk_forb_model_full, "Year", var="Drought"))
  # emmeans(fk_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(fk_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
#  ggplot(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="Forb_gm2"), aes(x=Drought, y=ANPP_gm2, col=Grazing)) +geom_point()+geom_smooth(method="lm",se=F)+facet_grid(Site~Year)

  # Save to writable tables
  fk_forb_anova_df <- data.frame(effect=row.names(fk_forb_anova), fk_forb_anova,fxn_type="Forb",  site="FK")
  fk_forb_emtrends <- data.frame(test(emtrends(fk_forb_model_full, "Year", var="Drought")), fxn_type="Forb", site="FK")
  
  # 2019
  fk_forb_2019_lme <- lme(log(ANPP_gm2+.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2019 & Site=="FK" & fxn_type=="Forb_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_forb_2019_lme, type="marginal")
  performance::r2(fk_forb_2019_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT - Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2020
  fk_forb_2020_lme <- lme(log(ANPP_gm2+.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2020 & Site=="FK" & fxn_type=="Forb_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(fk_forb_2020_lme, type="marginal")
  performance::r2(fk_forb_2020_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT.. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  
  ## AnnualGrass
  ##
  fk_annualgrass_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="AnnualGrass_gm2")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  fk_annualgrass_anova <- anova.lme(fk_annualgrass_model_full, type="marginal")
  plot(fk_annualgrass_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(fk_annualgrass_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="AnnualGrass_gm2")$ANPP_gm2+.1))
  
  emtrends(fk_annualgrass_model_full, "Year", var="Drought")
  test(emtrends(fk_annualgrass_model_full, "Year", var="Drought"))
  # emmeans(fk_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(fk_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  #  ggplot(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="FK" & fxn_type=="AnnualGrass_gm2"), aes(x=Drought, y=ANPP_gm2, col=Grazing)) +geom_point()+geom_smooth(method="lm",se=F)+facet_grid(Site~Year)
  
  # Save to writable tables
  fk_annualgrass_anova_df <- data.frame(effect=row.names(fk_annualgrass_anova), fk_annualgrass_anova, fxn_type="AnnualGrass", site="FK")
  fk_annualgrass_emtrends <- data.frame(test(emtrends(fk_annualgrass_model_full, "Year", var="Drought")), fxn_type="AnnualGrass", site="FK")
  
  # 2022
  fk_annualgrass_2022_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                          , data=filter(anpp_fxn_grps, Year==2022 & Site=="FK" & fxn_type=="AnnualGrass_gm2")
                          , random = ~1 |Block/Paddock
                          , na.action = na.omit)
  anova.lme(fk_annualgrass_2022_lme, type="marginal")
  performance::r2(fk_annualgrass_2022_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT - Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  # 2023
  fk_annualgrass_2023_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                          , data=filter(anpp_fxn_grps, Year==2023 & Site=="FK" & fxn_type=="AnnualGrass_gm2")
                          , random = ~1 |Block/Paddock
                          , na.action = na.omit)
  anova.lme(fk_annualgrass_2023_lme, type="marginal")
  performance::r2(fk_annualgrass_2023_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT.. Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  
  ### Thunder Basin
  ###
  
  ## C3p
  ##
  tb_c3p_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                           , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="C3P_gm2")
                           , random = ~1 |Block/Paddock/Plot
                           , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit)
  tb_c3p_anova <- anova.lme(tb_c3p_model_full, type="marginal")
  plot(tb_c3p_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(tb_c3p_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="C3P_gm2")$ANPP_gm2+0.1))
  
  emtrends(tb_c3p_model_full, "Year", var="Drought")
  test(emtrends(tb_c3p_model_full, "Year", var="Drought"))
  emmeans(tb_c3p_model_full, "Grazing", by="Year")
  pairs(emmeans(tb_c3p_model_full, "Grazing", by="Year")) ## Just a higher c3p in MMMMM than MLLMM in 2023... probably spurious
  
  # Save to writable tables
  tb_c3p_anova_df <- data.frame(effect=row.names(tb_c3p_anova), tb_c3p_anova, fxn_type="C3P", site="TB")
  tb_c3p_emtrends <- data.frame(test(emtrends(tb_c3p_model_full, "Year", var="Drought")), fxn_type="C3P", site="TB")
  
  ### Split by year to get R2 values for significant regressions 
  # 2021 (only significant year from emtrends output)
  tb_c3p_2021_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                         , data=filter(anpp_fxn_grps, Year==2021 & Site=="TB" & fxn_type=="C3P_gm2")
                         , random = ~1 |Block/Paddock
                         , na.action = na.omit)
  anova.lme(tb_c3p_2021_lme, type="marginal")
  performance::r2(tb_c3p_2021_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT Marginal R2 considers only the variance of the fixed effects, which is what I want

  ## C4p
  ##
  tb_c4p_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                           , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="C4P_gm2")
                           , random = ~1 |Block/Paddock/Plot
                           , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                           , control=lmeControl(returnObject=TRUE)
                           , na.action = na.omit)
  tb_c4p_anova <- anova.lme(tb_c4p_model_full, type="marginal")
  plot(tb_c4p_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(tb_c4p_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="C4P_gm2")$ANPP_gm2+.1))
  
  # no significant interactions
    # emtrends(tb_c4p_model_full, "Year", var="Drought")
  # test(emtrends(tb_c4p_model_full, "Year", var="Drought"))
  # emmeans(tb_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(tb_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  
  # Save to writable tables
  tb_c4p_anova_df <- data.frame(effect=row.names(tb_c4p_anova), tb_c4p_anova, fxn_type="C4P", site="TB")
#  tb_c4p_emtrends <- data.frame(test(emtrends(tb_c4p_model_full, "Year", var="Drought")), fxn_type="C4P", site="TB")
  

  ## Forb
  ##
  tb_forb_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                            , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="Forb_gm2")
                            , random = ~1 |Block/Paddock/Plot
                            , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                            , control=lmeControl(returnObject=TRUE)
                            , na.action = na.omit)
  tb_forb_anova <- anova.lme(tb_forb_model_full, type="marginal")
  plot(tb_forb_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(tb_forb_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="Forb_gm2")$ANPP_gm2+.1))
  
  # NO SIGNIFICANT INTERACTIONS SO STOP HERE
    # emtrends(tb_forb_model_full, "Year", var="Drought")
  # test(emtrends(tb_forb_model_full, "Year", var="Drought"))
  # # emmeans(tb_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(tb_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  #  ggplot(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="Forb_gm2"), aes(x=Drought, y=ANPP_gm2, col=Grazing)) +geom_point()+geom_smooth(method="lm",se=F)+facet_grid(Site~Year)
  
  # Save to writable tables
  tb_forb_anova_df <- data.frame(effect=row.names(tb_forb_anova), tb_forb_anova,fxn_type="Forb",  site="TB")
 # tb_forb_emtrends <- data.frame(test(emtrends(tb_forb_model_full, "Year", var="Drought")), fxn_type="Forb", site="TB")
  

  ## AnnualGrass
  ##
  tb_annualgrass_model_full <- lme(log(ANPP_gm2+.1) ~ as.factor(Year)*Drought+as.factor(Year)*Grazing+Drought*Grazing
                                   , data=filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="AnnualGrass_gm2")
                                   , random = ~1 |Block/Paddock/Plot
                                   , correlation=corAR1(form = ~1 |Block/Paddock/Plot)
                                   , control=lmeControl(returnObject=TRUE)
                                   , na.action = na.omit)
  tb_annualgrass_anova <- anova.lme(tb_annualgrass_model_full, type="marginal")
  plot(tb_annualgrass_model_full, type=c("p","smooth"), col.line=1)
  qqnorm(tb_annualgrass_model_full, abline = c(0,1)) ## qqplot
  hist(log(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="AnnualGrass_gm2")$ANPP_gm2+.1))
  
  emtrends(tb_annualgrass_model_full, "Year", var="Drought")
  test(emtrends(tb_annualgrass_model_full, "Year", var="Drought"))
  # emmeans(tb_anpp_model_full, "Grazing", by="Year")
  # pairs(emmeans(tb_anpp_model_full, "Grazing", by="Year")) ## Nothing comes out in the multiple comparison
  #  ggplot(filter(anpp_fxn_grps, Year %in% 2019:2023 & Site=="TB" & fxn_type=="AnnualGrass_gm2"), aes(x=Drought, y=ANPP_gm2, col=Grazing)) +geom_point()+geom_smooth(method="lm",se=F)+facet_grid(Site~Year)
  
  # Save to writable tables
  tb_annualgrass_anova_df <- data.frame(effect=row.names(tb_annualgrass_anova), tb_annualgrass_anova, fxn_type="AnnualGrass", site="TB")
  tb_annualgrass_emtrends <- data.frame(test(emtrends(tb_annualgrass_model_full, "Year", var="Drought")), fxn_type="AnnualGrass", site="TB")
  
  # 2021 (Only significant drought effect in this year from emtrends output)
  tb_annualgrass_2021_lme <- lme(log(ANPP_gm2+0.1) ~ Drought
                                 , data=filter(anpp_fxn_grps, Year==2021 & Site=="TB" & fxn_type=="AnnualGrass_gm2")
                                 , random = ~1 |Block/Paddock
                                 , na.action = na.omit)
  anova.lme(tb_annualgrass_2021_lme, type="marginal")
  performance::r2(tb_annualgrass_2021_lme) # CAN'T HAVE PLOT IN RANDOM STATEMENT - Marginal R2 considers only the variance of the fixed effects, which is what I want
  
  
  ### Write tables of model output
  anova_raw_fxn_df <- fk_c3p_anova_df %>% 
    bind_rows(fk_c4p_anova_df,
              fk_forb_anova_df,
              fk_annualgrass_anova_df,
              tb_c3p_anova_df,
              tb_c4p_anova_df,
              tb_forb_anova_df,
              tb_annualgrass_anova_df)
  emtrends_raw_fxn_df <- fk_c3p_emtrends %>% 
    bind_rows(fk_c4p_emtrends,
              fk_forb_emtrends,
              fk_annualgrass_emtrends,
              tb_c3p_emtrends,
         #     tb_c4p_emtrends, ##no significant interaction in main model
        #      tb_forb_emtrends,##no significant interaction in main model
              tb_annualgrass_emtrends)
  
  write.csv(anova_raw_fxn_df, file=paste0(table_dir,"fxn groups anpp lme ANCOVA output_both sites_2019-2023",Sys.Date(),".csv"), row.names=F)
  write.csv(emtrends_raw_fxn_df, file=paste0(table_dir,"fxn groups anpp emtrends_both sites_2019-2023",Sys.Date(),".csv"), row.names=F)
  
  
  
  
  
  
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






