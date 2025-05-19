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

### Read in cleaning script (check line numbers haven't changed in initial script)
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






