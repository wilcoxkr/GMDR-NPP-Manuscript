### Function to conduct type 3 anova's using lme
### Currently only works for models that have 2 or 3 independant variables, but could be coded to deal with 4 or 5 variables relatively easily
### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created: Jan. 25, 2021, Last Modified Jan 26, 2021



### example showing format for function inputs:

# anova_t3(IndVars=c('Time_Period','Grazing','Drought'),
#          DepVar='Soil_Moisture',
#          RndForm='~1 |Block/Paddock/Plot',
#          Data=subset(smoist_data_all, Site=="TB" & Year==2020)
# )  

anova_t3 <- function(IndVars=IndVars, DepVar=DepVar, RndForm=RndForm, Data=Data){
  anova_out <- {}
  IndVarMatrix <- matrix(nrow=length(IndVars),ncol=length(IndVars))
  IndVars2x <- c(IndVars,IndVars)
  
  for(REORDER in 1:length(IndVars)){
    IndVarMatrix[REORDER,] <- IndVars2x[REORDER:(length(IndVars)+(REORDER-1))]
  }
  rm(IndVars2x)
  
  for(RUN in 1:length(IndVars)){
    model_formula_temp <- paste0(DepVar,"~", paste0(IndVarMatrix[RUN,], collapse="*"))
    model_temp <- lme(as.formula(model_formula_temp)
                      , data=Data
                      , random = as.formula(RndForm)
                      , correlation=corCompSymm(form = as.formula(RndForm))
                      , control=lmeControl(returnObject=TRUE)
                      , na.action = na.omit)
    anova_out_temp <- anova(model_temp)
    
    if(length(IndVars)==3){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7:8,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[7,]
        )
      }
    } # End if vars==3 statement
    
    if(length(IndVars)==2){ ## Pulls model output for variable that is last
      if(RUN==1){ ### This currently works for 
        anova_partial_temp <- rbind(anova_out_temp[(length(IndVars)+1),],
                                    anova_out_temp[4,]
        )
      }
      if(RUN %in% 2:length(IndVars)){
        anova_partial_temp <- anova_out_temp[(length(IndVars)+1),]
        
      }
    } # End if vars == 2 statement
    anova_out <- rbind(anova_out, anova_partial_temp)
  } # End reorder loop
  return(anova_out)
} # End function
