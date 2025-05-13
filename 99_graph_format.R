################################################################################
##  graph_format.R: code for consistent graph format across all GxD publications
##
##  Author: Kimberly Komatsu
##  Date created: November 19, 2020
################################################################################

###setting the graph look
#theme set -- text size, background, and grid lines
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20))


#designating colors for drought and grazing treatments
droughtColor <- c('#feedde', '#fdbe85', '#fd8d3c', '#e6550d', '#a63603') #from 0 to 99 ##change this to be blue at 0 to red at 99
grazingColor <- c('#ABDEFF', '#469BEC', '#6D882B') #from HHMMM to MMMMM to MLLMM
droughtSymbol <- c('21','22','23','24','25') # Allows for fill
droughtSymbol <- c(21,22,23,24,25) # Allows for fill
grazingSymbol <- c('3','4','8') # Does not allow for fill

##come up with symbols for the treatments too (symbols for grazing, colors for drought?)


###want to know what it will look like?

# library(tidyverse)

#test data
# drought <- c(0,0,0,25,25,25,50,50,50,75,75,75,99,99,99)
# grazing<-c('HHMMM','MMMMM','MLLMM','HHMMM','MMMMM','MLLMM','HHMMM','MMMMM','MLLMM','HHMMM','MMMMM','MLLMM','HHMMM','MMMMM','MLLMM')
# drought_value <- c(0,0,0,25,25,25,50,50,50,75,75,75,99,99,99)
# grazing_value <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
# test_variance <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
# test_value <- c(runif(length(drought_value), min=1, max=5))
# data <- data.frame(drought, grazing, drought_value, grazing_value, test_value, test_variance)
# 
# #sample plot - color by grazing treatment
# ggplot(data=data, aes(x=as.factor(drought), y=test_value, fill=grazing)) +
#   geom_bar(stat='identity', position=position_dodge(), color='black') +
#   geom_errorbar(aes(ymin=test_value-test_variance, ymax=test_value+test_variance), position=position_dodge(0.9), width=0.2) +
#   xlab('Drought Treatment') + ylab('Rainfall Reduction (%)') + 
#   scale_fill_manual(values=grazingColor, name='Grazing\nTreatment')
# 
# #sample plot - color by drought treatment
# ggplot(data=data, aes(x=grazing, y=test_value, fill=as.factor(drought))) +
#   geom_bar(stat='identity', position=position_dodge(), color='black') +
#   geom_errorbar(aes(ymin=test_value-(test_variance/5), ymax=test_value+(test_variance/5)), position=position_dodge(0.9), width=0.2) +
#   xlab('Grazing Treatment') + ylab('Grazing Treatment Number') + 
#   scale_fill_manual(values=droughtColor, name='Drought\nTreatment')
