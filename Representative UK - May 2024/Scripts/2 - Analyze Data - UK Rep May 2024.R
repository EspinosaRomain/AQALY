#Import libraries
library(dplyr)
library(MASS)
library(survival)
library(ggplot2)
library(gridExtra)
library(doParallel) 
library(doSNOW)
library("plotrix")
library(tidyr)
library(survival)
library(splitstackshape)
library(msm)
library(plyr)
library(ggridges) 
library(hrbrthemes) 
library(logitr)
library(ggtext)
library(car)

#Set Working Directory
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#FUNCTION TO LOAD ELEMENTS#
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#Import Data
longdf=loadRData("TreatedData/UKRepMay2024.RData")
retainedDesign=loadRData("Design/Design.RData")

#import functions
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsUsedToCreateDesign.R"))
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsForMixedLogitEstimation.R"))
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsForAQALYscores.R"))
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsExportResultsTable.R"))

#Matrix to save results
ResultsMat=matrix(data=NA,nrow=16,ncol=2)

#Mixed logit estimation: Full sample
set.seed(123)
MLScaled_UKRepMay2024=mixedLogitEstimationScaledAllRandomButTime(longdf)
summary(MLScaled_UKRepMay2024)
save(MLScaled_UKRepMay2024, file = paste0("Output/EstimatesMLScaledUKRepMay2024_AQALY.RData"))
resSum=summary(MLScaled_UKRepMay2024)
resSum$coefTable[2:17,]
exportMixedLogEstimates("resultsUKRep",
                        betas_funct=resSum$coefTable[2:17,1],
                        se_funct=resSum$coefTable[2:17,2],
                        caption_funct="Estimated impact of welfare impairments on AQALY score. Mixed logit estimates, UK representative sample.",
                        label_funct="tab:resultsUKRep")

#Test Equality for Health and Behavior
linearHypothesis(MLScaled_UKRepMay2024,c("X3_5 = X4_5"),test="Chisq")

#Dividing the sample into two
longdf$healthAbove=ifelse(longdf$WeightsDomains.Health.>longdf$WeightsDomains.Behav.,1,0)
MLScaled_UKRepMay2024_HealthAbove=mixedLogitEstimationScaledAllRandomButTime(longdf[longdf$healthAbove==1,])
MLScaled_UKRepMay2024_HealthNOTAbove=mixedLogitEstimationScaledAllRandomButTime(longdf[longdf$healthAbove==0,])
MLScaled_UKRepMay2024_HealthAbove
MLScaled_UKRepMay2024_HealthNOTAbove

#Results - Estimates and SEs
ResultsMat[,1:2]=as.matrix(summary(MLScaled_UKRepMay2024)$coefTable[2:17,1:2])
round(ResultsMat[,1:2],4)
write.csv(ResultsMat, file="Output/ML_EstimatedMLScaledCoefsAndSE_UKRepMay2024.csv")

#Estimate specific welfare scores
AQALYExamples=matrix(data=NA,nrow=5,ncol=2)
AQALYExamples[1,]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"1111"))
AQALYExamples[2,]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"3333"))
AQALYExamples[3,]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"5555"))
AQALYExamples[4,]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"7777"))
AQALYExamples[5,]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"9999"))
AQALYExamples
write.csv(AQALYExamples, file="Output/AQALYExamplesMLScaled_UKRepMay2024.csv")

#Look at the distribution
allStates=expand.grid(dim1=1:9, dim2=1:9, dim3=1:9, dim4=1:9)
allStates

allStates=cbind(allStates,
                matrix(nrow=dim(allStates)[1],
                       data=unlist(lapply(1:dim(allStates)[1], 
                                          function(i){unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,paste0(allStates[i,1],allStates[i,2],allStates[i,3],allStates[i,4]))
                                                             )}
                                          ))
                       ,byrow=TRUE)
                )
colnames(allStates)=c("dim1","dim2","dim3","dim4","AQALY","SE_AQALY")
allStates=allStates[order(-allStates$AQALY),]
allScores=as.data.frame(allStates)
allScores$id=seq.int(nrow(allScores))
allScores$lb=allScores$AQALY-1.96*allScores$SE_AQALY
allScores$ub=allScores$AQALY+1.96*allScores$SE_AQALY
exportScores=allScores[,5:6]
exportScores$SE196=exportScores$SE_AQALY*1.96
write.csv(exportScores, file="Output/AQALYDistrib_UKRepMay2024.csv")

dim(allScores[allScores$AQALY<0,])[1]/dim(allScores)[1]

graphDistrib=ggplot(data=allScores, aes(x=id))+
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5)+
  geom_line(aes(y = AQALY), color = "#5F98D3", size=1, alpha=1)+
  geom_line(aes(y = lb), color = "#5F98D3", alpha=0)+
  geom_line(aes(y = ub), color = "#5F98D3", alpha=0)+
  geom_ribbon(aes(ymin=lb,ymax=ub), fill="#5F98D3", alpha=0.5)+theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=14,face="bold"))+
  xlab("Welfare states (from best to worst)")+
  ylab("AQALY estimate")
graphDistrib

ggsave(
  "Output/DistribAQALYScores.png",
  graphDistrib,
  width = 10,
  height = 8,
  dpi = 1200,
  bg='white'
)



#Heterogeneity Plot
set.seed(123)
nbDraws=1000
dfHet=data.frame(getRandomDrawPar(MLScaled_UKRepMay2024, ndraw_funct = nbDraws))[,2:17]

vecVar=c("Nutrition - Mild positive","Nutrition - Neutral","Nutrition - Mild nevative","Nutrition - Severe negative",
         "Environment - Mild positive","Environment - Neutral","Environment - Mild nevative","Environment - Severe negative",
         "Health - Mild positive","Health - Neutral","Health - Mild nevative","Health - Severe negative",
         "Behavior - Mild positive","Behavior - Neutral","Behavior - Mild nevative","Behavior - Severe negative")

dfForGraph=data.frame(matrix(unlist(lapply(1:nbDraws, function(i){c(vecVar[1],dfHet[i,1],1) })), ncol=3, byrow=TRUE))
for(j in 2:16){
  dfForGraph=rbind(dfForGraph,data.frame(matrix(unlist(lapply(1:nbDraws, function(i){c(vecVar[j],dfHet[i,j],j) })), ncol=3, byrow=TRUE)))  
}
head(dfForGraph)
tail(dfForGraph)
colnames(dfForGraph)=c("Var","Score","Order")
dfForGraph$Score=as.numeric(dfForGraph$Score)
dfForGraph$Order=as.integer(dfForGraph$Order)
for(var_funct in unique(dfForGraph$Var)){
  tmp_var=quantile(dfForGraph[dfForGraph$Var==var_funct,]$Score,probs=c(0.01,0.99), na.rm=TRUE)
  print(tmp_var)
  dfForGraph[dfForGraph$Var==var_funct,]$Score=ifelse(dfForGraph[dfForGraph$Var==var_funct,]$Score<tmp_var[1] | dfForGraph[dfForGraph$Var==var_funct,]$Score>tmp_var[2],NA,dfForGraph[dfForGraph$Var==var_funct,]$Score)
}
dfForGraph=dfForGraph[!is.na(dfForGraph$Score),]

#Graph about heterogeneity
graphHet=ggplot(dfForGraph, aes(x = Score, y = reorder(Var, -Order), height=stat(density), fill=stat(x)))+
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
  scale_y_discrete(labels = rev(vecVar))+
  expand_limits(y=c(0,17))+
  expand_limits(x=-1.5)+
  theme_ridges()+
  theme(
    plot.title = element_markdown(),
    legend.text = element_markdown(),
    plot.caption = element_markdown(),
    axis.text.y=element_text(size=11)
  )+
  geom_point(mapping=aes(x=ResultsMat[1,1], y=16), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[2,1], y=15), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[3,1], y=14), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[4,1], y=13), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[5,1], y=12), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[6,1], y=11), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[7,1], y=10), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[8,1], y=9), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[9,1], y=8), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[10,1], y=7), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[11,1], y=6), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[12,1], y=5), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[13,1], y=4), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[14,1], y=3), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[15,1], y=2), color="black", inherit.aes = FALSE, size=2)+
  geom_point(mapping=aes(x=ResultsMat[16,1], y=1), color="black", inherit.aes = FALSE, size=2)+
  annotate("text", x=ResultsMat[16,1], y=0.4, label="Average estimate")+
  theme(legend.position="none")+
  labs(title = "", y = "", x="",
       caption = "<span style='font-size:13pt'>**Impact on AQALY grade**</span><br/>
       <span style='font-size:10pt'>Heterogeneous effects from Scaled Mixed Logit Estimates</span><br/>
       <span style='font-size:10pt'>UK representative sample</span>")
graphHet

ggsave(
  "Output/HeterogeneousEffectScaled_RepresentativeUK.png",
  graphHet,
  width = 8,
  height = 8,
  dpi = 1200,
  bg='white'
)

  
