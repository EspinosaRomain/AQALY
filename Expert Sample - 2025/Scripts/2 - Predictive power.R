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
library(parallel)
library(pbmcapply)

#Set Working Directory
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#FUNCTION TO LOAD ELEMENTS#
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#Import Data
longdfExpert=loadRData("TreatedData/ExpertSample2025.RData")
retainedDesign=loadRData("Design/Design.RData")
longdfUKRep=loadRData("ImportedData/UKRepMay2024.RData")
MLScaled_UKRepMay2024=loadRData("ImportedData/EstimatesMLScaledUKRepMay2024_AQALY.RData")

#Merge datasets
tmpData=longdfExpert
tmpData$id=tmpData$id+max(longdfUKRep$id)
tmpData$obsId=tmpData$obsId+max(longdfUKRep$obsId)
longdfMerged=rbind(longdfUKRep,tmpData)

#import functions
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsUsedToCreateDesign.R"))
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsForMixedLogitEstimation.R"))
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsForAQALYscores.R"))
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsExportResultsTable.R"))

#Predicted power in the UK representative sample
resSum=summary(MLScaled_UKRepMay2024)
longdfUKRep$predictedProba=resSum$probabilities[,2]
cor(longdfUKRep$predictedProba,longdfUKRep$decision)

#Predicted power in the Expert representative sample
longdfExpert$MINUSyear=-longdfExpert$year
longdfExpertWithPredictions=longdfExpert
longdfExpertWithPredictions$predictedProba=predict(MLScaled_UKRepMay2024, newdata = longdfExpert, obsID = "obsId")[,2]
cor(longdfExpertWithPredictions$predictedProba,longdfExpertWithPredictions$decision)

#Mixed logit estimation: Full sample
listEstimates=list()
maxI=5
vecTotalWeight=rep(NA,maxI)
for(i in 1:maxI){
  longdfMerged$weightsSample=ifelse(longdfMerged$treatment=="Representative",1,i)
  set.seed(123)
  listEstimates[[i]]=mixedLogitEstimationScaledAllRandomButTimeWithWeights(longdfMerged,"weightsSample")
  vecTotalWeight[i]=sum(longdfMerged[longdfMerged$treatment=="ExpertSample",]$weightsSample)/sum(longdfMerged$weightsSample)
}
summary(listEstimates[[1]])$coefTable[2:17,]
tmp=listEstimates[[1]]

#Total weights in the sample:
vecTotalWeight

#Estimate specific welfare scores
AQALYExamplesWithWeights=matrix(data=NA,nrow=5,ncol=(maxI+1))
AQALYExamplesWithWeights[1,1]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"1111"))[1]
AQALYExamplesWithWeights[2,1]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"3333"))[1]
AQALYExamplesWithWeights[3,1]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"5555"))[1]
AQALYExamplesWithWeights[4,1]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"7777"))[1]
AQALYExamplesWithWeights[5,1]=unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,"9999"))[1]
for(i in 1:maxI){
  AQALYExamplesWithWeights[1,(i+1)]=unlist(AQALYscoreFullGradeScaleScaledML(listEstimates[[i]],"1111"))[1]
  AQALYExamplesWithWeights[2,(i+1)]=unlist(AQALYscoreFullGradeScaleScaledML(listEstimates[[i]],"3333"))[1]
  AQALYExamplesWithWeights[3,(i+1)]=unlist(AQALYscoreFullGradeScaleScaledML(listEstimates[[i]],"5555"))[1]
  AQALYExamplesWithWeights[4,(i+1)]=unlist(AQALYscoreFullGradeScaleScaledML(listEstimates[[i]],"7777"))[1]
  AQALYExamplesWithWeights[5,(i+1)]=unlist(AQALYscoreFullGradeScaleScaledML(listEstimates[[i]],"9999"))[1]
}
AQALYExamplesWithWeights
colnames(AQALYExamplesWithWeights)=c("Rep. Sample Only","Weight: 100%","Weight: 200%","Weight: 300%","Weight: 400%","Weight: 500%")
AQALYExamplesWithWeights
write.csv(AQALYExamplesWithWeights, file="Output/AQALYExamplesWithWeights.csv")


#Look at the distribution
allStates=expand.grid(dim1=1:9, dim2=1:9, dim3=1:9, dim4=1:9)
allStates

mydata=unlist(pbmclapply(1:dim(allStates)[1], 
                       function(i){c(unlist(AQALYscoreFullGradeScaleScaledML(MLScaled_UKRepMay2024,paste0(allStates[i,1],allStates[i,2],allStates[i,3],allStates[i,4]))),
                                     unlist(AQALYscoreFullGradeScaleScaledML(listEstimates[[5]],paste0(allStates[i,1],allStates[i,2],allStates[i,3],allStates[i,4]))))
                       }, mc.cores = detectCores()-1))

allStates=cbind(allStates,
                matrix(nrow=dim(allStates)[1],
                       data=mydata,
                       ,byrow=TRUE)
)

colnames(allStates)=c("dim1","dim2","dim3","dim4","AQALY","SE_AQALY","AQALY_500","SE_AQALY_500")
allStates=allStates[order(-allStates$AQALY),]
allScores=as.data.frame(allStates)
allScores$id=seq.int(nrow(allScores))
allScores=allScores[order(-allScores$AQALY_500),]
allScores$id_500=seq.int(nrow(allScores))

cor(allScores$AQALY,allScores$AQALY_500)

graphDistrib=ggplot(data=allScores)+
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=0.5)+
  geom_line(aes(x=id, y = AQALY), color = "#5F98D3", size=1, alpha=1)+
  geom_line(aes(x=id_500, y = AQALY_500), color = "#E77D73", size=1, alpha=0.7)+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=14,face="bold"))+
  xlab("Welfare states (from best to worst)")+
  ylab("AQALY estimate")+
  annotate("text", x = 6100, y = -0.3, label = "Representative \n sample", color="#5F98D3", fontface =2)+
  annotate("text", x = 5000, y = -0.6, label = "Weighted merged sample \n (Expert's relative weight: 5)", color = "#E77D73", fontface =2)
graphDistrib
#3281

saveRDS(allScores,"TreatedData/allScores_AQALY500.RDS")

ggsave(
   "Output/DistribAQALYScoresExperts.pdf",
   graphDistrib,
   width = 10,
   height = 8,
   dpi = 1200,
   bg='white'
)


