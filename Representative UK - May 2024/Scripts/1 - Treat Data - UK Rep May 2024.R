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

#Set Working Directory
setwd(dirname(dirname(rstudioapi::getSourceEditorContext()$path)))

#FUNCTION TO LOAD ELEMENTS#
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#Import Data
df=loadRData("RawData/UKRepMay2024_rawAnonymData.RData")

#Filter with attention checks
table(df$understandGreen)
df=df[df$understandGreen=="Factors that positively affect animal welfare",]

table(df$understandRed)
df=df[df$understandRed=="Factors that negatively affect animal welfare",]

table(df$attentionCheck1)
df=df[df$attentionCheck1=="The welfare of animals",]

table(df$attentionCheck2)
df=df[df$attentionCheck2=="Factors that explain ecosystem services",]

#import functions
source(paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/CommonFunctions/","FunctionsUsedToCreateDesign.R"))

#Rename IDs
df$id=1:dim(df)[1]

#Reconstruct of pairs displayed
df$pair1=df$t1
df$pair2=df$t29
df$pair3=df$t39
df$pair4=df$t49
df$pair5=df$t59
df$pair6=df$t69
df$pair7=df$t79
df$pair8=df$t89
df$pair9=df$t99
df$pair10=df$t109

#Distribution of pairs
displayedDf=data.frame(pair=1:200)
displayedDf$freq=0
for(i in 1:200){
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair1)[as.character(i)]),0,table(df$pair1)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair2)[as.character(i)]),0,table(df$pair2)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair3)[as.character(i)]),0,table(df$pair3)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair4)[as.character(i)]),0,table(df$pair4)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair5)[as.character(i)]),0,table(df$pair5)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair6)[as.character(i)]),0,table(df$pair6)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair7)[as.character(i)]),0,table(df$pair7)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair8)[as.character(i)]),0,table(df$pair8)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair9)[as.character(i)]),0,table(df$pair9)[as.character(i)])
  displayedDf[displayedDf$pair==i,]$freq=displayedDf[displayedDf$pair==i,]$freq+ifelse(is.na(table(df$pair10)[as.character(i)]),0,table(df$pair10)[as.character(i)])
}
displayedDf$relFreq=displayedDf$freq/sum(displayedDf$freq)
densityPairs=ggplot(data=displayedDf,mapping=aes(x=pair,y=relFreq))+
  geom_bar(stat="identity", fill="steelblue")+
  theme_minimal()+labs(x="Pair number",y="Relative Frequency (%)")
densityPairs

#Create variable to append with Pilot 3
df$treatment="Representative"

#Load the retained design
retainedDesign=loadRData("Design/Design.RData")

#Prepare dataset for estimation
listVars=c("Task1Choice","Task2Choice","Task3Choice","Task4Choice","Task5Choice","Task6Choice","Task7Choice","Task8Choice","Task9Choice","Task10Choice")
longdf=df[,c("id","treatment",listVars,"WeightsDomains.Nutr.","WeightsDomains.Env.","WeightsDomains.Health.","WeightsDomains.Behav.")] %>% 
  pivot_longer(cols=listVars,
               names_to="task",
               values_to="choice")
listVars=c("pair1","pair2","pair3","pair4","pair5","pair6","pair7","pair8","pair9","pair10")
longdf=cbind(longdf,
  df[,c(listVars)] %>% 
  pivot_longer(cols=listVars,
               names_to="pair",
               values_to="displayedPair"))
longdf$task=rep(1:10,max(longdf$id)) #I rename the task ID
longdf=cbind(longdf,retainedDesign[longdf$displayedPair,]) #I assign the card value for each situation
longdf=expandRows(longdf, count=2, count.is.col=FALSE) #I duplicate to have a row per card
longdf$card=rep(1:2,dim(longdf)[1]/2) #I create an indicator for the card
longdf$decision=ifelse(longdf$choice=="Situation 1 is better" & longdf$card==1,1,0)
longdf$decision=ifelse(longdf$choice=="Situation 2 is better" & longdf$card==2,1,longdf$decision)
longdf$state=ifelse(longdf$card==1,longdf$ordered_state1,longdf$ordered_state2)
longdf=subset(longdf, select = -c(choice, pair,ordered_state1,ordered_state2))
longdf=cbind(longdf,matrixOfXs=matrix(unlist(lapply(longdf$state,generateXrow)),ncol=17, byrow=TRUE))
colnames(longdf)[(dim(longdf)[2]-16):dim(longdf)[2]]=colnames=c("year", "X1_2","X1_3","X1_4","X1_5","X2_2","X2_3","X2_4","X2_5","X3_2","X3_3","X3_4","X3_5","X4_2","X4_3","X4_4","X4_5")

longdf$obsId=sort(rep(seq(1:(dim(longdf)[1]/2)),2)) #Identifiers for mixed logit estimation

save(longdf,file="TreatedData/UKRepMay2024.RData")

