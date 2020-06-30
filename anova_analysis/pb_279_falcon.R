#this R script analyzes how different stages of falcon-assembly & polishing and different versions of configuration files affect results of BUSCO assessment for genome-assemblies of pb_279(Leuge)

library(car)
library(ggplot2)
library(lme4)
library(lmerTest)
library(lmtest)

setwd('~/LUfungiProject/anova_analysis/')
falcon=read.csv("pb_279_falcon.csv")
str(falcon) #check 
falcon$step<-factor(falcon$step)
falcon$cfg_version<-factor(falcon$cfg_version)
complete<-falcon$Single+falcon$Duplicated
falcon<-cbind(falcon,complete)
str(falcon)

#Construct a model about the number of missing BUSCOs~the 2 factors
fit.Missing<-lm(Missing~step+cfg_version,data=falcon)
#Construct a model about the number of complete BUSCOs~the 2 factors
fit.complete<-lm(falcon$complete~falcon$step+falcon$cfg_version)
#Check assumptions
#Independent measurements
dwtest(fit.Missing)
dwtest(fit.complete)
#Assumption is fulfilled

#check the assessment results
anova(fit.Missing)
anova(fit.complete)
# $step is significant, whereas $cfg_version is insignicant 