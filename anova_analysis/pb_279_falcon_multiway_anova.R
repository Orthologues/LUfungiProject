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
str(falcon)