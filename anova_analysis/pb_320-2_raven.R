#this R script analyzes how different parameters affect results of BUSCO assessment for genome-assemblies of pb_320-2(Mysco)

library(car)
library(ggplot2)
library(lme4)
library(lmerTest)
library(lmtest)

setwd('~/LUfungiProject/anova_analysis/')
raven=read.csv("pb_320-2_raven.csv")
str(raven) #check 
raven$num_polish<-factor(raven$num_polish)
raven$parameter<-factor(raven$parameter)
complete<-raven$Single+raven$Duplicated
raven<-cbind(raven,complete)
str(raven)

#Construct a model about the number of missing BUSCOs~the 2 factors
fit.Missing<-lm(Missing~num_polish+parameter,data=raven)
#Construct a model about the number of complete BUSCOs~the 2 factors
fit.complete<-lm(raven$complete~raven$num_polish+raven$parameter)
#Construct a model about the number of single-copy complete BUSCOs~the 2 factors
fit.single<-lm(raven$Single~raven$num_polish+raven$parameter)
#Construct a model about the number of duplicated complete BUSCOs~the 2 factors
fit.duplicated<-lm(raven$Duplicated~raven$num_polish+raven$parameter)
#Construct a model about the value of N50~the 2 factors
fit.N50<-lm(raven$N50_Kbp~raven$num_polish+raven$parameter)
#Construct a model about the size of the largest contig~the 2 factors
fit.ctg<-lm(raven$Largest_ctg_Kbp~raven$num_polish+raven$parameter)
#Check assumptions
#Independent measurements
dwtest(fit.Missing)
dwtest(fit.complete)
dwtest(fit.single)
dwtest(fit.duplicated)
dwtest(fit.N50)
dwtest(fit.ctg)
par(mfrow=c(2,2))
plot(fit.Missing)
plot(fit.complete)
plot(fit.single)
plot(fit.duplicated)
plot(fit.N50)
plot(fit.ctg)
par(mfrow=c(1,1))
shapiro.test(fit.Missing$residuals)
shapiro.test(fit.complete$residuals)
shapiro.test(fit.single$residuals)
shapiro.test(fit.duplicated$residuals)
x <- fit.Missing$residuals
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
leveneTest(fit.Missing)
x <- fit.complete$residuals
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
leveneTest(fit.complete)
x <- fit.single$residuals
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
leveneTest(fit.single)
x <- fit.duplicated$residuals
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
leveneTest(fit.duplicated)
x <- fit.N50$residuals
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
leveneTest(fit.N50)
x <- fit.ctg$residuals
y <- pnorm(summary(x), mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE))
ks.test(x, y)
leveneTest(fit.ctg)
#Regarding BUSCO stats, despite leveneTest doesn't work and some Residuals vs Fitted plots doesn't look well here, other assumptions are all fulfilled
#Regarding N50 and largest contig size, each of them doesn't have its ksTest passed, so their residuals don't satisfy the assumption of normal distribution  

#check the assessment results
anova(fit.Missing)
summary(fit.Missing)
anova(fit.complete)
summary(fit.complete)
anova(fit.single)
summary(fit.single)
anova(fit.duplicated)
summary(fit.duplicated)
#these two models are useless though
anova(fit.N50) 
summary(fit.N50)
anova(fit.ctg)
summary(fit.ctg)