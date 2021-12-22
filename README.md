# GCaMP-Calcium-analysis
Code script for analyzing Calcium signaling in GCaMP cells with LMM

#GCaMP Calcium Analysis JMR

#load the library 
library(tidyverse)
library(plyr)
library(dplyr)
library(stringr)

#load the data 
filenames <- list.files(pattern="*.csv") #read all the .csv files 


#create a dataframe; Remove missing values
df<-purrr::map_df(filenames, read_csv, .id = 'filename') 
df<- na.omit(df) #remove missing values 

#change the names to lowercase
df$roiName<- tolower(df$roiName)
df$Spot<- tolower(df$Spot)
df$Condition<- tolower(df$Condition)
df$drug<- tolower(df$drug)

#Change nominal and date format
df$filename = factor(df$filename)
df$roiName = factor(df$roiName)
df$Animal = factor(df$Animal)
df$drug = factor(df$drug)
df$Spot = factor(df$Spot)
df$filename = factor(df$filename)
df$date = factor(df$date)
view(df)

#Modifying variables and cleaning data frame
#Short df and create a unique ROIname
df_short <- df[ ,c('date', 'drug', 'Condition', 'Spot', 'Animal', 'roiName', 'amplitude', 'halfWidth','peakAUC')]
#Replace after drugs for baseline
df_short$drug <- str_replace_all(df_short$drug, 'after drugs', 'baseline')
#adding unique Spot
df_short$Unique_Spot <-paste(df_short$Animal, df_short$Spot, sep= "_")
#adding unique ROI name
df_short$Unique_ROIname <-paste(df_short$Animal, df_short$Spot, df_short$roiName, sep= "_")
#adding Duration of event
df_short$Duration <-df_short$halfWidth*2
#create a new variable for process and soma
df_short$ROI<- ifelse(grepl(pattern = "p", df_short$roiName, ignore.case = T),"Process","Soma")
unique(df_short$drug)
unique(df_short$ROI)
unique(df_short$Animal)
df_short <- df_short[ ,c('date', 'Unique_Spot', 'Unique_ROIname', 'ROI', 'Condition', 'drug', 'amplitude', 'Duration','peakAUC')]
view(df_short)


#Mean of each trial by ROI 

# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(df_short, c("date", "Unique_Spot", "Unique_ROIname", "ROI", "Condition", "drug"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(Duration,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude))

#Remove negative amplitudes
ROI.means<- ROI.means%>%
  filter(amp_mean>0)

view(ROI.means)

#Converting to nominal factor the modified data frame

ROI.means$Unique_Spot= factor(ROI.means$Unique_Spot)
ROI.means$Unique_ROIname = factor(ROI.means$Unique_ROIname)
ROI.means$ROI = factor(ROI.means$ROI)
ROI.means$Condition = factor(ROI.means$Condition)
ROI.means$drug = factor(ROI.means$drug)


########################################################################################################
#STATISTICS

summary(ROI.means) 

#Amplitude

# explore the amplitude
library(plyr)
ddply(ROI.means, ~ drug * Condition, function(data) summary(data$amp_mean))
ddply(ROI.means, ~ drug * Condition, summarise, amp.mean=mean(amp_mean), amp.sd=sd(amp_mean))


# histograms for two factors
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "no stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "no stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "no stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "no stim",]$amp_mean)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "stim",]$amp_mean)


boxplot(amp_mean ~ drug * Condition, data=ROI.means, xlab="Drug.Condition", ylab="amplitude") # boxplots
with(ROI.means, interaction.plot(drug, Condition, amp_mean, ylim=c(0, max(ROI.means$amp_mean)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, amp_mean, ylim=c(0, max(ROI.means$amp_mean)))) # interaction?

#Testing for normality in the residuals
library(plyr)
m = aov(amp_mean ~ drug * Condition, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(amp_mean ~ drug * Condition, data=ROI.means, center=mean) # Levene's test
leveneTest(amp_mean ~ drug * Condition, data=ROI.means, center=median) # Brown-Forsythe test

#DO LMM Test before modification of the distribution of the data

# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"


# LMM with Spot as random effect for DRUG factor
m = lmer(amp_mean ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with Spot as random effect for Condition factor
m = lmer(amp_mean ~ (Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ Condition)), test=adjusted(type="holm"))


# LMM with Spot as random effect drug*Condition
m = lmer(amp_mean ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * Condition)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, Condition, amp_mean, ylim=c(0, max(ROI.means$amp_mean)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, amp_mean, ylim=c(0, max(ROI.means$amp_mean)))) # interaction?



# see if data seems Poisson-distributed (DATA DO NOT FIT IN THIS DISTRIBUTION TYPE)
#Scaling amp_mean  * 10
#ROI.means$amp_scale = abs(ROI.means$amp_mean)*10 
#View(ROI.means) # verify
#summary(ROI.means)
#library(fitdistrplus)
#fit = fitdist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "no stim",]$amp_mean, "pois", discrete=TRUE)
#gofstat(fit) # goodness-of-fit test

#MODIFICATION OF THE DATA WITH LOGNORMAL DISTRIBUTION
#creating log of amp_mean
ROI.means$logamp_mean = log(ROI.means$amp_mean) # log transform
View(ROI.means) # verify

hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "no stim",]$logamp_mean)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "stim",]$logamp_mean)

boxplot(logamp_mean ~ drug * Condition, data=ROI.means, xlab="Drug.Condition", ylab="amplitude") # boxplots
with(ROI.means, interaction.plot(drug, Condition, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction

#Testing for normality in the residuals
library(plyr)
m = aov(logamp_mean ~ drug * Condition, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(logamp_mean ~ drug * Condition, data=ROI.means, center=mean) # Levene's test
leveneTest(logamp_mean ~ drug * Condition, data=ROI.means, center=median) # Brown-Forsythe test

# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"


# LMM with Spot as random effect for DRUG factor
m = lmer(logamp_mean ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with Spot as random effect for Condition factor
m = lmer(logamp_mean ~ (Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ Condition)), test=adjusted(type="holm"))


# LMM with Spot as random effect drug*Condition
m = lmer(logamp_mean ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * Condition)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, Condition, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, amp_mean, ylim=c(0, max(ROI.means$logamp_mean)))) # interaction?

#####################################################################################################
#Stats nEvents

ROI.means$Unique_Spot= factor(ROI.means$Unique_Spot)
ROI.means$Unique_ROIname = factor(ROI.means$Unique_ROIname)
ROI.means$ROI = factor(ROI.means$ROI)
ROI.means$Condition = factor(ROI.means$Condition)
ROI.means$drug = factor(ROI.means$drug)

summary(ROI.means) 

# explore the nEvents
library(plyr)
ddply(ROI.means, ~ drug * Condition, function(data) summary(data$nEvents))
ddply(ROI.means, ~ drug * Condition, summarise, nEvents=mean(nEvents), nEvents=sd(nEvents))


# histograms for two factors
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "no stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "no stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "no stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "no stim",]$nEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "stim",]$nEvents)

boxplot(nEvents ~ drug * Condition, data=ROI.means, xlab="Drug.Condition", ylab="nEvents") # boxplots
with(ROI.means, interaction.plot(drug, Condition, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?

#Testing for normality in the residuals
library(plyr)
m = aov(nEvents ~ drug * Condition, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(nEvents ~ drug * Condition, data=ROI.means, center=mean) # Levene's test
leveneTest(nEvents ~ drug * Condition, data=ROI.means, center=median) # Brown-Forsythe test

# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"

#DOING LMM BEFORE TRANSFORMING THE DATA
# LMM with Spot as random effect for DRUG factor
m = lmer(nEvents ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with Spot as random effect for Condition factor
m = lmer(nEvents ~ (Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ Condition)), test=adjusted(type="holm"))

# LMM with Spot as random effect drug*Condition
m = lmer(nEvents ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * Condition)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, Condition, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction


# see if data seems Poisson-distributed (DATA FIT BUT THERE IS NOT MUCH OF A DIFFERENCE)
#library(fitdistrplus)
#fit = fitdist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "no stim",]$nEvents, "pois", discrete=TRUE)
#gofstat(fit) # goodness-of-fit test
#fit = fitdist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "stim",]$nEvents, "pois", discrete=TRUE)


# main GLMM test on nEvents with drug 
#m = glmer(nEvents ~ (drug) + (1|Unique_ROIname), data=ROI.means, family=poisson, nAGQ=1)
#Anova(m, type=3)
# perform post hoc pairwise comparisons
#library(multcomp) # for glht
#library(emmeans) # for emm
#summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# main GLMM test on nEvents with drug and Condition
#m = glmer(nEvents ~ (drug * Condition) + (1|Unique_ROIname), data=ROI.means, family=poisson, nAGQ=1)
#Anova(m, type=3)
# perform post hoc pairwise comparisons
#library(multcomp) # for glht
#library(emmeans) # for emm
#summary(glht(m, emm(pairwise ~ drug * Condition)), test=adjusted(type="holm"))
#with(ROI.means, interaction.plot(drug, Condition, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?
#with(ROI.means, interaction.plot(Condition, drug, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?


#Transforming to log
ROI.means$lognEvents = log(ROI.means$nEvents) # log transform
View(ROI.means) # verify

hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "no stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$Condition == "stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "no stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$Condition == "stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "no stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$Condition == "stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "no stim",]$lognEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$Condition == "stim",]$lognEvents)

boxplot(lognEvents ~ drug * Condition, data=ROI.means, xlab="Drug.Condition", ylab="lognEvents") # boxplots
with(ROI.means, interaction.plot(drug, Condition, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction

#Testing for normality in the residuals
library(plyr)
m = aov(lognEvents ~ drug * Condition, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(lognEvents ~ drug * Condition, data=ROI.means, center=mean) # Levene's test
leveneTest(lognEvents ~ drug * Condition, data=ROI.means, center=median) # Brown-Forsythe test

# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"


# LMM with Spot as random effect for DRUG factor
m = lmer(lognEvents ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with Spot as random effect for Condition factor
m = lmer(lognEvents ~ (Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ Condition)), test=adjusted(type="holm"))

# LMM with Spot as random effect drug*Condition
m = lmer(lognEvents ~ (drug * Condition)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * Condition)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, Condition, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction?
with(ROI.means, interaction.plot(Condition, drug, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction
##########################################################################################################

#Statistics  of Drug and ROI

#Stats nEvents

ROI.means$Unique_Spot= factor(ROI.means$Unique_Spot)
ROI.means$Unique_ROIname = factor(ROI.means$Unique_ROIname)
ROI.means$ROI = factor(ROI.means$ROI)
ROI.means$Condition = factor(ROI.means$Condition)
ROI.means$drug = factor(ROI.means$drug)

summary(ROI.means) 

# explore the nEvents
library(plyr)
ddply(ROI.means, ~ drug * ROI, function(data) summary(data$nEvents))
ddply(ROI.means, ~ drug * ROI, summarise, nEvents=mean(nEvents), nEvents=sd(nEvents))


# histograms for two factors
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$ROI == "Process",]$nEvents)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$ROI == "Soma",]$nEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$ROI == "Process",]$nEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$ROI == "Soma",]$nEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$ROI == "Process",]$nEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$ROI == "Soma",]$nEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$ROI == "Process",]$nEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$ROI == "Soma",]$nEvents)


boxplot(nEvents ~ drug * ROI, data=ROI.means, xlab="Drug.ROI", ylab="nEvents") # boxplots
with(ROI.means, interaction.plot(drug, ROI, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?
with(ROI.means, interaction.plot(ROI, drug, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?

#Testing for normality in the residuals
library(plyr)
m = aov(nEvents ~ drug * ROI, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(nEvents ~ drug * ROI, data=ROI.means, center=mean) # Levene's test
leveneTest(nEvents ~ drug * ROI, data=ROI.means, center=median) # Brown-Forsythe test

# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"

#DOING LMM before transforming the data
# LMM with Spot as random effect for DRUG factor
m = lmer(nEvents ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with Spot as random effect for ROI factor
m = lmer(nEvents ~ (ROI)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ ROI)), test=adjusted(type="holm"))

# LMM with Spot as random effect drug*ROI
m = lmer(nEvents ~ (drug * ROI)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * ROI)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, ROI, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?
with(ROI.means, interaction.plot(ROI, drug, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction



# see if data seems Poisson-distributed (THEY FIT BUT THERE IS NOT MUCH OF A DIFFERENCE)
#library(fitdistrplus)
#fit = fitdist(ROI.means[ROI.means$drug == "baseline" & ROI.means$ROI == "Process",]$nEvents, "pois", discrete=TRUE)
#gofstat(fit) # goodness-of-fit test
#fit = fitdist(ROI.means[ROI.means$drug == "baseline" & ROI.means$ROI == "Soma",]$nEvents, "pois", discrete=TRUE)
#gofstat(fit) # goodness-of-fit test


# main GLMM test on nEvents with drug 
#m = glmer(nEvents ~ (drug) + (1|Unique_ROIname), data=ROI.means, family=poisson, nAGQ=1)
#Anova(m, type=3)
# perform post hoc pairwise comparisons
#library(multcomp) # for glht
#library(emmeans) # for emm
#summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# main GLMM test on nEvents with ROI
#m = glmer(nEvents ~ (ROI) + (1|Unique_ROIname), data=ROI.means, family=poisson, nAGQ=1)
#Anova(m, type=3)
# perform post hoc pairwise comparisons
#library(multcomp) # for glht
#library(emmeans) # for emm
#summary(glht(m, emm(pairwise ~ ROI)), test=adjusted(type="holm"))

# main GLMM test on nEvents with drug and ROI
#m = glmer(nEvents ~ (drug * ROI) + (1|Unique_ROIname), data=ROI.means, family=poisson, nAGQ=1)
#Anova(m, type=3)
# perform post hoc pairwise comparisons
#library(multcomp) # for glht
#library(emmeans) # for emm
#summary(glht(m, emm(pairwise ~ drug * ROI)), test=adjusted(type="holm"))

#with(ROI.means, interaction.plot(drug, ROI, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?
#with(ROI.means, interaction.plot(ROI, drug, nEvents, ylim=c(0, max(ROI.means$nEvents)))) # interaction?


#LMM with lognEvents

hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$ROI == "Process",]$lognEvents)
hist(ROI.means[ROI.means$drug == "baseline" & ROI.means$ROI == "Soma",]$lognEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$ROI == "Process",]$lognEvents)
hist(ROI.means[ROI.means$drug == "peg" & ROI.means$ROI == "Soma",]$lognEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$ROI == "Process",]$lognEvents)
hist(ROI.means[ROI.means$drug == "nimodipine" & ROI.means$ROI == "Soma",]$lognEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$ROI == "Process",]$lognEvents)
hist(ROI.means[ROI.means$drug == "pyr3" & ROI.means$ROI == "Soma",]$lognEvents)

boxplot(lognEvents ~ drug * ROI, data=ROI.means, xlab="Drug.ROI", ylab="lognEvents") # boxplots
with(ROI.means, interaction.plot(drug, ROI, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction?
with(ROI.means, interaction.plot(ROI, drug, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction

#Testing for normality in the residuals
library(plyr)
m = aov(lognEvents ~ drug * ROI, data=ROI.means) # fit model
shapiro.test(residuals(m)) # test residuals
qqnorm(residuals(m)); qqline(residuals(m)) # plot residuals
# tests for homoscedasticity (homogeneity of variance)
library(car)
leveneTest(lognEvents ~ drug * ROI, data=ROI.means, center=mean) # Levene's test
leveneTest(lognEvents ~ drug * ROI, data=ROI.means, center=median) # Brown-Forsythe test

# libraries for LMMs 
library(lme4) # for lmer
library(lmerTest)
library(car) # for Anova

# set sum-to-zero contrasts for the Anova calls
contrasts(ROI.means$date) <- "contr.sum"
contrasts(ROI.means$Unique_Spot) <- "contr.sum"
contrasts(ROI.means$Unique_ROIname) <- "contr.sum"
contrasts(ROI.means$ROI) <- "contr.sum"
contrasts(ROI.means$Condition) <- "contr.sum"
contrasts(ROI.means$drug) <- "contr.sum"


# LMM with ROI as random effect for DRUG factor
m = lmer(lognEvents ~ (drug )  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug)), test=adjusted(type="holm"))

# LMM with ROI as random effect for ROI factor
m = lmer(lognEvents ~ (ROI)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ ROI)), test=adjusted(type="holm"))

# LMM with ROI as random effect drug*Condition
m = lmer(lognEvents ~ (drug * ROI)  + (1|Unique_ROIname), data=ROI.means)
Anova(m, type=3, test.statistic="F")
# perform post hoc pairwise comparisons
library(multcomp) # for glht
library(emmeans) # for emm
summary(glht(m, emm(pairwise ~ drug * ROI)), test=adjusted(type="holm"))

with(ROI.means, interaction.plot(drug, ROI, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction?
with(ROI.means, interaction.plot(ROI, drug, lognEvents, ylim=c(0, max(ROI.means$lognEvents)))) # interaction

