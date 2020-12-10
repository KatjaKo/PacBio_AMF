---
  title: "script_bayesian_stats"
author: "Katja Kozjek"
---

# Load libraries

library(phyloseq); packageVersion("phyloseq") 
library(vegan); packageVersion("vegan") 
library(ggplot2); packageVersion("ggplot2")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(ggpubr); packageVersion("ggpubr")
library(devtools); packageVersion("devtools")
library(ape); packageVersion("ape")
library(brms); packageVersion("brms") 
library(shinystan); packageVersion("shinystan") 
library(picante); packageVersion("picante")
library(ggvegan); packageVersion("ggvegan")
library(wesanderson); packageVersion("wesanderson")
  
# ALPHA diversity indices

#we calculate Observed and Shannon for each sample

alpha_div <- data.frame(estimate_richness(physeq_unite_Filtered_ASV, split=TRUE, measures = c("Observed", "Shannon")), sample_data(physeq_unite_Filtered_ASV))
plot_richness(physeq_unite_Filtered_ASV, x="treatment", color="farming_system", measures=c("Observed","Shannon")) + stat_boxplot(geom = "errorbar") + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Treatment") + ggtitle("AMF alpha diversity") + ylab("Alpha metrics")

write.csv(alpha_div, "results/alpha_div.csv")

## Prepare data for BRMS model

data_brms <- alpha_div
head(data_brms)
str(data_brms)

data_brms$farming_system <- as.factor(data_brms$farming_system)
data_brms$treatment <- as.factor(data_brms$treatment)
data_brms$time <- as.factor(data_brms$time)

# Negbinominal distribution for Observed

observed_brms <-brm(formula = (Observed)~time*treatment*farming_system+(1|block/plot),
                    data = data_brms,
                    family="negbinomial",
                    seed = 1041,
                    warmup = 2000,#The number of warmup iterations should 
                    #not be larger than iter and the default is iter/2.
                    iter = 6000, 
                    chains = 5,
                    control = list(adapt_delta = 0.999, stepsize=0.09, max_treedepth = 15))

pp_check(observed_brms)

#save an object to a file
saveRDS(observed_brms, "data/rds/observed_brms.rds")

#reload it so you don't have to run it again
observed_brms <- readRDS("data/rds/observed_brms.rds")

#check model fit 
yfit_observed <- data_brms$Observed
yrep_observed <- posterior_predict(observed_brms, draws = length(yfit_observed))

launch_shinystan(observed_brms)

#explore results
summary(observed_brms)
bayes_R2(observed_brms)

#marginal effects
marginal_effects(observed_brms)

## plot BRMS observed 
newdat <- expand.grid(
  time = factor(c("4w", "13w"),levels=levels(data_brms$time)), 
  farming_system = factor(c("BIODYN", "CONMIN"), levels=levels(data_brms$farming_system)),
  treatment = factor(c("C","R", "RC"),levels=levels(data_brms$treatment)))
head(newdat)

#use the function "fitted" from brms to get fitted means of the posterior for treatments
#It summarizes from the posterior and calculates meadian with 95% credible interval
#= we are 95% sure that the real estimate is within this interval

fit_observed = fitted(observed_brms, newdata = newdat, robust=T, re_formula = NA,summary = TRUE, probs = c(0.025, 0.975))  
fit_observed

# no colnames, I add them manually
fitdf_observed_plot = cbind(newdat, fit_observed)
fitdf_observed_plot

#change order of treatment in fitted values
fitdf_observed_plot$treatment<-factor(fitdf_observed_plot$treatment, levels=c("R", "RC", "C"))
data_brms$treatment<-factor(data_brms$treatment, levels=c("R", "RC", "C"))

#prepare and extract results table
fitdf_observed<-fitdf_observed_plot[, c("farming_system","treatment", "time", "Estimate","Q2.5","Q97.5")]

#make the final plot, without raw data
fitted_observed<-ggplot(data = fitdf_observed_plot, 
                        aes(x = farming_system,
                            y =Estimate,
                            color = treatment,
                            shape=treatment)) +
  geom_point(size=6,position=pos.dodge) +
  geom_errorbar(aes(ymin=Q2.5,
                    ymax=Q97.5),
                width=0,
                position=pos.dodge,
                size=1.2, show.legend = F) +
  scale_color_manual(values=wes_palette(n=3, name="Moonrise2"), name  ="") +
  scale_shape_manual(values = c(15, 17, 19),name  ="")+
  xlab("Farming system")+
  facet_grid(~time)+ 
  ylab("Observed ASV")+theme_bw()

print(fitted_observed)

#get location of the means including credible intervals from here:
fitdf_observed_plot<-as.data.frame(fitdf_observed_plot)

#if you want, you can also calculate differences between treatments based on the posterior
##calucate differences with CrI, you need the whole posterior
#the whole posterior distribution for me. I save it as dataframe (fit1)

fit_observed = as.data.frame(fitted(observed_brms, newdata=newdat,
                                    re_formula = NA, summary = FALSE))

#this is the whole posterior, not summarized as before

#add names
colnames<-as.character(interaction(newdat$farming_system, newdat$treatment, newdat$time))

colnames(fit_observed)<-colnames
names(fit_observed)

#calculate the mean differences between systems across the roofs for 4w
BIODYN_brms_obs_T1<-(fit_observed$BIODYN.C.4w+fit_observed$BIODYN.RC.4w+fit_observed$BIODYN.R.4w)/3
mean(BIODYN_brms_obs_T1)
sd(BIODYN_brms_obs_T1)
CONMIN_brms_obs_T1<-(fit_observed$CONMIN.C.4w+fit_observed$CONMIN.RC.4w+fit_observed$CONMIN.R.4w)/3
mean(CONMIN_brms_obs_T1)
sd(CONMIN_brms_obs_T1)

diff_brms_obs_T1<-(BIODYN_brms_obs_T1-CONMIN_brms_obs_T1)

#use the quantile function to get mean of that differences with CrI
round(quantile(diff_brms_obs_T1, probs=c(0.025, 0.5, 0.975)),2)
#this is a result
mean(diff_brms_obs_T1)

#calculate the mean differences between systems across the roofs for 13w
BIODYN_brms_obs_T3<-(fit_observed$BIODYN.C.13w+fit_observed$BIODYN.RC.13w+fit_observed$BIODYN.R.13w)/3
mean(BIODYN_brms_obs_T3)
sd(BIODYN_brms_obs_T3)
CONMIN_brms_obs_T3<-(fit_observed$CONMIN.C.13w+fit_observed$CONMIN.RC.13w+fit_observed$CONMIN.R.13w)/3
mean(CONMIN_brms_obs_T3)
sd(CONMIN_brms_obs_T3)

diff_brms_obs_T3<-(BIODYN_brms_obs_T3-CONMIN_brms_obs_T3)

#use the quantile function to get mean of that differences with CrI
round(quantile(diff_brms_obs_T3, probs=c(0.025, 0.5, 0.975)),2)
#this is a clear result, it includes zero so its likely that there is no difference
mean(diff_brms_obs_T3)

# Comparisons for Observed ASV

#calculate treatment effects
fitTest_observed<-fitted(observed_brms,
                         newdata = newdat,
                         summary = F,
                         re_formula = NA)

fitTest_observed<-as.data.frame(fitTest_observed)
colnames(fitTest_observed) <- interaction(newdat$time, newdat$farming_system, newdat$treatment)

contrastObserved<-cbind(
  BD.RRC.4w<-fitTest_observed$`4w.BIODYN.R`-fitTest_observed$`4w.BIODYN.RC`,
  BD.RC.4w<-fitTest_observed$`4w.BIODYN.R`-fitTest_observed$`4w.BIODYN.C`,
  BD.RCC.4w<-fitTest_observed$`4w.BIODYN.RC`-fitTest_observed$`4w.BIODYN.C`,
  CM.RRC.4w<-fitTest_observed$`4w.CONMIN.R`-fitTest_observed$`4w.CONMIN.RC`,
  CM.RC.4w<-fitTest_observed$`4w.CONMIN.R`-fitTest_observed$`4w.CONMIN.C`,
  CM.RCC.4w<-fitTest_observed$`4w.CONMIN.RC`-fitTest_observed$`4w.CONMIN.C`,
  BD.RRC.13w<-fitTest_observed$`13w.BIODYN.R`-fitTest_observed$`13w.BIODYN.RC`,
  BD.RC.13w<-fitTest_observed$`13w.BIODYN.R`-fitTest_observed$`13w.BIODYN.C`,
  BD.RCC.13w<-fitTest_observed$`13w.BIODYN.RC`-fitTest_observed$`13w.BIODYN.C`,
  CM.RRC.13w<-fitTest_observed$`13w.CONMIN.R`-fitTest_observed$`13w.CONMIN.RC`,
  CM.RC.13w<-fitTest_observed$`13w.CONMIN.R`-fitTest_observed$`13w.CONMIN.C`,
  CM.RCC.13w<-fitTest_observed$`13w.CONMIN.RC`-fitTest_observed$`13w.CONMIN.C`,
  BD.CM.4w<-(fitTest_observed$`4w.BIODYN.C`+fitTest_observed$`4w.BIODYN.R`+fitTest_observed$`4w.BIODYN.RC`)-
    (fitTest_observed$`4w.CONMIN.C`+fitTest_observed$`4w.CONMIN.R`+fitTest_observed$`4w.CONMIN.RC`),
  BD.CM.13w<-(fitTest_observed$`13w.BIODYN.C`+fitTest_observed$`13w.BIODYN.R`+fitTest_observed$`13w.BIODYN.RC`)-
    (fitTest_observed$`13w.CONMIN.C`+fitTest_observed$`13w.CONMIN.R`+fitTest_observed$`13w.CONMIN.RC`))

colnames2<-c("BD.RRC.4w", "BD.RC.4w", "BD.RCC.4w", "CM.RRC.4w", "CM.RC.4w","CM.RCC.4w", "BD.RRC.13w", "BD.RC.13w", "BD.RCC.13w", "CM.RRC.13w" ,"CM.RC.13w", "CM.RCC.13w", "BD.CM.4w", "BD.CM.13w")
colnames(contrastObserved)<-colnames2

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

#get the effect size data with 95% CrI
contrasts_observed <- t(round_df(plyr::colwise(quantile)(as.data.frame(contrastObserved), 
                                                         probs = c(0.5,0.025,0.975)),2))

# BRMS model for Shannon index 

shannon_brms <-brm(formula = (Shannon)~time*treatment*farming_system+(1|block/plot),
                   data = data_brms,
                   family="gaussian",
                   seed = 1041,
                   warmup = 2000,#The number of warmup iterations should 
                   #not be larger than iter and the default is iter/2.
                   iter = 6000, 
                   chains = 5,
                   control = list(adapt_delta = 0.999, max_treedepth = 15))

pp_check(shannon_brms)

#save an object to a file
saveRDS(shannon_brms, "data/rds/shannon_brms.rds")

#reload it so you don't have to run it again
shannon_brms <- readRDS("data/rds/shannon_brms.rds")

#check model fit 
yfit_shannon <- data_brms$Shannon
yrep_shannon <- posterior_predict(shannon_brms, draws = length(yfit_shannon))

launch_shinystan(shannon_brms)

#explore results
summary(shannon_brms)

#marginal effects
marginal_effects(shannon_brms)

## plot brms for Shannon index

fit_shannon = fitted(shannon_brms, newdata = newdat, robust=T, re_formula = NA,summary = TRUE,probs = c(0.025, 0.975))
fit_shannon

# no colnames, I add them manually
fitdf_shannon_plot = cbind(newdat, fit_shannon)
fitdf_shannon_plot

#change order of treatment in fitted values
fitdf_shannon_plot$treatment<-factor(fitdf_shannon_plot$treatment, levels=c("R", "RC", "C"))

#prepare and extract results table
fitdf_shannon<-fitdf_shannon_plot[, c("farming_system","treatment", "time", "Estimate","Q2.5","Q97.5")]

#make the final plot
fitted_shannon<-ggplot(data = fitdf_shannon_plot, 
                       aes(x = farming_system,
                           y =Estimate,
                           color = treatment,
                           shape=treatment)) +
  geom_point(size=6,position=pos.dodge) +
  geom_errorbar(aes(ymin=Q2.5,
                    ymax=Q97.5),
                width=0,
                position=pos.dodge,
                size=1.2, show.legend = F) +
  scale_color_manual(values=wes_palette(n=3, name="Moonrise2"), name  ="") +
  scale_shape_manual(values = c(15, 17, 19),name  ="")+
  xlab("Farming system")+
  facet_grid(~time) + 
  ylab("Shannon diversity index")+theme_bw()

print(fitted_shannon)

#get location of the means including credible intervals from here:
fitdf_shannon_plot<-as.data.frame(fitdf_shannon_plot)

#this is the whole posterior, not summarized as before
fit_shannon = as.data.frame(fitted(shannon_brms, newdata=newdat,
                                   re_formula = NA, summary = FALSE))

#add names
colnames<-as.character(interaction(newdat$farming_system, newdat$treatment, newdat$time))

colnames(fit_shannon)<-colnames
names(fit_shannon)

#calculate the mean differences between systems across the roofs at 4w
BIODYN_brms_shannon_T1<-(fit_shannon$BIODYN.C.4w+fit_shannon$BIODYN.RC.4w+fit_shannon$BIODYN.R.4w)/3
mean(BIODYN_brms_shannon_T1)
sd(BIODYN_brms_shannon_T1)
CONMIN_brms_shannon_T1<-(fit_shannon$CONMIN.C.4w+fit_shannon$CONMIN.RC.4w+fit_shannon$CONMIN.R.4w)/3
mean(CONMIN_brms_shannon_T1)
sd(CONMIN_brms_shannon_T1)

diff_brms_shannon_T1 <- (BIODYN_brms_shannon_T1-CONMIN_brms_shannon_T1)

#use the quantile function to get mean of that differences with CrI
round(quantile(diff_brms_shannon_T1, probs=c(0.025, 0.5, 0.975)),2)
#this is a result
mean(diff_brms_shannon_T1)


#calculate the mean differences between systems across the roofs at 13w
BIODYN_brms_shannon_T3<-(fit_shannon$BIODYN.C.13w+fit_shannon$BIODYN.RC.13w+fit_shannon$BIODYN.R.13w)/3
mean(BIODYN_brms_shannon_T3)
sd(BIODYN_brms_shannon_T3)
CONMIN_brms_shannon_T3<-(fit_shannon$CONMIN.C.13w+fit_shannon$CONMIN.RC.13w+fit_shannon$CONMIN.R.13w)/3
mean(CONMIN_brms_shannon_T3)
sd(CONMIN_brms_shannon_T3)

diff_brms_shannon_T3 <- (BIODYN_brms_shannon_T3-CONMIN_brms_shannon_T3)

#use the quantile function to get mean of that differences with CrI
round(quantile(diff_brms_shannon_T3, probs=c(0.025, 0.5, 0.975)),2)
#this is a result
mean(diff_brms_shannon_T3)

# Comparisons for Shannon

#calculate treatment effects
fitTest_shannon<-fitted(shannon_brms,
                        newdata = newdat,
                        summary = F,
                        re_formula = NA)

fitTest_shannon<-as.data.frame(fitTest_shannon)
colnames(fitTest_shannon) <- interaction(newdat$time, newdat$farming_system, newdat$treatment)

contrastShannon<-cbind(
  BD.RRC.4w<-fitTest_shannon$`4w.BIODYN.R`-fitTest_shannon$`4w.BIODYN.RC`,
  BD.RC.4w<-fitTest_shannon$`4w.BIODYN.R`-fitTest_shannon$`4w.BIODYN.C`,
  BD.RCC.4w<-fitTest_shannon$`4w.BIODYN.RC`-fitTest_shannon$`4w.BIODYN.C`,
  CM.RRC.4w<-fitTest_shannon$`4w.CONMIN.R`-fitTest_shannon$`4w.CONMIN.RC`,
  CM.RC.4w<-fitTest_shannon$`4w.CONMIN.R`-fitTest_shannon$`4w.CONMIN.C`,
  CM.RCC.4w<-fitTest_shannon$`4w.CONMIN.RC`-fitTest_shannon$`4w.CONMIN.C`,
  BD.RRC.13w<-fitTest_shannon$`13w.BIODYN.R`-fitTest_shannon$`13w.BIODYN.RC`,
  BD.RC.13w<-fitTest_shannon$`13w.BIODYN.R`-fitTest_shannon$`13w.BIODYN.C`,
  BD.RCC.13w<-fitTest_shannon$`13w.BIODYN.RC`-fitTest_shannon$`13w.BIODYN.C`,
  CM.RRC.13w<-fitTest_shannon$`13w.CONMIN.R`-fitTest_shannon$`13w.CONMIN.RC`,
  CM.RC.13w<-fitTest_shannon$`13w.CONMIN.R`-fitTest_shannon$`13w.CONMIN.C`,
  CM.RCC.13w<-fitTest_shannon$`13w.CONMIN.RC`-fitTest_shannon$`13w.CONMIN.C`,
  BD.CM.4w<-(fitTest_shannon$`4w.BIODYN.C`+fitTest_shannon$`4w.BIODYN.R`+fitTest_shannon$`4w.BIODYN.RC`)-
    (fitTest_shannon$`4w.CONMIN.C`+fitTest_shannon$`4w.CONMIN.R`+fitTest_shannon$`4w.CONMIN.RC`),
  BD.CM.13w<-(fitTest_shannon$`13w.BIODYN.C`+fitTest_shannon$`13w.BIODYN.R`+fitTest_shannon$`13w.BIODYN.RC`)-
    (fitTest_shannon$`13w.CONMIN.C`+fitTest_shannon$`13w.CONMIN.R`+fitTest_shannon$`13w.CONMIN.RC`))

colnames(contrastShannon)<-colnames2

#get the effect size data with 95% CrI
contrasts_shannon <- t(round_df(plyr::colwise(quantile)(as.data.frame(contrastShannon), 
                                                        probs = c(0.025, 0.5, 0.975)),2))

# add figure for water

View(water_content)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

water_summary <- data_summary(water_content, varname="water", 
                              groupnames=c("treatment", "farming_system", "time"))

water_summary$farming_system <- as.factor(water_summary$farming_system)
water_summary$time <- as.factor(water_summary$time)
#plot 4w before 13w
water_summary$time_order <- factor(water_summary$time, levels=c("4w", "13w"), labels = c("4 weeks", "13 weeks"))

water_summary$treatment<-factor(water_summary$treatment, levels=c("R", "RC", "C"))

p_water <- ggplot(water_summary, aes(farming_system, water)) +
  geom_errorbar(
    aes(ymin = water-sd, ymax = water+sd, color = treatment),
    position = position_dodge(0.3), width = 0.2
  ) +
  geom_point(aes(color = treatment), position = position_dodge(0.3), size=2) +
  facet_grid(.~time_order) + theme_bw() + xlab("Farming system") + ylab("WC / WHC") + 
  scale_color_manual(values = c("#B6854D", "#8D8680", "#0F0D0E")) + labs(color = "Drought treatment")

p_water

# BRMS model for water

brms_water <- brm(formula = (water)~time*treatment*farming_system+(1|block/plot),
                  data = water_content,
                  family="gaussian",
                  seed = 1041,
                  warmup = 2000,#The number of warmup iterations should 
                  #not be larger than iter and the default is iter/2.
                  iter = 6000, 
                  chains = 5,
                  control = list(adapt_delta = 0.999, max_treedepth = 15))

pp_check(brms_water)

#save an object to a file
saveRDS(brms_water, "data/rds/brms_water.rds")

#reload it so you don't have to run it again
brms_water <- readRDS("data/rds/brms_water.rds")

#check model fit 
yfit_water <- water_content$water
yrep_water<- posterior_predict(brms_water, draws = length(yfit_water))

launch_shinystan(brms_water)

#explore results
summary(brms_water)

#marginal effects
marginal_effects(brms_water)

## BRMS summary water

newdat_water <- expand.grid(
  time = factor(c("4w", "13w"),levels=levels(water_summary$time)), 
  farming_system = factor(c("BIODYN", "CONMIN"), levels=levels(water_summary$farming_system)),
  treatment = factor(c("C","R","RC"),levels=levels(water_summary$treatment)))
head(newdat_water)

fit_water = fitted(brms_water, newdata = newdat_water, robust=T, re_formula = NA,summary = TRUE, probs = c(0.025, 0.975))  
fit_water

# no colnames, I add them manually
fitdf_water_plot = cbind(newdat_water, fit_water)
fitdf_water_plot

#change order of treatment in fitted values
fitdf_water_plot$treatment<-factor(fitdf_water_plot$treatment, levels=c("R", "RC", "C"))
water_summary$treatment<-factor(water_summary$treatment, levels=c("R", "RC", "C"))

#prepare and extract results table
fitdf_water<-fitdf_water_plot[, c("farming_system","treatment", "time", "Estimate","Q2.5","Q97.5")]

#get location of the means including credible intervals from here:
fitdf_water_plot<-as.data.frame(fitdf_water_plot)

#this is the whole posterior, not summarized as before
fit_water = as.data.frame(fitted(brms_water, newdata=newdat_water,
                                 re_formula = NA, summary = FALSE))

#add names
colnames<-as.character(interaction(newdat_water$farming_system, newdat_water$treatment, newdat_water$time))

colnames(fit_water)<-colnames
names(fit_water)

#calculate the mean differences between systems across the roofs for 4w
BIODYN_brms_water_T1<-(fit_water$BIODYN.C.4w+fit_water$BIODYN.RC.4w+fit_water$BIODYN.R.4w)/3
mean(BIODYN_brms_water_T1)
sd(BIODYN_brms_obs_T1)
CONMIN_brms_water_T1<-(fit_water$CONMIN.C.4w+fit_water$CONMIN.RC.4w+fit_water$CONMIN.R.4w)/3
mean(CONMIN_brms_water_T1)
sd(CONMIN_brms_water_T1)

diff_brms_water_T1<-(BIODYN_brms_water_T1-CONMIN_brms_water_T1)

#use the quantile function to get mean of that differences with CrI
round(quantile(diff_brms_water_T1, probs=c(0.025, 0.5, 0.975)),2)
#this is a result
mean(diff_brms_water_T1)

#calculate the mean differences between systems across the roofs for 13w
BIODYN_brms_water_T3<-(fit_water$BIODYN.C.13w+fit_water$BIODYN.RC.13w+fit_water$BIODYN.R.13w)/3
mean(BIODYN_brms_water_T3)
sd(BIODYN_brms_water_T3)
CONMIN_brms_water_T3<-(fit_water$CONMIN.C.13w+fit_water$CONMIN.RC.13w+fit_water$CONMIN.R.13w)/3
mean(CONMIN_brms_water_T3)
sd(CONMIN_brms_water_T3)

diff_brms_water_T3<-(BIODYN_brms_water_T3-CONMIN_brms_water_T3)

#use the quantile function to get mean of that differences with CrI
round(quantile(diff_brms_water_T3, probs=c(0.025, 0.5, 0.975)),2)
#this is a clear result, it includes zero so its likely that there is no difference
mean(diff_brms_water_T3)

# Comparisons for soil moisture

#calculate treatment effects
fitTest<-fitted(brms_water,
                newdata = newdat_water,
                summary = F,
                re_formula = NA)

fitTest<-as.data.frame(fitTest)
colnames(fitTest) <- interaction(newdat_water$time, newdat_water$farming_system, newdat_water$treatment)

contrastWaterNew<-cbind(
  BD.RRC.4w<-fitTest$`4w.BIODYN.R`-fitTest$`4w.BIODYN.RC`,
  BD.RC.4w<-fitTest$`4w.BIODYN.R`-fitTest$`4w.BIODYN.C`,
  BD.RCC.4w<-fitTest$`4w.BIODYN.RC`-fitTest$`4w.BIODYN.C`,
  CM.RRC.4w<-fitTest$`4w.CONMIN.R`-fitTest$`4w.CONMIN.RC`,
  CM.RC.4w<-fitTest$`4w.CONMIN.R`-fitTest$`4w.CONMIN.C`,
  CM.RCC.4w<-fitTest$`4w.CONMIN.RC`-fitTest$`4w.CONMIN.C`,
  BD.RRC.13w<-fitTest$`13w.BIODYN.R`-fitTest$`13w.BIODYN.RC`,
  BD.RC.13w<-fitTest$`13w.BIODYN.R`-fitTest$`13w.BIODYN.C`,
  BD.RCC.13w<-fitTest$`13w.BIODYN.RC`-fitTest$`13w.BIODYN.C`,
  CM.RRC.13w<-fitTest$`13w.CONMIN.R`-fitTest$`13w.CONMIN.RC`,
  CM.RC.13w<-fitTest$`13w.CONMIN.R`-fitTest$`13w.CONMIN.C`,
  CM.RCC.13w<-fitTest$`13w.CONMIN.RC`-fitTest$`13w.CONMIN.C`,
  BD.CM.4w<-(fitTest$`4w.BIODYN.C`+fitTest$`4w.BIODYN.R`+fitTest$`4w.BIODYN.RC`)-
    (fitTest$`4w.CONMIN.C`+fitTest$`4w.CONMIN.R`+fitTest$`4w.CONMIN.RC`),
  BD.CM.13w<-(fitTest$`13w.BIODYN.C`+fitTest$`13w.BIODYN.R`+fitTest$`13w.BIODYN.RC`)-
    (fitTest$`13w.CONMIN.C`+fitTest$`13w.CONMIN.R`+fitTest$`13w.CONMIN.RC`))

colnames(contrastWaterNew)<-colnames2

#get the effect size data with 95% CrI
contrasts_water <- t(round_df(plyr::colwise(quantile)(as.data.frame(contrastWaterNew), 
                                                      probs = c(0.025, 0.5, 0.975)),2))
