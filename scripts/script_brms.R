

# ALPHA diversity indices

```{r alpha_div, echo=F}

#the alpha diversity metrics were calculated based on samples rarefied to the lowest number of reads, calcualted per sample
#alpha diversity estimates on rarefied data
#we calculate Observed and Shannon for each sample

alpha_div <- data.frame(estimate_richness(physeq_unite_Filtered_ASV, split=TRUE, measures = c("Observed", "Shannon")), sample_data(physeq_unite_Filtered_ASV))

alpha_div

plot_richness(physeq_unite_Filtered_ASV, x="treatment", color="farming_system", measures=c("Observed","Shannon")) + stat_boxplot(geom = "errorbar") + geom_boxplot() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("Treatment") + ggtitle("AMF alpha diversity") + ylab("Alpha metrics")

write.csv(alpha_div, "results/alpha_div.csv")

```
**Observed ASV richness and Shannon diversity index are calculated.**
  
  ## prepare data for BRMS model
  
  ```{r prepare data for BRMS model, echo=F}

#import data
head(alpha_div)

data_brms <- alpha_div
head(data_brms)
str(data_brms)

data_brms$farming_system <- as.factor(data_brms$farming_system)
data_brms$treatment <- as.factor(data_brms$treatment)
data_brms$time <- as.factor(data_brms$time)

```

+ negbinominal and gaussian BRMS model for observed richness 

```{r BRMS for observed richness, echo=F}

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

# Save an object to a file
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

```
**BRMS model for Observed richness is created.**
  
  ## plot BRMS observed richness
  
  ```{r plot brms for observed richness, echo=F}

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
write.table(fitdf_observed,"results/observed_brms.txt",sep="\t",quote=F)

#this is the dataframe i will use for plotting
pos.dodge <- position_dodge(0.7) # move them .7 to the left and right

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
ggsave("results/observed_brms.tiff", width = 18, height = 10, units = "cm")
ggsave("results/observed_brms.jpeg", width = 18, height = 10, units = "cm")
ggsave("results/observed_brms.pdf", width = 18, height = 10, units = "cm")

```


```{r, echo=F}
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

```

## BRMS model for Shannon index 

```{r BRMS for Shannon, echo=F}

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

# Save an object to a file
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

```
**BRMS model for Shannon index is created.**
  
  ## plot brms for Shannon index
  
  ```{r plot brms for Shannon, echo=F}

#use the function "fitted" from brms
#to get fitted means of the posterior for treatments
#It summarizes from the posterior and calculates meadian with 95% credible interval
#= we are 95% sure that the real estimate is within this interval

fit_shannon = fitted(shannon_brms, newdata = newdat, robust=T, re_formula = NA,summary = TRUE,probs = c(0.025, 0.975))
fit_shannon

# no colnames, I add them manually
fitdf_shannon_plot = cbind(newdat, fit_shannon)
fitdf_shannon_plot

#change order of treatment in fitted values
fitdf_shannon_plot$treatment<-factor(fitdf_shannon_plot$treatment, levels=c("R", "RC", "C"))

#prepare and extract results table
fitdf_shannon<-fitdf_shannon_plot[, c("farming_system","treatment", "time", "Estimate","Q2.5","Q97.5")]
write.table(fitdf_shannon,"results/shannon_brms.txt",sep="\t",quote=F)

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
ggsave("results/shannon_brms.tiff", width = 18, height = 10, units = "cm")
ggsave("results/shannon_brms.jpeg", width = 18, height = 10, units = "cm")
ggsave("results/shannon_brms.pdf", width = 18, height = 10, units = "cm")

```

```{r, echo=F}

#get location of the means including credible intervals from here:
fitdf_shannon_plot<-as.data.frame(fitdf_shannon_plot)

#if you want, you can also calculate differences between treatments based on the posterior
##calucate differences with CrI, you need the whole posterior
#the whole posterior distribution for me. I save it as dataframe (fit1)

fit_shannon = as.data.frame(fitted(shannon_brms, newdata=newdat,
                                   re_formula = NA, summary = FALSE))

#this is the whole posterior, not summarized as before

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

# add figure for water

```{r water content}

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
#plot 4w before 123w
water_summary$time_order <- factor(water_summary$time, levels=c("4w", "13w"), labels = c("4 weeks", "13 weeks"))

water_summary$treatment<-factor(water_summary$treatment, levels=c("R", "RC", "C"))

p_water <- ggplot(water_summary, aes(farming_system, water)) +
  geom_errorbar(
    aes(ymin = water-sd, ymax = water+sd, color = treatment),
    position = position_dodge(0.3), width = 0.2
  ) +
  geom_point(aes(color = treatment), position = position_dodge(0.3), size=2) +     facet_grid(.~time_order) + theme_bw() + xlab("Farming system") + ylab("WC / WHC") + scale_color_manual(values = c("#B6854D", "#8D8680", "#0F0D0E")) + labs(color = "Drought treatment")

p_water
ggsave("results/water.tiff", width = 20, height = 10, units = "cm", dpi=500)
ggsave("results/water.jpeg", width = 20, height = 10, units = "cm", dpi=500)
ggsave("results/water.pdf", width = 20, height = 10, units = "cm", dpi=500)

```

# BRMS model for water

```{r brms model for water}

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

# Save an object to a file
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

```

## plot BRMS water

```{r plot brms for water, echo=F}

newdat_water <- expand.grid(
  time = factor(c("4w", "13w"),levels=levels(water_summary$time)), 
  farming_system = factor(c("BIODYN", "CONMIN"), levels=levels(water_summary$farming_system)),
  treatment = factor(c("C","R","RC"),levels=levels(water_summary$treatment)))
head(newdat_water)

#use the function "fitted" from brms to get fitted means of the posterior for treatments
#It summarizes from the posterior and calculates meadian with 95% credible interval
#= we are 95% sure that the real estimate is within this interval

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
write.table(fitdf_water,"results/water_brms.txt",sep="\t",quote=F)

#this is the dataframe i will use for plotting
pos.dodge <- position_dodge(0.7) # move them .7 to the left and right

#make the final plot, without raw data
fitted_water<-ggplot(data = fitdf_water_plot, 
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

```

## BRMS summary water

```{r, echo=F}
#get location of the means including credible intervals from here:
fitdf_water_plot<-as.data.frame(fitdf_water_plot)

#if you want, you can also calculate differences between treatments based on the posterior
##calucate differences with CrI, you need the whole posterior
#the whole posterior distribution for me. I save it as dataframe (fit1)

fit_water = as.data.frame(fitted(brms_water, newdata=newdat_water,
                                 re_formula = NA, summary = FALSE))

#this is the whole posterior, not summarized as before

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

```