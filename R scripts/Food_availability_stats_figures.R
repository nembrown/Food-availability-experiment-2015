
# Read in packages --------------------------------------------------------
#### read in useful packages
library(here)
library(mgcv)
library(gratia)
library(tidyr)
library(bbmle) 
library(glmmTMB)
library(doBy)
library(plyr)
library(dplyr)
library(ggplot2) 
library(doBy)
library(grid)
library(lmtest)
library(fitdistrplus)
library(visreg)
library(lme4)
library(coefplot)
library(arm)
library(lmerTest)
library(boot)
library(MASS)
require(scales)
library(car)
library(knitr)
library(tidyverse)
library(kableExtra)
library(multcomp)
library(arm) ## for sim()
library(descr)  ## for LogRegR2
library(reshape2)
library(grid)
library(DHARMa)
library(gap)
library(qrnn)
library(colorspace)
library(cowplot)
library(devtools)
library(installr)


# dataframe management ----------------------------------------------------
#wd is top of project, use here() to find path

#data from inside the mesocosms:
food.exp.data.12.2019<-read.csv("C:Data//Mesocosm inventory data//food.exp.data.mesocosm.12.csv")

#data from the tiles:
food.exp.data.tile.all<-read.csv("C:Data//Mesocosm inventory data//food.exp.data.tile.all.csv")

#caprellid data:
food.caprellid.data<-read.csv("C:Data//Emily caprellid data.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )


#ordered and unordered factors
food.exp.data.12.2019$oFood.quality<-factor(food.exp.data.12.2019$Food.quality, levels=c("None", "Low", "High"), ordered=TRUE)
food.exp.data.12.2019$Food.quality<-factor(food.exp.data.12.2019$Food.quality, levels=c("None", "Low", "High"), ordered=FALSE)

#some variables need to be read over
food.exp.data.12.2019$caprellid.percent<-food.exp.data.tile.all$caprellid
food.exp.data.12.2019$hydroid<-food.exp.data.tile.all$hydroid
food.exp.data.12.2019$alive.bot<-food.exp.data.tile.all$alive.bot
food.exp.data.12.2019$formicula<-food.exp.data.tile.all$formicula
food.exp.data.12.2019$alive.mem<-food.exp.data.tile.all$alive.mem
food.exp.data.12.2019$didemnum<-food.exp.data.tile.all$didemnum
food.exp.data.12.2019$total<-food.exp.data.tile.all$total
food.exp.data.12.2019$bare<-food.exp.data.tile.all$bare
food.exp.data.12.2019$occupied.space<-100-food.exp.data.12.2019$bare
food.exp.data.12.2019$occupied.space.001<-0.01*(food.exp.data.12.2019$occupied.space)
food.exp.data.12.2019$everything.wet.weight<-food.exp.data.tile.all$everything.wet.weight
food.exp.data.12.2019$everything.wet.weight.per.1<-(food.exp.data.tile.all$everything.wet.weight)/food.exp.data.12.2019$occupied.space
food.exp.data.12.2019$Mussel.wet.weight<-food.exp.data.tile.all$Mussel.wet.weight
food.exp.data.12.2019$total_dry_biomass<-food.exp.data.tile.all$total_dry_biomass
food.exp.data.12.2019$total_dry_biomass_per1<-food.exp.data.tile.all$total_dry_biomass/food.exp.data.12.2019$occupied.space
food.exp.data.12.2019$hydroid_dry_biomass<-food.exp.data.tile.all$hydroid_dry_biomass
food.exp.data.12.2019$caprellid_dry_biomass<-food.exp.data.tile.all$caprellid_dry_biomass
food.exp.data.12.2019$tunicate_dry_biomass<-food.exp.data.tile.all$tunicate_dry_biomass
food.exp.data.12.2019$hydtobot<-(food.exp.data.12.2019$alive.bot)/(food.exp.data.12.2019$alive.bot+food.exp.data.12.2019$hydroid)
food.exp.data.12.2019$rest_dry_biomass<-food.exp.data.tile.all$rest_dry_biomass

#small negative biomass is within error of scale - change to zero
food.exp.data.12.2019$hydroid_dry_biomass[food.exp.data.12.2019$hydroid_dry_biomass<0]<-0
food.exp.data.12.2019$tunicate_dry_biomass[food.exp.data.12.2019$tunicate_dry_biomass<0]<-0
food.exp.data.12.2019$hydtobot_dry_biomass<-(food.exp.data.12.2019$tunicate_dry_biomass)/(food.exp.data.12.2019$tunicate_dry_biomass+food.exp.data.12.2019$hydroid_dry_biomass)


food.exp.data.12.2019$Mussel.wet.weight.per.1<-(food.exp.data.12.2019$Mussel.wet.weight)/(food.exp.data.12.2019$mussel_complete+1)

#making it a proportion instead of % cover
food.exp.data.12.2019$caprellid.percent.001<-(0.01*(food.exp.data.12.2019$caprellid.percent))+0.01
food.exp.data.12.2019$hydroid.001<-(0.01*(food.exp.data.12.2019$hydroid))+0.01
food.exp.data.12.2019$alive.bot.001<-(0.01*(food.exp.data.12.2019$alive.bot))+0.01
food.exp.data.12.2019$alive.mem.001<-(0.01*(food.exp.data.12.2019$alive.mem))+0.01
food.exp.data.12.2019$formicula.001<-(0.01*(food.exp.data.12.2019$formicula))+0.01
food.exp.data.12.2019$didemnum.001<-(0.01*(food.exp.data.12.2019$didemnum))+0.01

# need to have zscores for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
food.exp.data.12.2019_zscores<-food.exp.data.12.2019
food.exp.data.12.2019_zscores$hydrogen.concentration<-scale(food.exp.data.12.2019$hydrogen.concentration, center=TRUE, scale=TRUE)
food.exp.data.12.2019_zscores$av.pH<-scale(food.exp.data.12.2019$av.pH, center=TRUE, scale=TRUE)
food.exp.data.12.2019_zscores$min.10.pH<-scale(food.exp.data.12.2019$min.10.pH, center=TRUE, scale=TRUE)
food.exp.data.12.2019_zscores$Food.quality <- factor(food.exp.data.12.2019_zscores$Food.quality, levels = c("None", "Low","High" ))
food.exp.data.12.2019_zscores$Mesocosm <- as.factor(food.exp.data.tile.all$Mesocosm)
food.exp.data.12.2019_zscores$av.pH.unscaled <-food.exp.data.12.2019_zscores$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')
food.exp.data.12.2019_zscores$min.10.pH.unscaled <-food.exp.data.12.2019_zscores$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

food.caprellid.data_zscores<-food.caprellid.data
food.caprellid.data_zscores$min.10.pH<-scale(food.exp.data.12.2019$min.10.pH, center=TRUE, scale=TRUE)
food.caprellid.data_zscores$oFood.quality<-factor(food.caprellid.data_zscores$Food.quality, levels=c("None", "Low", "High"), ordered=TRUE)
food.caprellid.data_zscores$Food.quality<-factor(food.caprellid.data_zscores$Food.quality, levels=c("None", "Low", "High"), ordered=FALSE)


# Visualizing histograms of pH to use as continuous vs. discrete variable
# food.exp.data.12.2019_none<-food.exp.data.12.2019 %>% filter(Food.quality=="None")
# food.exp.data.12.2019_high<-food.exp.data.12.2019 %>% filter(Food.quality=="High")
# food.exp.data.12.2019_low<-food.exp.data.12.2019 %>% filter(Food.quality=="Low")
# 
# hist(food.exp.data.12.2019_none$min.10.pH, breaks=10)#10
# hist(food.exp.data.12.2019_high$min.10.pH, breaks=5)#5
# hist(food.exp.data.12.2019_low$min.10.pH, breaks=5)#5

# ggplot(food.exp.data.12.2019, aes(x=min.10.pH))+geom_density()+theme_classic()


# Notes on contrasts ------------------------------------------------------

#Notes on ordered factors: the factor food.quality is ordered
#options(contrasts = c("contr.sum", "contr.poly"))
#The contrast function, contr.sum(), gives orthogonal contrasts where you compare every level to the overall mean.


# GAM information and resources---------------------------------------------------------------------
#Gavin Simpson's gratia - used to visualize gams

#library(gratia)

#We have a varying coefficient model aka ANCOVA, so use "by"

#We're going to use gam(...,select=TRUE), automatic model selection via null space penalization

##Food.quality is an ordered factor - meaning "none" is the base level and the other two relate to that....
# see https://www.fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/



# Plotting settings -------------------------------------------------------

colorset2 = c("High"="#F8A02E" ,"Low"="#439E5F","None"= "#666666")
theme_set(theme_classic(base_size = 6))
theme_update(plot.margin = unit(c(0,0,0,0), "cm"))



# GAM beta hydroid / gam.beta.hydroid.12 --------------------------------------------------------

#distributions - binomial or beta with various families

gam.binomial.hydroid.12<- gam(formula = cbind(hydroid, 100-hydroid)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")
gam.beta.hydroid.12<- gam(hydroid.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydroid.12.1<- gam(hydroid.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydroid.12.2<- gam(hydroid.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.hydroid.12.3<- gam(hydroid.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.beta.hydroid.12, gam.beta.hydroid.12.1, gam.beta.hydroid.12.2, gam.beta.hydroid.12.3, gam.binomial.hydroid.12)
#cauchit is the best

plot(gam.beta.hydroid.12.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.hydroid.12.3)
#qq_plot(gam.beta.hydroid.12.3, method = 'simulate')
#k_check(gam.beta.hydroid.12.3)

gam.beta.hydroid.12.3.unordered<- gam(hydroid.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.beta.hydroid.12.3)
summary(gam.beta.hydroid.12.3.unordered)
# need an unordered and an unordered model to get estimates for overall Food.quality effect


### plotting based on gam confidence intervals
#using code from https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/

fam.gam.hydroid <- family(gam.beta.hydroid.12.3)
fam.gam.hydroid 
str(fam.gam.hydroid )
ilink.gam.hydroid <- fam.gam.hydroid$linkinv
ilink.gam.hydroid

food.exp.data.12.2019_zscores$min.10.pH.unscaled <-food.exp.data.12.2019_zscores$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')
head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)
            
mod.hydroid<-gam.beta.hydroid.12.3
ndata.hydroid <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydroid <- add_column(ndata.hydroid, fit = predict(mod.hydroid, newdata = ndata.hydroid, type = 'response'))
ndata.hydroid <- bind_cols(ndata.hydroid, setNames(as_tibble(predict(mod.hydroid, ndata.hydroid, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform
ndata.hydroid <- mutate(ndata.hydroid,
                fit_resp  = ilink.gam.hydroid(fit_link),
                right_upr = ilink.gam.hydroid(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.hydroid(fit_link - (2 * se_link)))

#make sure pH is unscaled in the plot
ndata.hydroid$min.10.pH.unscaled<-ndata.hydroid$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


plt.gam.hydroid <- ggplot(ndata.hydroid, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = hydroid.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Obelia")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.hydroid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.gam.hydroid 
ggsave("C:Data//Graphs March 2020//hydroid_pred.png")



# GAM beta botryllus / gam.beta.alive.bot.12.2 ----------------------------------------------------

#binomial first
gam.binomial.alive.bot.12<- gam(formula = cbind(alive.bot, 100-alive.bot)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, REML=TRUE)

#beta next
gam.beta.alive.bot.12<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.beta.alive.bot.12.1<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.beta.alive.bot.12.2<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.beta.alive.bot.12.3<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.beta.alive.bot.12, gam.beta.alive.bot.12.1, gam.beta.alive.bot.12.2, gam.beta.alive.bot.12.3, gam.binomial.alive.bot.12)
#cauchit is the best

plot(gam.beta.alive.bot.12.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.alive.bot.12.3)
#qq_plot(gam.beta.alive.bot.12.3, method = 'simulate')
#k_check(gam.beta.alive.bot.12.3)
summary(gam.beta.alive.bot.12.3)

gam.beta.alive.bot.12.3.unordered<- gam(alive.bot.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)


fam.gam.alive.bot <- family(gam.beta.alive.bot.12.3)
fam.gam.alive.bot
ilink.gam.alive.bot<- fam.gam.alive.bot$linkinv
ilink.gam.alive.bot


mod.alive.bot<-gam.beta.alive.bot.12.3
ndata.alive.bot <- with(food.exp.data.12.2019_zscores, 
                        data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                        length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the mod.alive.botel for the new data
ndata.alive.bot <- add_column(ndata.alive.bot, fit = predict(mod.alive.bot, newdata = ndata.alive.bot, type = 'response'))


ndata.alive.bot <- bind_cols(ndata.alive.bot, setNames(as_tibble(predict(mod.alive.bot, ndata.alive.bot, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.alive.bot <- mutate(ndata.alive.bot,
                fit_resp  = ilink.gam.alive.bot(fit_link),
                right_upr = ilink.gam.alive.bot(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.alive.bot(fit_link - (2 * se_link)))


food.exp.data.12.2019_zscores$min.10.pH.unscaled<-food.exp.data.12.2019_zscores$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

ndata.alive.bot$min.10.pH.unscaled<-ndata.alive.bot$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.alive.bot <- ggplot(ndata.alive.bot, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = alive.bot.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Botryllus")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.alive.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.alive.bot
ggsave("C:Data//Graphs March 2020//alive.bot_pred.png")



# GAM negbin caprellid / gam.nb.caprellid.12 -----------------------------------------------------------

#variety of potential distributions
gam.nb.caprellid.12<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = negbin(nbinom12.caprellid$estimate[[1]]), select=TRUE, method="REML")
gam.nb.caprellid.12.1<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.caprellid.12<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = poisson, select=TRUE, method="REML")
gam.lm.caprellid.12<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = gaussian, select=TRUE, method="REML")
gam.log.lm.caprellid.12<- gam(log(total.caprellids+1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = gaussian, select=TRUE, method="REML")


AICtab(gam.lm.caprellid.12, gam.log.lm.caprellid.12, gam.nb.caprellid.12.1,gam.nb.caprellid.12, gam.poisson.caprellid.12)
#gam.log.lm.caprellid.12 by far the best

#appraise(gam.log.lm.caprellid.12)
#qq_plot(gam.log.lm.caprellid.12, method = 'simulate')
plot(gam.log.lm.caprellid.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.log.lm.caprellid.12)
summary(gam.log.lm.caprellid.12.unordered)

gam.log.lm.caprellid.12.unordered<- gam(log(total.caprellids+1) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, select=TRUE, method="REML")

fam.gam.caprellid <- family(gam.log.lm.caprellid.12)
fam.gam.caprellid
ilink.gam.caprellid<- fam.gam.caprellid$linkinv
ilink.gam.caprellid

mod.caprellid<-gam.log.lm.caprellid.12
ndata.caprellid <- with(food.caprellid.data_zscores, 
                        data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                        length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.caprellid <- add_column(ndata.caprellid, fit = predict(mod.caprellid, newdata = ndata.caprellid, type = 'response'))

predict(mod.caprellid, newdata = ndata.caprellid, type = 'response')
ndata.caprellid <- bind_cols(ndata.caprellid, setNames(as_tibble(predict(mod.caprellid, ndata.caprellid, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid <- mutate(ndata.caprellid,
                          fit_resp  = ilink.gam.caprellid(fit_link),
                          right_upr = ilink.gam.caprellid(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.caprellid(fit_link - (2 * se_link)))


food.caprellid.data_zscores$min.10.pH.unscaled<-food.caprellid.data_zscores$min.10.pH * attr(food.caprellid.data_zscores$min.10.pH, 'scaled:scale') + attr(food.caprellid.data_zscores$min.10.pH, 'scaled:center')
ndata.caprellid$min.10.pH.unscaled<-ndata.caprellid$min.10.pH * attr(food.caprellid.data_zscores$min.10.pH, 'scaled:scale') + attr(food.caprellid.data_zscores$min.10.pH, 'scaled:center')




# plot 
plt.caprellid <- ggplot(ndata.caprellid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = log(total.caprellids+1), shape=CO2, colour=oFood.quality), data = food.caprellid.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Caprella")~ "abundance"), textstyle("(Log # of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid

ggsave("C:Data//Graphs March 2020//caprellid_pred.png")


# GAM caprellid percent -----------------------------------------------------------

gam.binomial.caprellid.percent.12<- gam(formula = cbind(caprellid.percent, 100-caprellid.percent)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.beta.caprellid.percent.12<- gam(caprellid.percent.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.caprellid.percent.12.1<- gam(caprellid.percent.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.caprellid.percent.12.2<- gam(caprellid.percent.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.caprellid.percent.12.3<- gam(caprellid.percent.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.caprellid.percent.12, gam.beta.caprellid.percent.12.1, gam.beta.caprellid.percent.12.3, gam.beta.caprellid.percent.12.2,gam.binomial.caprellid.percent.12)
#12,12.1, 12.2 are tied ... go with simplest logit


plot(gam.beta.caprellid.percent.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

#appraise(gam.beta.caprellid.percent.12)
#qq_plot(gam.beta.caprellid.percent.12, method = 'simulate')
#does not look great
#k_check(gam.beta.caprellid.percent.12)
summary(gam.beta.caprellid.percent.12.unordered)

gam.beta.caprellid.percent.12.unordered<- gam(caprellid.percent.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.caprellid.percent <- family(gam.beta.caprellid.percent.12)
fam.gam.caprellid.percent
str(fam.gam.caprellid.percent)
ilink.gam.caprellid.percent<- fam.gam.caprellid.percent$linkinv
ilink.gam.caprellid.percent


mod.caprellid.percent<-gam.beta.caprellid.percent.12
ndata.caprellid.percent <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.caprellid.percentel for the new data
ndata.caprellid.percent <- add_column(ndata.caprellid.percent, fit = predict(mod.caprellid.percent, newdata = ndata.caprellid.percent, type = 'response'))


ndata.caprellid.percent <- bind_cols(ndata.caprellid.percent, setNames(as_tibble(predict(mod.caprellid.percent, ndata.caprellid.percent, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid.percent <- mutate(ndata.caprellid.percent,
                          fit_resp  = ilink.gam.caprellid.percent(fit_link),
                          right_upr = ilink.gam.caprellid.percent(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.caprellid.percent(fit_link - (2 * se_link)))


ndata.caprellid.percent$min.10.pH.unscaled<-ndata.caprellid.percent$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.caprellid.percent <- ggplot(ndata.caprellid.percent, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = caprellid.percent.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Caprella")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid.percent,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid.percent
ggsave("C:Data//Graphs March 2020//caprellid.percent_pred.png")



# GAM beta formicula / gam.beta.formicula.12 -----------------------------------------------------------

gam.binomial.formicula.12<- gam(formula = cbind(formicula, 100-formicula)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")
gam.beta.formicula.12<- gam(formicula.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.formicula.12.1<- gam(formicula.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.formicula.12.2<- gam(formicula.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.formicula.12.3<- gam(formicula.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.formicula.12, gam.beta.formicula.12.1, gam.beta.formicula.12.2,gam.binomial.formicula.12, gam.beta.formicula.12.3)
#12, 12.1, 12.2 best - go with simplest logit


plot(gam.beta.formicula.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.formicula.12)
#qq_plot(gam.beta.formicula.12, method = 'simulate')
#k_check(gam.beta.formicula.12)
summary(gam.beta.formicula.12)

gam.beta.formicula.12.unordered<- gam(formicula.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.formicula <- family(gam.beta.formicula.12)
fam.gam.formicula
str(fam.gam.formicula)
ilink.gam.formicula<- fam.gam.formicula$linkinv
ilink.gam.formicula


mod.formicula<-gam.beta.formicula.12
ndata.formicula <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.formiculael for the new data
ndata.formicula <- add_column(ndata.formicula, fit = predict(mod.formicula, newdata = ndata.formicula, type = 'response'))


ndata.formicula <- bind_cols(ndata.formicula, setNames(as_tibble(predict(mod.formicula, ndata.formicula, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.formicula <- mutate(ndata.formicula,
                          fit_resp  = ilink.gam.formicula(fit_link),
                          right_upr = ilink.gam.formicula(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.formicula(fit_link - (2 * se_link)))


ndata.formicula$min.10.pH.unscaled<-ndata.formicula$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.formicula <- ggplot(ndata.formicula, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = formicula.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Folliculina")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.formicula,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.formicula
ggsave("C:Data//Graphs March 2020//formicula_pred.png")



# GAM beta membranipora / gam.beta.alive.mem.12 --------------------------------------------------------

#binomial first
gam.binomial.alive.mem.12<- gam(formula = cbind(alive.mem, 100-alive.mem)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.beta.alive.mem.12<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.alive.mem.12.1<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.alive.mem.12.2<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.alive.mem.12.3<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.beta.alive.mem.12, gam.beta.alive.mem.12.1, gam.beta.alive.mem.12.2, gam.binomial.alive.mem.12, gam.beta.alive.mem.12.3)
#12, 12.1, 12.2, 12.3 all equal, go with logit


plot(gam.beta.alive.mem.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.alive.mem.12)
#qq_plot(gam.beta.alive.mem.12, method = 'simulate')
#k_check(gam.beta.alive.mem.12)
summary(gam.beta.alive.mem.12)
vis.gam(gam.beta.alive.mem.12)

gam.beta.alive.mem.12.unordered<- gam(alive.mem.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.alive.mem <- family(gam.beta.alive.mem.12)
fam.gam.alive.mem
ilink.gam.alive.mem<- fam.gam.alive.mem$linkinv
ilink.gam.alive.mem


mod.alive.mem<-gam.beta.alive.mem.12
ndata.alive.mem <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.alive.memel for the new data
ndata.alive.mem <- add_column(ndata.alive.mem, fit = predict(mod.alive.mem, newdata = ndata.alive.mem, type = 'response'))


ndata.alive.mem <- bind_cols(ndata.alive.mem, setNames(as_tibble(predict(mod.alive.mem, ndata.alive.mem, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.alive.mem <- mutate(ndata.alive.mem,
                          fit_resp  = ilink.gam.alive.mem(fit_link),
                          right_upr = ilink.gam.alive.mem(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.alive.mem(fit_link - (2 * se_link)))

ndata.alive.mem$min.10.pH.unscaled<-ndata.alive.mem$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.alive.mem <- ggplot(ndata.alive.mem, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = alive.mem.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Membranipora")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.alive.mem,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.alive.mem
ggsave("C:Data//Graphs March 2020//alive.mem_pred.png")


# GAM beta didemnum / gam.beta.didemnum.12 ------------------------------------------------------------

#binomial first
gam.binomial.didemnum.12<- gam(formula = cbind(didemnum, 100-didemnum)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

#beta next
gam.beta.didemnum.12<- gam(didemnum.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.didemnum.12.1<- gam(didemnum.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.didemnum.12.2<- gam(didemnum.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.didemnum.12.3<- gam(didemnum.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.didemnum.12, gam.beta.didemnum.12.1, gam.beta.didemnum.12.2,  gam.binomial.didemnum.12, gam.beta.didemnum.12.3)
#logit, all the betas are equal go with logit

plot(gam.beta.didemnum.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.didemnum.12)
#qq_plot(gam.beta.didemnum.12, method = 'simulate')
#not very good
#k_check(gam.beta.didemnum.12)
summary(gam.beta.didemnum.12)
vis.gam(gam.beta.didemnum.12)

#appraise doesn't fit that well .... but I think there's just not enough data

gam.beta.didemnum.12.unordered<- gam(didemnum.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")

fam.gam.didemnum <- family(gam.beta.didemnum.12)
fam.gam.didemnum
str(fam.gam.didemnum)
ilink.gam.didemnum<- fam.gam.didemnum$linkinv
ilink.gam.didemnum


mod.didemnum<-gam.beta.didemnum.12
ndata.didemnum <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.didemnumel for the new data
ndata.didemnum <- add_column(ndata.didemnum, fit = predict(mod.didemnum, newdata = ndata.didemnum, type = 'response'))


ndata.didemnum <- bind_cols(ndata.didemnum, setNames(as_tibble(predict(mod.didemnum, ndata.didemnum, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.didemnum <- mutate(ndata.didemnum,
                          fit_resp  = ilink.gam.didemnum(fit_link),
                          right_upr = ilink.gam.didemnum(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.didemnum(fit_link - (2 * se_link)))

ndata.didemnum$min.10.pH.unscaled<-ndata.didemnum$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.didemnum <- ggplot(ndata.didemnum, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = didemnum.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Didemnum")~ "abundance"), textstyle("(proportion cover)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.didemnum,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.didemnum
ggsave("C:Data//Graphs March 2020//didemnum_pred.png")



# GAM negbin mussel incomplete / gam.nb.mussel_complete.12 ---------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$mussel_complete, "Poisson")
qqp(food.exp.data.12.2019_zscores$mussel_complete, "pois", lambda=poisson.12$estimate[[1]])
#estimating lambda

nbinom12.mussel <- fitdistr(food.exp.data.12.2019_zscores$mussel_complete, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$mussel_complete, "nbinom", size = nbinom12.mussel$estimate[[1]], mu = nbinom12.mussel$estimate[[2]])
#estimating theta

gam.nb.mussel_complete.12<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")

gam.nb.mussel_complete.12.1<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(link="log"), select=TRUE, method="REML")
gam.nb.mussel_complete.12.2<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(link="sqrt"), select=TRUE, method="REML")
gam.nb.mussel_complete.12.3<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")

gam.poisson.mussel_complete.12<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson(link="log"), select=TRUE, method="REML")
gam.poisson.mussel_complete.12.1<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson(link="identity"), select=TRUE, method="REML")
#gam.poisson.mussel_complete.12.2<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson(link="sqrt"), select=TRUE, method="REML")
#won't run

AICtab(gam.nb.mussel_complete.12, glm.nb.mussel_complete.12.hydrogen, gam.nb.mussel_complete.12.1, gam.nb.mussel_complete.12.2, gam.nb.mussel_complete.12.3, gam.poisson.mussel_complete.12.1, gam.poisson.mussel_complete.12)
#gam with theta estaimted from data is best
#gam.nb.mussel_complete.12 


plot(gam.nb.mussel_complete.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.nb.mussel_complete.12)
#qq_plot(gam.nb.mussel_complete.12, method = 'simulate')
#looks good!
#k_check(gam.nb.mussel_complete.12)
summary(gam.nb.mussel_complete.12)


gam.nb.mussel_complete.12.unordered<- gam(mussel_complete ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")
summary(gam.nb.mussel_complete.12.unordered)


fam.gam.mussel_complete <- family(gam.nb.mussel_complete.12)
fam.gam.mussel_complete
ilink.gam.mussel_complete<- fam.gam.mussel_complete$linkinv
ilink.gam.mussel_complete

mod.mussel_complete<-gam.nb.mussel_complete.12
ndata.mussel_complete <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.mussel_complete <- add_column(ndata.mussel_complete, fit = predict(mod.mussel_complete, newdata = ndata.mussel_complete, type = 'response'))

predict(mod.mussel_complete, newdata = ndata.mussel_complete, type = 'response')
ndata.mussel_complete <- bind_cols(ndata.mussel_complete, setNames(as_tibble(predict(mod.mussel_complete, ndata.mussel_complete, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.mussel_complete <- mutate(ndata.mussel_complete,
                          fit_resp  = ilink.gam.mussel_complete(fit_link),
                          right_upr = ilink.gam.mussel_complete(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.mussel_complete(fit_link - (2 * se_link)))

ndata.mussel_complete$min.10.pH.unscaled<-ndata.mussel_complete$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.mussel_complete <- ggplot(ndata.mussel_complete, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = mussel_complete, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Mytilus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.mussel_complete,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.mussel_complete
ggsave("C:Data//Graphs March 2020//mussel_complete_pred.png")


# GAM negbin barnacles / gam.nb.num.barn.alive.12 -----------------------------------------------------------


nbinom12.barn.alive <- fitdistr(food.exp.data.12.2019_zscores$num.barn.alive, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.barn.alive, "nbinom", size = nbinom12.barn.alive$estimate[[1]], mu = nbinom12.mussel$estimate[[2]])

#negative binomial 
gam.nb.num.barn.alive.12<- gam(num.barn.alive ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.barnacles$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.barn.alive.12.1<- gam(num.barn.alive ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.barn.alive.12<- gam(num.barn.alive ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.nb.num.barn.alive.12, glm.nb.num.barn.alive.12.hydrogen, gam.nb.num.barn.alive.12.1, gam.poisson.num.barn.alive.12)

plot(gam.nb.num.barn.alive.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.nb.num.barn.alive.12)
#qq_plot(gam.nb.num.barn.alive.12, method = 'simulate')
#looks really good
#k_check(gam.nb.num.barn.alive.12)
summary(gam.nb.num.barn.alive.12)


#a few outside the area
#appraise a bit funnelly

gam.nb.num.barn.alive.12.unordered<- gam(num.barn.alive ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.barnacles$estimate[[1]]), select=TRUE, method="REML")


fam.gam.num.barn.alive <- family(gam.nb.num.barn.alive.12)
fam.gam.num.barn.alive
str(fam.gam.num.barn.alive)
ilink.gam.num.barn.alive<- fam.gam.num.barn.alive$linkinv
ilink.gam.num.barn.alive


mod.num.barn.alive<-gam.nb.num.barn.alive.12
ndata.num.barn.alive <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.num.barn.alive


## add the fitted values by predicting from the model for the new data
ndata.num.barn.alive <- add_column(ndata.num.barn.alive, fit = predict(mod.num.barn.alive, newdata = ndata.num.barn.alive, type = 'response'))

predict(mod.num.barn.alive, newdata = ndata.num.barn.alive, type = 'response')
ndata.num.barn.alive <- bind_cols(ndata.num.barn.alive, setNames(as_tibble(predict(mod.num.barn.alive, ndata.num.barn.alive, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.barn.alive <- mutate(ndata.num.barn.alive,
                                  fit_resp  = ilink.gam.num.barn.alive(fit_link),
                                  right_upr = ilink.gam.num.barn.alive(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.num.barn.alive(fit_link - (2 * se_link)))


ndata.num.barn.alive$min.10.pH.unscaled<-ndata.num.barn.alive$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.num.barn.alive <- ggplot(ndata.num.barn.alive, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = num.barn.alive, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Balanus")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.num.barn.alive,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.num.barn.alive
ggsave("C:Data//Graphs March 2020//num.barn.alive_pred.png")



# GAM negbin disporella / gam.nb.disporella.12 ----------------------------------------------------------

nbinom12.disporella <- fitdistr(food.exp.data.12.2019_zscores$disporella, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$disporella, "nbinom", size = nbinom12.disporella$estimate[[1]], mu = nbinom12.disporella$estimate[[2]])
#getting theta

gam.nb.disporella.12<- gam(disporella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")
gam.nb.disporella.12.1<- gam(disporella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.disporella.12<- gam(disporella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson(), select=TRUE, method="REML")

AICtab(gam.nb.disporella.12.1,gam.nb.disporella.12,gam.poisson.disporella.12)
#used estimated theta

plot(gam.nb.disporella.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.nb.disporella.12)
#not bad!
#qq_plot(gam.nb.disporella.12, method = 'simulate')
#k_check(gam.nb.disporella.12)
summary(gam.nb.disporella.12)

#a few outside the area
#appraise a bit funnelly

gam.nb.disporella.12.unordered<- gam(disporella ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")


fam.gam.disporella <- family(gam.nb.disporella.12)
fam.gam.disporella
str(fam.gam.disporella)
ilink.gam.disporella<- fam.gam.disporella$linkinv
ilink.gam.disporella


mod.disporella<-gam.nb.disporella.12
ndata.disporella <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                       length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.disporella

str(food.exp.data.12.2019_zscores)
str(ndata.disporella)

## add the fitted values by predicting from the model for the new data
ndata.disporella <- add_column(ndata.disporella, fit = predict(mod.disporella, newdata = ndata.disporella, type = 'response'))

predict(mod.disporella, newdata = ndata.disporella, type = 'response')
ndata.disporella <- bind_cols(ndata.disporella, setNames(as_tibble(predict(mod.disporella, ndata.disporella, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.disporella <- mutate(ndata.disporella,
                               fit_resp  = ilink.gam.disporella(fit_link),
                               right_upr = ilink.gam.disporella(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.disporella(fit_link - (2 * se_link)))

ndata.disporella$min.10.pH.unscaled<-ndata.disporella$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.disporella <- ggplot(ndata.disporella, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = disporella, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Disporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.disporella,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.disporella
ggsave("C:Data//Graphs March 2020//disporella_pred.png")


# GAM negbin schizo / gam.nb.schizo.12 --------------------------------------------------------------

nbinom12.schizo <- fitdistr(food.exp.data.12.2019_zscores$schizo, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$schizo, "nbinom", size = nbinom12.schizo$estimate[[1]], mu = nbinom12.schizo$estimate[[2]])
#getting theta

gam.nb.schizo.12<- gam(schizo ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.schizo$estimate[[1]]), select=TRUE, method="REML")
gam.nb.schizo.12.1<- gam(schizo ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.schizo.12<- gam(schizo ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.schizo.12, gam.nb.schizo.12.1, gam.poisson.schizo.12)


plot(gam.nb.schizo.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.nb.schizo.12)
#looks pretty good
#qq_plot(gam.nb.schizo.12, method = 'simulate')
#k_check(gam.nb.schizo.12)
summary(gam.nb.schizo.12)

gam.nb.schizo.12.unordered<- gam(schizo ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.schizo$estimate[[1]]), select=TRUE, method="REML")
summary(gam.nb.schizo.12.unordered)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.schizo <- family(gam.nb.schizo.12)
fam.gam.schizo
str(fam.gam.schizo)
ilink.gam.schizo<- fam.gam.schizo$linkinv
ilink.gam.schizo

mod.schizo<-gam.nb.schizo.12
ndata.schizo <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                   length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.schizo

str(food.exp.data.12.2019_zscores)
str(ndata.schizo)

## add the fitted values by predicting from the model for the new data
ndata.schizo <- add_column(ndata.schizo, fit = predict(mod.schizo, newdata = ndata.schizo, type = 'response'))

predict(mod.schizo, newdata = ndata.schizo, type = 'response')
ndata.schizo <- bind_cols(ndata.schizo, setNames(as_tibble(predict(mod.schizo, ndata.schizo, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))

## create the interval and backtransform

ndata.schizo <- mutate(ndata.schizo,
                           fit_resp  = ilink.gam.schizo(fit_link),
                           right_upr = ilink.gam.schizo(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.schizo(fit_link - (2 * se_link)))


ndata.schizo$min.10.pH.unscaled<-ndata.schizo$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 
plt.schizo <- ggplot(ndata.schizo, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = schizo, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Schizoporella")~ "abundance"), textstyle("(# of colonies)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.schizo,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.schizo
ggsave("C:Data//Graphs March 2020//schizo_pred.png")


# GAM poisson num nudi / gam.poisson.num.nudi.12  ------------------------------------------------------------
nbinom12.num.nudi <- fitdistr(food.exp.data.12.2019_zscores$num.nudi, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.nudi, "nbinom", size = nbinom12.num.nudi$estimate[[1]], mu = nbinom12.num.nudi$estimate[[2]])
#theta

gam.nb.num.nudi.12.1<- gam(num.nudi ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.nb.num.nudi.12<- gam(num.nudi ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.nudi$estimate[[1]]), select=TRUE, method="REML")
gam.poisson.num.nudi.12<- gam(num.nudi ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.nudi.12, gam.nb.num.nudi.12.1, gam.poisson.num.nudi.12)
###poisson is the best fit by 2 dAIC


#appraise(gam.poisson.num.nudi.12)
#qq_plot(gam.poisson.num.nudi.12, method = 'simulate')
#looks quite good!
plot(gam.poisson.num.nudi.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.poisson.num.nudi.12)
summary(gam.poisson.num.nudi.12)
#resids a bit funny but same in neg bin

gam.poisson.num.nudi.12.unordered<- gam(num.nudi ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")
want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.num.nudi <- family(gam.poisson.num.nudi.12)
ilink.gam.num.nudi<- fam.gam.num.nudi$linkinv

mod.num.nudi<-gam.poisson.num.nudi.12
ndata.num.nudi <- with(food.exp.data.12.2019_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.num.nudi <- add_column(ndata.num.nudi, fit = predict(mod.num.nudi, newdata = ndata.num.nudi, type = 'response'))

predict(mod.num.nudi, newdata = ndata.num.nudi, type = 'response')
ndata.num.nudi <- bind_cols(ndata.num.nudi, setNames(as_tibble(predict(mod.num.nudi, ndata.num.nudi, se.fit = TRUE)[1:2]),
                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.nudi <- mutate(ndata.num.nudi,
                       fit_resp  = ilink.gam.num.nudi(fit_link),
                       right_upr = ilink.gam.num.nudi(fit_link + (2 * se_link)),
                       right_lwr = ilink.gam.num.nudi(fit_link - (2 * se_link)))


ndata.num.nudi$min.10.pH.unscaled<-ndata.num.nudi$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.nudi <- ggplot(ndata.num.nudi, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = num.nudi, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Hermissenda")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0, "cm"), legend.margin=margin(0, 0.05, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.num.nudi
ggsave("C:Data//Graphs March 2020//num.nudi_pred.png")


# GAM nb() serpulids / gam.nb.num.serpulid.12.1 -----------------------------------------------------------

nbinom12.num.serpulid <- fitdistr(food.exp.data.12.2019_zscores$num.serpulid, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.serpulid, "nbinom", size = nbinom12.num.serpulid$estimate[[1]], mu = nbinom12.num.serpulid$estimate[[2]])
#theta

#negative binomial first
gam.nb.num.serpulid.12<- gam(num.serpulid ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.serpulid.12.1<- gam(num.serpulid ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.serpulid.12<- gam(num.serpulid ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.serpulid.12, gam.nb.num.serpulid.12.1, gam.poisson.num.serpulid.12)

##gam.nb.num.serpulid.12.1 is best

#appraise(gam.nb.num.serpulid.12.1)
#qq_plot(gam.nb.num.serpulid.12.1, method = 'simulate')
#looks good
plot(gam.nb.num.serpulid.12.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.nb.num.serpulid.12.1)
summary(gam.nb.num.serpulid.12.1)

#residuals a bit weird .... but they are the same in glm as in gam (doesn't improve)
#I think because of the zeros? 

gam.nb.num.serpulid.12.1.unordered<- gam(num.serpulid ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.num.serpulid <- family(gam.nb.num.serpulid.12.1)
fam.gam.num.serpulid
str(fam.gam.num.serpulid)
ilink.gam.num.serpulid<- fam.gam.num.serpulid$linkinv
ilink.gam.num.serpulid

mod.num.serpulid<-gam.nb.num.serpulid.12.1
ndata.num.serpulid <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.num.serpulid <- add_column(ndata.num.serpulid, fit = predict(mod.num.serpulid, newdata = ndata.num.serpulid, type = 'response'))

predict(mod.num.serpulid, newdata = ndata.num.serpulid, type = 'response')
ndata.num.serpulid <- bind_cols(ndata.num.serpulid, setNames(as_tibble(predict(mod.num.serpulid, ndata.num.serpulid, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.serpulid <- mutate(ndata.num.serpulid,
                         fit_resp  = ilink.gam.num.serpulid(fit_link),
                         right_upr = ilink.gam.num.serpulid(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.num.serpulid(fit_link - (2 * se_link)))


ndata.num.serpulid$min.10.pH.unscaled<-ndata.num.serpulid$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 
plt.num.serpulid <- ggplot(ndata.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = num.serpulid, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.num.serpulid
ggsave("C:Data//Graphs March 2020//num.serpulid_pred.png")

# colorset2 = c("High"="#F8A02E" ,"Low"="#439E5F","None"= "#666666")
# colorset_none = c("High"="#FFFFFF" ,"Low"="#FFFFFF","None"= "#666666")
# colorset_low = c("High"="#FFFFFF" ,"Low"="#439E5F","None"= "#666666")
# 
# ### Plot for powerpoint:
# plt.num.serpulid <- ggplot(ndata.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
#   
#   geom_line(aes(colour=oFood.quality)) +
#   geom_point(aes(y = num.serpulid, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
#   xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
#   scale_color_manual(values=colorset_none)+
#   scale_fill_manual(values=colorset_none)+
#   scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
#   geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
#   theme(legend.position='none')
# plt.num.serpulid
# ggsave("C:Data//Graphs March 2020//num.serpulid_pred_none.png")
# 
# plt.num.serpulid <- ggplot(ndata.num.serpulid, aes(x = min.10.pH.unscaled, y = fit)) + 
#   
#   geom_line(aes(colour=oFood.quality)) +
#   geom_point(aes(y = num.serpulid, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
#   xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Serpulid abundance"), textstyle("(# of individuals)")))))+  
#   scale_color_manual(values=colorset_low)+
#   scale_fill_manual(values=colorset_low)+
#   scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
#   geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
#   theme(legend.position='none')
# plt.num.serpulid
# ggsave("C:Data//Graphs March 2020//num.serpulid_pred_low.png")
# 
# GAM negbin orange sponge / gam.nb.orange_sponge.12 -------------------------------------------------------

nbinom12.orange_sponge <- fitdistr(food.exp.data.12.2019_zscores$orange_sponge, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$orange_sponge, "nbinom", size = nbinom12.orange_sponge$estimate[[1]], mu = nbinom12.orange_sponge$estimate[[2]])
#theta

gam.nb.orange_sponge.12<- gam(orange_sponge ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.orange_sponge$estimate[[1]]), select=TRUE, method="REML")
gam.nb.orange_sponge.12.1<- gam(orange_sponge ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.orange_sponge.12<- gam(orange_sponge ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.orange_sponge.12, gam.nb.orange_sponge.12.1,  gam.poisson.orange_sponge.12)

###.12

#appraise(gam.nb.orange_sponge.12)
#a bit funnel-y
#qq_plot(gam.nb.orange_sponge.12, method = 'simulate')
plot(gam.nb.orange_sponge.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.nb.orange_sponge.12)
summary(gam.nb.orange_sponge.12)

gam.nb.orange_sponge.12.unordered<- gam(orange_sponge ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.orange_sponge$estimate[[1]]), select=TRUE, method="REML")

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)
fam.gam.orange_sponge <- family(gam.nb.orange_sponge.12)
fam.gam.orange_sponge
str(fam.gam.orange_sponge)
ilink.gam.orange_sponge<- fam.gam.orange_sponge$linkinv
ilink.gam.orange_sponge

mod.orange_sponge<-gam.nb.orange_sponge.12
ndata.orange_sponge <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                     length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.orange_sponge <- add_column(ndata.orange_sponge, fit = predict(mod.orange_sponge, newdata = ndata.orange_sponge, type = 'response'))

predict(mod.orange_sponge, newdata = ndata.orange_sponge, type = 'response')
ndata.orange_sponge <- bind_cols(ndata.orange_sponge, setNames(as_tibble(predict(mod.orange_sponge, ndata.orange_sponge, se.fit = TRUE)[1:2]),
                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.orange_sponge <- mutate(ndata.orange_sponge,
                             fit_resp  = ilink.gam.orange_sponge(fit_link),
                             right_upr = ilink.gam.orange_sponge(fit_link + (2 * se_link)),
                             right_lwr = ilink.gam.orange_sponge(fit_link - (2 * se_link)))


ndata.orange_sponge$min.10.pH.unscaled<-ndata.orange_sponge$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 
plt.orange_sponge <- ggplot(ndata.orange_sponge, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = orange_sponge, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle("Sponge abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.orange_sponge,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.orange_sponge
ggsave("C:Data//Graphs March 2020//orange_sponge_pred.png")



# GAM negbin corella / gam.nb.num.corella.12 -------------------------------------------------------------

nbinom12.num.corella <- fitdistr(food.exp.data.12.2019_zscores$num.corella, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.corella, "nbinom", size = nbinom12.num.corella$estimate[[1]], mu = nbinom12.num.corella$estimate[[2]])
#extracting theta

gam.nb.num.corella.12<- gam(num.corella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.corella.12.1<- gam(num.corella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.num.corella.12<- gam(num.corella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.corella.12, gam.nb.num.corella.12.1, gam.poisson.num.corella.12)
#.12

#appraise(gam.nb.num.corella.12)
#looks pretty good - slight pattern
#qq_plot(gam.nb.num.corella.12, method = 'simulate')
plot(gam.nb.num.corella.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.nb.num.corella.12)
summary(gam.nb.num.corella.12)


gam.nb.num.corella.12.unordered<- gam(num.corella ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)
fam.gam.num.corella <- family(gam.nb.num.corella.12)
ilink.gam.num.corella<- fam.gam.num.corella$linkinv


mod.num.corella<-gam.nb.num.corella.12
ndata.num.corella <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                      length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.num.corella <- add_column(ndata.num.corella, fit = predict(mod.num.corella, newdata = ndata.num.corella, type = 'response'))

predict(mod.num.corella, newdata = ndata.num.corella, type = 'response')
ndata.num.corella <- bind_cols(ndata.num.corella, setNames(as_tibble(predict(mod.num.corella, ndata.num.corella, se.fit = TRUE)[1:2]),
                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.num.corella <- mutate(ndata.num.corella,
                              fit_resp  = ilink.gam.num.corella(fit_link),
                              right_upr = ilink.gam.num.corella(fit_link + (2 * se_link)),
                              right_lwr = ilink.gam.num.corella(fit_link - (2 * se_link)))


ndata.num.corella$min.10.pH.unscaled<-ndata.num.corella$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

par(lheight=0.2) 
# plot 
plt.num.corella <- ggplot(ndata.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = num.corella, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(NA,atop(textstyle(italic("Corella")~ "abundance"), textstyle("(# of individuals)")))))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.num.corella
ggsave("C:Data//Graphs March 2020//num.corella_pred.png")



# GAM poisson clam / gam.poisson.clam.12 ----------------------------------------------------------------

nbinom12.clam <- fitdistr(food.exp.data.12.2019_zscores$clam, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$clam, "nbinom", size = nbinom12.clam$estimate[[1]], mu = nbinom12.clam$estimate[[2]])
#theta

#negative binomial first
gam.nb.clam.12<- gam(clam ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.clam$estimate[[1]]), select=TRUE, method="REML")
gam.nb.clam.12.1<- gam(clam ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.clam.12<- gam(clam ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.clam.12, gam.nb.clam.12.1, gam.poisson.clam.12)

###gam poisson is best

#appraise(gam.poisson.clam.12)
#qq_plot(gam.poisson.clam.12, method = 'simulate')
plot(gam.poisson.clam.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#k_check(gam.poisson.clam.12)
summary(gam.poisson.clam.12)

#residuals a bit patterny as well 
gam.poisson.clam.12.unordered<- gam(clam ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")
summary(gam.poisson.clam.12.unordered)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.clam <- family(gam.poisson.clam.12)
fam.gam.clam
str(fam.gam.clam)
ilink.gam.clam<- fam.gam.clam$linkinv
ilink.gam.clam

mod.clam<-gam.poisson.clam.12
ndata.clam <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                    length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.clam <- add_column(ndata.clam, fit = predict(mod.clam, newdata = ndata.clam, type = 'response'))

predict(mod.clam, newdata = ndata.clam, type = 'response')
ndata.clam <- bind_cols(ndata.clam, setNames(as_tibble(predict(mod.clam, ndata.clam, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.clam <- mutate(ndata.clam,
                            fit_resp  = ilink.gam.clam(fit_link),
                            right_upr = ilink.gam.clam(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.clam(fit_link - (2 * se_link)))


ndata.clam$min.10.pH.unscaled<-ndata.clam$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 
plt.clam <- ggplot(ndata.clam, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = clam, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Clam abundance\n(# of individuals)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.clam,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.clam
ggsave("C:Data//Graphs March 2020//clam_pred.png")


# Fig 2 plot generation ---------------------------------------------------
library(patchwork)
fig.2<-wrap_plots(plt.gam.hydroid,plt.alive.bot,plt.formicula,plt.caprellid.percent,plt.alive.mem,plt.didemnum,
          plt.mussel_complete,plt.num.barn.alive,plt.disporella,plt.schizo,plt.num.nudi,plt.num.serpulid,
          plt.orange_sponge,plt.num.corella,plt.clam, ncol=5)+
          plot_annotation(tag_levels = 'a')

fig.2

ggplot2::ggsave(plot=fig.2, "C:Data//For submission//For resubmission//RESUB2//First look//Fig2_resized.pdf", width=8.75, height=4.75, units="in")

#ggplot2::ggsave("C:Data//For submission//For resubmission//Fig2.png", width=65, height=35, units="cm")


# Pulling model results to a table ----------------------------------------

hydroid.gam<-summary(gam.beta.hydroid.12.3)
botryllus.gam<-summary(gam.beta.alive.bot.12.3)
caprellid.gam<-summary(gam.beta.caprellid.percent.12)
formicula.gam<-summary(gam.beta.formicula.12)
alive.mem.gam<-summary(gam.beta.alive.mem.12)
didemnum.gam<-summary(gam.beta.didemnum.12)
mussel_complete.gam<-summary(gam.nb.mussel_complete.12)
alive.barn.gam<-summary(gam.nb.num.barn.alive.12)
disporella.gam<-summary(gam.nb.disporella.12)
schizo.gam<-summary(gam.nb.schizo.12)
num.nudi.gam<-summary(gam.poisson.num.nudi.12)
num.serpulid.gam<-summary(gam.nb.num.serpulid.12.1)
clam.gam<-summary(gam.poisson.clam.12)
corella.gam<-summary(gam.nb.num.corella.12)
orange_sponge.gam<-summary(gam.nb.orange_sponge.12)

hydroid.gam.unordered<-summary(gam.beta.hydroid.12.3.unordered)
botryllus.gam.unordered<-summary(gam.beta.alive.bot.12.3.unordered)
caprellid.gam.unordered<-summary(gam.beta.caprellid.percent.12.unordered)
formicula.gam.unordered<-summary(gam.beta.formicula.12.unordered)
alive.mem.gam.unordered<-summary(gam.beta.alive.mem.12.unordered)
didemnum.gam.unordered<-summary(gam.beta.didemnum.12.unordered)
mussel_complete.gam.unordered<-summary(gam.nb.mussel_complete.12.unordered)
alive.barn.gam.unordered<-summary(gam.nb.num.barn.alive.12.unordered)
disporella.gam.unordered<-summary(gam.nb.disporella.12.unordered)
schizo.gam.unordered<-summary(gam.nb.schizo.12.unordered)
num.nudi.gam.unordered<-summary(gam.poisson.num.nudi.12.unordered)
num.serpulid.gam.unordered<-summary(gam.nb.num.serpulid.12.1.unordered)
clam.gam.unordered<-summary(gam.poisson.clam.12.unordered)
corella.gam.unordered<-summary(gam.nb.num.corella.12.unordered)
orange_sponge.gam.unordered<-summary(gam.nb.orange_sponge.12.unordered)

hydroid.gam.p.table<-as.data.frame(hydroid.gam.unordered$p.table)
hydroid.gam.s.table<-as.data.frame(hydroid.gam$s.table)

botryllus.gam.p.table<-as.data.frame(botryllus.gam.unordered$p.table)
botryllus.gam.s.table<-as.data.frame(botryllus.gam$s.table)

caprellid.gam.p.table<-as.data.frame(caprellid.gam.unordered$p.table)
caprellid.gam.s.table<-as.data.frame(caprellid.gam$s.table)

formicula.gam.p.table<-as.data.frame(formicula.gam.unordered$p.table)
formicula.gam.s.table<-as.data.frame(formicula.gam$s.table)

alive.mem.gam.p.table<-as.data.frame(alive.mem.gam.unordered$p.table)
alive.mem.gam.s.table<-as.data.frame(alive.mem.gam$s.table)

didemnum.gam.p.table<-as.data.frame(didemnum.gam.unordered$p.table)
didemnum.gam.s.table<-as.data.frame(didemnum.gam$s.table)

mussel_complete.gam.p.table<-as.data.frame(mussel_complete.gam.unordered$p.table)
mussel_complete.gam.s.table<-as.data.frame(mussel_complete.gam$s.table)

alive.barn.gam.p.table<-as.data.frame(alive.barn.gam.unordered$p.table)
alive.barn.gam.s.table<-as.data.frame(alive.barn.gam$s.table)

disporella.gam.p.table<-as.data.frame(disporella.gam.unordered$p.table)
disporella.gam.s.table<-as.data.frame(disporella.gam$s.table)

schizo.gam.p.table<-as.data.frame(schizo.gam.unordered$p.table)
schizo.gam.s.table<-as.data.frame(schizo.gam$s.table)

num.nudi.gam.p.table<-as.data.frame(num.nudi.gam.unordered$p.table)
num.nudi.gam.s.table<-as.data.frame(num.nudi.gam$s.table)

num.serpulid.gam.p.table<-as.data.frame(num.serpulid.gam.unordered$p.table)
num.serpulid.gam.s.table<-as.data.frame(num.serpulid.gam$s.table)

orange_sponge.gam.p.table<-as.data.frame(orange_sponge.gam.unordered$p.table)
orange_sponge.gam.s.table<-as.data.frame(orange_sponge.gam$s.table)

corella.gam.p.table<-as.data.frame(corella.gam.unordered$p.table)
corella.gam.s.table<-as.data.frame(corella.gam$s.table)

clam.gam.p.table<-as.data.frame(clam.gam.unordered$p.table)
clam.gam.s.table<-as.data.frame(clam.gam$s.table)



head(clam.gam.p.table)
#### Building the stats table
ptable<-rbind(hydroid.gam.p.table, 
               botryllus.gam.p.table, 
               caprellid.gam.p.table,
               formicula.gam.p.table,
               alive.mem.gam.p.table,
               didemnum.gam.p.table,
               mussel_complete.gam.p.table,
               alive.barn.gam.p.table,
               disporella.gam.p.table,
               schizo.gam.p.table,
               num.nudi.gam.p.table,
               num.serpulid.gam.p.table,
               orange_sponge.gam.p.table,
               corella.gam.p.table,
               clam.gam.p.table)


colnames(ptable) <- c("Estimate", "SE", "z", "p")
ptable$Factor<-rep(c("Intercept", "Low quality food", "High quality food"))



#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia", 1, 3) %>%
  group_rows("Botryllus", 4, 6) %>% 
  group_rows("Caprella",7, 9) %>% 
  group_rows("Folliculina", 10, 12) %>% 
  group_rows("Membranipora", 13,15) %>% 
  group_rows("Didemnum", 16, 18) %>% 
  group_rows("Mussels", 19, 21) %>% 
  group_rows("Barnacles", 22, 24) %>% 
  group_rows("Disporella", 25, 27) %>% 
  group_rows("Schizoporella", 28, 30) %>% 
  group_rows("Hermissenda", 31, 33) %>% 
  group_rows("Serpulid", 34, 36) %>% 
  group_rows("Sponge", 37, 39) %>% 
  group_rows("Corella", 40, 42) %>% 
  group_rows("Clams", 43, 45) %>% 
save_kable(file = "C:Data//For submission//ptable.html", self_contained = T)


### s table
stable<-rbind(hydroid.gam.s.table, 
              botryllus.gam.s.table, 
              caprellid.gam.s.table,
              formicula.gam.s.table,
              alive.mem.gam.s.table,
              didemnum.gam.s.table,
              mussel_complete.gam.s.table,
              alive.barn.gam.s.table,
              disporella.gam.s.table,
              schizo.gam.s.table,
              num.nudi.gam.s.table,
              num.serpulid.gam.s.table,
              orange_sponge.gam.s.table,
              corella.gam.s.table,
              clam.gam.s.table)


colnames(stable) <- c("Estimated_df", "Reference_df", "Chi_squared", "p_smooth")
stable$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Low quality food", "smooth min.10.pH * High quality food"))

#stable$species<-rep(c("hydroid", 
#                      "botryllus", 
#                       "caprellid",
#                     "formicula",
#                     "alive.mem",
#                    "didemnum",
  #                    "mussel_complete",
  #                   "alive.barn",
  #                   "disporella",
  #                   "schizo",
  #                   "num.nudi",
  #                   "num.serpulid",
  #                   "orange_sponge",
  #                   "corella",
  #                   "clam"), each=3)


stable %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia", 1, 3) %>%
  group_rows("Botryllus", 4, 6) %>% 
  group_rows("Caprella",7, 9) %>% 
  group_rows("Folliculina", 10, 12) %>% 
  group_rows("Membranipora", 13,15) %>% 
  group_rows("Didemnum", 16, 18) %>% 
  group_rows("Mussels", 19, 21) %>% 
  group_rows("Barnacles", 22, 24) %>% 
  group_rows("Disporella", 25, 27) %>% 
  group_rows("Schizoporella", 28, 30) %>% 
  group_rows("Hermissenda", 31, 33) %>% 
  group_rows("Serpulid", 34, 36) %>% 
  group_rows("Sponge", 37, 39) %>% 
  group_rows("Corella", 40, 42) %>% 
  group_rows("Clams", 43, 45) %>% 
  save_kable(file = "C:Data//For submission//stable.html", self_contained = T)

  
pstable<-cbind(ptable, stable)

pstable %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth, Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia, beta", 1, 3) %>%
  group_rows("Botryllus, beta", 4, 6) %>% 
  group_rows("Caprella, beta",7, 9) %>% 
  group_rows("Folliculina, beta", 10, 12) %>% 
  group_rows("Membranipora, beta", 13,15) %>% 
  group_rows("Didemnum, beta", 16, 18) %>% 
  group_rows("Mussels, negative binomial", 19, 21) %>% 
  group_rows("Barnacles, negative binomial", 22, 24) %>% 
  group_rows("Disporella, negative binomial", 25, 27) %>% 
  group_rows("Schizoporella, negative binomial", 28, 30) %>% 
  group_rows("Hermissenda, negative binomial", 31, 33) %>% 
  group_rows("Serpulid, negative binomial", 34, 36) %>% 
  group_rows("Sponge, negative binomial", 37, 39) %>% 
  group_rows("Corella, negative binomial", 40, 42) %>% 
  group_rows("Clams, poisson", 43, 45) %>% 
  save_kable(file = "C:Data//For submission//For resubmission//RESUB2//First look//pstable.html", self_contained = T)


# Richness ----------------------------------------------------------------
gam.nb.richness.12.1<- gam(richness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.richness.12<- gam(richness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.richness.12.1, gam.poisson.richness.12)

plot(gam.poisson.richness.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.poisson.richness.12)
#okay but qq plot not the best on ends
#qq_plot(gam.poisson.richness.12, method = 'simulate')
#k_check(gam.poisson.richness.12)
summary(gam.poisson.richness.12)
#a few outside the QQ plot on both ends
#low p value for k - but NS and edf is not super close to k-index

gam.poisson.richness.12.unordered<- gam(richness ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")
fam.gam.richness <- family(gam.poisson.richness.12)
ilink.gam.richness<- fam.gam.richness$linkinv
ilink.gam.richness


mod.richness<-gam.poisson.richness.12
ndata.richness <- with(food.exp.data.12.2019_zscores, 
                       data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                       length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.richness <- add_column(ndata.richness, fit = predict(mod.richness, newdata = ndata.richness, type = 'response'))

predict(mod.richness, newdata = ndata.richness, type = 'response')
ndata.richness <- bind_cols(ndata.richness, setNames(as_tibble(predict(mod.richness, ndata.richness, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.richness <- mutate(ndata.richness,
                               fit_resp  = ilink.gam.richness(fit_link),
                               right_upr = ilink.gam.richness(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.richness(fit_link - (2 * se_link)))


ndata.richness$min.10.pH.unscaled<-ndata.richness$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.richness <- ggplot(ndata.richness, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = richness, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species richness")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.richness,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,20)
plt.richness
ggsave("C:Data//Graphs March 2020//richness_pred.png")







# Evenness ----------------------------------------------------------------

gam.lm.evenness.12<- gam(evenness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.evenness.12.1<- gam(evenness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

AICtab( gam.lm.evenness.12, gam.gamma.evenness.12.1)

plot(gam.lm.evenness.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.evenness.12)
#looks very good
#qq_plot(gam.lm.evenness.12, method = 'simulate')
#k_check(gam.lm.evenness.12)
summary(gam.lm.evenness.12)

gam.lm.evenness.12.unordered<- gam(evenness ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")

fam.gam.evenness <- family(gam.lm.evenness.12)
ilink.gam.evenness<- fam.gam.evenness$linkinv
mod.evenness<-gam.lm.evenness.12
ndata.evenness <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.evenness <- add_column(ndata.evenness, fit = predict(mod.evenness, newdata = ndata.evenness, type = 'response'))

predict(mod.evenness, newdata = ndata.evenness, type = 'response')
ndata.evenness <- bind_cols(ndata.evenness, setNames(as_tibble(predict(mod.evenness, ndata.evenness, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.evenness <- mutate(ndata.evenness,
                         fit_resp  = ilink.gam.evenness(fit_link),
                         right_upr = ilink.gam.evenness(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.evenness(fit_link - (2 * se_link)))


ndata.evenness$min.10.pH.unscaled<-ndata.evenness$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.evenness <- ggplot(ndata.evenness, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = evenness, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species evenness")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.evenness
ggsave("C:Data//Graphs March 2020//evenness_pred.png")


# Occupied space ----------------------------------------------------------

gam.lm.occupied.space.12<- gam(occupied.space ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.occupied.space.12.1<- gam(occupied.space*0.01 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.binomial.occupied.space.12<- gam(formula = cbind(occupied.space, 100-occupied.space)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

gam.beta.occupied.space.12<- gam(occupied.space.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.occupied.space.12.1<- gam(occupied.space.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.occupied.space.12.2<- gam(occupied.space.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.occupied.space.12.3<- gam(occupied.space.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")

AICtab(gam.lm.occupied.space.12, gam.beta.occupied.space.12.3, gam.gamma.occupied.space.12.1, gam.beta.occupied.space.12, gam.beta.occupied.space.12.1, gam.beta.occupied.space.12.2, gam.binomial.occupied.space.12)
#beta cauchit is best 

plot(gam.beta.occupied.space.12.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.occupied.space.12.3)
#a bit blocky
#qq_plot(gam.beta.occupied.space.12.3, method = 'simulate')
#k_check(gam.beta.occupied.space.12.3)
#k gettinga bit low but ns
summary(gam.beta.occupied.space.12.3)
gam.beta.occupied.space.12.3.unordered<- gam(occupied.space.001~ s(min.10.pH, k=15)+ Food.quality + s(min.10.pH, by=oFood.quality, k=15), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.beta.occupied.space.12.3.unordered)

fam.gam.occupied.space <- family(gam.beta.occupied.space.12.3)
fam.gam.occupied.space
str(fam.gam.occupied.space)
ilink.gam.occupied.space<- fam.gam.occupied.space$linkinv
ilink.gam.occupied.space


mod.occupied.space<-gam.beta.occupied.space.12.3
ndata.occupied.space <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                 length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.occupied.space <- add_column(ndata.occupied.space, fit = predict(mod.occupied.space, newdata = ndata.occupied.space, type = 'response'))

predict(mod.occupied.space, newdata = ndata.occupied.space, type = 'response')
ndata.occupied.space <- bind_cols(ndata.occupied.space, setNames(as_tibble(predict(mod.occupied.space, ndata.occupied.space, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.occupied.space <- mutate(ndata.occupied.space,
                         fit_resp  = ilink.gam.occupied.space(fit_link),
                         right_upr = ilink.gam.occupied.space(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.occupied.space(fit_link - (2 * se_link)))


ndata.occupied.space$min.10.pH.unscaled<-ndata.occupied.space$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.occupied.space <- ggplot(ndata.occupied.space, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = occupied.space.001, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.occupied.space
ggsave("C:Data//Graphs March 2020//occupied.space_pred.png")



# Everything wet weight ---------------------------------------------------

gam.lm.everything.wet.weight.12<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.everything.wet.weight.12<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.everything.wet.weight.12.1, gam.lm.log.everything.wet.weight.12, gam.tweedie.everything.wet.weight.12.1,  gam.lm.everything.wet.weight.12, gam.gamma.everything.wet.weight.12.1)

#gam.lm.log.everything.wet.weight.12 is best 


plot(gam.lm.log.everything.wet.weight.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.log.everything.wet.weight.12)
#Looks quite good
#qq_plot(gam.lm.log.everything.wet.weight.12, method = 'simulate')
#k_check(gam.lm.log.everything.wet.weight.12)
summary(gam.lm.log.everything.wet.weight.12)

gam.lm.log.everything.wet.weight.12.unordered<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")

fam.gam.everything.wet.weight <- family(gam.lm.log.everything.wet.weight.12)
fam.gam.everything.wet.weight
ilink.gam.everything.wet.weight<- fam.gam.everything.wet.weight$linkinv
ilink.gam.everything.wet.weight


mod.everything.wet.weight<-gam.lm.log.everything.wet.weight.12
ndata.everything.wet.weight <- with(food.exp.data.12.2019_zscores, 
                                    data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                    length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.everything.wet.weight <- add_column(ndata.everything.wet.weight, fit = predict(mod.everything.wet.weight, newdata = ndata.everything.wet.weight, type = 'response'))

predict(mod.everything.wet.weight, newdata = ndata.everything.wet.weight, type = 'response')
ndata.everything.wet.weight <- bind_cols(ndata.everything.wet.weight, setNames(as_tibble(predict(mod.everything.wet.weight, ndata.everything.wet.weight, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.everything.wet.weight <- mutate(ndata.everything.wet.weight,
                         fit_resp  = ilink.gam.everything.wet.weight(fit_link),
                         right_upr = ilink.gam.everything.wet.weight(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.everything.wet.weight(fit_link - (2 * se_link)))


ndata.everything.wet.weight$min.10.pH.unscaled<-ndata.everything.wet.weight$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.everything.wet.weight <- ggplot(ndata.everything.wet.weight, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = log(everything.wet.weight), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total wet biomass per mesocosm\n(g, log scale)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.everything.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.everything.wet.weight
ggsave("C:Data//Graphs March 2020//everything.wet.weight_pred.png")


# Everything wet weight per 1 ---------------------------------------------

gam.lm.everything.wet.weight.per.1.12<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.everything.wet.weight.per.1.12<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.everything.wet.weight.per.1.12.1, gam.lm.log.everything.wet.weight.per.1.12, gam.tweedie.everything.wet.weight.per.1.12.1,gam.lm.everything.wet.weight.per.1.12, gam.gamma.everything.wet.weight.per.1.12.1)

#gam.lm.log.everything.wet.weight.per.1.12 is best by far


plot(gam.lm.log.everything.wet.weight.per.1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.log.everything.wet.weight.per.1.12)
#good but mabe a bit funnelly
#qq_plot(gam.lm.log.everything.wet.weight.per.1.12, method = 'simulate')
#k_check(gam.lm.log.everything.wet.weight.per.1.12)
summary(gam.lm.log.everything.wet.weight.per.1.12)

gam.lm.log.everything.wet.weight.per.1.12.unordered<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")

fam.gam.everything.wet.weight.per.1 <- family(gam.lm.log.everything.wet.weight.per.1.12)
fam.gam.everything.wet.weight.per.1
ilink.gam.everything.wet.weight.per.1<- fam.gam.everything.wet.weight.per.1$linkinv
ilink.gam.everything.wet.weight.per.1


mod.everything.wet.weight.per.1<-gam.lm.log.everything.wet.weight.per.1.12
ndata.everything.wet.weight.per.1 <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.everything.wet.weight.per.1 <- add_column(ndata.everything.wet.weight.per.1, fit = predict(mod.everything.wet.weight.per.1, newdata = ndata.everything.wet.weight.per.1, type = 'response'))

predict(mod.everything.wet.weight.per.1, newdata = ndata.everything.wet.weight.per.1, type = 'response')
ndata.everything.wet.weight.per.1 <- bind_cols(ndata.everything.wet.weight.per.1, setNames(as_tibble(predict(mod.everything.wet.weight.per.1, ndata.everything.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.everything.wet.weight.per.1 <- mutate(ndata.everything.wet.weight.per.1,
                                      fit_resp  = ilink.gam.everything.wet.weight.per.1(fit_link),
                                      right_upr = ilink.gam.everything.wet.weight.per.1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.everything.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.everything.wet.weight.per.1$min.10.pH.unscaled<-ndata.everything.wet.weight.per.1$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.everything.wet.weight.per.1 <- ggplot(ndata.everything.wet.weight.per.1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = log(everything.wet.weight.per.1), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Log Biomass per 1 % cover \n(wet weight)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.everything.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.everything.wet.weight.per.1
ggsave("C:Data//Graphs March 2020//everything.wet.weight.per.1_pred.png")




# Dry weight --------------------------------------------------------------
gam.lm.total_dry_biomass.12<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.total_dry_biomass.12<- gam(log(total_dry_biomass) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.total_dry_biomass.12.1, gam.lm.log.total_dry_biomass.12, gam.tweedie.total_dry_biomass.12.1, gam.lm.total_dry_biomass.12, gam.gamma.total_dry_biomass.12.1)

#gam.lm.loglink is best but lm is only 0.1 so going with simplest.... 


plot(gam.lm.total_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.total_dry_biomass.12 )
#Looks v good
#qq_plot(gam.lm.total_dry_biomass.12 , method = 'simulate')
#k_check(gam.lm.total_dry_biomass.12 )
summary(gam.lm.total_dry_biomass.12 )

gam.lm.total_dry_biomass.12.unordered<- gam(total_dry_biomass ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
summary(gam.lm.total_dry_biomass.12.unordered)

fam.gam.total_dry_biomass <- family(gam.lm.total_dry_biomass.12)
fam.gam.total_dry_biomass
str(fam.gam.total_dry_biomass)
ilink.gam.total_dry_biomass<- fam.gam.total_dry_biomass$linkinv
ilink.gam.total_dry_biomass


mod.total_dry_biomass<-gam.lm.total_dry_biomass.12
ndata.total_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.total_dry_biomass <- add_column(ndata.total_dry_biomass, fit = predict(mod.total_dry_biomass, newdata = ndata.total_dry_biomass, type = 'response'))

predict(mod.total_dry_biomass, newdata = ndata.total_dry_biomass, type = 'response')
ndata.total_dry_biomass <- bind_cols(ndata.total_dry_biomass, setNames(as_tibble(predict(mod.total_dry_biomass, ndata.total_dry_biomass, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.total_dry_biomass <- mutate(ndata.total_dry_biomass,
                                      fit_resp  = ilink.gam.total_dry_biomass(fit_link),
                                      right_upr = ilink.gam.total_dry_biomass(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.total_dry_biomass(fit_link - (2 * se_link)))


ndata.total_dry_biomass$min.10.pH.unscaled<-ndata.total_dry_biomass$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.total_dry_biomass <- ggplot(ndata.total_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (total_dry_biomass), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.total_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.total_dry_biomass
ggplot2::ggsave("C:Data//Graphs March 2020//total_dry_biomass_pred.png")


# Dry weight per 1 % cover --------------------------------------------------------------

gam.lm.total_dry_biomass_per1.12<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.total_dry_biomass_per1.12.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.total_dry_biomass_per1.12<- gam(log(total_dry_biomass_per1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.total_dry_biomass_per1.12<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.total_dry_biomass_per1.12.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.total_dry_biomass_per1.12.1, gam.lm.log.total_dry_biomass_per1.12, gam.tweedie.total_dry_biomass_per1.12, gam.lm.total_dry_biomass_per1.12, gam.gamma.total_dry_biomass_per1.12.1)

#tweedie is best by 1.7

plot(gam.tweedie.total_dry_biomass_per1.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.tweedie.total_dry_biomass_per1.12 )
#not available for tweedie
#qq_plot(gam.tweedie.total_dry_biomass_per1.12 , method = 'simulate')
#looks good
#k_check(gam.tweedie.total_dry_biomass_per1.12 )
summary(gam.tweedie.total_dry_biomass_per1.12 )

gam.tweedie.total_dry_biomass_per1.12.unordered<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=Food.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
summary(gam.tweedie.total_dry_biomass_per1.12.unordered)

fam.gam.total_dry_biomass_per1 <- family(gam.tweedie.total_dry_biomass_per1.12)
fam.gam.total_dry_biomass_per1
str(fam.gam.total_dry_biomass_per1)
ilink.gam.total_dry_biomass_per1<- fam.gam.total_dry_biomass_per1$linkinv
ilink.gam.total_dry_biomass_per1


mod.total_dry_biomass_per1<-gam.tweedie.total_dry_biomass_per1.12
ndata.total_dry_biomass_per1 <- with(food.exp.data.12.2019_zscores, 
                                     data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                     length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.total_dry_biomass_per1 <- add_column(ndata.total_dry_biomass_per1, fit = predict(mod.total_dry_biomass_per1, newdata = ndata.total_dry_biomass_per1, type = 'response'))

predict(mod.total_dry_biomass_per1, newdata = ndata.total_dry_biomass_per1, type = 'response')
ndata.total_dry_biomass_per1 <- bind_cols(ndata.total_dry_biomass_per1, setNames(as_tibble(predict(mod.total_dry_biomass_per1, ndata.total_dry_biomass_per1, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.total_dry_biomass_per1 <- mutate(ndata.total_dry_biomass_per1,
                                  fit_resp  = ilink.gam.total_dry_biomass_per1(fit_link),
                                  right_upr = ilink.gam.total_dry_biomass_per1(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.total_dry_biomass_per1(fit_link - (2 * se_link)))


ndata.total_dry_biomass_per1$min.10.pH.unscaled<-ndata.total_dry_biomass_per1$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.total_dry_biomass_per1 <- ggplot(ndata.total_dry_biomass_per1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (total_dry_biomass_per1), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.total_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.total_dry_biomass_per1
ggplot2::ggsave("C:Data//Graphs March 2020//total_dry_biomass_per1_pred.png")

# hydroid biomass ---------------------------------------------------------

gam.lm.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.hydroid_dry_biomass.12<- gam(log(hydroid_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.loglink.hydroid_dry_biomass.12, gam.lm.log.hydroid_dry_biomass.12, gam.tweedie.hydroid_dry_biomass.12,  gam.lm.hydroid_dry_biomass.12, gam.gamma.hydroid_dry_biomass.12)

#gamma is best 
plot(gam.gamma.hydroid_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.hydroid_dry_biomass.12 )
#looks good, maybe slightly funnelly
#qq_plot(gam.gamma.hydroid_dry_biomass.12 , method = 'simulate')
#k_check(gam.gamma.hydroid_dry_biomass.12 )
summary(gam.gamma.hydroid_dry_biomass.12 )

gam.gamma.hydroid_dry_biomass.12.unordered<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
fam.gam.hydroid_dry_biomass <- family(gam.gamma.hydroid_dry_biomass.12)
fam.gam.hydroid_dry_biomass
str(fam.gam.hydroid_dry_biomass)
ilink.gam.hydroid_dry_biomass<- fam.gam.hydroid_dry_biomass$linkinv
ilink.gam.hydroid_dry_biomass


mod.hydroid_dry_biomass<-gam.gamma.hydroid_dry_biomass.12
ndata.hydroid_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydroid_dry_biomass <- add_column(ndata.hydroid_dry_biomass, fit = predict(mod.hydroid_dry_biomass, newdata = ndata.hydroid_dry_biomass, type = 'response'))

predict(mod.hydroid_dry_biomass, newdata = ndata.hydroid_dry_biomass, type = 'response')
ndata.hydroid_dry_biomass <- bind_cols(ndata.hydroid_dry_biomass, setNames(as_tibble(predict(mod.hydroid_dry_biomass, ndata.hydroid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.hydroid_dry_biomass <- mutate(ndata.hydroid_dry_biomass,
                                  fit_resp  = ilink.gam.hydroid_dry_biomass(fit_link),
                                  right_upr = ilink.gam.hydroid_dry_biomass(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.hydroid_dry_biomass(fit_link - (2 * se_link)))


ndata.hydroid_dry_biomass$min.10.pH.unscaled<-ndata.hydroid_dry_biomass$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.hydroid_dry_biomass <- ggplot(ndata.hydroid_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (hydroid_dry_biomass+0.01), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Obelia") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.hydroid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.hydroid_dry_biomass
ggplot2::ggsave("C:Data//Graphs March 2020//hydroid_dry_biomass_pred.png")




# botryllus biomass -------------------------------------------------------
food.exp.data.12.2019_zscores$tunicate_dry_biomass[food.exp.data.12.2019_zscores$tunicate_dry_biomass<0]<-0

gam.lm.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.tunicate_dry_biomass.12<- gam(log(tunicate_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.tunicate_dry_biomass.12, gam.lm.log.tunicate_dry_biomass.12, gam.tweedie.tunicate_dry_biomass.12,gam.lm.tunicate_dry_biomass.12, gam.gamma.tunicate_dry_biomass.12)

#gamma is the best

plot(gam.gamma.tunicate_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.tunicate_dry_biomass.12 )
#a bit patterny
#qq_plot(gam.gamma.tunicate_dry_biomass.12 , method = 'simulate')
#k_check(gam.gamma.tunicate_dry_biomass.12 )
summary(gam.gamma.tunicate_dry_biomass.12 )

gam.gamma.tunicate_dry_biomass.12.unordered<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
fam.gam.tunicate_dry_biomass <- family(gam.gamma.tunicate_dry_biomass.12)
fam.gam.tunicate_dry_biomass
str(fam.gam.tunicate_dry_biomass)
ilink.gam.tunicate_dry_biomass<- fam.gam.tunicate_dry_biomass$linkinv
ilink.gam.tunicate_dry_biomass


mod.tunicate_dry_biomass<-gam.gamma.tunicate_dry_biomass.12
ndata.tunicate_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.tunicate_dry_biomass <- add_column(ndata.tunicate_dry_biomass, fit = predict(mod.tunicate_dry_biomass, newdata = ndata.tunicate_dry_biomass, type = 'response'))

predict(mod.tunicate_dry_biomass, newdata = ndata.tunicate_dry_biomass, type = 'response')
ndata.tunicate_dry_biomass <- bind_cols(ndata.tunicate_dry_biomass, setNames(as_tibble(predict(mod.tunicate_dry_biomass, ndata.tunicate_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.tunicate_dry_biomass <- mutate(ndata.tunicate_dry_biomass,
                                    fit_resp  = ilink.gam.tunicate_dry_biomass(fit_link),
                                    right_upr = ilink.gam.tunicate_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.tunicate_dry_biomass(fit_link - (2 * se_link)))


ndata.tunicate_dry_biomass$min.10.pH.unscaled<-ndata.tunicate_dry_biomass$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.tunicate_dry_biomass <- ggplot(ndata.tunicate_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (tunicate_dry_biomass+0.01), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.tunicate_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.tunicate_dry_biomass
ggplot2::ggsave("C:Data//Graphs March 2020//tunicate_dry_biomass_pred.png")



# caprellid biomass -------------------------------------------------------
gam.lm.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.caprellid_dry_biomass.12<- gam(log(caprellid_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.caprellid_dry_biomass.12, gam.lm.log.caprellid_dry_biomass.12, gam.tweedie.caprellid_dry_biomass.12, gam.lm.caprellid_dry_biomass.12, gam.gamma.caprellid_dry_biomass.12)
#gamma by 1.2

plot(gam.gamma.caprellid_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.caprellid_dry_biomass.12 )
# look good
#qq_plot(gam.gamma.caprellid_dry_biomass.12 , method = 'simulate')
#k_check(gam.gamma.caprellid_dry_biomass.12 )
summary(gam.gamma.caprellid_dry_biomass.12 )
gam.gamma.caprellid_dry_biomass.12.unordered<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.caprellid_dry_biomass <- family(gam.gamma.caprellid_dry_biomass.12)
ilink.gam.caprellid_dry_biomass<- fam.gam.caprellid_dry_biomass$linkinv
mod.caprellid_dry_biomass<-gam.gamma.caprellid_dry_biomass.12
ndata.caprellid_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.caprellid_dry_biomass <- add_column(ndata.caprellid_dry_biomass, fit = predict(mod.caprellid_dry_biomass, newdata = ndata.caprellid_dry_biomass, type = 'response'))

predict(mod.caprellid_dry_biomass, newdata = ndata.caprellid_dry_biomass, type = 'response')
ndata.caprellid_dry_biomass <- bind_cols(ndata.caprellid_dry_biomass, setNames(as_tibble(predict(mod.caprellid_dry_biomass, ndata.caprellid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid_dry_biomass <- mutate(ndata.caprellid_dry_biomass,
                                    fit_resp  = ilink.gam.caprellid_dry_biomass(fit_link),
                                    right_upr = ilink.gam.caprellid_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.caprellid_dry_biomass(fit_link - (2 * se_link)))


ndata.caprellid_dry_biomass$min.10.pH.unscaled<-ndata.caprellid_dry_biomass$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.caprellid_dry_biomass <- ggplot(ndata.caprellid_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (caprellid_dry_biomass+0.01), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.caprellid_dry_biomass
ggplot2::ggsave("C:Data//Graphs March 2020//caprellid_dry_biomass_pred.png")



# caprellid biomass per invididual -------------------------------------------------------

food.exp.data.12.2019_zscores$caprellid_dry_biomass_per1<-food.exp.data.12.2019_zscores$caprellid_dry_biomass/(food.caprellid.data_zscores$total.caprellids)

gam.lm.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.caprellid_dry_biomass_per1.12<- gam(log(caprellid_dry_biomass_per1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.caprellid_dry_biomass_per1.12, gam.lm.log.caprellid_dry_biomass_per1.12, gam.tweedie.caprellid_dry_biomass_per1.12,  gam.lm.caprellid_dry_biomass_per1.12, gam.gamma.caprellid_dry_biomass_per1.12)

#gamma is the best by 2.2
plot(gam.gamma.caprellid_dry_biomass_per1.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.caprellid_dry_biomass_per1.12 )
#a bit funnelly and qq right tail
#qq_plot(gam.gamma.caprellid_dry_biomass_per1.12 , method = 'simulate')
#k_check(gam.gamma.caprellid_dry_biomass_per1.12 )
summary(gam.gamma.caprellid_dry_biomass_per1.12 )
gam.gamma.caprellid_dry_biomass_per1.12.unordered<- gam(caprellid_dry_biomass_per1+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.caprellid_dry_biomass_per1 <- family(gam.gamma.caprellid_dry_biomass_per1.12)
fam.gam.caprellid_dry_biomass_per1
str(fam.gam.caprellid_dry_biomass_per1)
ilink.gam.caprellid_dry_biomass_per1<- fam.gam.caprellid_dry_biomass_per1$linkinv
ilink.gam.caprellid_dry_biomass_per1


mod.caprellid_dry_biomass_per1<-gam.gamma.caprellid_dry_biomass_per1.12
ndata.caprellid_dry_biomass_per1 <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.caprellid_dry_biomass_per1 <- add_column(ndata.caprellid_dry_biomass_per1, fit = predict(mod.caprellid_dry_biomass_per1, newdata = ndata.caprellid_dry_biomass_per1, type = 'response'))

predict(mod.caprellid_dry_biomass_per1, newdata = ndata.caprellid_dry_biomass_per1, type = 'response')
ndata.caprellid_dry_biomass_per1 <- bind_cols(ndata.caprellid_dry_biomass_per1, setNames(as_tibble(predict(mod.caprellid_dry_biomass_per1, ndata.caprellid_dry_biomass_per1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid_dry_biomass_per1 <- mutate(ndata.caprellid_dry_biomass_per1,
                                      fit_resp  = ilink.gam.caprellid_dry_biomass_per1(fit_link),
                                      right_upr = ilink.gam.caprellid_dry_biomass_per1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.caprellid_dry_biomass_per1(fit_link - (2 * se_link)))


ndata.caprellid_dry_biomass_per1$min.10.pH.unscaled<-ndata.caprellid_dry_biomass_per1$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.caprellid_dry_biomass_per1 <- ggplot(ndata.caprellid_dry_biomass_per1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (caprellid_dry_biomass_per1), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.caprellid_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')+ylim(0,0.012)
plt.caprellid_dry_biomass_per1
ggplot2::ggsave("C:Data//Graphs March 2020//caprellid_dry_biomass_per1_pred.png")


# rest biomass ------------------------------------------------------------

#kcheck was significant so increase k from 10 to 11
gam.lm.rest_dry_biomass.12<- gam(rest_dry_biomass ~ s(min.10.pH, k=12)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=12),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=12),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.rest_dry_biomass.12<- gam(log(rest_dry_biomass+0.1) ~ s(min.10.pH, k=12)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=12),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=12),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH, k=12)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=12),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.rest_dry_biomass.12, gam.lm.log.rest_dry_biomass.12, gam.tweedie.rest_dry_biomass.12, gam.lm.rest_dry_biomass.12, gam.gamma.rest_dry_biomass.12)

#tweedie the best 


plot(gam.tweedie.rest_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#qq_plot(gam.tweedie.rest_dry_biomass.12 , method = 'simulate')
#k_check(gam.tweedie.rest_dry_biomass.12 )
summary(gam.tweedie.rest_dry_biomass.12)

gam.tweedie.rest_dry_biomass.12.unordered<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
summary(gam.tweedie.rest_dry_biomass.12.unordered)

fam.gam.rest_dry_biomass <- family(gam.tweedie.rest_dry_biomass.12)
fam.gam.rest_dry_biomass
str(fam.gam.rest_dry_biomass)
ilink.gam.rest_dry_biomass<- fam.gam.rest_dry_biomass$linkinv
ilink.gam.rest_dry_biomass


mod.rest_dry_biomass<-gam.tweedie.rest_dry_biomass.12
ndata.rest_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                            length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))
## add the fitted values by predicting from the model for the new data
ndata.rest_dry_biomass <- add_column(ndata.rest_dry_biomass, fit = predict(mod.rest_dry_biomass, newdata = ndata.rest_dry_biomass, type = 'response'))

predict(mod.rest_dry_biomass, newdata = ndata.rest_dry_biomass, type = 'response')
ndata.rest_dry_biomass <- bind_cols(ndata.rest_dry_biomass, setNames(as_tibble(predict(mod.rest_dry_biomass, ndata.rest_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.rest_dry_biomass <- mutate(ndata.rest_dry_biomass,
                                    fit_resp  = ilink.gam.rest_dry_biomass(fit_link),
                                    right_upr = ilink.gam.rest_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.rest_dry_biomass(fit_link - (2 * se_link)))


ndata.rest_dry_biomass$min.10.pH.unscaled<-ndata.rest_dry_biomass$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.rest_dry_biomass <- ggplot(ndata.rest_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = (rest_dry_biomass+0.01), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression("Remaining dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.rest_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.rest_dry_biomass
ggplot2::ggsave("C:Data//Graphs March 2020//rest_dry_biomass_pred.png")




# Mussel wet weight -------------------------------------------------------

gam.lm.Mussel.wet.weight.12<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.Mussel.wet.weight.12<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.Mussel.wet.weight.12.1, gam.lm.log.Mussel.wet.weight.12, gam.tweedie.Mussel.wet.weight.12.1,gam.lm.Mussel.wet.weight.12, gam.gamma.Mussel.wet.weight.12.1)

#gam.lm.log.Mussel.wet.weight.12 is best 

plot(gam.lm.log.Mussel.wet.weight.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.log.Mussel.wet.weight.12)
#look pretty good, esp qq
#qq_plot(gam.lm.log.Mussel.wet.weight.12, method = 'simulate')
#k_check(gam.lm.log.Mussel.wet.weight.12)
summary(gam.lm.log.Mussel.wet.weight.12)

#residuals many in a straight line but otherwise good 

gam.lm.log.Mussel.wet.weight.12.unordered<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
fam.gam.Mussel.wet.weight <- family(gam.lm.log.Mussel.wet.weight.12)
fam.gam.Mussel.wet.weight
ilink.gam.Mussel.wet.weight<- fam.gam.Mussel.wet.weight$linkinv
ilink.gam.Mussel.wet.weight


mod.Mussel.wet.weight<-gam.lm.log.Mussel.wet.weight.12
ndata.Mussel.wet.weight <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.Mussel.wet.weight <- add_column(ndata.Mussel.wet.weight, fit = predict(mod.Mussel.wet.weight, newdata = ndata.Mussel.wet.weight, type = 'response'))

predict(mod.Mussel.wet.weight, newdata = ndata.Mussel.wet.weight, type = 'response')
ndata.Mussel.wet.weight <- bind_cols(ndata.Mussel.wet.weight, setNames(as_tibble(predict(mod.Mussel.wet.weight, ndata.Mussel.wet.weight, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.Mussel.wet.weight <- mutate(ndata.Mussel.wet.weight,
                                      fit_resp  = ilink.gam.Mussel.wet.weight(fit_link),
                                      right_upr = ilink.gam.Mussel.wet.weight(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.Mussel.wet.weight(fit_link - (2 * se_link)))


ndata.Mussel.wet.weight$min.10.pH.unscaled<-ndata.Mussel.wet.weight$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.Mussel.wet.weight <- ggplot(ndata.Mussel.wet.weight, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = log(Mussel.wet.weight+0.1), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Mytilus") ~"wet weight (g, log scale)"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.Mussel.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.Mussel.wet.weight
ggsave("C:Data//Graphs March 2020//Mussel.wet.weight_pred.png")






# Mussel wet weight per individual mussel ---------------------------------
#big outlier .... 5 mussels weigh 105.5 g?? so 17 g per mussel? 

food.exp.data.12.2019_zscores$Mussel.wet.weight[food.exp.data.12.2019_zscores$Mesocosm==8]<-14.0
food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1<-(food.exp.data.12.2019_zscores$Mussel.wet.weight)/(food.exp.data.12.2019_zscores$mussel_complete+1)

gam.lm.Mussel.wet.weight.per.1.12<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.Mussel.wet.weight.per.1.12<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.Mussel.wet.weight.per.1.12<- gam(log(Mussel.wet.weight.per.1+0.01) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.Mussel.wet.weight.per.1.12.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.Mussel.wet.weight.per.1.12.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")

AICtab(gam.loglink.Mussel.wet.weight.per.1.12.1, gam.lm.log.Mussel.wet.weight.per.1.12, gam.tweedie.Mussel.wet.weight.per.1.12.1,  gam.lm.Mussel.wet.weight.per.1.12, gam.gamma.Mussel.wet.weight.per.1.12)

#gamma is best 
plot(gam.gamma.Mussel.wet.weight.per.1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.gamma.Mussel.wet.weight.per.1.12)
# good but a bit funnelly
#qq_plot(gam.gamma.Mussel.wet.weight.per.1.12, method = 'simulate')
#k_check(gam.gamma.Mussel.wet.weight.per.1.12)
summary(gam.gamma.Mussel.wet.weight.per.1.12)

gam.gamma.Mussel.wet.weight.per.1.12.unordered<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.Mussel.wet.weight.per.1 <- family(gam.gamma.Mussel.wet.weight.per.1.12)
fam.gam.Mussel.wet.weight.per.1
str(fam.gam.Mussel.wet.weight.per.1)
ilink.gam.Mussel.wet.weight.per.1<- fam.gam.Mussel.wet.weight.per.1$linkinv
ilink.gam.Mussel.wet.weight.per.1


mod.Mussel.wet.weight.per.1<-gam.gamma.Mussel.wet.weight.per.1.12
ndata.Mussel.wet.weight.per.1 <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.Mussel.wet.weight.per.1 <- add_column(ndata.Mussel.wet.weight.per.1, fit = predict(mod.Mussel.wet.weight.per.1, newdata = ndata.Mussel.wet.weight.per.1, type = 'response'))

predict(mod.Mussel.wet.weight.per.1, newdata = ndata.Mussel.wet.weight.per.1, type = 'response')
ndata.Mussel.wet.weight.per.1 <- bind_cols(ndata.Mussel.wet.weight.per.1, setNames(as_tibble(predict(mod.Mussel.wet.weight.per.1, ndata.Mussel.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.Mussel.wet.weight.per.1 <- mutate(ndata.Mussel.wet.weight.per.1,
                                  fit_resp  = ilink.gam.Mussel.wet.weight.per.1(fit_link),
                                  right_upr = ilink.gam.Mussel.wet.weight.per.1(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.Mussel.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.Mussel.wet.weight.per.1$min.10.pH.unscaled<-ndata.Mussel.wet.weight.per.1$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.Mussel.wet.weight.per.1 <- ggplot(ndata.Mussel.wet.weight.per.1, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = Mussel.wet.weight.per.1, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Biomass per mussel \n(wet weight)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.Mussel.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.Mussel.wet.weight.per.1
ggsave("C:Data//Graphs March 2020//Mussel.wet.weight.per.1_pred.png")




# Hydtobot ----------------------------------------------------------------
food.exp.data.12.2019_zscores$hydtobot[food.exp.data.12.2019_zscores$hydtobot==1]<-0.99
food.exp.data.12.2019_zscores$hydtobot[food.exp.data.12.2019_zscores$hydtobot==0]<-0.01

gam.lm.hydtobot.12<- gam(log(hydtobot+0.01)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.log.lm.hydtobot.12<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.tweedie.hydtobot.12<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family =tw(), select=TRUE, method="REML")

gam.beta.hydtobot.12<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydtobot.12.1<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydtobot.12.2<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")


AICtab( gam.tweedie.hydtobot.12 ,gam.beta.hydtobot.12, gam.lm.hydtobot.12, gam.log.lm.hydtobot.12, gam.beta.hydtobot.12.1, gam.beta.hydtobot.12.2)
### beta logit is the best 

plot(gam.beta.hydtobot.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.hydtobot.12)
#qq a bit wiggly
#qq_plot(gam.beta.hydtobot.12, method = 'simulate')
#k_check(gam.beta.hydtobot.12)
gam.check(gam.tweedie.hydtobot.12)
summary(gam.beta.hydtobot.12)
#not the best qq plot
#resids are in rows

gam.beta.hydtobot.12.unordered<- gam(hydtobot~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
fam.gam.hydtobot <- family(gam.beta.hydtobot.12)
fam.gam.hydtobot 

ilink.gam.hydtobot <- fam.gam.hydtobot$linkinv
ilink.gam.hydtobot


food.exp.data.12.2019_zscores$min.10.pH.unscaled <-food.exp.data.12.2019_zscores$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')
head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

mod.hydtobot<-gam.beta.hydtobot.12
ndata.hydtobot <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                                               length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydtobot <- add_column(ndata.hydtobot, fit = predict(mod.hydtobot, newdata = ndata.hydtobot, type = 'response'))


ndata.hydtobot <- bind_cols(ndata.hydtobot, setNames(as_tibble(predict(mod.hydtobot, ndata.hydtobot, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.hydtobot <- mutate(ndata.hydtobot,
                        fit_resp  = ilink.gam.hydtobot(fit_link),
                        right_upr = ilink.gam.hydtobot(fit_link + (2 * se_link)),
                        right_lwr = ilink.gam.hydtobot(fit_link - (2 * se_link)))




ndata.hydtobot$min.10.pH.unscaled<-ndata.hydtobot$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.gam.hydtobot <- ggplot(ndata.hydtobot, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = hydtobot, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "cover ratio"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.hydtobot,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.2, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=6), legend.title = element_text(size=7))+ 
  geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)+coord_cartesian(ylim = c(0, 1)) 
plt.gam.hydtobot 


# Hydtobot by weight ------------------------------------------------------
food.exp.data.12.2019_zscores$hydtobot_dry_biomass[food.exp.data.12.2019_zscores$hydtobot_dry_biomass==1]<-0.99
food.exp.data.12.2019_zscores$hydtobot_dry_biomass[food.exp.data.12.2019_zscores$hydtobot_dry_biomass==0]<-0.01

gam.lm.hydtobot_dry_biomass.12<- gam(log(hydtobot_dry_biomass+0.01)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.log.lm.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.tweedie.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family =tw(), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.12.1<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.12.2<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")

AICtab(gam.tweedie.hydtobot_dry_biomass.12 ,gam.beta.hydtobot_dry_biomass.12, gam.lm.hydtobot_dry_biomass.12, gam.log.lm.hydtobot_dry_biomass.12, gam.beta.hydtobot_dry_biomass.12.1, gam.beta.hydtobot_dry_biomass.12.2, lm.hydtobot_dry_biomass, lm.log.hydtobot_dry_biomass)
### logit is the best 

plot(gam.beta.hydtobot_dry_biomass.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.beta.hydtobot_dry_biomass.12)
#look ok, qq is bloacky
#qq_plot(gam.beta.hydtobot_dry_biomass.12, method = 'simulate')
#k_check(gam.beta.hydtobot_dry_biomass.12)
summary(gam.beta.hydtobot_dry_biomass.12)
gam.beta.hydtobot_dry_biomass.12.unordered<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")

fam.gam.hydtobot_dry_biomass <- family(gam.beta.hydtobot_dry_biomass.12)
fam.gam.hydtobot_dry_biomass 
str(fam.gam.hydtobot_dry_biomass )
ilink.gam.hydtobot_dry_biomass <- fam.gam.hydtobot_dry_biomass$linkinv
ilink.gam.hydtobot_dry_biomass

food.exp.data.12.2019_zscores$min.10.pH.unscaled <-food.exp.data.12.2019_zscores$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')
head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

mod.hydtobot_dry_biomass<-gam.beta.hydtobot_dry_biomass.12
ndata.hydtobot_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH= seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.hydtobot_dry_biomass <- add_column(ndata.hydtobot_dry_biomass, fit = predict(mod.hydtobot_dry_biomass, newdata = ndata.hydtobot_dry_biomass, type = 'response'))


ndata.hydtobot_dry_biomass <- bind_cols(ndata.hydtobot_dry_biomass, setNames(as_tibble(predict(mod.hydtobot_dry_biomass, ndata.hydtobot_dry_biomass, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.hydtobot_dry_biomass <- mutate(ndata.hydtobot_dry_biomass,
                         fit_resp  = ilink.gam.hydtobot_dry_biomass(fit_link),
                         right_upr = ilink.gam.hydtobot_dry_biomass(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.hydtobot_dry_biomass(fit_link - (2 * se_link)))




ndata.hydtobot_dry_biomass$min.10.pH.unscaled<-ndata.hydtobot_dry_biomass$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 

plt.gam.hydtobot_dry_biomass <- ggplot(ndata.hydtobot_dry_biomass, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y = hydtobot_dry_biomass, shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "biomass ratio"))+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.hydtobot_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)#+coord_cartesian(ylim = c(0, 1)) 
plt.gam.hydtobot_dry_biomass 


#mesocosm #8 didn't have either tunicates or hydroids weight - 0% hydroids, 2% tunicates


# CAP1 --------------------------------------------------------------------
#negative values so can't do gamma
gam.lm.CAP1.12<- gam(CAP1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.loglink.CAP1.12.1<- gam(CAP1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( gam.loglink.CAP1.12.1, gam.lm.CAP1.12)
#gam.lm.CAP1

plot(gam.lm.CAP1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.CAP1.12)
#look good
#qq_plot(gam.lm.CAP1.12, method = 'simulate')
#k_check(gam.lm.CAP1.12)
summary(gam.lm.CAP1.12)
gam.lm.CAP1.12.unordered<- gam(CAP1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
summary(gam.lm.CAP1.12.unordered)

fam.gam.CAP1 <- family(gam.lm.CAP1.12)
fam.gam.CAP1
str(fam.gam.CAP1)
ilink.gam.CAP1<- fam.gam.CAP1$linkinv
ilink.gam.CAP1


mod.CAP1<-gam.lm.CAP1.12
ndata.CAP1 <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.CAP1 <- add_column(ndata.CAP1, fit = predict(mod.CAP1, newdata = ndata.CAP1, type = 'response'))

predict(mod.CAP1, newdata = ndata.CAP1, type = 'response')
ndata.CAP1 <- bind_cols(ndata.CAP1, setNames(as_tibble(predict(mod.CAP1, ndata.CAP1, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.CAP1 <- mutate(ndata.CAP1,
                                      fit_resp  = ilink.gam.CAP1(fit_link),
                                      right_upr = ilink.gam.CAP1(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.CAP1(fit_link - (2 * se_link)))


ndata.CAP1$min.10.pH.unscaled<-ndata.CAP1$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.CAP1 <- ggplot(ndata.CAP1, aes(x = min.10.pH.unscaled, y = fit)) + 
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y =(CAP1), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Partial-dbRDA axis 1\n(36% of constrained variation)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='bottom', legend.box='horizontal', legend.spacing=unit(0.1, "cm"), legend.margin=margin(0, 0, 0, 0, "cm"), legend.key.size = unit(0, "cm"), legend.text = element_text(size=3), legend.title = element_text(size=4))
plt.CAP1
ggsave("C:Data//Graphs March 2020//CAP1_pred.png")




# Distances ---------------------------------------------------------------

#k check was significant so increased k from 10 to 11

gam.lm.distances.12<- gam(distcentroid ~ s(min.10.pH, k=11)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=11),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.loglink.distances.12.1<- gam(distcentroid~ s(min.10.pH, k=11)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=11),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")
gam.gamma.distances.12<- gam(distcentroid~ s(min.10.pH, k=11)+ oFood.quality + s(min.10.pH, by=oFood.quality, k=11),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

AICtab(gam.loglink.distances.12.1,  gam.lm.distances.12, gam.gamma.distances.12)
#gam.lm.distances although both are equal


plot(gam.lm.distances.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
#appraise(gam.lm.distances.12)
#looks good
#qq_plot(gam.lm.distances.12, method = 'simulate')
#k_check(gam.lm.distances.12)
#good now
summary(gam.lm.distances.12)

gam.lm.distances.12.unordered<- gam(distcentroid~ s(min.10.pH, k=11)+ Food.quality + s(min.10.pH, by=oFood.quality, k=11),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
summary(gam.lm.distances.12.unordered)

fam.gam.distances <- family(gam.lm.distances.12)
fam.gam.distances
str(fam.gam.distances)
ilink.gam.distances<- fam.gam.distances$linkinv
ilink.gam.distances


mod.distances<-gam.lm.distances.12
ndata.distances <- with(food.exp.data.12.2019_zscores, tibble(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                             length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.distances <- add_column(ndata.distances, fit = predict(mod.distances, newdata = ndata.distances, type = 'response'))

predict(mod.distances, newdata = ndata.distances, type = 'response')
ndata.distances <- bind_cols(ndata.distances, setNames(as_tibble(predict(mod.distances, ndata.distances, se.fit = TRUE)[1:2]),
                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.distances <- mutate(ndata.distances,
                     fit_resp  = ilink.gam.distances(fit_link),
                     right_upr = ilink.gam.distances(fit_link + (2 * se_link)),
                     right_lwr = ilink.gam.distances(fit_link - (2 * se_link)))


ndata.distances$min.10.pH.unscaled<-ndata.distances$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')


# plot 
plt.distances <- ggplot(ndata.distances, aes(x = min.10.pH.unscaled, y = fit)) + 
  
  geom_line(aes(colour=oFood.quality)) +
  geom_point(aes(y =(distcentroid), shape=CO2, colour=oFood.quality), data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Heterogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset2, guide = guide_legend(title="Food supplement", title.position = "top"))+
  scale_fill_manual(values=colorset2, guide = FALSE)+
  scale_shape_manual(values=c(19,17), labels=c("Ambient", "Low pH"), guide = guide_legend(title="pH treatment", title.position = "top"))+
  geom_ribbon(data = ndata.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position='none')
plt.distances
ggsave("C:Data//Graphs March 2020//distances_pred.png")

# Community plotting ------------------------------------------------------
library(cowplot)

#### revised community fig
fig.4.community<-wrap_plots( plt.occupied.space,plt.total_dry_biomass,
                            plt.richness, plt.evenness,
                            plt.CAP1, plt.distances, ncol=2)+
                            plot_annotation(tag_levels = 'a')

fig.4.community

ggplot2::ggsave(plot=fig.4.community, "C:Data//For submission//For resubmission//RESUB2//First look//Fig4.community.pdf", width=3, height=5, units="in")



#hyd tobot figs
fig.s3.hydtobot<-wrap_plots(plt.gam.hydtobot, plt.gam.hydtobot_dry_biomass, ncol=2) + plot_annotation(tag_levels = 'a')

fig.s3.hydtobot

ggplot2::ggsave("C:Data//For submission//For resubmission//RESUB2//First look//Fig.S3.hydotobot.png", width=6, height=3, units="in", dpi=600)



# Community level tables --------------------------------------------------

richness.gam<- summary(gam.poisson.richness.12)
richness.gam.unordered<- summary(gam.poisson.richness.12.unordered)

evenness.gam<-summary(gam.lm.evenness.12)
evenness.gam.unordered<-summary(gam.lm.evenness.12.unordered)


occupied.space.gam<- summary(gam.beta.occupied.space.12.3)
occupied.space.gam.unordered<- summary(gam.beta.occupied.space.12.3.unordered)

distances.gam<- summary(gam.lm.distances.12)
distances.gam.unordered<- summary(gam.lm.distances.12.unordered)


CAP1.gam <- summary(gam.lm.CAP1.12)
CAP1.gam.unordered <- summary(gam.lm.CAP1.12.unordered)

#dry biomass
total_dry_biomass.gam <- summary(gam.lm.total_dry_biomass.12)
total_dry_biomass.gam.unordered <- summary(gam.lm.total_dry_biomass.12.unordered)

hydroid_dry_biomass.gam.unordered <- summary(gam.gamma.hydroid_dry_biomass.12.unordered) 
hydroid_dry_biomass.gam<- summary(gam.gamma.hydroid_dry_biomass.12) 

caprellid_dry_biomass.gam.unordered <- summary(gam.gamma.caprellid_dry_biomass.12.unordered)
caprellid_dry_biomass.gam <- summary(gam.gamma.caprellid_dry_biomass.12)

tunicate_dry_biomass.gam <- summary(gam.gamma.tunicate_dry_biomass.12)
tunicate_dry_biomass.gam.unordered <- summary(gam.gamma.tunicate_dry_biomass.12.unordered)

rest_dry_biomass.gam.unordered <- summary(gam.tweedie.rest_dry_biomass.12.unordered)
rest_dry_biomass.gam <- summary(gam.tweedie.rest_dry_biomass.12)

hydtobot_dry_biomass.gam.unordered<-summary(gam.beta.hydtobot_dry_biomass.12.unordered)
hydtobot_dry_biomass.gam<-summary(gam.beta.hydtobot_dry_biomass.12)

#wet biomass
everything.wet.weight.gam <-summary(gam.lm.log.everything.wet.weight.12)
everything.wet.weight.gam.unordered <- summary(gam.lm.log.everything.wet.weight.12.unordered)

everything.wet.weight.per.1.gam <-summary(gam.lm.log.everything.wet.weight.per.1.12)
everything.wet.weight.per.1.gam.unordered <- summary(gam.lm.log.everything.wet.weight.per.1.12.unordered)


Mussel.wet.weight.gam <- summary(gam.lm.log.Mussel.wet.weight.12)
Mussel.wet.weight.gam.unordered <-summary(gam.lm.log.Mussel.wet.weight.12.unordered)

Mussel.wet.weight.per.1.gam <- summary(gam.gamma.Mussel.wet.weight.per.1.12)
Mussel.wet.weight.per.1.gam.unordered <-summary(gam.gamma.Mussel.wet.weight.per.1.12.unordered)


#competition metric 
hydtobot.gam <- summary(gam.beta.hydtobot.12)
hydtobot.gam.unordered <- summary(gam.beta.hydtobot.12.unordered)

#ptable building
richness.gam.p.table<-as.data.frame(richness.gam.unordered$p.table)
richness.gam.s.table<-as.data.frame(richness.gam$s.table)

evenness.gam.p.table<-as.data.frame(evenness.gam.unordered$p.table)
evenness.gam.s.table<-as.data.frame(evenness.gam$s.table)

occupied.space.gam.p.table<-as.data.frame(occupied.space.gam.unordered$p.table)
occupied.space.gam.s.table<-as.data.frame(occupied.space.gam$s.table)

distances.gam.p.table<-as.data.frame(distances.gam.unordered$p.table)
distances.gam.s.table<-as.data.frame(distances.gam$s.table)

CAP1.gam.p.table<-as.data.frame(CAP1.gam.unordered$p.table)
CAP1.gam.s.table<-as.data.frame(CAP1.gam$s.table)

total_dry_biomass.gam.p.table<-as.data.frame(total_dry_biomass.gam.unordered$p.table)
total_dry_biomass.gam.s.table<-as.data.frame(total_dry_biomass.gam$s.table)

hydroid_dry_biomass.gam.p.table<-as.data.frame(hydroid_dry_biomass.gam.unordered$p.table)
hydroid_dry_biomass.gam.s.table<-as.data.frame(hydroid_dry_biomass.gam$s.table)

caprellid_dry_biomass.gam.p.table<-as.data.frame(caprellid_dry_biomass.gam.unordered$p.table)
caprellid_dry_biomass.gam.s.table<-as.data.frame(caprellid_dry_biomass.gam$s.table)

tunicate_dry_biomass.gam.p.table<-as.data.frame(tunicate_dry_biomass.gam.unordered$p.table)
tunicate_dry_biomass.gam.s.table<-as.data.frame(tunicate_dry_biomass.gam$s.table)

rest_dry_biomass.gam.p.table<-as.data.frame(rest_dry_biomass.gam.unordered$p.table)
rest_dry_biomass.gam.s.table<-as.data.frame(rest_dry_biomass.gam$s.table)

everything.wet.weight.gam.p.table<-as.data.frame(everything.wet.weight.gam.unordered$p.table)
everything.wet.weight.gam.s.table<-as.data.frame(everything.wet.weight.gam$s.table)

everything.wet.weight.per.1.gam.p.table<-as.data.frame(everything.wet.weight.per.1.gam.unordered$p.table)
everything.wet.weight.per.1.gam.s.table<-as.data.frame(everything.wet.weight.per.1.gam$s.table)

Mussel.wet.weight.gam.p.table<-as.data.frame(Mussel.wet.weight.gam.unordered$p.table)
Mussel.wet.weight.gam.s.table<-as.data.frame(Mussel.wet.weight.gam$s.table)

Mussel.wet.weight.per.1.gam.p.table<-as.data.frame(Mussel.wet.weight.per.1.gam.unordered$p.table)
Mussel.wet.weight.per.1.gam.s.table<-as.data.frame(Mussel.wet.weight.per.1.gam$s.table)

hydtobot.gam.p.table<-as.data.frame(hydtobot.gam.unordered$p.table)
hydtobot.gam.s.table<-as.data.frame(hydtobot.gam$s.table)

hydtobot_dry_biomass.gam.p.table<-as.data.frame(hydtobot_dry_biomass.gam.unordered$p.table)
hydtobot_dry_biomass.gam.s.table<-as.data.frame(hydtobot_dry_biomass.gam$s.table)

hydtobot_dry_biomass.gam.p.table
hydtobot_dry_biomass.gam.s.table

#richness.gam.p.table and  hydtobot.gam.p.table, is with z value 
colnames(richness.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot_dry_biomass.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(occupied.space.gam.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")


#### Building the stats table
ptable.community.t<-rbind(richness.gam.p.table,
              evenness.gam.p.table,
              occupied.space.gam.p.table,
              total_dry_biomass.gam.p.table,
              CAP1.gam.p.table,
              distances.gam.p.table,
              hydtobot.gam.p.table, 
              hydtobot_dry_biomass.gam.p.table,
              everything.wet.weight.gam.p.table
              )


colnames(ptable.community.t) <- c("Estimate", "SE", "t", "p")
ptable.community.t$Factor<-rep(c("Intercept", "Low quality food", "High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.community.t %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (z)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (z)", 7,9) %>% 
  group_rows("Total dry biomass, normal", 10,12) %>% 
  group_rows("Partial dbRDA (1st axis), normal",13,15) %>% 
  group_rows("Heterogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (z)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (z)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:Data//For submission//ptable.community.t.html", self_contained = T)


#again hydtobot and richness
#richness.gam.p.table and  hydtobot.gam.p.table, is with Chisq
colnames(richness.gam.s.table) <- c("edf", "Ref.df", "F", "p-value")
colnames(hydtobot.gam.s.table) <- c("edf", "Ref.df",  "F", "p-value")
colnames(hydtobot_dry_biomass.gam.s.table) <- c("edf", "Ref.df",  "F", "p-value")

colnames(occupied.space.gam.s.table) <- c("edf", "Ref.df",  "F", "p-value")


### s table
stable.community.f<-rbind(richness.gam.s.table,
                          evenness.gam.s.table,
                          occupied.space.gam.s.table,
                          total_dry_biomass.gam.s.table,
                          CAP1.gam.s.table,
                          distances.gam.s.table,
                          hydtobot.gam.s.table, 
                          hydtobot_dry_biomass.gam.s.table,
                          everything.wet.weight.gam.s.table
)


colnames(stable.community.f) <- c("Estimated_df", "Reference_df", "F", "p_smooth")
stable.community.f$Smooth_terms<-rep(c("smooth min.10.pH", "smooth min.10.pH * Low quality food", "smooth min.10.pH * High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

stable.community.f %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth) %>% 
  kable(escape=F, digits=4) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (Chi-square)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (Chi-square)", 7,9) %>% 
  group_rows("Total dry biomass, normal", 10,12) %>% 
  group_rows("Partial dbRDA (1st axis), normal",13,15) %>% 
  group_rows("Heterogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (Chi-square)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi-square)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:Data//For submission//stable.community.f.html", self_contained = T)

pstable.community<-cbind(ptable.community.t, stable.community.f)


pstable.community %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.051, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth, Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (Chi-square, z)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (Chi-square, z)", 7,9) %>% 
  group_rows("Total dry biomass, normal", 10,12) %>% 
  group_rows("Partial dbRDA (1st axis), normal",13,15) %>% 
  group_rows("Heterogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (Chi-square, z)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi-square, z)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  save_kable(file = "C:Data//For submission//For resubmission//RESUB2//First look//pstable.community.html ", self_contained = T)
