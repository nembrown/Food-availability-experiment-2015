
# read in useful packages
setwd("C:/Users/norahbrown/Dropbox/Projects/Summer 2015 food availability experiments/Data")





library(tidyr)
library(bbmle) 
library(glmmTMB)
library(doBy)
library(plyr)
require(dplyr)
library(ggplot2) 
library(doBy)
library(grid)
library(glmmADMB)
library(betareg)
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
require(reshape2)
library(ggplot2)
library(grid)
library(DHARMa)
library(gap)
library(qrnn)
library(mgcv)
library(colorspace)
library(gratia)
library(cowplot)



# dataframe management ----------------------------------------------------


head(food.exp.data.12.2019)
food.exp.data.12.2019<-read.csv("C:Mesocosm inventory data//food.exp.data.mesocosm.12.csv")
food.exp.data.tile.all<-read.csv("C:Mesocosm inventory data//food.exp.data.tile.all.csv")
food.caprellid.data<-read.csv("c:Emily caprellid data.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )

food.exp.data.12.2019$oFood.quality<-factor(food.exp.data.12.2019$Food.quality, levels=c("None", "Low", "High"), ordered=TRUE)
food.exp.data.12.2019$Food.quality<-factor(food.exp.data.12.2019$Food.quality, levels=c("None", "Low", "High"), ordered=FALSE)

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

food.exp.data.12.2019$hydroid_dry_biomass[food.exp.data.12.2019$hydroid_dry_biomass<0]<-0
food.exp.data.12.2019$tunicate_dry_biomass[food.exp.data.12.2019$tunicate_dry_biomass<0]<-0
food.exp.data.12.2019$hydtobot_dry_biomass<-(food.exp.data.12.2019$tunicate_dry_biomass)/(food.exp.data.12.2019$tunicate_dry_biomass+food.exp.data.12.2019$hydroid_dry_biomass)


food.exp.data.12.2019$CAP1<-model.meso.bray.scores$CAP1
food.exp.data.12.2019$distances<-bd.bray$distances

food.exp.data.12.2019$Mussel.wet.weight.per.1<-(food.exp.data.12.2019$Mussel.wet.weight)/(food.exp.data.12.2019$mussel_complete+1)


food.exp.data.12.2019$hydroid.001<-0.01*(food.exp.data.12.2019$hydroid)
food.exp.data.12.2019$alive.bot.001<-0.01*(food.exp.data.12.2019$alive.bot)
food.exp.data.12.2019$alive.mem.001<-0.01*(food.exp.data.12.2019$alive.mem)
food.exp.data.12.2019$formicula.001<-0.01*(food.exp.data.12.2019$formicula)
food.exp.data.12.2019$didemnum.001<-0.01*(food.exp.data.12.2019$didemnum+0.01)

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


#attributes(food.exp.data.12.2019_zscores$min.10.pH)
#You can use the attributes to unscale:
  

# Notes on contrasts ------------------------------------------------------


#Notes on ordered factors: my factor food.quality is ordered
#use this one for overall model
#options(contrasts = c("contr.sum", "contr.poly"))
#The contrast function, contr.sum(), gives orthogonal contrasts where you compare every level to the overall mean.

#use this one for the "summary" part to partition low vs high food compared to none - none has to be the first level ... 
#options(contrasts = c("contr.treatment", "contr.poly"))



# #visualizing data -------------------------------------------------------
hist(food.exp.data.12.2019$hydroid)


#Binomial: (doesn't work for hydroid) 
fitBinom=fitdist(data=food.exp.data.12.2019$hydroid*0.01, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$hydroid, "binom", size = 1, prob = fitBinom$estimate[[1]])

#log-normal
qqp(food.exp.data.12.2019$hydroid, "lnorm")
#not great at lower end

#normal
qqp(food.exp.data.12.2019$hydroid, "norm")
#not great at lower end

#Beta distribution
beta.12<-fitdist(0.01*(food.exp.data.12.2019$hydroid+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$hydroid+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#not great at the upper end

#beta-binomial? 

?fitdist
?qqp
# Model building ----------------------------------------------------------


#Binomial
glm.binomial.hydroid.12.hydrogen<- glm(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial(link="logit"))
summary(glm.binomial.hydroid.12.hydrogen)


m1<-glm.binomial.hydroid.12.hydrogen
op<-par(mfrow=c(2,2))
plot(m1)
#7 or 8 points with high leverage

#improved qq plot
qq.gam(m1,pch=16)
#doesn't look good

#check for overdsispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m1)
#definitely overdispersed
#the variance is greater than the mean 


#okay so move on from binomial

#options are 1. alternative distributions: quasi binomial, beta-binomial
#Keep in mind that once you switch to quasi-likelihood you will either have to eschew inferential methods 
#such as the likelihood ratio test, profile confidence intervals, AIC, etc., or make more heroic assumptions 
#to compute "quasi-" analogs of all of the above (such as QAIC).
#not too many options for beta-binomial... 

#2. observation-level random effects
#3. GAM


# #quasibinomial: ---------------------------------------------------------


quasi.hydroid.1<- glm(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = quasibinomial)
summary(quasi.hydroid.1)


quasi.hydroid.2<- glm(formula = cbind(hydroid, 100-hydroid)~ min.10.pH+Food.quality, data = food.exp.data.12.2019_zscores, family = quasibinomial)
summary(quasi.hydroid.2)

anova(quasi.hydroid.1, quasi.hydroid.2, test="F")

#But the assumptions are off here according to Bolker? 


# GLMM: Observation-level random effect added -----------------------------

theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))

#let's try observation-level random effects
#glmmadmb is no longer the best one - Bolker suggests glmmTMB, using Template Model Builder

?glmmTMB


glmm.binomial.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality+(1|Mesocosm), data = food.exp.data.12.2019_zscores, family = binomial)
summary(glmm.binomial.hydroid.12)

glmm.betabinomial.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality+(1|Mesocosm), data = food.exp.data.12.2019_zscores, family = betabinomial)
summary(glmm.betabinomial.hydroid.12)

glmm.betabinomial.nr.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = betabinomial)
summary(glmm.betabinomial.nr.hydroid.12)

glmm.betabinomial.reml.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = betabinomial, REML=TRUE)
summary(glmm.betabinomial.reml.hydroid.12)

AICtab(glmm.binomial.hydroid.12, glmm.betabinomial.hydroid.12, glmm.betabinomial.nr.hydroid.12, glmm.betabinomial.reml.hydroid.12)
#betabinomial best 


#glmmTMB has the capability to simulate from a fitted model. These simulations resample random effects from their estimated distribution.
simo=simulate(glmm.betabinomial.hydroid.12, seed=1)
Simdat=food.exp.data.12.2019_zscores
Simdat$hydroid=simo[[1]]

Simdat=transform(Simdat,
                 type="simulated")
food.exp.data.12.2019_zscores$type = "observed"  
Dat=rbind(food.exp.data.12.2019_zscores, Simdat) 

#Then we can plot the simulated data against the observed data to check if they are similar.

ggplot(Dat,  aes(hydroid, colour=type))+geom_density()+facet_grid(food.exp.data.12.2019_zscores$Food.quality)

#For low food not very similar

#glmm.betabinomial.nr.hydroid.12 -> not great, but glmm.betabinomial.hydroid.12 a bit better? 

#before evaulating inference - need to check diagnostics - search for glmmtmb diagnostics? 

#ggplot(data = NULL) + geom_point(aes(y = residuals(glmm.betabinomial.nr.hydroid.12,type = "pearson", scaled = TRUE), x = fitted(glmm.betabinomial.nr.hydroid.12)))
#cant' yet do fitted objects from betabinomial

#dharma

simulationOutput <- simulateResiduals(fittedModel = glmm.betabinomial.hydroid.12)
plotSimulatedResiduals(simulationOutput = simulationOutput)
#nr. not a perfect fit though - lines should fit  but one is curved and the other two are flat. QQ plot residuals are pretty good... 
#next one down is not great either but a bit better? 

testZeroInflation(simulationOutput)
plotSimulatedResiduals(simulationOutput = simulationOutput, quantreg = T)
#so our model is not great still.... 



augDat <- data.frame(food.exp.data.12.2019_zscores,resid=residuals(glmm.betabinomial.nr.hydroid.12,type="pearson"),
                     fitted=fitted(glmm.betabinomial.nr.hydroid.12))
ggplot(augDat,aes(x=Food.quality,y=resid))+geom_boxplot()+coord_flip()
#doesn't work for glmmTMB



# GAM ---------------------------------------------------------------------
#Gavin Simpson's gratia
library(devtools)
#devtools::install_github('gavinsimpson/gratia')
library(gratia)

set.seed(20)
dat <- gamSim(1, n = 400, dist = "normal", scale = 2, verbose = FALSE)

head(dat)

mod <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat, method = "REML")
summary(mod)

k.check(mod)
#function for checking whether sufficient numbers of basis functions were used in each smooth in the model.
#k-index reported. The further below 1 this is, the more likely it is that there is missed pattern left in the residuals.


draw(mod)
#In gratia, we use the draw() generic to produce ggplot2-like plots from objects. 
#To visualize the four estimated smooth functions in the GAM mod, we would use draw(mod)


evaluate_smooth(mod, "x1")
# if you wanted to extract most of the data used to build the plot 


appraise(mod)
#diagnostics

qq_plot(mod, method = 'simulate')
#the data simulation procedure (method = 'simulate') to generate reference quantiles, 
#which typically have better performance for GLM-like models (Augustin et al., 2012).

#We have a varying coefficient model aka ANCOVA, so use "by"

#Let smoothness selection do all the work by adding a penalty on the null space of each smooth. gam(...,select=TRUE) does this.

#This approach leaves the original smoothing penalty unchanged, but constructs an additional penalty for each smooth, 
#which penalizes only functions in the null space of the original penalty (the 'completely smooth' functions). 
#Hence, if all the smoothing parameters for a term tend to infinity, the term will be selected out of the model. 
#This latter approach is more expensive computationally, but has the advantage that it can be applied automatically to any smooth term. 
#The select argument to gam turns on this method.

#In fact, as implemented, both approaches operate by eigen-decomposiong the original penalty matrix. 
#A new penalty is created on the null space: it is the matrix with the same eigenvectors as the original penalty, 
#but with the originally postive egienvalues set to zero, and the originally zero eigenvalues set to something positive. 
#The first approach just addes a multiple of this penalty to the original penalty, where the multiple is chosen so that the new penalty 
#can not dominate the original. The second approach treats the new penalty as an extra penalty, with its own smoothing parameter.

#Of course, as with all model selection methods, some care must be take to ensure that the automatic selection is sensible, 
#and a decision about the effective degrees of freedom at which to declare a term 'negligible' has to be made


## I do in fact have ordered factor - meaning "none" is the base level and the other two relate to that....
#https://www.fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/
#The first s(x) in the model is the smooth effect of x on the reference level of the ordered factor of. 
#The second smoother, s(x, by = of) is the set of L???1 difference smooths, which model the smooth differences 
#between the reference level smoother and those of the individual levels (excluding the reference one)
#the intercept is then an estimate of the mean response for the reference level of the factor, and the remaining 
#model coefficients estimate the differences between the mean response of the reference level and that of the other factor levels.

#"oFood.quality" is ordered "Food.quality" is unordered

#select=TRUE is automatic model selection via null space penalization



# GAM beta hydroid / gam.beta.hydroid.12 --------------------------------------------------------


gam.binomial.hydroid.12<- gam(formula = cbind(hydroid, 100-hydroid)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")


food.exp.data.12.2019_zscores$hydroid.001<-0.01*(food.exp.data.12.2019_zscores$hydroid+0.01)

gam.beta.hydroid.12<- gam(hydroid.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydroid.12.1<- gam(hydroid.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydroid.12.2<- gam(hydroid.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.hydroid.12.3<- gam(hydroid.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(gam.beta.hydroid.12, gam.beta.hydroid.12.1, gam.beta.hydroid.12.2, gam.binomial.hydroid.12, glmm.betabinomial.nr.hydroid.12)
#logit REML and cloglog same ... just logit.... 

plot(gam.beta.hydroid.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.hydroid.12)
qq_plot(gam.beta.hydroid.12, method = 'simulate')
k.check(gam.beta.hydroid.12)

summary(gam.beta.hydroid.12)

gam.beta.hydroid.12.unordered<- gam(hydroid.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.beta.hydroid.12)
summary(gam.beta.hydroid.12.unordered)


#Confidence intervals
#https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
#install.packages("colorspace", repos = "http://R-Forge.R-project.org")
#library(colorspace)
#library(ggplot2)




fam.gam.hydroid <- family(gam.beta.hydroid.12)
fam.gam.hydroid 
str(fam.gam.hydroid )
ilink.gam.hydroid <- fam.gam.hydroid$linkinv
ilink.gam.hydroid


str(ndata.hydroid)
str(food.exp.data.12.2019_zscores)

attributes(food.exp.data.12.2019_zscores$min.10.pH)
attributes(ndata.hydroid$min.10.pH)


food.exp.data.12.2019_zscores$min.10.pH.unscaled <-food.exp.data.12.2019_zscores$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')
head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)
            
mod.hydroid<-gam.beta.hydroid.12
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




ndata.hydroid$min.10.pH.unscaled<-ndata.hydroid$min.10.pH * attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$min.10.pH, 'scaled:center')

# plot 



colorset2 = c("High"="#F8A02E" ,"Low"="#439E5F","None"= "#666666")

plt.gam.hydroid <- ggplot(ndata.hydroid, aes(x = min.10.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = hydroid.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Obelia")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.hydroid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.gam.hydroid 
ggsave("C:Graphs May 2019//hydroid_pred.png")

head(ndata.hydroid)


# GAM beta botryllus / gam.beta.alive.bot.12.2 ----------------------------------------------------

#binomial first
gam.binomial.alive.bot.12<- gam(formula = cbind(alive.bot, 100-alive.bot)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, REML=TRUE)

food.exp.data.12.2019_zscores$alive.bot.001<-0.01*(food.exp.data.12.2019_zscores$alive.bot+0.01)

#beta next
gam.beta.alive.bot.12<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.beta.alive.bot.12.1<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.beta.alive.bot.12.2<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.beta.alive.bot.12.3<- gam(alive.bot.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.beta.alive.bot.12, gam.beta.alive.bot.12.1, gam.beta.alive.bot.12.2, gam.binomial.alive.bot.12)
#clog log  is the best
#cauchit doesn't work well with REML

plot(gam.beta.alive.bot.12.2, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.alive.bot.12.2)
qq_plot(gam.beta.alive.bot.12.2, method = 'simulate')
k.check(gam.beta.alive.bot.12.2)
summary(gam.beta.alive.bot.12.2)

gam.beta.alive.bot.12.2.unordered<- gam(alive.bot.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)


fam.gam.alive.bot <- family(gam.beta.alive.bot.12.2)
fam.gam.alive.bot
str(fam.gam.alive.bot)
ilink.gam.alive.bot<- fam.gam.alive.bot$linkinv
ilink.gam.alive.bot


mod.alive.bot<-gam.beta.alive.bot.12.2
ndata.alive.bot <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = alive.bot.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Botryllus")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.alive.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.alive.bot
ggsave("C:Graphs May 2019//alive.bot_pred.png")



# GAM negbin caprellid / gam.nb.caprellid.12 -----------------------------------------------------------
poisson.12<-fitdistr(food.caprellid.data_zscores$total.caprellids, "Poisson")
qqp(food.caprellid.data_zscores$total.caprellids, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.caprellid <- fitdistr(food.caprellid.data_zscores$total.caprellids, "Negative Binomial")
qqp(food.caprellid.data_zscores$total.caprellids, "nbinom", size = nbinom12.caprellid$estimate[[1]], mu = nbinom12.caprellid$estimate[[2]])


glm.poisson.total.caprellids.12.hydrogen<- glm(total.caprellids~ min.10.pH*oFood.quality, data = food.caprellid.data_zscores, family="poisson")
plot(glm.poisson.total.caprellids.12.hydrogen)
overdisp_fun(glm.poisson.total.caprellids.12.hydrogen)
#definitely overdispersed --> neg binomial way better

glm.nb.total.caprellids.12.hydrogen<- glm.nb(total.caprellids~ min.10.pH*oFood.quality, data = food.caprellid.data_zscores)
plot(glm.nb.total.caprellids.12.hydrogen)
#This one is actually okay..... 

gam.nb.caprellid.12<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = negbin(nbinom12.caprellid$estimate[[1]]), select=TRUE, method="REML")
gam.nb.caprellid.12.1<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.caprellid.12<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = poisson, select=TRUE, method="REML")
gam.lm.caprellid.12<- gam(total.caprellids ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = gaussian, select=TRUE, method="REML")
gam.log.lm.caprellid.12<- gam(log(total.caprellids+1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = gaussian, select=TRUE, method="REML")


AICtab(gam.lm.caprellid.12, gam.log.lm.caprellid.12, gam.nb.caprellid.12.1,gam.nb.caprellid.12, gam.poisson.caprellid.12, glm.nb.total.caprellids.12.hydrogen, glm.poisson.total.caprellids.12.hydrogen)
#gam is better - need to make sure it's reading the correct theta ... 
#also make sure min.10.pH is reading as a matrix? 

appraise(gam.log.lm.caprellid.12)
qq_plot(gam.log.lm.caprellid.122, method = 'simulate')
plot(gam.log.lm.caprellid.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.log.lm.caprellid.12)
summary(gam.log.lm.caprellid.12.unordered)

gam.log.lm.caprellid.12.unordered<- gam(log(total.caprellids+1) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.caprellid.data_zscores, select=TRUE, method="REML")

fam.gam.caprellid <- family(gam.log.lm.caprellid.12)
fam.gam.caprellid
str(fam.gam.caprellid)
ilink.gam.caprellid<- fam.gam.caprellid$linkinv
ilink.gam.caprellid

mod.caprellid<-gam.log.lm.caprellid.12
ndata.caprellid <- with(food.caprellid.data_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.caprellid

str(food.caprellid.data_zscores)
str(ndata.caprellid)

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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(total.caprellids+1), shape=CO2, colour=oFood.quality), size=3, data = food.caprellid.data_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Caprella")~ "abundance", paste("(Log # of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.caprellid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.caprellid

ggsave("C:Graphs May 2019//caprellid_pred.png")


# GAM caprellid percent -----------------------------------------------------------
fitBinom=fitdist(data=food.exp.data.12.2019$caprellid.percent, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$caprellid.percent, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$caprellid.percent+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$caprellid.percent+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#good



glm.binomial.caprellid.percent.12.hydrogen<-glm(formula = cbind(caprellid.percent, 100-caprellid.percent)~ (min.10.pH)*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)

gam.binomial.caprellid.percent.12<- gam(formula = cbind(caprellid.percent, 100-caprellid.percent)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$caprellid.percent.001<-0.01*(food.exp.data.12.2019_zscores$caprellid.percent+0.01)

#beta next
gam.beta.caprellid.percent.12<- gam(caprellid.percent.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.caprellid.percent.12.1<- gam(caprellid.percent.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.caprellid.percent.12.2<- gam(caprellid.percent.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.caprellid.percent.12.3<- gam(caprellid.percent.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.caprellid.percent.12, gam.beta.caprellid.percent.12.1, gam.beta.caprellid.percent.12.2,gam.binomial.caprellid.percent.12, glm.binomial.caprellid.percent.12.hydrogen)
#cloglog is the best along with logit - go with simpler logit


plot(gam.beta.caprellid.percent.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

appraise(gam.beta.caprellid.percent.12)
qq_plot(gam.beta.caprellid.percent.12, method = 'simulate')
k.check(gam.beta.caprellid.percent.12)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = caprellid.percent.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Caprella")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.caprellid.percent,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.caprellid.percent
ggsave("C:Graphs May 2019//caprellid.percent_pred.png")



# GAM beta formicula / gam.beta.formicula.12 -----------------------------------------------------------
str(food.exp.data.12.2019_zscores)
par(op)

fitBinom=fitdist(data=food.exp.data.12.2019$formicula, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$formicula, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$formicula+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$formicula+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#not great at the upper end


#binomial first

glm.binomial.formicula.12.hydrogen<-glm(formula = cbind(formicula, 100-formicula)~ (min.10.pH)*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)

gam.binomial.formicula.12<- gam(formula = cbind(formicula, 100-formicula)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$formicula.001<-0.01*(food.exp.data.12.2019_zscores$formicula+0.01)

#beta next
gam.beta.formicula.12<- gam(formicula.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.formicula.12.1<- gam(formicula.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.formicula.12.2<- gam(formicula.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.formicula.12.3<- gam(formicula.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.formicula.12, gam.beta.formicula.12.1, gam.beta.formicula.12.2,gam.binomial.formicula.12, glm.binomial.formicula.12.hydrogen)
#cloglog is the best along with logit - go with simpler logit


plot(gam.beta.formicula.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

appraise(gam.beta.formicula.12)
qq_plot(gam.beta.formicula.12, method = 'simulate')
k.check(gam.beta.formicula.12)
summary(gam.beta.formicula.12.4)

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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = formicula.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Folliculina")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.formicula,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.formicula
ggsave("C:Graphs May 2019//formicula_pred.png")



# GAM beta membranipora / gam.beta.alive.mem.12 --------------------------------------------------------

fitBinom=fitdist(data=food.exp.data.12.2019$alive.mem, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$alive.mem, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$alive.mem+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$alive.mem+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#good

####GLM betareg??? 


glm.binomial.alive.mem.12.hydrogen<-glm(formula = cbind(alive.mem, 100-alive.mem)~ (min.10.pH)*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)

#binomial first
gam.binomial.alive.mem.12<- gam(formula = cbind(alive.mem, 100-alive.mem)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$alive.mem.001<-0.01*(food.exp.data.12.2019_zscores$alive.mem+0.01)

#beta next
gam.beta.alive.mem.12<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.alive.mem.12.1<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.alive.mem.12.2<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.alive.mem.12.3<- gam(alive.mem.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.beta.alive.mem.12, gam.beta.alive.mem.12.1, gam.beta.alive.mem.12.2, gam.binomial.alive.mem.12, glm.binomial.alive.mem.12.hydrogen)
#logit is of the best (but 0 between it and cloglog so simpler is better)


plot(gam.beta.alive.mem.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.alive.mem.12)
qq_plot(gam.beta.alive.mem.12, method = 'simulate')
k.check(gam.beta.alive.mem.12)
summary(gam.beta.alive.mem.12)
vis.gam(gam.beta.alive.mem.12)

gam.beta.alive.mem.12.unordered<- gam(alive.mem.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.alive.mem <- family(gam.beta.alive.mem.12)
fam.gam.alive.mem
str(fam.gam.alive.mem)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = alive.mem.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Membranipora")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.alive.mem,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.alive.mem
ggsave("C:Graphs May 2019//alive.mem_pred.png")


# GAM beta didemnum / gam.beta.didemnum.12 ------------------------------------------------------------


str(food.exp.data.12.2019_zscores)
par(op)
fitBinom=fitdist(data=food.exp.data.12.2019$didemnum, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$didemnum, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$didemnum+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$didemnum+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#not great at all


#binomial first
gam.binomial.didemnum.12<- gam(formula = cbind(didemnum, 100-didemnum)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$didemnum.001<-0.01*(food.exp.data.12.2019_zscores$didemnum+0.01)

#beta next
gam.beta.didemnum.12<- gam(didemnum.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.didemnum.12.1<- gam(didemnum.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.didemnum.12.2<- gam(didemnum.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.didemnum.12.3<- gam(didemnum.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.didemnum.12, gam.beta.didemnum.12.1, gam.beta.didemnum.12.2,  gam.binomial.didemnum.12)
#logit is the best

plot(gam.beta.didemnum.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.didemnum.12)
qq_plot(gam.beta.didemnum.12, method = 'simulate')
k.check(gam.beta.didemnum.12)
summary(gam.beta.didemnum.12)
vis.gam(gam.beta.didemnum.12)

#appraise doesn't fit that well .... *** what to do? 

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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = didemnum.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Didemnum")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.didemnum,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.didemnum
ggsave("C:Graphs May 2019//didemnum_pred.png")



# GAM negbin mussel incomplete / gam.nb.mussel_complete.12 ---------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$mussel_complete, "Poisson")
qqp(food.exp.data.12.2019_zscores$mussel_complete, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.mussel <- fitdistr(food.exp.data.12.2019_zscores$mussel_complete, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$mussel_complete, "nbinom", size = nbinom12.mussel$estimate[[1]], mu = nbinom12.mussel$estimate[[2]])
#better than poisson



glm.poisson.mussel_complete.12<- glm(mussel_complete~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")

glm.nb.mussel_complete.12.hydrogen<- glm.nb(mussel_complete~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

AICtab(glm.poisson.mussel_complete.12, glm.nb.mussel_complete.12.hydrogen)


plot(glm.nb.mussel_complete.12.hydrogen)
plot(glm.poisson.mussel_complete.12.hydrogen)
#This one is actually okay..... 

overdisp_fun(glm.poisson.mussel_complete.12)
#very overdispersed


#negative binomial first
gam.nb.mussel_complete.12<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")

gam.nb.mussel_complete.12.1<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(link="log"), select=TRUE, method="REML")
gam.nb.mussel_complete.12.2<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(link="sqrt"), select=TRUE, method="REML")
gam.nb.mussel_complete.12.3<- gam(mussel_complete ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")



nb(link="sqrt")
#negbin with theta estimated is better

AICtab(gam.nb.mussel_complete.12, glm.nb.mussel_complete.12.hydrogen, gam.nb.mussel_complete.12.1, gam.nb.mussel_complete.12.2, gam.nb.mussel_complete.12.3)
#gam is better


plot(gam.nb.mussel_complete.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.nb.mussel_complete.12)
qq_plot(gam.nb.mussel_complete.12, method = 'simulate')
k.check(gam.nb.mussel_complete.12)
summary(gam.nb.mussel_complete.12)


gam.nb.mussel_complete.12.unordered<- gam(mussel_complete ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")
summary(gam.nb.mussel_complete.12.unordered)

#a few outside the area

fam.gam.mussel_complete <- family(gam.nb.mussel_complete.12)
fam.gam.mussel_complete
str(fam.gam.mussel_complete)
ilink.gam.mussel_complete<- fam.gam.mussel_complete$linkinv
ilink.gam.mussel_complete

mod.mussel_complete<-gam.nb.mussel_complete.12
ndata.mussel_complete <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.mussel_complete

str(food.exp.data.12.2019_zscores)
str(ndata.mussel_complete)

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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = mussel_complete, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Mytilus")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.mussel_complete,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.mussel_complete
ggsave("C:Graphs May 2019//mussel_complete_pred.png")


# GAM negbin barnacles / gam.nb.num.barn.alive.12 -----------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.barn.alive, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.barn.alive, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.barnacles <- fitdistr(food.exp.data.12.2019_zscores$num.barn.alive, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.barn.alive, "nbinom", size = nbinom12.barnacles$estimate[[1]], mu = nbinom12.barnacles$estimate[[2]])
#better than poisson


glm.nb.num.barn.alive.12.hydrogen<- glm.nb(num.barn.alive~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.binomial.num.barn.alive.12.hydrogen)

glm.poisson.num.barn.alive.12.hydrogen<- glm(num.barn.alive~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.barn.alive.12.hydrogen)
#yes overdispersed

#negative binomial 
gam.nb.num.barn.alive.12<- gam(num.barn.alive ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.barnacles$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.barn.alive.12.1<- gam(num.barn.alive ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")

AICtab(gam.nb.num.barn.alive.12, glm.nb.num.barn.alive.12.hydrogen, gam.nb.num.barn.alive.12.1)

plot(gam.nb.num.barn.alive.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.nb.num.barn.alive.12)
qq_plot(gam.nb.num.barn.alive.12, method = 'simulate')
k.check(gam.nb.num.barn.alive.12)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.barn.alive, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Balanus")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.num.barn.alive,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.num.barn.alive
ggsave("C:Graphs May 2019//num.barn.alive_pred.png")



# GAM negbin disporella / gam.nb.disporella.12 ----------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$disporella, "Poisson")
qqp(food.exp.data.12.2019_zscores$disporella, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.disporella <- fitdistr(food.exp.data.12.2019_zscores$disporella, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$disporella, "nbinom", size = nbinom12.disporella$estimate[[1]], mu = nbinom12.disporella$estimate[[2]])
#better than poisson
#strill not great



glm.binomial.disporella.12.hydrogen<- glm.nb(disporella~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.binomial.disporella.12.hydrogen)
#not great

glm.poisson.disporella.12.hydrogen<- glm(disporella~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.disporella.12.hydrogen)
#poisson is overdispersed

#negative binomial
gam.nb.disporella.12<- gam(disporella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")
gam.nb.disporella.12.1<- gam(disporella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")



AICtab(gam.nb.disporella.12.1,gam.nb.disporella.12, glm.binomial.disporella.12.hydrogen)

plot(gam.nb.disporella.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.nb.disporella.12)
qq_plot(gam.nb.disporella.12, method = 'simulate')
k.check(gam.nb.disporella.12)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = disporella, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Disporella")~ "abundance", paste("(# of colonies)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.disporella,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.disporella
ggsave("C:Graphs May 2019//disporella_pred.png")


# GAM negbin schizo / gam.nb.schizo.12 --------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$schizo, "Poisson")
qqp(food.exp.data.12.2019_zscores$schizo, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.schizo <- fitdistr(food.exp.data.12.2019_zscores$schizo, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$schizo, "nbinom", size = nbinom12.schizo$estimate[[1]], mu = nbinom12.schizo$estimate[[2]])
#strill not great



glm.nb.schizo.12.hydrogen<- glm.nb(schizo~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.nb.schizo.12.hydrogen)


glm.poisson.schizo.12.hydrogen<- glm(schizo~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.schizo.12.hydrogen)
#poisson is not overdispersed

?gam
#negative binomial first
gam.nb.schizo.12<- gam(schizo ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.schizo$estimate[[1]]), select=TRUE, method="REML")
gam.nb.schizo.12.1<- gam(schizo ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.poisson.schizo.12<- gam(schizo ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.schizo.12, gam.nb.schizo.12.1, glm.nb.schizo.12.hydrogen, glm.poisson.schizo.12.hydrogen, gam.poisson.schizo.12)


plot(gam.nb.schizo.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.nb.schizo.12)
qq_plot(gam.nb.schizo.12, method = 'simulate')
k.check(gam.nb.schizo.12)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = schizo, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Schizoporella")~ "abundance", paste("(# of colonies)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.schizo,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.schizo
ggsave("C:Graphs May 2019//schizo_pred.png")

help(summary.gam)


# GAM poisson num nudi / gam.poisson.num.nudi.12  ------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.nudi, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.nudi, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.num.nudi <- fitdistr(food.exp.data.12.2019_zscores$num.nudi, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.nudi, "nbinom", size = nbinom12.num.nudi$estimate[[1]], mu = nbinom12.num.nudi$estimate[[2]])
#strill not great

glm.nb.num.nudi.12.hydrogen<- glm.nb(num.nudi~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019_zscores)
plot(glm.nb.num.nudi.12.hydrogen)



glm.poisson.num.nudi.12.hydrogen<- glm(num.nudi~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.nudi.12.hydrogen)
#poisson is not overdispersed

?gam
#negative binomial first
gam.nb.num.nudi.12.1<- gam(num.nudi ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.nb.num.nudi.12<- gam(num.nudi ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.nudi$estimate[[1]]), select=TRUE, method="REML")


gam.poisson.num.nudi.12<- gam(num.nudi ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.nudi.12, gam.nb.num.nudi.12.1, glm.nb.num.nudi.12.hydrogen, glm.poisson.num.nudi.12.hydrogen, gam.poisson.num.nudi.12)

###poisson is the best fit by 2 dAIC


appraise(gam.poisson.num.nudi.12)
qq_plot(gam.poisson.num.nudi.12, method = 'simulate')
plot(gam.poisson.num.nudi.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.poisson.num.nudi.12)
summary(gam.poisson.num.nudi.12)

#resids a bit funny but same in neg bin
gam.poisson.num.nudi.12.unordered<- gam(num.nudi ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.num.nudi <- family(gam.poisson.num.nudi.12)
fam.gam.num.nudi
str(fam.gam.num.nudi)
ilink.gam.num.nudi<- fam.gam.num.nudi$linkinv
ilink.gam.num.nudi

mod.num.nudi<-gam.poisson.num.nudi.12
ndata.num.nudi <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
                                                                               length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.num.nudi

str(food.exp.data.12.2019_zscores)
str(ndata.num.nudi)

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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.nudi, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Hermissenda")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.num.nudi
ggsave("C:Graphs May 2019//num.nudi_pred.png")


# GAM nb() serpulids / gam.nb.num.serpulid.12.1 -----------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.serpulid, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.serpulid, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.num.serpulid <- fitdistr(food.exp.data.12.2019_zscores$num.serpulid, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.serpulid, "nbinom", size = nbinom12.num.serpulid$estimate[[1]], mu = nbinom12.num.serpulid$estimate[[2]])
#better, not great

glm.nb.num.serpulid.12<- glm.nb(num.serpulid~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.nb.num.serpulid.12)
summary(glm.nb.num.serpulid.12)


glm.poisson.num.serpulid.12.hydrogen<- glm(num.serpulid~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.serpulid.12.hydrogen)
#poisson is definitely overdispersed

#negative binomial first
gam.nb.num.serpulid.12<- gam(num.serpulid ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.serpulid.12.1<- gam(num.serpulid ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.poisson.num.serpulid.12<- gam(num.serpulid ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.serpulid.12, gam.nb.num.serpulid.12.1, glm.nb.num.serpulid.12, glm.poisson.num.serpulid.12.hydrogen, gam.poisson.num.serpulid.12)

###glm is best but same as GAM nb() --> i.e. not negbin... 
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.nb.num.serpulid.12.1)
qq_plot(gam.nb.num.serpulid.12.1, method = 'simulate')
plot(gam.nb.num.serpulid.12.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.nb.num.serpulid.12.1)
summary(gam.nb.num.serpulid.12.1)

#residuals a bit weird .... but they are the same in glm as in gam??? whic
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.serpulid, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop("Serpulid abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.num.serpulid
ggsave("C:Graphs May 2019//num.serpulid_pred.png")


# GAM negbin orange sponge / gam.nb.orange_sponge.12 -------------------------------------------------------



poisson.12<-fitdistr(food.exp.data.12.2019_zscores$orange_sponge, "Poisson")
qqp(food.exp.data.12.2019_zscores$orange_sponge, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.orange_sponge <- fitdistr(food.exp.data.12.2019_zscores$orange_sponge, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$orange_sponge, "nbinom", size = nbinom12.orange_sponge$estimate[[1]], mu = nbinom12.orange_sponge$estimate[[2]])
#better, not great

glm.nb.orange_sponge.12<- glm.nb(orange_sponge~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

glm.poisson.orange_sponge.12.hydrogen<- glm(orange_sponge~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.orange_sponge.12.hydrogen)
#poisson is definitely overdispersed

#negative binomial first
gam.nb.orange_sponge.12<- gam(orange_sponge ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.orange_sponge$estimate[[1]]), select=TRUE, method="REML")
gam.nb.orange_sponge.12.1<- gam(orange_sponge ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")

gam.poisson.orange_sponge.12<- gam(orange_sponge ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.orange_sponge.12, gam.nb.orange_sponge.12.1, glm.nb.orange_sponge.12, glm.poisson.orange_sponge.12.hydrogen, gam.poisson.orange_sponge.12)

###gam negbinis best but same as GAM nb() 
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.nb.orange_sponge.12)
qq_plot(gam.nb.orange_sponge.12, method = 'simulate')
plot(gam.nb.orange_sponge.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.nb.orange_sponge.12)
summary(gam.nb.orange_sponge.12)

#residuals a bit weird .... but they are the same in glm as in gam??? whic
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = orange_sponge, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop("Sponge abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.orange_sponge,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.orange_sponge
ggsave("C:Graphs May 2019//orange_sponge_pred.png")



# GAM negbin corella / gam.nb.num.corella.12 -------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.corella, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.corella, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.num.corella <- fitdistr(food.exp.data.12.2019_zscores$num.corella, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.corella, "nbinom", size = nbinom12.num.corella$estimate[[1]], mu = nbinom12.num.corella$estimate[[2]])
#better, not great

glm.nb.num.corella.12<- glm.nb(num.corella~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

glm.poisson.num.corella.12.hydrogen<- glm(num.corella~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.corella.12.hydrogen)
#poisson is definitely overdispersed

#negative binomial first
gam.nb.num.corella.12<- gam(num.corella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.nb.num.corella.12.1<- gam(num.corella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.poisson.num.corella.12<- gam(num.corella ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.num.corella.12, gam.nb.num.corella.12.1, glm.nb.num.corella.12, glm.poisson.num.corella.12.hydrogen, gam.poisson.num.corella.12)

###gam negbinis best
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.nb.num.corella.12)
qq_plot(gam.nb.num.corella.12, method = 'simulate')
plot(gam.nb.num.corella.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.nb.num.corella.12)
summary(gam.nb.num.corella.12)

#residuals a bit patterny? 
gam.nb.num.corella.12.unordered<- gam(num.corella ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.num.corella <- family(gam.nb.num.corella.12)
fam.gam.num.corella
str(fam.gam.num.corella)
ilink.gam.num.corella<- fam.gam.num.corella$linkinv
ilink.gam.num.corella

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

# plot 
plt.num.corella <- ggplot(ndata.num.corella, aes(x = min.10.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.corella, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(atop(italic("Corella")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.num.corella
ggsave("C:Graphs May 2019//num.corella_pred.png")



# GAM poisson clam / gam.poisson.clam.12 ----------------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$clam, "Poisson")
qqp(food.exp.data.12.2019_zscores$clam, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.clam <- fitdistr(food.exp.data.12.2019_zscores$clam, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$clam, "nbinom", size = nbinom12.clam$estimate[[1]], mu = nbinom12.clam$estimate[[2]])
#better, not great

glm.nb.clam.12<- glm.nb(clam~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

glm.poisson.clam.12.hydrogen<- glm(clam~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.clam.12.hydrogen)
#poisson is not overdispersed

#negative binomial first
gam.nb.clam.12<- gam(clam ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.clam$estimate[[1]]), select=TRUE, method="REML")
gam.nb.clam.12.1<- gam(clam ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.poisson.clam.12<- gam(clam ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.nb.clam.12, gam.nb.clam.12.1, glm.nb.clam.12, glm.poisson.clam.12.hydrogen, gam.poisson.clam.12)

###gam poisson is best
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.poisson.clam.12)
qq_plot(gam.poisson.clam.12, method = 'simulate')
plot(gam.poisson.clam.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.poisson.clam.12)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = clam, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Clam abundance\n(# of individuals)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.clam,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.clam
ggsave("C:Graphs May 2019//clam_pred.png")


### legend plot
par(xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    oma = c(2, 2, 0, 0), # move plot to the right and up
    mgp = c(2, 1, 0) # move axis labels closer to axis
) 

plt.legend <- ggplot(ndata.clam, aes(x = min.10.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = clam, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Clam abundance\n(# of individuals)")+  
  scale_color_manual(values=colorset2, name = "Food quality", labels = c("High-quality", "Low quality", "None"))+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17), name = "pH")+
  guides(color = guide_legend(order = 2), shape = guide_legend(order = 1), fill = FALSE)+
  geom_ribbon(data = ndata.clam,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)
plt.legend 

legend_food <- get_legend(plt.legend + theme(legend.position="bottom"))

# Fig 1 plot generation ---------------------------------------------------


library(cowplot)
fig.2<-plot_grid(plt.gam.hydroid,plt.alive.bot,plt.formicula,plt.caprellid.percent,plt.alive.mem,plt.didemnum,
          plt.mussel_complete,plt.num.barn.alive,plt.disporella,plt.schizo,plt.num.nudi,plt.num.serpulid,
          plt.orange_sponge,plt.num.corella,plt.clam,legend_food, ncol=5, rel_heights = c(1,1,1,.2), align='v', axis = 'l', 
          labels=c('(a)', '(b)','(c)', '(d)', '(e)', '(f)', '(g)', 
                   '(h)', '(i)', '(j)','(k)','(l)','(m)','(n)','(o)', ''))

fig.2

ggplot2::ggsave("C:For submission//Fig2.png", width=65, height=35, units="cm")



# Pulling model results to a table ----------------------------------------


#edf is low 0.0124
#You can often find spline terms with EDF <1 that are not linear. 
#In that case it seems that the null space (the linear, perfectly smooth bits)
#has been shrunk as well as the range space (the wiggly bits) but the smoothness parameter(s) 
#for the wiggly bit allow some small amount of wiggliness and the shrinkage to the null space drags 
#the EDF below 1. This is fine; the model is estimating a slightly wiggly function but uncertainty in 
#that estimate probably means that a linear function will reside in the 95% confidence interval for the smooth.
#https://stats.stackexchange.com/questions/179591/gam-with-low-e-d-f-estimated-degrees-of-freedom-value-in-main-effect-not-inte

hydroid.gam<-summary(gam.beta.hydroid.12)
botryllus.gam<-summary(gam.beta.alive.bot.12.2)
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

hydroid.gam.unordered<-summary(gam.beta.hydroid.12.unordered)
botryllus.gam.unordered<-summary(gam.beta.alive.bot.12.2.unordered)
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
save_kable(file = "C:For submission//ptable.html", self_contained = T)


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
  save_kable(file = "C:For submission//stable.html", self_contained = T)

  
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
  save_kable(file = "C:For submission//pstable.html", self_contained = T)


# Richness ----------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$richness, "Poisson")
qqp(food.exp.data.12.2019_zscores$richness, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.richness<- fitdistr(food.exp.data.12.2019_zscores$richness, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$richness, "nbinom", size = nbinom12.richness$estimate[[1]], mu = nbinom12.richness$estimate[[2]])
#worse than poiisson


glm.nb.richness.12.hydrogen<- glm.nb(richness~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.binomial.richness.12.hydrogen)

glm.poisson.richness.12.hydrogen<- glm(richness~ min.10.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.richness.12.hydrogen)
#not overdispersed

#negative binomial 
gam.nb.richness.12<- gam(richness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.richness$estimate[[1]]), select=TRUE, method="REML")
gam.nb.richness.12.1<- gam(richness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.poisson.richness.12<- gam(richness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")



AICtab(gam.nb.richness.12, glm.nb.richness.12.hydrogen, gam.nb.richness.12.1, gam.poisson.richness.12)

plot(gam.poisson.richness.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.poisson.richness.12)
qq_plot(gam.poisson.richness.12, method = 'simulate')
k.check(gam.poisson.richness.12)
summary(gam.poisson.richness.12)


#a few outside the QQ plot on both ends... doesn't look great 


gam.poisson.richness.12.unordered<- gam(richness ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")


fam.gam.richness <- family(gam.poisson.richness.12)
fam.gam.richness
str(fam.gam.richness)
ilink.gam.richness<- fam.gam.richness$linkinv
ilink.gam.richness


mod.richness<-gam.poisson.richness.12
ndata.richness <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = richness, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species richness")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.richness,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ylim(0,20)
plt.richness
ggsave("C:Graphs May 2019//richness_pred.png")







# Evenness ----------------------------------------------------------------
gamma.12.evenness<-fitdistr(food.exp.data.12.2019_zscores$evenness+0.01, "gamma")
qqp(food.exp.data.12.2019$evenness, "gamma", shape = gamma.12.evenness$estimate[[1]], rate = gamma.12.evenness$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$evenness)
#normal works

lm.evenness<-lm(evenness ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
glm.evenness<-glm(evenness ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.evenness.12<- gam(evenness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.evenness.12.1<- gam(evenness ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")



AICtab(lm.evenness, glm.evenness, gam.lm.evenness.12, gam.gamma.evenness.12.1)

plot(gam.lm.evenness.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.evenness.12)
qq_plot(gam.lm.evenness.12, method = 'simulate')
k.check(gam.lm.evenness.12)
summary(gam.lm.evenness.12)





gam.lm.evenness.12.unordered<- gam(evenness ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.evenness <- family(gam.lm.evenness.12)
fam.gam.evenness
str(fam.gam.evenness)
ilink.gam.evenness<- fam.gam.evenness$linkinv
ilink.gam.evenness


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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = evenness, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Species evenness")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.evenness
ggsave("C:Graphs May 2019//evenness_pred.png")





# Occupied space ----------------------------------------------------------

food.exp.data.12.2019_zscores$occupied.space<-100-food.exp.data.12.2019_zscores$bare

beta.12<-fitdist(0.01*(food.exp.data.12.2019_zscores$occupied.space), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019_zscores$occupied.space), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#good! 

gamma.12.occupied.space<-fitdistr((food.exp.data.12.2019_zscores$occupied.space*0.01), "gamma")
qqp(food.exp.data.12.2019_zscores$occupied.space, "gamma", shape = gamma.12.occupied.space$estimate[[1]], rate = gamma.12.occupied.space$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$occupied.space)
#normal is okay not great 

lm.occupied.space<-lm(occupied.space ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
glm.occupied.space<-glm(occupied.space*0.01 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.occupied.space.12<- gam(occupied.space ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.occupied.space.12.1<- gam(occupied.space*0.01 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

gam.binomial.occupied.space.12<- gam(formula = cbind(occupied.space, 100-occupied.space)~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")


food.exp.data.12.2019_zscores$occupied.space.001<-0.01*(food.exp.data.12.2019_zscores$occupied.space)

gam.beta.occupied.space.12<- gam(occupied.space.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.occupied.space.12.1<- gam(occupied.space.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.occupied.space.12.2<- gam(occupied.space.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.occupied.space.12.3<- gam(occupied.space.001~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(lm.occupied.space, glm.occupied.space, gam.lm.occupied.space.12, gam.gamma.occupied.space.12.1, gam.beta.occupied.space.12, gam.beta.occupied.space.12.1, gam.beta.occupied.space.12.2, gam.binomial.occupied.space.12)


#beta is best 

plot(gam.beta.occupied.space.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.occupied.space.12)
qq_plot(gam.beta.occupied.space.12, method = 'simulate')
k.check(gam.beta.occupied.space.12)
summary(gam.beta.occupied.space.12)





gam.beta.occupied.space.12.unordered<- gam(occupied.space.001~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.beta.occupied.space.12.unordered)

fam.gam.occupied.space <- family(gam.beta.occupied.space.12)
fam.gam.occupied.space
str(fam.gam.occupied.space)
ilink.gam.occupied.space<- fam.gam.occupied.space$linkinv
ilink.gam.occupied.space


mod.occupied.space<-gam.beta.occupied.space.12
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = occupied.space.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.occupied.space
ggsave("C:Graphs May 2019//occupied.space_pred.png")



# Everything wet weight ---------------------------------------------------

gamma.12.everything.wet.weight<-fitdistr(food.exp.data.12.2019_zscores$everything.wet.weight+0.01, "gamma")
qqp(food.exp.data.12.2019$everything.wet.weight, "gamma", shape = gamma.12.everything.wet.weight$estimate[[1]], rate = gamma.12.everything.wet.weight$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$everything.wet.weight)

qqp(log(food.exp.data.12.2019_zscores$everything.wet.weight))
#log normal the best

qqp(food.exp.data.12.2019_zscores$everything.wet.weight, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.everything.wet.weight<-lm(everything.wet.weight ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.everything.wet.weight<-lm(log(everything.wet.weight) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.everything.wet.weight<-glm(everything.wet.weight ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.loglink.everything.wet.weight<-glm(everything.wet.weight ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.lm.everything.wet.weight.12<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.everything.wet.weight.12<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.everything.wet.weight.12.1, gam.lm.log.everything.wet.weight.12, gam.tweedie.everything.wet.weight.12.1, lm.log.everything.wet.weight, glm.loglink.everything.wet.weight, lm.everything.wet.weight, glm.everything.wet.weight, gam.lm.everything.wet.weight.12, gam.gamma.everything.wet.weight.12.1)

#gam.lm.log.everything.wet.weight.12 is best 


plot(gam.lm.log.everything.wet.weight.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.log.everything.wet.weight.12)
qq_plot(gam.lm.log.everything.wet.weight.12, method = 'simulate')
k.check(gam.lm.log.everything.wet.weight.12)
summary(gam.lm.log.everything.wet.weight.12)





gam.lm.log.everything.wet.weight.12.unordered<- gam(log(everything.wet.weight) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.everything.wet.weight <- family(gam.lm.log.everything.wet.weight.12)
fam.gam.everything.wet.weight
str(fam.gam.everything.wet.weight)
ilink.gam.everything.wet.weight<- fam.gam.everything.wet.weight$linkinv
ilink.gam.everything.wet.weight


mod.everything.wet.weight<-gam.lm.log.everything.wet.weight.12
ndata.everything.wet.weight <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(everything.wet.weight), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total wet biomass per mesocosm\n(g, log scale)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.everything.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.everything.wet.weight
ggsave("C:Graphs May 2019//everything.wet.weight_pred.png")


# Everything wet weight per 1 ---------------------------------------------

gamma.12.everything.wet.weight.per.1<-fitdistr(food.exp.data.12.2019_zscores$everything.wet.weight.per.1+0.01, "gamma")
qqp(food.exp.data.12.2019$everything.wet.weight.per.1, "gamma", shape = gamma.12.everything.wet.weight.per.1$estimate[[1]], rate = gamma.12.everything.wet.weight.per.1$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$everything.wet.weight.per.1)

qqp(log(food.exp.data.12.2019_zscores$everything.wet.weight.per.1))
#log normal the best

qqp(food.exp.data.12.2019_zscores$everything.wet.weight.per.1, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.everything.wet.weight.per.1<-lm(everything.wet.weight.per.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.everything.wet.weight.per.1<-lm(log(everything.wet.weight.per.1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.everything.wet.weight.per.1<-glm(everything.wet.weight.per.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.loglink.everything.wet.weight.per.1<-glm(everything.wet.weight.per.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.lm.everything.wet.weight.per.1.12<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.everything.wet.weight.per.1.12<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.everything.wet.weight.per.1.12.1, gam.lm.log.everything.wet.weight.per.1.12, gam.tweedie.everything.wet.weight.per.1.12.1, lm.log.everything.wet.weight.per.1, glm.loglink.everything.wet.weight.per.1, lm.everything.wet.weight.per.1, glm.everything.wet.weight.per.1, gam.lm.everything.wet.weight.per.1.12, gam.gamma.everything.wet.weight.per.1.12.1)

#gam.lm.log.everything.wet.weight.per.1.12 is best 


plot(gam.lm.log.everything.wet.weight.per.1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.log.everything.wet.weight.per.1.12)
qq_plot(gam.lm.log.everything.wet.weight.per.1.12, method = 'simulate')
k.check(gam.lm.log.everything.wet.weight.per.1.12)
summary(gam.lm.log.everything.wet.weight.per.1.12)





gam.lm.log.everything.wet.weight.per.1.12.unordered<- gam(log(everything.wet.weight.per.1) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.everything.wet.weight.per.1 <- family(gam.lm.log.everything.wet.weight.per.1.12)
fam.gam.everything.wet.weight.per.1
str(fam.gam.everything.wet.weight.per.1)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(everything.wet.weight.per.1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Log Biomass per 1 % cover \n(wet weight)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.everything.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.everything.wet.weight.per.1
ggsave("C:Graphs May 2019//everything.wet.weight.per.1_pred.png")




# Dry weight --------------------------------------------------------------


gamma.12.total_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$total_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$total_dry_biomass, "gamma", shape = gamma.12.total_dry_biomass$estimate[[1]], rate = gamma.12.total_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$total_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$total_dry_biomass))


qqp(food.exp.data.12.2019_zscores$total_dry_biomass, "lnorm")




lm.total_dry_biomass<-lm(total_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.total_dry_biomass<-lm(log(total_dry_biomass) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.total_dry_biomass<-glm(total_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.loglink.total_dry_biomass<-glm(total_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.lm.total_dry_biomass.12<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.total_dry_biomass.12<- gam(log(total_dry_biomass) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.total_dry_biomass.12.1, gam.lm.log.total_dry_biomass.12, gam.tweedie.total_dry_biomass.12.1, lm.log.total_dry_biomass, glm.loglink.total_dry_biomass, lm.total_dry_biomass, glm.total_dry_biomass, gam.lm.total_dry_biomass.12, gam.gamma.total_dry_biomass.12.1)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.lm.total_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.total_dry_biomass.12 )
qq_plot(gam.lm.total_dry_biomass.12 , method = 'simulate')
k.check(gam.lm.total_dry_biomass.12 )
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (total_dry_biomass), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.total_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.total_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//total_dry_biomass_pred.png")


# Dry weight per 1 % cover --------------------------------------------------------------


gamma.12.total_dry_biomass_per1<-fitdistr(food.exp.data.12.2019_zscores$total_dry_biomass_per1+0.01, "gamma")
qqp(food.exp.data.12.2019$total_dry_biomass_per1, "gamma", shape = gamma.12.total_dry_biomass_per1$estimate[[1]], rate = gamma.12.total_dry_biomass_per1$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$total_dry_biomass_per1)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$total_dry_biomass_per1))


qqp(food.exp.data.12.2019_zscores$total_dry_biomass_per1, "lnorm")




lm.total_dry_biomass_per1<-lm(total_dry_biomass_per1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.total_dry_biomass_per1<-lm(log(total_dry_biomass_per1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.total_dry_biomass_per1<-glm(total_dry_biomass_per1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.loglink.total_dry_biomass_per1<-glm(total_dry_biomass_per1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.lm.total_dry_biomass_per1.12<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.total_dry_biomass_per1.12.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.total_dry_biomass_per1.12<- gam(log(total_dry_biomass_per1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.total_dry_biomass_per1.12<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.total_dry_biomass_per1.12.1<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.total_dry_biomass_per1.12.1, gam.lm.log.total_dry_biomass_per1.12, gam.tweedie.total_dry_biomass_per1.12, lm.log.total_dry_biomass_per1, glm.loglink.total_dry_biomass_per1, lm.total_dry_biomass_per1, glm.total_dry_biomass_per1, gam.lm.total_dry_biomass_per1.12, gam.gamma.total_dry_biomass_per1.12.1)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.tweedie.total_dry_biomass_per1.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.tweedie.total_dry_biomass_per1.12 )
qq_plot(gam.tweedie.total_dry_biomass_per1.12 , method = 'simulate')
k.check(gam.tweedie.total_dry_biomass_per1.12 )
summary(gam.tweedie.total_dry_biomass_per1.12 )





gam.tweedie.total_dry_biomass_per1.12.unordered<- gam(total_dry_biomass_per1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=Food.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
summary(gam.tweedie.total_dry_biomass_per1.12.unordered)

fam.gam.total_dry_biomass_per1 <- family(gam.tweedie.total_dry_biomass_per1.12)
fam.gam.total_dry_biomass_per1
str(fam.gam.total_dry_biomass_per1)
ilink.gam.total_dry_biomass_per1<- fam.gam.total_dry_biomass_per1$linkinv
ilink.gam.total_dry_biomass_per1


mod.total_dry_biomass_per1<-gam.tweedie.total_dry_biomass_per1.12
ndata.total_dry_biomass_per1 <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (total_dry_biomass_per1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.total_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.total_dry_biomass_per1
ggplot2::ggsave("C:Graphs May 2019//total_dry_biomass_per1_pred.png")

# hydroid biomass ---------------------------------------------------------

hist(food.exp.data.12.2019_zscores$hydroid_dry_biomass)

gamma.12.hydroid_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$hydroid_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$hydroid_dry_biomass, "gamma", shape = gamma.12.hydroid_dry_biomass$estimate[[1]], rate = gamma.12.hydroid_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$hydroid_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$hydroid_dry_biomass))


qqp(food.exp.data.12.2019_zscores$hydroid_dry_biomass, "lnorm")




lm.hydroid_dry_biomass<-lm(hydroid_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.hydroid_dry_biomass<-lm(log(hydroid_dry_biomass+0.1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.hydroid_dry_biomass<-glm(hydroid_dry_biomass+0.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.hydroid_dry_biomass.12<- gam(log(hydroid_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.hydroid_dry_biomass.12, gam.lm.log.hydroid_dry_biomass.12, gam.tweedie.hydroid_dry_biomass.12, lm.log.hydroid_dry_biomass, lm.hydroid_dry_biomass, glm.hydroid_dry_biomass, gam.lm.hydroid_dry_biomass.12, gam.gamma.hydroid_dry_biomass.12)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.gamma.hydroid_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.gamma.hydroid_dry_biomass.12 )
qq_plot(gam.gamma.hydroid_dry_biomass.12 , method = 'simulate')
k.check(gam.gamma.hydroid_dry_biomass.12 )
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (hydroid_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Obelia") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.hydroid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.hydroid_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//hydroid_dry_biomass_pred.png")




# botryllus biomass -------------------------------------------------------
food.exp.data.12.2019_zscores$tunicate_dry_biomass[food.exp.data.12.2019_zscores$tunicate_dry_biomass<0]<-0


gamma.12.tunicate_dry_biomass<-fitdistr((food.exp.data.12.2019_zscores$tunicate_dry_biomass+0.01), "gamma")
qqp(food.exp.data.12.2019_zscores$tunicate_dry_biomass, "gamma", shape = gamma.12.tunicate_dry_biomass$estimate[[1]], rate = gamma.12.tunicate_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$tunicate_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$tunicate_dry_biomass))


qqp(food.exp.data.12.2019_zscores$tunicate_dry_biomass, "lnorm")




lm.tunicate_dry_biomass<-lm(tunicate_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.tunicate_dry_biomass<-lm(log(tunicate_dry_biomass+0.1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.tunicate_dry_biomass<-glm(tunicate_dry_biomass+0.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.tunicate_dry_biomass.12<- gam(log(tunicate_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.tunicate_dry_biomass.12, gam.lm.log.tunicate_dry_biomass.12, gam.tweedie.tunicate_dry_biomass.12, lm.log.tunicate_dry_biomass, lm.tunicate_dry_biomass, glm.tunicate_dry_biomass, gam.lm.tunicate_dry_biomass.12, gam.gamma.tunicate_dry_biomass.12)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.gamma.tunicate_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.gamma.tunicate_dry_biomass.12 )
qq_plot(gam.gamma.tunicate_dry_biomass.12 , method = 'simulate')
k.check(gam.gamma.tunicate_dry_biomass.12 )
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (tunicate_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.tunicate_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.tunicate_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//tunicate_dry_biomass_pred.png")



# caprellid biomass -------------------------------------------------------


gamma.12.caprellid_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$caprellid_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$caprellid_dry_biomass, "gamma", shape = gamma.12.caprellid_dry_biomass$estimate[[1]], rate = gamma.12.caprellid_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$caprellid_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$caprellid_dry_biomass))


qqp(food.exp.data.12.2019_zscores$caprellid_dry_biomass, "lnorm")




lm.caprellid_dry_biomass<-lm(caprellid_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.caprellid_dry_biomass<-lm(log(caprellid_dry_biomass+0.1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.caprellid_dry_biomass<-glm(caprellid_dry_biomass+0.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.caprellid_dry_biomass.12<- gam(log(caprellid_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.caprellid_dry_biomass.12, gam.lm.log.caprellid_dry_biomass.12, gam.tweedie.caprellid_dry_biomass.12, lm.log.caprellid_dry_biomass, lm.caprellid_dry_biomass, glm.caprellid_dry_biomass, gam.lm.caprellid_dry_biomass.12, gam.gamma.caprellid_dry_biomass.12)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.gamma.caprellid_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.gamma.caprellid_dry_biomass.12 )
qq_plot(gam.gamma.caprellid_dry_biomass.12 , method = 'simulate')
k.check(gam.gamma.caprellid_dry_biomass.12 )
summary(gam.gamma.caprellid_dry_biomass.12 )

summary(gam.gamma.caprellid_dry_biomass.12.unordered)




gam.gamma.caprellid_dry_biomass.12.unordered<- gam(caprellid_dry_biomass+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.caprellid_dry_biomass <- family(gam.gamma.caprellid_dry_biomass.12)
fam.gam.caprellid_dry_biomass
str(fam.gam.caprellid_dry_biomass)
ilink.gam.caprellid_dry_biomass<- fam.gam.caprellid_dry_biomass$linkinv
ilink.gam.caprellid_dry_biomass


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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (caprellid_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.caprellid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.caprellid_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//caprellid_dry_biomass_pred.png")



# caprellid biomass per invididual -------------------------------------------------------

head(food.exp.data.12.2019_zscores)
head(food.caprellid.data_zscores)
food.exp.data.12.2019_zscores$caprellid_dry_biomass_per1<-food.exp.data.12.2019_zscores$caprellid_dry_biomass/(food.caprellid.data_zscores$total.caprellids)

gamma.12.caprellid_dry_biomass_per1<-fitdistr(food.exp.data.12.2019_zscores$caprellid_dry_biomass_per1, "gamma")
qqp(food.exp.data.12.2019$caprellid_dry_biomass_per1, "gamma", shape = gamma.12.caprellid_dry_biomass_per1$estimate[[1]], rate = gamma.12.caprellid_dry_biomass_per1$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$caprellid_dry_biomass_per1)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$caprellid_dry_biomass_per1))


qqp(food.exp.data.12.2019_zscores$caprellid_dry_biomass_per1, "lnorm")



lm.caprellid_dry_biomass_per1<-lm(caprellid_dry_biomass_per1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.caprellid_dry_biomass_per1<-lm(log(caprellid_dry_biomass_per1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.caprellid_dry_biomass_per1<-glm(caprellid_dry_biomass_per1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.caprellid_dry_biomass_per1.12<- gam(log(caprellid_dry_biomass_per1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.caprellid_dry_biomass_per1.12<- gam(caprellid_dry_biomass_per1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.caprellid_dry_biomass_per1.12, gam.lm.log.caprellid_dry_biomass_per1.12, gam.tweedie.caprellid_dry_biomass_per1.12, lm.log.caprellid_dry_biomass_per1, lm.caprellid_dry_biomass_per1, glm.caprellid_dry_biomass_per1, gam.lm.caprellid_dry_biomass_per1.12, gam.gamma.caprellid_dry_biomass_per1.12)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.gamma.caprellid_dry_biomass_per1.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.gamma.caprellid_dry_biomass_per1.12 )
qq_plot(gam.gamma.caprellid_dry_biomass_per1.12 , method = 'simulate')
k.check(gam.gamma.caprellid_dry_biomass_per1.12 )
summary(gam.gamma.caprellid_dry_biomass_per1.12 )

summary(gam.gamma.caprellid_dry_biomass_per1.12.unordered)




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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (caprellid_dry_biomass_per1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.caprellid_dry_biomass_per1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ylim(0,0.012)
plt.caprellid_dry_biomass_per1
ggplot2::ggsave("C:Graphs May 2019//caprellid_dry_biomass_per1_pred.png")


# rest biomass ------------------------------------------------------------


gamma.12.rest_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$rest_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$rest_dry_biomass, "gamma", shape = gamma.12.rest_dry_biomass$estimate[[1]], rate = gamma.12.rest_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$rest_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$rest_dry_biomass))


qqp(food.exp.data.12.2019_zscores$rest_dry_biomass, "lnorm")




lm.rest_dry_biomass<-lm(rest_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.rest_dry_biomass<-lm(log(rest_dry_biomass+0.1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.rest_dry_biomass<-glm(rest_dry_biomass+0.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.lm.rest_dry_biomass.12<- gam(rest_dry_biomass ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.rest_dry_biomass.12<- gam(log(rest_dry_biomass+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.rest_dry_biomass.12, gam.lm.log.rest_dry_biomass.12, gam.tweedie.rest_dry_biomass.12, lm.log.rest_dry_biomass, lm.rest_dry_biomass, glm.rest_dry_biomass, gam.lm.rest_dry_biomass.12, gam.gamma.rest_dry_biomass.12)

#gam.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.gamma.rest_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.gamma.rest_dry_biomass.12 )
qq_plot(gam.gamma.rest_dry_biomass.12 , method = 'simulate')
k.check(gam.gamma.rest_dry_biomass.12 )
summary(gam.gamma.rest_dry_biomass.12 )

###k check is significant




gam.gamma.rest_dry_biomass.12.unordered<- gam(rest_dry_biomass+0.1 ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
summary(gam.gamma.rest_dry_biomass.12.unordered)

fam.gam.rest_dry_biomass <- family(gam.gamma.rest_dry_biomass.12)
fam.gam.rest_dry_biomass
str(fam.gam.rest_dry_biomass)
ilink.gam.rest_dry_biomass<- fam.gam.rest_dry_biomass$linkinv
ilink.gam.rest_dry_biomass


mod.rest_dry_biomass<-gam.gamma.rest_dry_biomass.12
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (rest_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression("Remaining dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.rest_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.rest_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//rest_dry_biomass_pred.png")




# Mussel wet weight -------------------------------------------------------


gamma.12.Mussel.wet.weight<-fitdistr(food.exp.data.12.2019_zscores$Mussel.wet.weight+0.01, "gamma")
qqp(food.exp.data.12.2019$Mussel.wet.weight, "gamma", shape = gamma.12.Mussel.wet.weight$estimate[[1]], rate = gamma.12.Mussel.wet.weight$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight)

qqp(log(food.exp.data.12.2019_zscores$Mussel.wet.weight+1))
#log normal the best

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight, "lnorm")


#none of these look good 
View(food.exp.data.12.2019_zscores)


lm.Mussel.wet.weight<-lm(Mussel.wet.weight ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.Mussel.wet.weight<-lm(log(Mussel.wet.weight+0.1) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

gam.lm.Mussel.wet.weight.12<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight+0.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.Mussel.wet.weight.12<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.Mussel.wet.weight.12.1, gam.lm.log.Mussel.wet.weight.12, gam.tweedie.Mussel.wet.weight.12.1, lm.log.Mussel.wet.weight, lm.Mussel.wet.weight, gam.lm.Mussel.wet.weight.12, gam.gamma.Mussel.wet.weight.12.1)

#gam.lm.log.Mussel.wet.weight.12 is best 


plot(gam.lm.log.Mussel.wet.weight.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.log.Mussel.wet.weight.12)
qq_plot(gam.lm.log.Mussel.wet.weight.12, method = 'simulate')
k.check(gam.lm.log.Mussel.wet.weight.12)
summary(gam.lm.log.Mussel.wet.weight.12)

#residuals many in a straight line but otherwise good 



gam.lm.log.Mussel.wet.weight.12.unordered<- gam(log(Mussel.wet.weight+0.1) ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.Mussel.wet.weight <- family(gam.lm.log.Mussel.wet.weight.12)
fam.gam.Mussel.wet.weight
str(fam.gam.Mussel.wet.weight)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(Mussel.wet.weight+0.1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Mytilus") ~"wet weight (g, log scale)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.Mussel.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.Mussel.wet.weight
ggsave("C:Graphs May 2019//Mussel.wet.weight_pred.png")






# Mussel wet weight per individual mussel ---------------------------------
#big outlier .... 5 mussels weigh 105.5 g?? so 17 g per mussel? 
View(food.exp.data.12.2019_zscores)

hist(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1)

food.exp.data.12.2019_zscores$Mussel.wet.weight[food.exp.data.12.2019_zscores$Mesocosm==8]<-14.0

food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1<-(food.exp.data.12.2019_zscores$Mussel.wet.weight)/(food.exp.data.12.2019_zscores$mussel_complete+1)

gamma.12.Mussel.wet.weight.per.1<-fitdistr(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1+0.01, "gamma")
qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1, "gamma", shape = gamma.12.Mussel.wet.weight.per.1$estimate[[1]], rate = gamma.12.Mussel.wet.weight.per.1$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1)

qqp(log(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1+1))
#log normal the best

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1, "lnorm")


#none of these look good 
head(food.exp.data.12.2019_zscores)


lm.Mussel.wet.weight.per.1<-lm(Mussel.wet.weight.per.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.Mussel.wet.weight.per.1<-lm(log(Mussel.wet.weight.per.1+0.01) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

#glm.Mussel.wet.weight.per.1<-glm(Mussel.wet.weight.per.1+0.01 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
#glm.loglink.Mussel.wet.weight.per.1<-glm(Mussel.wet.weight.per.1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.lm.Mussel.wet.weight.per.1.12<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.gamma.Mussel.wet.weight.per.1.12<- gam(Mussel.wet.weight.per.1+0.01 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.lm.log.Mussel.wet.weight.per.1.12<- gam(log(Mussel.wet.weight.per.1+0.01) ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.tweedie.Mussel.wet.weight.per.1.12.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.loglink.Mussel.wet.weight.per.1.12.1<- gam(Mussel.wet.weight.per.1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.loglink.Mussel.wet.weight.per.1.12.1, gam.lm.log.Mussel.wet.weight.per.1.12, gam.tweedie.Mussel.wet.weight.per.1.12.1, lm.log.Mussel.wet.weight.per.1, lm.Mussel.wet.weight.per.1, gam.lm.Mussel.wet.weight.per.1.12, gam.gamma.Mussel.wet.weight.per.1.12)

#gamma is best but produces very weird plkot so going with loglink


plot(gam.gamma.Mussel.wet.weight.per.1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.gamma.Mussel.wet.weight.per.1.12)
qq_plot(gam.gamma.Mussel.wet.weight.per.1.12, method = 'simulate')
k.check(gam.gamma.Mussel.wet.weight.per.1.12)
summary(gam.gamma.Mussel.wet.weight.per.1.12)

#don't look great... 


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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = Mussel.wet.weight.per.1, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Biomass per mussel \n(wet weight)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.Mussel.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.Mussel.wet.weight.per.1
ggsave("C:Graphs May 2019//Mussel.wet.weight.per.1_pred.png")




# Hydtobot ----------------------------------------------------------------
food.exp.data.12.2019_zscores$hydtobot[food.exp.data.12.2019_zscores$hydtobot==1]<-0.99
food.exp.data.12.2019_zscores$hydtobot[food.exp.data.12.2019_zscores$hydtobot==0]<-0.01



beta.12<-fitdist(0.01*(food.exp.data.12.2019$hydtobot), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$hydtobot), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

lm.hydtobot<-lm(hydtobot ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.hydtobot<-lm(log(hydtobot+0.01) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.loglink.hydtobot<-glm(hydtobot ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))


gam.lm.hydtobot.12<- gam(log(hydtobot+0.01)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.log.lm.hydtobot.12<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


gam.tweedie.hydtobot.12<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family =tw(), select=TRUE, method="REML")

gam.beta.hydtobot.12<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydtobot.12.1<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydtobot.12.2<- gam(hydtobot~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(gam.tweedie.hydtobot.12 ,gam.beta.hydtobot.12, gam.lm.hydtobot.12, gam.log.lm.hydtobot.12, gam.beta.hydtobot.12.1, gam.beta.hydtobot.12.2, lm.hydtobot, lm.log.hydtobot)
#logit REML and cloglog same ... just logit.... 

### logit is the best 

plot(gam.beta.hydtobot.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.hydtobot.12)
qq_plot(gam.beta.hydtobot.12, method = 'simulate')
k.check(gam.beta.hydtobot.12)
gam.check(gam.tweedie.hydtobot.12)
summary(gam.beta.hydtobot.12)


#not the best qq plot
#resids are in rows

gam.beta.hydtobot.12.unordered<- gam(hydtobot~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")

fam.gam.hydtobot <- family(gam.beta.hydtobot.12)
fam.gam.hydtobot 
str(fam.gam.hydtobot )
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = hydtobot, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "cover ratio"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.hydtobot,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)+coord_cartesian(ylim = c(0, 1)) 
plt.gam.hydtobot 







# Hydtobot by weight ------------------------------------------------------
food.exp.data.12.2019_zscores$hydtobot_dry_biomass[food.exp.data.12.2019_zscores$hydtobot_dry_biomass==1]<-0.99
food.exp.data.12.2019_zscores$hydtobot_dry_biomass[food.exp.data.12.2019_zscores$hydtobot_dry_biomass==0]<-0.01

View(food.exp.data.12.2019_zscores)


beta.12<-fitdist(0.01*(food.exp.data.12.2019$hydtobot_dry_biomass), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$hydtobot_dry_biomass), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

lm.hydtobot_dry_biomass<-lm(hydtobot_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.hydtobot_dry_biomass<-lm(log(hydtobot_dry_biomass+0.01) ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.loglink.hydtobot_dry_biomass<-glm(hydtobot_dry_biomass ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))


gam.lm.hydtobot_dry_biomass.12<- gam(log(hydtobot_dry_biomass+0.01)~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.log.lm.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


gam.tweedie.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family =tw(), select=TRUE, method="REML")

gam.beta.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.12.1<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.hydtobot_dry_biomass.12.2<- gam(hydtobot_dry_biomass~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(gam.tweedie.hydtobot_dry_biomass.12 ,gam.beta.hydtobot_dry_biomass.12, gam.lm.hydtobot_dry_biomass.12, gam.log.lm.hydtobot_dry_biomass.12, gam.beta.hydtobot_dry_biomass.12.1, gam.beta.hydtobot_dry_biomass.12.2, lm.hydtobot_dry_biomass, lm.log.hydtobot_dry_biomass)
#logit REML and cloglog same ... just logit.... 

### logit is the best 

plot(gam.beta.hydtobot_dry_biomass.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.beta.hydtobot_dry_biomass.12)
qq_plot(gam.beta.hydtobot_dry_biomass.12, method = 'simulate')
k.check(gam.beta.hydtobot_dry_biomass.12)
summary(gam.beta.hydtobot_dry_biomass.12)


#not the best qq plot
#resids are in rows

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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = hydtobot_dry_biomass, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "biomass ratio"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.hydtobot_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)#+coord_cartesian(ylim = c(0, 1)) 
plt.gam.hydtobot_dry_biomass 


#mesocosm #8 didn't have either tunicates or hydroids weight - 0% hydroids, 2% tunicates


# CAP1 --------------------------------------------------------------------

qqp(food.exp.data.12.2019_zscores$CAP1)


qqp(food.exp.data.12.2019_zscores$CAP1, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.CAP1<-lm(CAP1 ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

gam.lm.CAP1.12<- gam(CAP1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.loglink.CAP1.12.1<- gam(CAP1 ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( lm.CAP1,  gam.lm.CAP1.12)

#gam.lm.CAP1


plot(gam.lm.CAP1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.CAP1.12)
qq_plot(gam.lm.CAP1.12, method = 'simulate')
k.check(gam.lm.CAP1.12)
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y =(CAP1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Partial-dbRDA axis 1\n(36% of constrained variation)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.CAP1
ggsave("C:Graphs May 2019//CAP1_pred.png")




# Distances ---------------------------------------------------------------

qqp(food.exp.data.12.2019_zscores$distances)


qqp(food.exp.data.12.2019_zscores$distances, "lnorm")

head(food.exp.data.12.2019_zscores)

str(food.exp.data.12.2019_zscores)

lm.distances<-lm(distances ~ min.10.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

gam.lm.distances.12<- gam(distances ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.loglink.distances.12.1<- gam(distances ~ s(min.10.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( lm.distances,  gam.lm.distances.12)

#gam.lm.distances


plot(gam.lm.distances.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.lm.distances.12)
qq_plot(gam.lm.distances.12, method = 'simulate')
k.check(gam.lm.distances.12)
summary(gam.lm.distances.12)





gam.lm.distances.12.unordered<- gam(distances ~ s(min.10.pH)+ Food.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
summary(gam.lm.distances.12.unordered)

fam.gam.distances <- family(gam.lm.distances.12)
fam.gam.distances
str(fam.gam.distances)
ilink.gam.distances<- fam.gam.distances$linkinv
ilink.gam.distances


mod.distances<-gam.lm.distances.12
ndata.distances <- with(food.exp.data.12.2019_zscores, data_frame(min.10.pH = seq(min(min.10.pH), max(min.10.pH),
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
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y =(distances), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Homogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.distances
ggsave("C:Graphs May 2019//distances_pred.png")




# Community plotting ------------------------------------------------------
library(cowplot)
fig.biomass<-plot_grid(plt.total_dry_biomass,plt.everything.wet.weight, plt.hydroid_dry_biomass,
                 plt.tunicate_dry_biomass, plt.Mussel.wet.weight,ncol=5, align='v', 
                 labels=c('(a)', '(b)','(c)', '(d)', '(e)', 
                          label_size=12))

fig.biomass

ggplot2::ggsave("C:For submission//Fig.biomass.png", width=65, height=10, units="cm")

#### revised community fig
fig.3.community<-plot_grid( plt.occupied.space,plt.total_dry_biomass,
                            plt.richness, plt.evenness,
                            plt.CAP1, plt.distances,legend_food, ncol=2, rel_heights = c(1,1,1,.2),
                            align='v', axis='l',
                           labels=c('(a)', '(b)','(c)', '(d)', '(e)', '(f)'))

fig.3.community

ggplot2::ggsave("C:For submission//Fig3.community.png", width=30, height=40, units="cm")



#hyd tobot figs
fig.s1.hydtobot<-plot_grid(plt.gam.hydtobot, plt.gam.hydtobot_dry_biomass, legend_food, 
                           ncol=2, rel_heights = c(1,.2),align='v',axis = 'l', 
                           labels=c('(a)', '(b)', ''))

fig.s1.hydtobot

ggplot2::ggsave("C:For submission//Supplemental//Fig.S1.hydotobot.png", width=26, height=12, units="cm")




# Community level tables --------------------------------------------------


richness.gam<- summary(gam.poisson.richness.12)
richness.gam.unordered<- summary(gam.poisson.richness.12.unordered)

evenness.gam<-summary(gam.lm.evenness.12)
evenness.gam.unordered<-summary(gam.lm.evenness.12.unordered)


occupied.space.gam<- summary(gam.beta.occupied.space.12)
occupied.space.gam.unordered<- summary(gam.beta.occupied.space.12.unordered)

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

rest_dry_biomass.gam.unordered <- summary(gam.gamma.rest_dry_biomass.12.unordered)
rest_dry_biomass.gam <- summary(gam.gamma.rest_dry_biomass.12)

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
  group_rows("Homogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (z)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (z)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:For submission//ptable.community.t.html", self_contained = T)


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
  group_rows("Homogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (Chi-square)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi-square)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:For submission//stable.community.f.html", self_contained = T)

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
  group_rows("Homogeneity of dispersions, normal", 16,18) %>% 
  group_rows("Botryllus to Obelia dominance ratio by space, beta (Chi-square, z)", 19,21) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi-square, z)", 22,24) %>% 
  group_rows("Total wet biomass, normal (log)", 25,27) %>% 
  
  save_kable(file = "C:For submission//pstable.community.html", self_contained = T)








# Old code ----------------------------------------------------------------



summary(glm.binomial.hydroid.12.hydrogen)
profile(glm.binomial.hydroid.12.hydrogen)

resid.ssq <- sum(residuals(glm.binomial.num.nudi.12.hydrogen,type="pearson")^2)  ## sum of squares of Pearson resids
resid.df <- nrow(food.exp.data.12.2019_zscores)-length(coef(glm.binomial.hydroid.12.hydrogen))        ## estimated resid df (N-p)
resid.ssq/resid.df                                ## ratio should be approx 1

plot(glm.binomial.num.nudi.12.hydrogen)

library(descr)
LogRegR2(glm.binomial.num.nudi.12.hydrogen)

library(mgcv)
qq.gam(glm.binomial.num.nudi.12.hydrogen,pch=16)



y.transf.betareg <- function(y){  
  n.obs <- sum(!is.na(y))  
  (y * (n.obs - 1) + 0.5) / n.obs  
} 

beta.hydroid.12.1<- betareg(y.transf.betareg(hydroid/100) ~ min.10.pH*Food.quality, data =  food.exp.data.12.2019_zscores)
summary(beta.hydroid.12.1)
plot(beta.hydroid.12.1)
coeftest(beta.hydroid.12.1)





plot(beta.hydroid.12.1)
plot(beta.hydroid.12.1, which = 6, type = "deviance", sub.caption = "")
plot(beta.hydroid.12.1, which = 4, type = "deviance", sub.caption = "")
plot(beta.hydroid.12.1, which = 2, type = "deviance", sub.caption = "")
plot(beta.hydroid.12.1, which = 5, type = "deviance", sub.caption = "")


beta.hydroid.12.2<- betareg(y.transf.betareg(hydroid/100) ~ min.10.pH*Food.quality, data =  food.exp.data.12.2019_zscores, link="loglog")
summary(beta.hydroid.12.2)
plot(beta.hydroid.12.2)
coeftest(beta.hydroid.12.2)

AIC(beta.hydroid.12.2, beta.hydroid.12.1, beta.hydroid.12.3)

sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"), function(x) logLik(update(beta.hydroid.12.3, link = x)))

beta.hydroid.12.3<- betareg(y.transf.betareg(hydroid/100) ~ hydrogen.concentration*Food.quality, data =  food.exp.data.12.2019_zscores, link="cauchit")


coeftest(beta.hydroid.12.3)


Anova(beta.hydroid.12.3, "3", test.statistic = "Chisq")


beta.hydroid.12.1<- betareg(y.transf.betareg(hydroid/100) ~ hydrogen.concentration*Food.quality , data =  food.exp.data.12.2019_zscores)

sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"), function(x) logLik(update(beta.hydroid.12.1, link = x)))

beta.hydroid.12.2<- betareg(y.transf.betareg(hydroid/100) ~ min.10.pH*Food.quality , data =  food.exp.data.12.2019_zscores)

sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"), function(x) logLik(update(beta.hydroid.12.2, link = x)))

beta.hydroid.12.3<- betareg(y.transf.betareg(hydroid/100) ~ min.10.pH*Food.quality , data =  food.exp.data.12.2019_zscores, link="cauchit")
beta.hydroid.12.4<- betareg(y.transf.betareg(hydroid/100) ~ hydrogen.concentration*Food.quality , data =  food.exp.data.12.2019_zscores,  link="cauchit")


AIC(beta.hydroid.12.2, beta.hydroid.12.1, beta.hydroid.12.3, beta.hydroid.12.4)

#best one: 
beta.hydroid.12.3<- betareg(y.transf.betareg(hydroid/100) ~ min.10.pH*Food.quality , data =  food.exp.data.12.2019_zscores, link="cauchit")



# old code

gamma.12<-fitdistr(food.exp.data.12.2019$stand.biomass, "gamma")
qqp(food.exp.data.12.2019$stand.biomass, "gamma", shape = gamma.12$estimate[[1]], rate = gamma.12$estimate[[2]])

qqp(food.exp.data.12.2019$stand.biomass, "norm")


glm.Gamma.stand.biomass.12<- glm(formula = stand.biomass~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = Gamma)
plot(glm.Gamma.stand.biomass.12)
summary(glm.Gamma.stand.biomass.12)
Anova(glm.Gamma.stand.biomass.12, "3")


#####bottohyd
##this for ANOVA (eat treat to mean of others )
options(contrasts = c("contr.sum", "contr.poly"))


###this for Summary (each treatment to ref level)
options(contrasts = c("contr.treatment", "contr.poly"))
 
qqp(food.exp.data.12.2019$bottohyd_dry, "norm")


lm.bottohyd_dry<-lm(bottohyd_dry~CO2*Food.quality, data=food.exp.data.12.2019)
plot(lm.bottohyd_dry)
summary(lm.bottohyd_dry)
Anova(lm.bottohyd_dry, type="3")

lm.bottohyd_dry<-lm(bottohyd_dry~hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.bottohyd_dry)
summary(lm.bottohyd_dry)
Anova(lm.bottohyd_dry, type="3")
### hydroid_dry







######## 
glm.gamma.hydroid_dry_biomass_per1.12<-glm(formula = (hydroid_dry_biomass_per1+0.0001) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "Gamma")
summary(glm.gamma.hydroid_dry_biomass_per1.12)
plot(glm.gamma.hydroid_dry_biomass_per1.12)
Anova(glm.gamma.hydroid_dry_biomass_per1.12, type="3")
####

glm.gamma.tunicate_dry_biomass_per1.12<-glm(formula = (tunicate_dry_biomass_per1_all+0.0001) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "Gamma")
summary(glm.gamma.tunicate_dry_biomass_per1.12)
plot(glm.gamma.tunicate_dry_biomass_per1.12)
Anova(glm.gamma.tunicate_dry_biomass_per1.12, type="3")

###
glm.gamma.rest_dry_biomass_per1.12<-glm(formula = (rest_dry_biomass_per1+0.0001) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "Gamma")
summary(glm.gamma.rest_dry_biomass_per1.12)
plot(glm.gamma.rest_dry_biomass_per1.12)
Anova(glm.gamma.rest_dry_biomass_per1.12, type="3")


###
lm.caprellid_dry_biomass_per1<-lm(caprellid_dry_biomass_per1~hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.caprellid_dry_biomass_per1)
summary(lm.caprellid_dry_biomass_per1)
Anova(lm.caprellid_dry_biomass_per1, type="3")



options(contrasts = c("contr.sum", "contr.poly"))


####
lm.total_dry_biomass<-lm(total_dry_biomass~hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.total_dry_biomass)
summary(lm.total_dry_biomass)
Anova(lm.total_dry_biomass, type="3")

lm.total_dry_biomass.av<-lm(total_dry_biomass~hydrogen.concentration.av*Food.quality, data=food.exp.data.12.2019)
plot(lm.total_dry_biomass.av)
summary(lm.total_dry_biomass.av)
Anova(lm.total_dry_biomass.av, type="3")

lm.total_dry_biomass_per1<-lm(total_dry_biomass_per1~hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.total_dry_biomass_per1)
summary(lm.total_dry_biomass_per1)
Anova(lm.total_dry_biomass_per1, type="3")

glm.gamma.total_dry_biomass_per1.12<-glm(formula = (total_dry_biomass_per1) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "Gamma")
summary(glm.gamma.total_dry_biomass_per1.12)
plot(glm.gamma.total_dry_biomass_per1.12)
Anova(glm.gamma.total_dry_biomass_per1.12, type="3")

#hydroid binomial

#### by hydrogen concentration

contrasts(food.exp.data.12.2019$Food.quality)=contr.poly(3) 

options(contrasts = c("contr.treatment", "contr.poly"))
options(contrasts = c("contr.sum", "contr.poly"))


fitBinom=fitdist(data=food.exp.data.12.2019$hydroid, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$hydroid, "binom", size = 100, prob = fitBinom$estimate[[1]])


food.exp.data.12.2019$hydrogen.concentration<- 10^(-food.exp.data.12.2019$min.10.pH)

glm.binomial.hydroid.12.hydrogen<- glm(formula = cbind(hydroid, 100-hydroid)~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.hydroid.12.hydrogen)
summary(glm.binomial.hydroid.12.hydrogen)
profile(glm.binomial.hydroid.12.hydrogen)
Anova(glm.binomial.hydroid.12.hydrogen, "3", test.statistic = "LR")

1-(1596.4/1711.3)

food.exp.data.12.2019$hydrogen.concentration.av<- 10^(-food.exp.data.12.2019$av.pH)
glm.binomial.hydroid.12.hydrogen.av<- glm(formula = cbind(hydroid, 100-hydroid)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.hydroid.12.hydrogen.av)
summary(glm.binomial.hydroid.12.hydrogen.av)
profile(glm.binomial.hydroid.12.hydrogen.av)
Anova(glm.binomial.hydroid.12.hydrogen.av, "3", test.statistic = "LR")


########
#formicula binomial
head(food.exp.data.12.2019)

fitBinom=fitdist(data=food.exp.data.12.2019$formicula, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$formicula, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.formicula.12<- glm(formula = cbind(formicula, 100-formicula)~ hydrogen.concentration*Food.quality*Mussel.wet.weight, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.formicula.12)
summary(glm.binomial.formicula.12)
Anova(glm.binomial.formicula.12, "3")

glm.binomial.formicula.12.av<- glm(formula = cbind(formicula, 100-formicula)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.formicula.12.av)
summary(glm.binomial.formicula.12.av)
Anova(glm.binomial.formicula.12.av, "3")


###
#nudi.combined binomial
fitBinom=fitdist(data=food.exp.data.12.2019$nudi.combined, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$nudi.combined, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.nudi.combined.12<- glm(formula = cbind(nudi.combined, 100-nudi.combined)~ hydrogen.concentration*Food.quality*Mussel.wet.weight, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.nudi.combined.12)

summary(glm.binomial.nudi.combined.12)
Anova(glm.binomial.nudi.combined.12, "3")

glm.binomial.nudi.combined.12.av<- glm(formula = cbind(nudi.combined, 100-nudi.combined)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.nudi.combined.12.av)
summary(glm.binomial.nudi.combined.12.av)
Anova(glm.binomial.nudi.combined.12.av, "3")

### formicula beta
food.exp.data.12.2019$formicula.2<-(food.exp.data.12.2019$formicula*0.01)+0.01
y.transf.betareg <- function(y){  
  n.obs <- sum(!is.na(y))  
  (y * (n.obs - 1) + 0.5) / n.obs  
} 

y.transf.betareg(food.exp.data.12.2019$formicula/food.exp.data.12.2019$total) 

beta.12<-fitdist(0.01*(food.exp.data.12.2019$formicula+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$formicula+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

beta.12<-fitdist(y.transf.betareg(food.exp.data.12.2019$formicula/food.exp.data.12.2019$total) , "beta", start=NULL)
qqp(y.transf.betareg(food.exp.data.12.2019$formicula/food.exp.data.12.2019$total) , "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])


beta.formicula.12.1<- betareg(y.transf.betareg(food.exp.data.12.2019$formicula/food.exp.data.12.2019$total)  ~ min.10.pH*Food.quality, data =  food.exp.data.12.2019)
summary(beta.formicula.12.1)
plot(beta.formicula.12.1)
Anova(beta.formicula.12.1)

## I think binomial is better fitat lower distrib and the model looks better ... 

#alive.bot binomial
fitBinom=fitdist(data=food.exp.data.12.2019$alive.bot, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$alive.bot, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.alive.bot.12<- glm(formula = cbind(alive.bot, 100-alive.bot)~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.alive.bot.12)
summary(glm.binomial.alive.bot.12)
Anova(glm.binomial.alive.bot.12, "3")

glm.binomial.alive.bot.12.av<- glm(formula = cbind(alive.bot, 100-alive.bot)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.alive.bot.12.av)
summary(glm.binomial.alive.bot.12.av)
Anova(glm.binomial.alive.bot.12.av, "3")


#beta is also a good fit ... 
beta.12<-fitdist(0.01*(food.exp.data.12.2019$alive.bot+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$alive.bot+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])





#caprellid binomial
fitBinom=fitdist(data=food.exp.data.12.2019$caprellid, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$caprellid, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.caprellid.12<- glm(formula = cbind(caprellid, 100-caprellid)~ hydrogen.concentration*Food.quality*Mussel.wet.weight, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.caprellid.12)
summary(glm.binomial.caprellid.12)
Anova(glm.binomial.caprellid.12, "3")

glm.binomial.caprellid.12.av<- glm(formula = cbind(caprellid, 100-caprellid)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.caprellid.12.av)
summary(glm.binomial.caprellid.12.av)
Anova(glm.binomial.caprellid.12.av, "3")

head(food.exp.data.12.2019)
#alive.mem binomial
fitBinom=fitdist(data=food.exp.data.12.2019$alive.mem, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$alive.mem, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.alive.mem.12<- glm(formula = cbind(alive.mem, 100-alive.mem)~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.alive.mem.12)
summary(glm.binomial.alive.mem.12)
Anova(glm.binomial.alive.mem.12, "3")

glm.binomial.alive.mem.12.av<- glm(formula = cbind(alive.mem, 100-alive.mem)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.alive.mem.12.av)
summary(glm.binomial.alive.mem.12.av)
Anova(glm.binomial.alive.mem.12.av, "3")

#didemnum binomial
fitBinom=fitdist(data=food.exp.data.12.2019$didemnum, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$didemnum, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.didemnum.12<- glm(formula = cbind(didemnum, 100-didemnum)~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.didemnum.12)
summary(glm.binomial.didemnum.12)
Anova(glm.binomial.didemnum.12, "3")

glm.binomial.didemnum.12.av<- glm(formula = cbind(didemnum, 100-didemnum)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.didemnum.12.av)
summary(glm.binomial.didemnum.12.av)
Anova(glm.binomial.didemnum.12.av, "3")

#bowerbankia binomial
fitBinom=fitdist(data=food.exp.data.12.2019$bowerbankia, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$bowerbankia, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.bowerbankia.12<- glm(formula = cbind(bowerbankia, 100-bowerbankia)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.bowerbankia.12)
summary(glm.binomial.bowerbankia.12)
Anova(glm.binomial.bowerbankia.12, "3")

head(food.exp.data.12.2019)


### num.barn.alive
poisson.12<-fitdistr(food.exp.data.12.2019$num.barn.alive, "Poisson")
qqp(food.exp.data.12.2019$num.barn.alive, "pois", poisson.12$estimate)

glm.poisson.num.barn.alive.12<-glm(formula = (num.barn.alive) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.barn.alive.12)
plot(glm.poisson.num.barn.alive.12)
Anova(glm.poisson.num.barn.alive.12, type="III")

glm.poisson.num.barn.alive.12.av<-glm(formula = (num.barn.alive) ~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.barn.alive.12.av)
plot(glm.poisson.num.barn.alive.12.av)
Anova(glm.poisson.num.barn.alive.12.av, type="III")



### num.disporella
poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella, "Poisson")
qqp(food.exp.data.12.2019$num.disporella, "pois", poisson.12$estimate)

glm.poisson.num.disporella.12<-glm(formula = (num.disporella) ~hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.12)
plot(glm.poisson.num.disporella.12)
Anova(glm.poisson.num.disporella.12, type="III")

glm.poisson.num.disporella.12.av<-glm(formula = (num.disporella) ~hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.12.av)
plot(glm.poisson.num.disporella.12.av)
Anova(glm.poisson.num.disporella.12.av, type="III")



### num.serpulid
poisson.12<-fitdistr(food.exp.data.12.2019$num.serpulid, "Poisson")
qqp(food.exp.data.12.2019$num.serpulid, "pois", poisson.12$estimate)

glm.poisson.num.serpulid.12<-glm(formula = (num.serpulid) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.serpulid.12)
plot(glm.poisson.num.serpulid.12)
Anova(glm.poisson.num.serpulid.12, type = "III")

glm.poisson.num.serpulid.12.av<-glm(formula = (num.serpulid) ~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.serpulid.12.av)
plot(glm.poisson.num.serpulid.12.av)
Anova(glm.poisson.num.serpulid.12.av, type = "III")


head(food.exp.data.12.2019)

### num.nudi
poisson.12<-fitdistr(food.exp.data.12.2019$num.nudi, "Poisson")
qqp(food.exp.data.12.2019$num.nudi, "pois", poisson.12$estimate)

fitBinom=fitdist(data=food.exp.data.12.2019$num.nudi, dist="binom", fix.arg=list(size=1), start=list(prob=0.01))
qqp(food.exp.data.12.2019$num.nudi, "binom", size = 1, prob = fitBinom$estimate[[1]])


glm.poisson.num.nudi.12<-glm(formula = (num.nudi) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.nudi.12)
plot(glm.poisson.num.nudi.12)
anova(glm.poisson.num.nudi.12, test = "Chi")



glm.binomial.num.nudi.12<-glm(formula = (num.nudi) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "binomial"(link="cloglog"))
summary(glm.binomial.num.nudi.12)
plot(glm.binomial.num.nudi.12)
anova(glm.binomial.num.nudi.12, test = "Chi")

### num.nudi.eggs
poisson.12<-fitdistr(food.exp.data.12.2019$num.nudi, "Poisson")
qqp(food.exp.data.12.2019$num.nudi, "pois", lambda=poisson.12$estimate[[1]])

nbinom12 <- fitdistr(food.exp.data.12.2019$num.nudi, "Negative Binomial")
qqp(food.exp.data.12.2019$num.nudi, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])


glm.neg.binom.num.nudi.eggs.12<-glm.nb(formula = (num.nudi.eggs) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019)
summary(glm.neg.binom.num.nudi.eggs.12)
plot(glm.neg.binom.num.nudi.eggs.12)
anova(glm.neg.binom.num.nudi.eggs.12, test = "Chi")

glm.poisson.num.nudi.eggs.12<-glm(formula = (num.nudi.eggs) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family="poisson")
summary(glm.poisson.num.nudi.eggs.12)
plot(glm.poisson.num.nudi.eggs.12)
anova(glm.poisson.num.nudi.eggs.12, test = "Chi")

### num.corella
poisson.12<-fitdistr(food.exp.data.12.2019$num.corella, "Poisson")
qqp(food.exp.data.12.2019$num.corella, "pois", poisson.12$estimate)

nbinom12 <- fitdistr(food.exp.data.12.2019$num.corella, "Negative Binomial")
qqp(food.exp.data.12.2019$num.corella, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])

food.exp.data.12.2019$Food.quality<-as.factor(food.exp.data.12.2019$Food.quality)
food.exp.data.12.2019$Food.quality<-relevel(food.exp.data.12.2019$Food.quality, "None")

### richness
library(fitdistrplus)
poisson.12<-fitdistr(food.exp.data.12.2019$richness, "Poisson")
qqp(food.exp.data.12.2019$richness, "pois", poisson.12$estimate)


glm.poisson.richness.12<-glm(formula = (richness) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.richness.12)
plot(glm.poisson.richness.12)
Anova(glm.poisson.richness.12, type="3")

glm.poisson.richness.12.av<-glm(formula = (richness) ~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.richness.12.av)
plot(glm.poisson.richness.12.av)
Anova(glm.poisson.richness.12.av, type="3")

glm.poisson.richness.mesotile.12<-glm(formula = (richness.mesotile) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.richness.mesotile.12)
plot(glm.poisson.richness.mesotile.12)
Anova(glm.poisson.richness.mesotile.12, type="3")

glm.poisson.richness.mesotile.12.av<-glm(formula = (richness.mesotile) ~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.richness.mesotile.12.av)
plot(glm.poisson.richness.mesotile.12.av)
Anova(glm.poisson.richness.mesotile.12.av, type="3")

#use this one for the "summary" part to partition low vs high food. 
options(contrasts = c("contr.treatment", "contr.poly"))

#use this one for overall model
options(contrasts = c("contr.sum", "contr.poly"))



gamma.12<-fitdistr(food.exp.data.12.2019$eveness + 0.01, "gamma")
qqp(food.exp.data.12.2019$eveness, "gamma", shape = gamma.12$estimate[[1]], rate = gamma.12$estimate[[2]])
qqp(food.exp.data.12.2019$eveness, "norm")

#lm is better

lm.eveness<-lm(eveness ~ hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.eveness)
summary(lm.eveness)
Anova(lm.eveness, type="3" )

lm.eveness.av<-lm(eveness ~ hydrogen.concentration.av*Food.quality, data=food.exp.data.12.2019)
plot(lm.eveness.av)
summary(lm.eveness.av)
Anova(lm.eveness.av, type="3" )

head(food.exp.data.12.2019)
#occupied.space binomial
fitBinom=fitdist(data=food.exp.data.12.2019$occupied.space, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$occupied.space, "binom", size = 100, prob = fitBinom$estimate[[1]])

glm.binomial.occupied.space.12<- glm(formula = cbind(occupied.space, 100-occupied.space)~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.occupied.space.12)
summary(glm.binomial.occupied.space.12)
Anova(glm.binomial.occupied.space.12, "3")

glm.binomial.occupied.space.12.av<- glm(formula = cbind(occupied.space, 100-occupied.space)~ hydrogen.concentration.av*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.occupied.space.12.av)
summary(glm.binomial.occupied.space.12.av)
Anova(glm.binomial.occupied.space.12.av, "3")

head(glm.binomial.occupied.space.12)
fitBinom=fitdist(data=food.exp.data.12.2019$second.occ, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$second.occ, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$second.occ+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$second.occ+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

food.exp.data.12.2019$second.occ

food.exp.data.12.2019$second.occ[food.exp.data.12.2019$second.occ==0]<-0.0001


beta.second.occ.12.1<- betareg((second.occ) ~ hydrogen.concentration*Food.quality, data =  food.exp.data.12.2019)
summary(beta.second.occ.12.1)
plot(beta.second.occ.12)
plot(beta.second.occ.12.1, which = 6, type = "deviance", sub.caption = "")
plot(beta.second.occ.12.1, which = 4, type = "deviance", sub.caption = "")
plot(beta.second.occ.12.1, which = 1, type = "deviance", sub.caption = "")
Anova(beta.second.occ.12.1, test="Chi")



glm.binomial.second.occ.12<- glm(formula = cbind(second.occ, occupied.space-second.occ)~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
plot(glm.binomial.second.occ.12)
summary(glm.binomial.second.occ.12)
Anova(glm.binomial.second.occ.12, "3")


cbind(food.exp.data.12.2019$occupied.space, 100-food.exp.data.12.2019$occupied.space)

require(pscl)

glm.poisson.num.corella.12<-glm(formula = (num.corella) ~hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.corella.12)
plot(glm.poisson.num.corella.12)
anova(glm.poisson.num.corella.12, test = "Chi")


glm.poisson.num.corella.12<-zeroinfl(num.corella ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, dist = "negbin")
summary(glm.poisson.num.corella.12)
plot(glm.poisson.num.corella.12)
anova(glm.poisson.num.corella.12, test = "Chi")

glm.neg.binom.num.corella.12<-glm.nb(formula = (num.corella) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019)
summary(glm.neg.binom.num.corella.12)
plot(glm.neg.binom.num.corella.12)
anova(glm.neg.binom.num.corella.12, test = "Chi")


### whatabout pairwise differences between the levels?? 
### is there to testing whether there is any heterogeneity in the means of the levels of the predictor is what ANOVS is doing ... 

library(fitdistrplus)
library
nbinom12 <- fitdistr(food.exp.data.12.2019$hydroid, "Negative Binomial")
qqp(food.exp.data.12.2019$hydroid, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])

glm.neg.binomial.hydroid.12<- glm.nb(hydroid~ min.10.pH*Food.quality, data = food.exp.data.12.2019)
plot(glm.neg.binomial.hydroid.12)
summary(glm.neg.binomial.hydroid.12)
Anova(glm.neg.binomial.hydroid.12, type=3)

library(mgcv)
library(car)
influencePlot(glm.binomial.hydroid.12.1.hydrogen)


qqp(food.exp.data.12.2019$hydroid, "norm")
hydroid.lm<-lm(hydroid~min.10.pH*Food.quality, data = food.exp.data.12.2019)
summary(hydroid.lm)
anova(hydroid.lm, test = "Chi")

####### try a gam
gam.hydroid.12<- gam(0.01*hydroid ~ s(min.10.pH, fx=FALSE, k=-1, bs="cr")*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)
plot(gam.hydroid.12)
summary(gam.hydroid.12)


anova(gam.hydroid.12, glm.binomial.hydroid.12.hydrogen)

#### quasibinomial
glm.quasi.hydroid.12<- glm(0.01*hydroid ~ min.10.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = quasibinomial(link="probit"))
plot(glm.quasi.hydroid.12)
summary(glm.quasi.hydroid.12)
drop1(glm.quasi.hydroid.12, test="F")
AIC(glm.quasi.hydroid.12)

simulationOutput.beta <- simulateResiduals(fittedModel = glm.binomial.hydroid.12.hydrogen)
plotSimulatedResiduals(simulationOutput = simulationOutput.beta)
testZeroInflation(simulationOutput.beta)
testUniformity(simulationOutput = simulationOutput.beta)
testDispersion(simulationOutput = simulationOutput.beta)

#hydroid
head(food.exp.data.12.2019_zscores)

glm.overdispersed.hydroid<- glmer(formula = cbind(hydroid, 100-hydroid)~ min.10.pH*Food.quality + (1|Mesocosm), data = food.exp.data.12.2019_zscores, family =binomial("logit"))
plot(glm.overdispersed.hydroid)
hydroid.Anova<- Anova(glm.binomial.hydroid.12.hydrogen, "3", test.statistic = "Chi")
summary(glm.overdispersed.hydroid)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

cmod_lme4_L<-glm.overdispersed.hydroid
p1 <- plot(cmod_lme4_L,id=0.05,idLabels=~.obs)
p2 <- plot(cmod_lme4_L,ylim=c(-1.5,1),type=c("p","smooth"))
grid.arrange(p1,p2,nrow=1)

plot(cmod_lme4_L,Food.quality~resid(.,type="pearson"),xlim=c(-1.5,1))

dotplot(ranef(cmod_lme4_L,condVar=TRUE))






########
##### Beta
logit.hyd<-fitdistr(food.exp.data.12.2019$hydroid, "logistic")

# beta regression ---------------------------------------------------------


qqp(food.exp.data.12.2019$hydroid, "logit", location=logit.hyd$estimate[[1]], scale=logit.hyd$estimate[[2]])
?qqp
beta.12<-fitdist(0.01*(food.exp.data.12.2019$hydroid+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$hydroid+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])




beta.12<-fitdist(food.exp.data.12.2019$hydroid.2, "beta", start=NULL)
qqp(food.exp.data.12.2019$hydroid.2, "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

beta.12<-fitdist(food.exp.data.12.2019$hydroid.3, "beta", start=NULL)
qqp(food.exp.data.12.2019$hydroid.3, "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

food.exp.data.12.2019$hydroid.2<-0.01*(food.exp.data.12.2019$hydroid+0.01)
food.exp.data.12.2019$hydroid.3<-(food.exp.data.12.2019$hydroid/food.exp.data.12.2019$occupied.space)+0.01

beta.hydroid.12.1<- betareg(hydroid.2 ~ min.10.pH*Food.quality, data =  food.exp.data.12.2019)
summary(beta.hydroid.12.1)
plot(beta.hydroid.12)
plot(beta.hydroid.12, which = 6, type = "deviance", sub.caption = "")
plot(beta.hydroid.12, which = 4, type = "deviance", sub.caption = "")
plot(beta.hydroid.12, which = 1, type = "deviance", sub.caption = "")

gy_logit4 <- update(beta.hydroid.12.1, subset = -3)

coef(beta.hydroid.12, model = "precision")
coef(gy_logit4, model = "precision")


### beta regression with min.10.pH as an additional regressor for the precision parameterto account for heteroskedasticit
beta.hydroid.12.2<- betareg(hydroid.2 ~ min.10.pH*Food.quality|min.10.pH, data =  food.exp.data.12.2019)
beta.hydroid.12.3<- betareg(hydroid.2 ~ min.10.pH*Food.quality|Food.quality, data =  food.exp.data.12.2019)
beta.hydroid.12.4<- betareg(hydroid.2 ~ min.10.pH*Food.quality|Food.quality, data =  food.exp.data.12.2019, link="probit")
beta.hydroid.12.5<- betareg(hydroid.2 ~ min.10.pH*Food.quality|Food.quality, data =  food.exp.data.12.2019, link="loglog")


lrtest(beta.hydroid.12.1, beta.hydroid.12.2)
lrtest(beta.hydroid.12.1, beta.hydroid.12.5)


beta.hydroid.12.space<- betareg(hydroid.3 ~ min.10.pH*Food.quality, data =  food.exp.data.12.2019)
plot(beta.hydroid.12.space, which = 1:5, type = "deviance", sub.caption = "")
summary(beta.hydroid.12.space)
cooks.distance(beta.hydroid.12.space)
AIC(beta.hydroid.12.1)
AIC(beta.hydroid.12.6)


beta.hydroid.12.6<- betareg(hydroid.2 ~ min.10.pH*Food.quality, data =  food.exp.data.12.2019)

devtools::install_github("hilaryparker/explainr")
library(explainr)

summary(beta.hydroid.12.5)
plot(beta.hydroid.12)
plot(beta.hydroid.12, which = 5, type = "deviance", sub.caption = "")


Anova(beta.hydroid.12.space, type=3)
fitBinom=fitdist(data=scorebinom, dist="binom", fix.arg=list(size=8), start=list(prob=0.3))

binom.12<-fitdist(food.exp.data.12.2019$hydroid.2, "binom", fix.arg=list(size=10), start=list(prob=0.3))
qqp(rbinom(food.exp.data$num.prop.disporella.erect), size=0.1)

?fitdist



?betareg
library(betareg)
library(lmtest)

# test of constant precision based on beta distribution:
cp = betareg(min.10.pH~hydroid, food.exp.data.12.2019)
vp = betareg(x~f|f, df)
lrtest(cp, vp)

glm.quasi.hydroid.12<- glm(hydroid.2 ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = quasibinomial(link="probit"))
summary(glm.quasi.hydroid.12)


#### beta binomial from https://cran.r-project.org/web/packages/bbmle/vignettes/mle2.pdf
install.packages("bbmle")
install.packages("emdbook")
library("bbmle")
load(system.file("vignetteData","orob1.rda",package="bbmle"))
summary(orob1)
head(orob1)

set.seed(1001)
library(emdbook)
x1 <- rbetabinom(n=1000,prob=0.1,size=50,theta=10)

mtmp <- function(prob,size,theta) {
  -sum(dbetabinom(x1,prob,size,theta,log=TRUE))
}
m0 <- mle2(mtmp,start=list(prob=0.2,theta=9),data=list(size=50))


ML1 <- function(prob1,prob2,prob3,theta,food.exp.data.12.2019) {
  prob <- c(prob1,prob2,prob3)[as.numeric(food.exp.data.12.2019$hydroid)]
  size <- food.exp.data.12.2019$occupied.space
  -sum(dbetabinom(food.exp.data.12.2019$min.10.pH,prob,size,theta,log=TRUE))
}

m1 <- mle2(ML1,start=list(prob1=0.5,prob2=0.5,prob3=0.5,theta=1),
            data=list(food.exp.data.12.2019=food.exp.data.12.2019))




ML1 <- function(prob1,prob2,prob3,theta,x) {
  prob <- c(prob1,prob2,prob3)[as.numeric(x$dilution)]
  size <- x$n
  -sum(dbetabinom(x$m,prob,size,theta,log=TRUE))
}







#Disporella percent
head(food.exp.data.12.2019)
food.exp.data.12.2019$disporella<-(food.exp.data.12.2019$disporella)*0.01
glm.binomial.disporella.12<- glm(formula = cbind(disporella, 1.0001-disporella)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.disporella.12)
summary(glm.binomial.disporella.12)
anova(glm.binomial.disporella.12, test = "Chi")


### Number Disporella
nbinom12 <- fitdistr(food.exp.data.12.2019$num.disporella, "Negative Binomial")
qqp(food.exp.data.12.2019$num.disporella, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])


poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella, "Poisson")
qqp(food.exp.data.12.2019$num.disporella, "pois", poisson.12$estimate)

glm.poisson.num.disporella.12<-glm(formula = (num.disporella) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.12)
plot(glm.poisson.num.disporella.12)
anova(glm.poisson.num.disporella.12, test = "Chi")


glm.nb.num.disporella.12<-glm.nb(formula = (num.disporella) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019)
summary(glm.nb.num.disporella.12)
anova(glm.nb.num.disporella.12, test = "Chi")


visreg(glm.poisson.num.disporella.12, xvar = "min.10.pH", by= "Food.quality", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))


### Proportion erect. all
food.exp.data.12.2019$num.prop.disporella.erect<-(food.exp.data.12.2019$num.disporella.erect.all)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.erect.12<- glm(formula = cbind(num.prop.disporella.erect, 1-num.prop.disporella.erect)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.12)

summary(glm.binomial.num.prop.disporella.erect.12)
anova(glm.binomial.num.prop.disporella.erect.12, test = "Chi")



### Proportion erect.knobs
food.exp.data.12.2019$num.prop.disporella.erect.knobs<-(food.exp.data.12.2019$num.disporella.erect.knobs)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.erect.knobs.12<- glm(formula = cbind(num.prop.disporella.erect.knobs, 1-num.prop.disporella.erect.knobs)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.knobs.12)

summary(glm.binomial.num.prop.disporella.erect.knobs.12)
anova(glm.binomial.num.prop.disporella.erect.knobs.12, test = "Chi")

### Proportion flat. all
food.exp.data.12.2019$num.prop.disporella.flat<-(food.exp.data.12.2019$num.disporella.flat)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.flat.12<- glm(formula = cbind(num.prop.disporella.flat, 1-num.prop.disporella.flat)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.flat.12)

summary(glm.binomial.num.prop.disporella.flat.12)
anova(glm.binomial.num.prop.disporella.flat.12, test = "Chi")

### Proportion ruffles. all
food.exp.data.12.2019$num.prop.disporella.ruffles<-(food.exp.data.12.2019$num.disporella.erect.ruffles)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.ruffles.12<- glm(formula = cbind(num.prop.disporella.ruffles, 1-num.prop.disporella.ruffles)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.ruffles.12)

summary(glm.binomial.num.prop.disporella.ruffles.12)
anova(glm.binomial.num.prop.disporella.ruffles.12, test = "Chi")


### Proportion fan. all
food.exp.data.12.2019$num.prop.disporella.fan<-(food.exp.data.12.2019$num.disporella.erect.fan)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.fan.12<- glm(formula = cbind(num.prop.disporella.fan, 1-num.prop.disporella.fan)~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.fan.12)

summary(glm.binomial.num.prop.disporella.fan.12)
anova(glm.binomial.num.prop.disporella.fan.12, test = "Chi")

#### Num. erect fan

poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella.erect.fan, "Poisson")
qqp(food.exp.data.12.2019$num.disporella.erect.fan, "pois", poisson.12$estimate)

glm.poisson.num.disporella.erect.fan.12<-glm(formula = (num.disporella.erect.fan) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.erect.fan.12)
plot(glm.poisson.num.disporella.erect.fan.12)
anova(glm.poisson.num.disporella.erect.fan.12, test = "Chi")

#### Num ruffled
poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella.erect.ruffles, "Poisson")
qqp(food.exp.data.12.2019$num.disporella.erect.ruffles, "pois", poisson.12$estimate)

glm.poisson.num.disporella.erect.ruffles.12<-glm(formula = (num.disporella.erect.ruffles) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.erect.ruffles.12)
plot(glm.poisson.num.disporella.erect.ruffles.12)
anova(glm.poisson.num.disporella.erect.ruffles.12, test = "Chi")


###Num knobs
poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella.erect.knobs, "Poisson")
qqp(food.exp.data.12.2019$num.disporella.erect.knobs, "pois", poisson.12$estimate)

glm.poisson.num.disporella.erect.knobs.12<-glm(formula = (num.disporella.erect.knobs) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.erect.knobs.12)
plot(glm.poisson.num.disporella.erect.knobs.12)
anova(glm.poisson.num.disporella.erect.knobs.12, test = "Chi")

### Num didemnum
poisson.12<-fitdistr(food.exp.data.12.2019$num.didemnum, "Poisson")
qqp(food.exp.data.12.2019$num.didemnum, "pois", poisson.12$estimate)

glm.poisson.num.didemnum.12<-glm(formula = (num.didemnum) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.didemnum.12)
plot(glm.poisson.num.didemnum.12)
anova(glm.poisson.num.didemnum.12, test = "Chi")



####################### NORMAL
library(lme4)


qqp(log(food.exp.data.12.2019$richness+1), "norm")

qqp(food.exp.data.12.2019$richness, "lnorm")

qqp(sqrt(food.exp.data.12.2019$richness), "norm")

lm.hydroid<-lm(formula = log(hydroid+1) ~ min.10.pH*Food.quality, data = food.exp.data.12.2019)
summary(lm.hydroid)


summary(lm.richness)
plot(lm.richness)
anova(lm.richness, test = "F")


lognormal.12<-fitdistr(food.exp.data.12.2019$richness+0.01, "lognormal")
qqp(food.exp.data.12.2019$richness, "lnorm", lognormal.12$estimate)


glm.gamma.richness.12<- glm(formula = richness ~ min.10.pH*Food.quality, data = food.exp.data.12.2019, family = "lognormal")


qqp(food.exp.data.12.2019$richness, "lnorm")
?glm
?qqp

visreg(lm.richness, xvar = "min.10.pH", by= "Food.quality", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))




###might need gamlss pacage to do exponential functions
