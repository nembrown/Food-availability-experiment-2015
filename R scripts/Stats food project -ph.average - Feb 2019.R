# read in useful packages
setwd("C:/Users/Norah/Dropbox/Projects/Summer 2015 food availability experiments/Data")




#library(viridis)
library(bbmle) 
library(glmmTMB)
library(doBy)
library(plyr)
library(dplyr)
library(ggplot2) 
library(doBy)
library(grid)
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
library(scales)
library(car)
library(knitr)
library(tidyverse)
library(kableExtra)
library(multcomp)
library(arm) ## for sim()
library(descr)  ## for LogRegR2
library(reshape2)
library(colorspace)
library(grid)
library(DHARMa)
library(gap)
library(qrnn)
library(mgcv)
library(magrittr)
devtools::install_github("haozhu233/kableExtra")


# dataframe management ----------------------------------------------------

head(food.exp.data.12.2019)

food.exp.data.12.2019<-read.csv("C:Mesocosm inventory data//food.exp.data.mesocosm.12.csv")
food.exp.data.tile.all<-read.csv("C:Mesocosm inventory data//food.exp.data.tile.all.csv")
food.caprellid.data<-read.csv("c:Emily caprellid data.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )

food.exp.data.12.2019$oFood.quality<-factor(food.exp.data.12.2019$Food.quality, levels=c("None", "Low", "High"), ordered=TRUE)
food.exp.data.12.2019$Food.quality<-factor(food.exp.data.12.2019$Food.quality, levels=c("None", "Low", "High"), ordered=FALSE)

food.exp.data.12.2019$occupied.space<-100-food.exp.data.12.2019$bare
food.exp.data.12.2019$occupied.space.001<-0.01*(food.exp.data.12.2019$occupied.space)

food.exp.data.12.2019$hydroid<-food.exp.data.tile.all$hydroid
food.exp.data.12.2019$alive.bot<-food.exp.data.tile.all$alive.bot
food.exp.data.12.2019$formicula<-food.exp.data.tile.all$formicula
food.exp.data.12.2019$alive.mem<-food.exp.data.tile.all$alive.mem
food.exp.data.12.2019$didemnum<-food.exp.data.tile.all$didemnum
food.exp.data.12.2019$total<-food.exp.data.tile.all$total
food.exp.data.12.2019$bare<-food.exp.data.tile.all$bare
food.exp.data.12.2019$everything.wet.weight<-food.exp.data.tile.all$everything.wet.weight
food.exp.data.12.2019$everything.wet.weight.per.1<-(food.exp.data.tile.all$everything.wet.weight)/food.exp.data.12.2019$occupied.space
food.exp.data.12.2019$Mussel.wet.weight<-food.exp.data.tile.all$Mussel.wet.weight
food.exp.data.12.2019$total_dry_biomass<-food.exp.data.tile.all$total_dry_biomass
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
food.caprellid.data_zscores$av.pH<-scale(food.exp.data.12.2019$av.pH, center=TRUE, scale=TRUE)
food.caprellid.data_zscores$oFood.quality<-factor(food.caprellid.data_zscores$Food.quality, levels=c("None", "Low", "High"), ordered=TRUE)
food.caprellid.data_zscores$Food.quality<-factor(food.caprellid.data_zscores$Food.quality, levels=c("None", "Low", "High"), ordered=FALSE)


#attributes(food.exp.data.12.2019_zscores$av.pH)
#You can use the attributes to unscale:
  

# Notes on contrasts ------------------------------------------------------


#Notes on ordered factors: my factor food.quality is ordered
#use this one for overall model
options(contrasts = c("contr.sum", "contr.poly"))
#The contrast function, contr.sum(), gives orthogonal contrasts where you compare every level to the overall mean.

#use this one for the "summary" part to partition low vs high food compared to none - none has to be the first level ... 
options(contrasts = c("contr.treatment", "contr.poly"))



# #visualizing data -------------------------------------------------------

cbbPalette.all<- c( "#F8766D", "#00BA38", "#619CFF", "#F8766D", "#00BA38", "#619CFF")

plot.hydroid.12.hydrogen<- ggplot(food.exp.data.tile.all, aes(x=av.pH, y=(hydroid/total), colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality, weight=total), alpha=0.15,size=1.5, method="glm", method.args = list(family = "binomial")) +scale_shape_manual(values=c(19,17))
plot.hydroid.12.hydrogen<- plot.hydroid.12.hydrogen + theme_bw() +  xlab(expression("Average pH")) + ylab(expression(atop(italic("Obelia")~"abundance","(proportion cover)"))) + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))
plot.hydroid.12.hydrogen<- plot.hydroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.hydroid.12.hydrogen<- plot.hydroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.hydroid.12.hydrogen

#Binomial: (doesn't work for hydroid) 
fitBinom=fitdist(data=food.exp.data.12.2019$hydroid, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$hydroid, "binom", size = 100, prob = fitBinom$estimate[[1]])

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



# Model building ----------------------------------------------------------


#Binomial
glm.binomial.hydroid.12.hydrogen<- glm(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial(link="logit"))
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


quasi.hydroid.1<- glm(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = quasibinomial)
summary(quasi.hydroid.1)


quasi.hydroid.2<- glm(formula = cbind(hydroid, 100-hydroid)~ av.pH+Food.quality, data = food.exp.data.12.2019_zscores, family = quasibinomial)
summary(quasi.hydroid.2)

anova(quasi.hydroid.1, quasi.hydroid.2, test="F")

#But the assumptions are off here according to Bolker? 


# GLMM: Observation-level random effect added -----------------------------

theme_set(theme_bw()+theme(panel.spacing=grid::unit(0,"lines")))

#let's try observation-level random effects
#glmmadmb is no longer the best one - Bolker suggests glmmTMB, using Template Model Builder

?glmmTMB


glmm.av.binomial.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality+(1|Mesocosm), data = food.exp.data.12.2019_zscores, family = binomial)
summary(glmm.av.binomial.hydroid.12)

glmm.av.betabinomial.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality+(1|Mesocosm), data = food.exp.data.12.2019_zscores, family = betabinomial)
summary(glmm.av.betabinomial.hydroid.12)

glmm.av.betabinomial.nr.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = betabinomial)
summary(glmm.av.betabinomial.nr.hydroid.12)

glmm.av.betabinomial.reml.hydroid.12<- glmmTMB(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = betabinomial, REML=TRUE)
summary(glmm.av.betabinomial.reml.hydroid.12)

AICtab(glmm.av.binomial.hydroid.12, glmm.av.betabinomial.hydroid.12, glmm.av.betabinomial.nr.hydroid.12, glmm.av.betabinomial.reml.hydroid.12)
#betabinomial best 


#glmmTMB has the capability to simulate from a fitted model. These simulations resample random effects from their estimated distribution.
simo=simulate(glmm.av.betabinomial.hydroid.12, seed=1)
Simdat=food.exp.data.12.2019_zscores
Simdat$hydroid=simo[[1]]

Simdat=transform(Simdat,
                 type="simulated")
food.exp.data.12.2019_zscores$type = "observed"  
Dat=rbind(food.exp.data.12.2019_zscores, Simdat) 

#Then we can plot the simulated data against the observed data to check if they are similar.

ggplot(Dat,  aes(hydroid, colour=type))+geom_density()+facet_grid(food.exp.data.12.2019_zscores$Food.quality)

#For low food not very similar

#glmm.av.betabinomial.nr.hydroid.12 -> not great, but glmm.av.betabinomial.hydroid.12 a bit better? 

#before evaulating inference - need to check diagnostics - search for glmmtmb diagnostics? 

#ggplot(data = NULL) + geom_point(aes(y = residuals(glmm.av.betabinomial.nr.hydroid.12,type = "pearson", scaled = TRUE), x = fitted(glmm.av.betabinomial.nr.hydroid.12)))
#cant' yet do fitted objects from betabinomial

#dharma

simulationOutput <- simulateResiduals(fittedModel = glmm.av.betabinomial.hydroid.12)
plotSimulatedResiduals(simulationOutput = simulationOutput)
#nr. not a perfect fit though - lines should fit  but one is curved and the other two are flat. QQ plot residuals are pretty good... 
#next one down is not great either but a bit better? 

testZeroInflation(simulationOutput)
plotSimulatedResiduals(simulationOutput = simulationOutput, quantreg = T)
#so our model is not great still.... 



augDat <- data.frame(food.exp.data.12.2019_zscores,resid=residuals(glmm.av.betabinomial.nr.hydroid.12,type="pearson"),
                     fitted=fitted(glmm.av.betabinomial.nr.hydroid.12))
ggplot(augDat,aes(x=Food.quality,y=resid))+geom_boxplot()+coord_flip()
#doesn't work for glmmTMB



# GAM ---------------------------------------------------------------------
#Gavin Simpson's gratia
library(devtools)
devtools::install_github('gavinsimpson/gratia')
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



# GAM beta hydroid / gam.av.beta.hydroid.12.3 --------------------------------------------------------


gam.av.binomial.hydroid.12<- gam(formula = cbind(hydroid, 100-hydroid)~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")


food.exp.data.12.2019_zscores$hydroid.001<-0.01*(food.exp.data.12.2019_zscores$hydroid+0.01)

gam.av.beta.hydroid.12<- gam(hydroid.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.hydroid.12.1<- gam(hydroid.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.hydroid.12.2<- gam(hydroid.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.av.beta.hydroid.12.3<- gam(hydroid.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(gam.av.beta.hydroid.12, gam.av.beta.hydroid.12.1, gam.av.beta.hydroid.12.3, gam.av.beta.hydroid.12.2, gam.av.binomial.hydroid.12)
#cauchit the best

plot(gam.av.beta.hydroid.12.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.hydroid.12.3)
qq_plot(gam.av.beta.hydroid.12.3, method = 'simulate')
k.check(gam.av.beta.hydroid.12.3)

summary(gam.av.beta.hydroid.12.3)

gam.av.beta.hydroid.12.3.unordered<- gam(hydroid.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")
summary(gam.av.beta.hydroid.12.3.unordered)


#Confidence intervals
#https://www.fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/

fam.gam.av.hydroid <- family(gam.av.beta.hydroid.12.3)
fam.gam.av.hydroid 
str(fam.gam.av.hydroid )
ilink.gam.av.hydroid <- fam.gam.av.hydroid$linkinv
ilink.gam.av.hydroid


head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)
            
mod.av.hydroid<-gam.av.beta.hydroid.12.3
ndata.av.hydroid <- with(food.exp.data.12.2019_zscores, data_frame(av.pH= seq(min(av.pH), max(av.pH),
                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.hydroid <- add_column(ndata.av.hydroid, fit = predict(mod.av.hydroid, newdata = ndata.av.hydroid, type = 'response'))


ndata.av.hydroid <- bind_cols(ndata.av.hydroid, setNames(as_tibble(predict(mod.av.hydroid, ndata.av.hydroid, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.hydroid <- mutate(ndata.av.hydroid,
                fit_resp  = ilink.gam.av.hydroid(fit_link),
                right_upr = ilink.gam.av.hydroid(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.av.hydroid(fit_link - (2 * se_link)))




ndata.av.hydroid$av.pH.unscaled<-ndata.av.hydroid$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

colorset2 = c("High"="#F8A02E" ,"Low"="#439E5F","None"= "#666666")

# plot 

plt.av.gam.av.hydroid <- ggplot(ndata.av.hydroid, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = hydroid.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Obelia")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.hydroid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.gam.av.hydroid 
ggsave("C:Graphs May 2019//hydroid_pred_av.png")






# GAM beta botryllus / gam.av.beta.alive.bot.12.3 ----------------------------------------------------

#binomial first
gam.av.binomial.alive.bot.12<- gam(formula = cbind(alive.bot, 100-alive.bot)~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, REML=TRUE)

food.exp.data.12.2019_zscores$alive.bot.001<-0.01*(food.exp.data.12.2019_zscores$alive.bot+0.01)

#beta next
gam.av.beta.alive.bot.12<- gam(alive.bot.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, REML=TRUE)
gam.av.beta.alive.bot.12.1<- gam(alive.bot.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, REML=TRUE)
gam.av.beta.alive.bot.12.2<- gam(alive.bot.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, REML=TRUE)
gam.av.beta.alive.bot.12.3<- gam(alive.bot.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


AICtab(gam.av.beta.alive.bot.12, gam.av.beta.alive.bot.12.3, gam.av.beta.alive.bot.12.1, gam.av.beta.alive.bot.12.2, gam.av.binomial.alive.bot.12)
#cauchit the best

plot(gam.av.beta.alive.bot.12.3, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.alive.bot.12.3)
qq_plot(gam.av.beta.alive.bot.12.3, method = 'simulate')
k.check(gam.av.beta.alive.bot.12.3)
summary(gam.av.beta.alive.bot.12.3)

gam.av.beta.alive.bot.12.3.unordered<- gam(alive.bot.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, REML=TRUE)


fam.gam.av.alive.bot <- family(gam.av.beta.alive.bot.12.3)
fam.gam.av.alive.bot
str(fam.gam.av.alive.bot)
ilink.gam.av.alive.bot<- fam.gam.av.alive.bot$linkinv
ilink.gam.av.alive.bot


mod.av.alive.bot<-gam.av.beta.alive.bot.12.3
ndata.av.alive.bot <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                        length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.av.alive.botel for the new data
ndata.av.alive.bot <- add_column(ndata.av.alive.bot, fit = predict(mod.av.alive.bot, newdata = ndata.av.alive.bot, type = 'response'))


ndata.av.alive.bot <- bind_cols(ndata.av.alive.bot, setNames(as_tibble(predict(mod.av.alive.bot, ndata.av.alive.bot, se.fit = TRUE)[1:2]),
                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.alive.bot <- mutate(ndata.av.alive.bot,
                fit_resp  = ilink.gam.av.alive.bot(fit_link),
                right_upr = ilink.gam.av.alive.bot(fit_link + (2 * se_link)),
                right_lwr = ilink.gam.av.alive.bot(fit_link - (2 * se_link)))


food.exp.data.12.2019_zscores$av.pH.unscaled<-food.exp.data.12.2019_zscores$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

ndata.av.alive.bot$av.pH.unscaled<-ndata.av.alive.bot$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.alive.bot <- ggplot(ndata.av.alive.bot, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = alive.bot.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Botryllus")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.alive.bot,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.alive.bot
ggsave("C:Graphs May 2019//alive.bot_pred_av.png")



# GAM negbin caprellid / gam.av.nb.caprellid.12 -----------------------------------------------------------
poisson.12<-fitdistr(food.caprellid.data_zscores$total.caprellids, "Poisson")
qqp(food.caprellid.data_zscores$total.caprellids, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.caprellid <- fitdistr(food.caprellid.data_zscores$total.caprellids, "Negative Binomial")
qqp(food.caprellid.data_zscores$total.caprellids, "nbinom", size = nbinom12.caprellid$estimate[[1]], mu = nbinom12.caprellid$estimate[[2]])


glm.poisson.total.caprellids.12.hydrogen<- glm(total.caprellids~ av.pH*oFood.quality, data = food.caprellid.data_zscores, family="poisson")
#plot(glm.poisson.total.caprellids.12.hydrogen)
overdisp_fun(glm.poisson.total.caprellids.12.hydrogen)
#definitely overdispersed --> neg binomial way better

glm.nb.total.caprellids.12.hydrogen<- glm.nb(total.caprellids~ av.pH*oFood.quality, data = food.caprellid.data_zscores)
#plot(glm.nb.total.caprellids.12.hydrogen)
#This one is actually okay..... 

gam.av.nb.caprellid.12<- gam(total.caprellids ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = negbin(nbinom12.caprellid$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.caprellid.12.1<- gam(total.caprellids ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = nb(), select=TRUE, method="REML")


AICtab(gam.av.nb.caprellid.12.1,gam.av.nb.caprellid.12, glm.nb.total.caprellids.12.hydrogen, glm.poisson.total.caprellids.12.hydrogen)
#gam is better - need to make sure it's reading the correct theta ... 
#also make sure av.pH is reading as a matrix? 

appraise(gam.av.nb.caprellid.12)
qq_plot(gam.av.nb.caprellid.12, method = 'simulate')
plot(gam.av.nb.caprellid.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.av.nb.caprellid.12)
summary(gam.av.nb.caprellid.12)

gam.av.nb.caprellid.12.unordered<- gam(total.caprellids ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.caprellid.data_zscores, family = negbin(nbinom12.caprellid$estimate[[1]]), select=TRUE, method="REML")


fam.gam.av.caprellid <- family(gam.av.nb.caprellid.12)
fam.gam.av.caprellid
str(fam.gam.av.caprellid)
ilink.gam.av.caprellid<- fam.gam.av.caprellid$linkinv
ilink.gam.av.caprellid

mod.av.caprellid<-gam.av.nb.caprellid.12
ndata.av.caprellid <- with(food.caprellid.data_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.av.caprellid

str(food.caprellid.data_zscores)
str(ndata.av.caprellid)

## add the fitted values by predicting from the model for the new data
ndata.av.caprellid <- add_column(ndata.av.caprellid, fit = predict(mod.av.caprellid, newdata = ndata.av.caprellid, type = 'response'))

predict(mod.av.caprellid, newdata = ndata.av.caprellid, type = 'response')
ndata.av.caprellid <- bind_cols(ndata.av.caprellid, setNames(as_tibble(predict(mod.av.caprellid, ndata.av.caprellid, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.caprellid <- mutate(ndata.av.caprellid,
                          fit_resp  = ilink.gam.av.caprellid(fit_link),
                          right_upr = ilink.gam.av.caprellid(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.av.caprellid(fit_link - (2 * se_link)))


food.caprellid.data_zscores$av.pH.unscaled<-food.caprellid.data_zscores$av.pH * attr(food.caprellid.data_zscores$av.pH, 'scaled:scale') + attr(food.caprellid.data_zscores$av.pH, 'scaled:center')
ndata.av.caprellid$av.pH.unscaled<-ndata.av.caprellid$av.pH * attr(food.caprellid.data_zscores$av.pH, 'scaled:scale') + attr(food.caprellid.data_zscores$av.pH, 'scaled:center')

attributes(food.caprellid.data_zscores)
#unscaled is not working!!!! WHYYY



# plot 

plt.av.caprellid <- ggplot(ndata.av.caprellid, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = total.caprellids, shape=CO2, colour=oFood.quality), size=3, data = food.caprellid.data_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Caprella")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.caprellid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.caprellid
ggsave("C:Graphs May 2019//caprellid_pred_av.png")



# GAM caprellid percent -----------------------------------------------------------

glm.binomial.caprellid.av.percent.12.hydrogen<-glm(formula = cbind(caprellid.percent, 100-caprellid.percent)~ (av.pH)*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)

gam.binomial.caprellid.av.percent.12<- gam(formula = cbind(caprellid.percent, 100-caprellid.percent)~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$caprellid.percent.001<-0.01*(food.exp.data.12.2019_zscores$caprellid.percent+0.01)

#beta next
gam.beta.caprellid.av.percent.12<- gam(caprellid.percent.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.beta.caprellid.av.percent.12.1<- gam(caprellid.percent.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.beta.caprellid.av.percent.12.2<- gam(caprellid.percent.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.beta.caprellid.av.percent.12.3<- gam(caprellid.percent.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.beta.caprellid.av.percent.12, gam.beta.caprellid.av.percent.12.1, gam.beta.caprellid.av.percent.12.2,gam.binomial.caprellid.av.percent.12, glm.binomial.caprellid.av.percent.12.hydrogen)
#cloglog is the best along with logit - go with simpler logit


plot(gam.beta.caprellid.av.percent.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

appraise(gam.beta.caprellid.av.percent.12)
qq_plot(gam.beta.caprellid.av.percent.12, method = 'simulate')
k.check(gam.beta.caprellid.av.percent.12)
summary(gam.beta.caprellid.av.percent.12.unordered)

gam.beta.caprellid.av.percent.12.unordered<- gam(caprellid.percent.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.caprellid.av.percent <- family(gam.beta.caprellid.av.percent.12)
fam.gam.caprellid.av.percent
str(fam.gam.caprellid.av.percent)
ilink.gam.caprellid.av.percent<- fam.gam.caprellid.av.percent$linkinv
ilink.gam.caprellid.av.percent


mod.caprellid.av.percent<-gam.beta.caprellid.av.percent.12
ndata.caprellid.av.percent <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.caprellid.av.percentel for the new data
ndata.caprellid.av.percent <- add_column(ndata.caprellid.av.percent, fit = predict(mod.caprellid.av.percent, newdata = ndata.caprellid.av.percent, type = 'response'))


ndata.caprellid.av.percent <- bind_cols(ndata.caprellid.av.percent, setNames(as_tibble(predict(mod.caprellid.av.percent, ndata.caprellid.av.percent, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.caprellid.av.percent <- mutate(ndata.caprellid.av.percent,
                                  fit_resp  = ilink.gam.caprellid.av.percent(fit_link),
                                  right_upr = ilink.gam.caprellid.av.percent(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.caprellid.av.percent(fit_link - (2 * se_link)))


ndata.caprellid.av.percent$av.pH.unscaled<-ndata.caprellid.av.percent$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.caprellid.av.percent <- ggplot(ndata.caprellid.av.percent, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = caprellid.percent.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Caprella")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.caprellid.av.percent,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.caprellid.av.percent
ggsave("C:Graphs May 2019//caprellid.av.percent_pred.png")


# GAM beta formicula / gam.av.beta.formicula.12 -----------------------------------------------------------
str(food.exp.data.12.2019_zscores)
par(op)

fitBinom=fitdist(data=food.exp.data.12.2019$formicula, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$formicula, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$formicula+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$formicula+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#not great at the upper end


#binomial first

glm.binomial.formicula.12.hydrogen<-glm(formula = cbind(formicula, 100-formicula)~ (av.pH)*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)

gam.av.binomial.formicula.12<- gam(formula = cbind(formicula, 100-formicula)~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$formicula.001<-0.01*(food.exp.data.12.2019_zscores$formicula+0.01)

#beta next
gam.av.beta.formicula.12<- gam(formicula.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.formicula.12.1<- gam(formicula.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.formicula.12.2<- gam(formicula.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.av.beta.formicula.12.3<- gam(formicula.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.av.beta.formicula.12, gam.av.beta.formicula.12.3, gam.av.beta.formicula.12.1, gam.av.beta.formicula.12.2,gam.av.binomial.formicula.12, glm.binomial.formicula.12.hydrogen)
#cloglog is the best along with logit - go with simpler logit


plot(gam.av.beta.formicula.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)

appraise(gam.av.beta.formicula.12)
qq_plot(gam.av.beta.formicula.12, method = 'simulate')
k.check(gam.av.beta.formicula.12)
summary(gam.av.beta.formicula.12)


gam.av.beta.formicula.12.unordered<- gam(formicula.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")



fam.gam.av.formicula <- family(gam.av.beta.formicula.12)
fam.gam.av.formicula
str(fam.gam.av.formicula)
ilink.gam.av.formicula<- fam.gam.av.formicula$linkinv
ilink.gam.av.formicula


mod.av.formicula<-gam.av.beta.formicula.12
ndata.av.formicula <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.av.formiculael for the new data
ndata.av.formicula <- add_column(ndata.av.formicula, fit = predict(mod.av.formicula, newdata = ndata.av.formicula, type = 'response'))


ndata.av.formicula <- bind_cols(ndata.av.formicula, setNames(as_tibble(predict(mod.av.formicula, ndata.av.formicula, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.formicula <- mutate(ndata.av.formicula,
                          fit_resp  = ilink.gam.av.formicula(fit_link),
                          right_upr = ilink.gam.av.formicula(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.av.formicula(fit_link - (2 * se_link)))


ndata.av.formicula$av.pH.unscaled<-ndata.av.formicula$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.formicula <- ggplot(ndata.av.formicula, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = formicula.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Folliculina")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.formicula,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.formicula
ggsave("C:Graphs May 2019//formicula_pred_av.png")



# GAM beta membranipora / gam.av.beta.alive.mem.12 --------------------------------------------------------

fitBinom=fitdist(data=food.exp.data.12.2019$alive.mem, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$alive.mem, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$alive.mem+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$alive.mem+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#good



glm.binomial.alive.mem.12.hydrogen<-glm(formula = cbind(alive.mem, 100-alive.mem)~ (av.pH)*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)

#binomial first
gam.av.binomial.alive.mem.12<- gam(formula = cbind(alive.mem, 100-alive.mem)~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$alive.mem.001<-0.01*(food.exp.data.12.2019_zscores$alive.mem+0.01)

#beta next
gam.av.beta.alive.mem.12<- gam(alive.mem.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.alive.mem.12.1<- gam(alive.mem.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.alive.mem.12.2<- gam(alive.mem.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.av.beta.alive.mem.12.3<- gam(alive.mem.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab( gam.av.beta.alive.mem.12, gam.av.beta.alive.mem.12.1, gam.av.beta.alive.mem.12.2, gam.av.binomial.alive.mem.12, glm.binomial.alive.mem.12.hydrogen)
#logit is of the best (but 0 between it and cloglog so simpler is better)


plot(gam.av.beta.alive.mem.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.alive.mem.12)
qq_plot(gam.av.beta.alive.mem.12, method = 'simulate')
k.check(gam.av.beta.alive.mem.12)
summary(gam.av.beta.alive.mem.12)

gam.av.beta.alive.mem.12.unordered<- gam(alive.mem.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.av.alive.mem <- family(gam.av.beta.alive.mem.12)
fam.gam.av.alive.mem
str(fam.gam.av.alive.mem)
ilink.gam.av.alive.mem<- fam.gam.av.alive.mem$linkinv
ilink.gam.av.alive.mem


mod.av.alive.mem<-gam.av.beta.alive.mem.12
ndata.av.alive.mem <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.av.alive.memel for the new data
ndata.av.alive.mem <- add_column(ndata.av.alive.mem, fit = predict(mod.av.alive.mem, newdata = ndata.av.alive.mem, type = 'response'))


ndata.av.alive.mem <- bind_cols(ndata.av.alive.mem, setNames(as_tibble(predict(mod.av.alive.mem, ndata.av.alive.mem, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.alive.mem <- mutate(ndata.av.alive.mem,
                          fit_resp  = ilink.gam.av.alive.mem(fit_link),
                          right_upr = ilink.gam.av.alive.mem(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.av.alive.mem(fit_link - (2 * se_link)))

ndata.av.alive.mem$av.pH.unscaled<-ndata.av.alive.mem$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.alive.mem <- ggplot(ndata.av.alive.mem, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = alive.mem.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Membranipora")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.alive.mem,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.alive.mem
ggsave("C:Graphs May 2019//alive.mem_pred_av.png")


# GAM beta didemnum / gam.av.beta.didemnum.12 ------------------------------------------------------------


str(food.exp.data.12.2019_zscores)
par(op)
fitBinom=fitdist(data=food.exp.data.12.2019$didemnum, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$didemnum, "binom", size = 100, prob = fitBinom$estimate[[1]])

beta.12<-fitdist(0.01*(food.exp.data.12.2019$didemnum+0.01), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$didemnum+0.01), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#not great at all


#binomial first
gam.av.binomial.didemnum.12<- gam(formula = cbind(didemnum, 100-didemnum)~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")

food.exp.data.12.2019_zscores$didemnum.001<-0.01*(food.exp.data.12.2019_zscores$didemnum+0.01)

#beta next
gam.av.beta.didemnum.12<- gam(didemnum.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.didemnum.12.1<- gam(didemnum.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.didemnum.12.2<- gam(didemnum.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.av.beta.didemnum.12.3<- gam(didemnum.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")


AICtab(gam.av.beta.didemnum.12, gam.av.beta.didemnum.12.1, gam.av.beta.didemnum.12.2,  gam.av.binomial.didemnum.12)
#logit is the best

plot(gam.av.beta.didemnum.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.didemnum.12)
qq_plot(gam.av.beta.didemnum.12, method = 'simulate')
k.check(gam.av.beta.didemnum.12)
summary(gam.av.beta.didemnum.12)

#appraise doesn't fit that well .... *** what to do? 

gam.av.beta.didemnum.12.unordered<- gam(didemnum.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")


fam.gam.av.didemnum <- family(gam.av.beta.didemnum.12)
fam.gam.av.didemnum
str(fam.gam.av.didemnum)
ilink.gam.av.didemnum<- fam.gam.av.didemnum$linkinv
ilink.gam.av.didemnum


mod.av.didemnum<-gam.av.beta.didemnum.12
ndata.av.didemnum <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the mod.av.didemnumel for the new data
ndata.av.didemnum <- add_column(ndata.av.didemnum, fit = predict(mod.av.didemnum, newdata = ndata.av.didemnum, type = 'response'))


ndata.av.didemnum <- bind_cols(ndata.av.didemnum, setNames(as_tibble(predict(mod.av.didemnum, ndata.av.didemnum, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.didemnum <- mutate(ndata.av.didemnum,
                          fit_resp  = ilink.gam.av.didemnum(fit_link),
                          right_upr = ilink.gam.av.didemnum(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.av.didemnum(fit_link - (2 * se_link)))

ndata.av.didemnum$av.pH.unscaled<-ndata.av.didemnum$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.didemnum <- ggplot(ndata.av.didemnum, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = didemnum.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Didemnum")~ "abundance", paste("(proportion cover)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.didemnum,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.didemnum
ggsave("C:Graphs May 2019//didemnum_pred_av.png")



# GAM negbin mussel incomplete / gam.av.nb.mussel_complete.12 ---------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$mussel_complete, "Poisson")
qqp(food.exp.data.12.2019_zscores$mussel_complete, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.mussel <- fitdistr(food.exp.data.12.2019_zscores$mussel_complete, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$mussel_complete, "nbinom", size = nbinom12.mussel$estimate[[1]], mu = nbinom12.mussel$estimate[[2]])
#better than poisson



glm.av.poisson.mussel_complete.12<- glm(mussel_complete~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")

glm.av.nb.mussel_complete.12.hydrogen<- glm.nb(mussel_complete~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

AICtab(glm.av.poisson.mussel_complete.12, glm.av.nb.mussel_complete.12.hydrogen)


overdisp_fun(glm.av.poisson.mussel_complete.12)
#very overdispersed


#negative binomial first
gam.av.nb.mussel_complete.12<- gam(mussel_complete ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.mussel_complete.12.1<- gam(mussel_complete ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")

gam.av.nb.mussel_complete.12.2<- gam(mussel_complete ~ s(av.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(link="log"), select=TRUE, method="REML")
gam.av.nb.mussel_complete.12.3<- gam(mussel_complete ~ s(av.pH)+ oFood.quality + s(min.10.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(link="sqrt"), select=TRUE, method="REML")



nb(link="sqrt")
#negbin with theta estimated is better

AICtab(gam.av.nb.mussel_complete.12.2, gam.av.nb.mussel_complete.12.3, gam.av.nb.mussel_complete.12, glm.nb.mussel_complete.12.hydrogen, gam.av.nb.mussel_complete.12.1)
#gam is better


plot(gam.av.nb.mussel_complete.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.nb.mussel_complete.12)
qq_plot(gam.av.nb.mussel_complete.12, method = 'simulate')
k.check(gam.av.nb.mussel_complete.12)
summary(gam.av.nb.mussel_complete.12)


gam.av.nb.mussel_complete.12.unordered<- gam(mussel_complete ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.mussel$estimate[[1]]), select=TRUE, method="REML")


fam.gam.av.mussel_complete <- family(gam.av.nb.mussel_complete.12)
fam.gam.av.mussel_complete
str(fam.gam.av.mussel_complete)
ilink.gam.av.mussel_complete<- fam.gam.av.mussel_complete$linkinv
ilink.gam.av.mussel_complete

mod.av.mussel_complete<-gam.av.nb.mussel_complete.12
ndata.av.mussel_complete <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.av.mussel_complete

str(food.exp.data.12.2019_zscores)
str(ndata.av.mussel_complete)

## add the fitted values by predicting from the model for the new data
ndata.av.mussel_complete <- add_column(ndata.av.mussel_complete, fit = predict(mod.av.mussel_complete, newdata = ndata.av.mussel_complete, type = 'response'))

predict(mod.av.mussel_complete, newdata = ndata.av.mussel_complete, type = 'response')
ndata.av.mussel_complete <- bind_cols(ndata.av.mussel_complete, setNames(as_tibble(predict(mod.av.mussel_complete, ndata.av.mussel_complete, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.mussel_complete <- mutate(ndata.av.mussel_complete,
                          fit_resp  = ilink.gam.av.mussel_complete(fit_link),
                          right_upr = ilink.gam.av.mussel_complete(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.av.mussel_complete(fit_link - (2 * se_link)))

ndata.av.mussel_complete$av.pH.unscaled<-ndata.av.mussel_complete$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.mussel_complete <- ggplot(ndata.av.mussel_complete, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = mussel_complete, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Mytilus")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.mussel_complete,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.mussel_complete
ggsave("C:Graphs May 2019//mussel_complete_pred_av.png")


# GAM negbin barnacles / gam.av.nb.num.barn.alive.12 -----------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.barn.alive, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.barn.alive, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.barnacles <- fitdistr(food.exp.data.12.2019_zscores$num.barn.alive, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.barn.alive, "nbinom", size = nbinom12.barnacles$estimate[[1]], mu = nbinom12.barnacles$estimate[[2]])
#better than poisson


glm.nb.num.barn.alive.12.hydrogen<- glm.nb(num.barn.alive~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.binomial.num.barn.alive.12.hydrogen)

glm.poisson.num.barn.alive.12.hydrogen<- glm(num.barn.alive~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.barn.alive.12.hydrogen)
#yes overdispersed

#negative binomial 
gam.av.nb.num.barn.alive.12<- gam(num.barn.alive ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.barnacles$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.num.barn.alive.12.1<- gam(num.barn.alive ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")

AICtab(gam.av.nb.num.barn.alive.12, glm.nb.num.barn.alive.12.hydrogen, gam.av.nb.num.barn.alive.12.1)

plot(gam.av.nb.num.barn.alive.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.nb.num.barn.alive.12)
qq_plot(gam.av.nb.num.barn.alive.12, method = 'simulate')
k.check(gam.av.nb.num.barn.alive.12)
summary(gam.av.nb.num.barn.alive.12)


#a few outside the area
#appraise a bit funnelly

gam.av.nb.num.barn.alive.12.unordered<- gam(num.barn.alive ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.barnacles$estimate[[1]]), select=TRUE, method="REML")


fam.gam.av.num.barn.alive <- family(gam.av.nb.num.barn.alive.12)
fam.gam.av.num.barn.alive
str(fam.gam.av.num.barn.alive)
ilink.gam.av.num.barn.alive<- fam.gam.av.num.barn.alive$linkinv
ilink.gam.av.num.barn.alive


mod.av.num.barn.alive<-gam.av.nb.num.barn.alive.12
ndata.av.num.barn.alive <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.av.num.barn.alive


## add the fitted values by predicting from the model for the new data
ndata.av.num.barn.alive <- add_column(ndata.av.num.barn.alive, fit = predict(mod.av.num.barn.alive, newdata = ndata.av.num.barn.alive, type = 'response'))

predict(mod.av.num.barn.alive, newdata = ndata.av.num.barn.alive, type = 'response')
ndata.av.num.barn.alive <- bind_cols(ndata.av.num.barn.alive, setNames(as_tibble(predict(mod.av.num.barn.alive, ndata.av.num.barn.alive, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.num.barn.alive <- mutate(ndata.av.num.barn.alive,
                                  fit_resp  = ilink.gam.av.num.barn.alive(fit_link),
                                  right_upr = ilink.gam.av.num.barn.alive(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.av.num.barn.alive(fit_link - (2 * se_link)))


ndata.av.num.barn.alive$av.pH.unscaled<-ndata.av.num.barn.alive$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.num.barn.alive <- ggplot(ndata.av.num.barn.alive, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.barn.alive, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Balanus")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.num.barn.alive,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.num.barn.alive
ggsave("C:Graphs May 2019//num.barn.alive_pred_av.png")



# GAM negbin disporella / gam.av.nb.disporella.12 ----------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$disporella, "Poisson")
qqp(food.exp.data.12.2019_zscores$disporella, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.disporella <- fitdistr(food.exp.data.12.2019_zscores$disporella, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$disporella, "nbinom", size = nbinom12.disporella$estimate[[1]], mu = nbinom12.disporella$estimate[[2]])
#better than poisson
#strill not great



glm.binomial.disporella.12.hydrogen<- glm.nb(disporella~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.binomial.disporella.12.hydrogen)
#not great

glm.poisson.disporella.12.hydrogen<- glm(disporella~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.disporella.12.hydrogen)
#poisson is overdispersed

#negative binomial
gam.av.nb.disporella.12<- gam(disporella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.disporella.12.1<- gam(disporella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")



AICtab(gam.av.nb.disporella.12.1,gam.av.nb.disporella.12, glm.binomial.disporella.12.hydrogen)

plot(gam.av.nb.disporella.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.nb.disporella.12)
qq_plot(gam.av.nb.disporella.12, method = 'simulate')
k.check(gam.av.nb.disporella.12)
summary(gam.av.nb.disporella.12)

#a few outside the area
#appraise a bit funnelly

gam.av.nb.disporella.12.unordered<- gam(disporella ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.disporella$estimate[[1]]), select=TRUE, method="REML")


fam.gam.av.disporella <- family(gam.av.nb.disporella.12)
fam.gam.av.disporella
str(fam.gam.av.disporella)
ilink.gam.av.disporella<- fam.gam.av.disporella$linkinv
ilink.gam.av.disporella


mod.av.disporella<-gam.av.nb.disporella.12
ndata.av.disporella <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                       length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.av.disporella

str(food.exp.data.12.2019_zscores)
str(ndata.av.disporella)

## add the fitted values by predicting from the model for the new data
ndata.av.disporella <- add_column(ndata.av.disporella, fit = predict(mod.av.disporella, newdata = ndata.av.disporella, type = 'response'))

predict(mod.av.disporella, newdata = ndata.av.disporella, type = 'response')
ndata.av.disporella <- bind_cols(ndata.av.disporella, setNames(as_tibble(predict(mod.av.disporella, ndata.av.disporella, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.disporella <- mutate(ndata.av.disporella,
                               fit_resp  = ilink.gam.av.disporella(fit_link),
                               right_upr = ilink.gam.av.disporella(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.av.disporella(fit_link - (2 * se_link)))

ndata.av.disporella$av.pH.unscaled<-ndata.av.disporella$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.disporella <- ggplot(ndata.av.disporella, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = disporella, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Disporella")~ "abundance", paste("(# of colonies)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.disporella,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.disporella
ggsave("C:Graphs May 2019//disporella_pred_av.png")


# GAM negbin schizo / gam.av.nb.schizo.12 --------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$schizo, "Poisson")
qqp(food.exp.data.12.2019_zscores$schizo, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.schizo <- fitdistr(food.exp.data.12.2019_zscores$schizo, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$schizo, "nbinom", size = nbinom12.schizo$estimate[[1]], mu = nbinom12.schizo$estimate[[2]])
#strill not great



glm.nb.schizo.12.hydrogen<- glm.nb(schizo~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.nb.schizo.12.hydrogen)


glm.poisson.schizo.12.hydrogen<- glm(schizo~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.schizo.12.hydrogen)
#poisson is not overdispersed

#negative binomial first
gam.av.nb.schizo.12<- gam(schizo ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.schizo$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.schizo.12.1<- gam(schizo ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.av.poisson.schizo.12<- gam(schizo ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.av.nb.schizo.12, gam.av.nb.schizo.12.1, glm.nb.schizo.12.hydrogen, glm.poisson.schizo.12.hydrogen, gam.av.poisson.schizo.12)


plot(gam.av.nb.schizo.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.nb.schizo.12)
qq_plot(gam.av.nb.schizo.12, method = 'simulate')
k.check(gam.av.nb.schizo.12)
summary(gam.av.nb.schizo.12)

##### kindex pvalue is signifiant! ???? 

gam.av.nb.schizo.12.unordered<- gam(schizo ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.schizo$estimate[[1]]), select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.av.schizo <- family(gam.av.nb.schizo.12)
fam.gam.av.schizo
str(fam.gam.av.schizo)
ilink.gam.av.schizo<- fam.gam.av.schizo$linkinv
ilink.gam.av.schizo

mod.av.schizo<-gam.av.nb.schizo.12
ndata.av.schizo <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                   length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.av.schizo

str(food.exp.data.12.2019_zscores)
str(ndata.av.schizo)

## add the fitted values by predicting from the model for the new data
ndata.av.schizo <- add_column(ndata.av.schizo, fit = predict(mod.av.schizo, newdata = ndata.av.schizo, type = 'response'))

predict(mod.av.schizo, newdata = ndata.av.schizo, type = 'response')
ndata.av.schizo <- bind_cols(ndata.av.schizo, setNames(as_tibble(predict(mod.av.schizo, ndata.av.schizo, se.fit = TRUE)[1:2]),
                                                         c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.schizo <- mutate(ndata.av.schizo,
                           fit_resp  = ilink.gam.av.schizo(fit_link),
                           right_upr = ilink.gam.av.schizo(fit_link + (2 * se_link)),
                           right_lwr = ilink.gam.av.schizo(fit_link - (2 * se_link)))


ndata.av.schizo$av.pH.unscaled<-ndata.av.schizo$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 
plt.av.schizo <- ggplot(ndata.av.schizo, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = schizo, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Schizoporella")~ "abundance", paste("(# of colonies)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.schizo,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.schizo
ggsave("C:Graphs May 2019//schizo_pred_av.png")



# GAM negbin num nudi / gam.av.nb.num.nudi.12 ------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.nudi, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.nudi, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.num.nudi <- fitdistr(food.exp.data.12.2019_zscores$num.nudi, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.nudi, "nbinom", size = nbinom12.num.nudi$estimate[[1]], mu = nbinom12.num.nudi$estimate[[2]])
#strill not great

glm.nb.num.nudi.12.hydrogen<- glm.nb(num.nudi~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019_zscores)
plot(glm.nb.num.nudi.12.hydrogen)



glm.poisson.num.nudi.12.hydrogen<- glm(num.nudi~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.nudi.12.hydrogen)
#poisson is not overdispersed


#negative binomial first
gam.av.nb.num.nudi.12.1<- gam(num.nudi ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.av.nb.num.nudi.12<- gam(num.nudi ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.nudi$estimate[[1]]), select=TRUE, method="REML")


gam.av.poisson.num.nudi.12<- gam(num.nudi ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.av.nb.num.nudi.12, gam.av.nb.num.nudi.12.1, glm.nb.num.nudi.12.hydrogen, glm.poisson.num.nudi.12.hydrogen, gam.av.poisson.num.nudi.12)

###poisson is the best ---> but not gam.... 
#####gam.av.nb is the best gam.... 
#so if I were to go with gam... 
#<1 AIC


appraise(gam.av.nb.num.nudi.12)
qq_plot(gam.av.nb.num.nudi.12, method = 'simulate')
plot(gam.av.nb.num.nudi.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.av.nb.num.nudi.12)
summary(gam.av.nb.num.nudi.12)

#K CHECK IS BAD!!!!! 
gam.av.nb.num.nudi.12.unordered<- gam(num.nudi ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.nudi$estimate[[1]]), select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.av.num.nudi <- family(gam.av.nb.num.nudi.12)
fam.gam.av.num.nudi
str(fam.gam.av.num.nudi)
ilink.gam.av.num.nudi<- fam.gam.av.num.nudi$linkinv
ilink.gam.av.num.nudi

mod.av.num.nudi<-gam.av.nb.num.nudi.12
ndata.av.num.nudi <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                               length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

ndata.av.num.nudi

str(food.exp.data.12.2019_zscores)
str(ndata.av.num.nudi)

## add the fitted values by predicting from the model for the new data
ndata.av.num.nudi <- add_column(ndata.av.num.nudi, fit = predict(mod.av.num.nudi, newdata = ndata.av.num.nudi, type = 'response'))

predict(mod.av.num.nudi, newdata = ndata.av.num.nudi, type = 'response')
ndata.av.num.nudi <- bind_cols(ndata.av.num.nudi, setNames(as_tibble(predict(mod.av.num.nudi, ndata.av.num.nudi, se.fit = TRUE)[1:2]),
                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.num.nudi <- mutate(ndata.av.num.nudi,
                       fit_resp  = ilink.gam.av.num.nudi(fit_link),
                       right_upr = ilink.gam.av.num.nudi(fit_link + (2 * se_link)),
                       right_lwr = ilink.gam.av.num.nudi(fit_link - (2 * se_link)))


ndata.av.num.nudi$av.pH.unscaled<-ndata.av.num.nudi$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 
plt.av.num.nudi <- ggplot(ndata.av.num.nudi, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.nudi, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Hermissenda")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.num.nudi,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.num.nudi
ggsave("C:Graphs May 2019//num.nudi_pred_av.png")


# GAM nb() serpulids / gam.av.nb.num.serpulid.12.1 -----------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.serpulid, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.serpulid, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.num.serpulid <- fitdistr(food.exp.data.12.2019_zscores$num.serpulid, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.serpulid, "nbinom", size = nbinom12.num.serpulid$estimate[[1]], mu = nbinom12.num.serpulid$estimate[[2]])
#better, not great

glm.nb.num.serpulid.12<- glm.nb(num.serpulid~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.nb.num.serpulid.12)
summary(glm.nb.num.serpulid.12)


glm.poisson.num.serpulid.12.hydrogen<- glm(num.serpulid~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.serpulid.12.hydrogen)
#poisson is definitely overdispersed

#negative binomial first
gam.av.nb.num.serpulid.12<- gam(num.serpulid ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.serpulid$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.num.serpulid.12.1<- gam(num.serpulid ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.av.poisson.num.serpulid.12<- gam(num.serpulid ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.av.nb.num.serpulid.12, gam.av.nb.num.serpulid.12.1, glm.nb.num.serpulid.12, glm.poisson.num.serpulid.12.hydrogen, gam.av.poisson.num.serpulid.12)

###gam.av.nb.num.serpulid.12.1 is best  
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

### h

appraise(gam.av.nb.num.serpulid.12.1)
qq_plot(gam.av.nb.num.serpulid.12.1, method = 'simulate')
plot(gam.av.nb.num.serpulid.12.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.av.nb.num.serpulid.12.1)
summary(gam.av.nb.num.serpulid.12.1)

#residuals a bit weird .... but they are the same in glm as in gam??? whic

gam.av.nb.num.serpulid.12.1.unordered<- gam(num.serpulid ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.av.num.serpulid <- family(gam.av.nb.num.serpulid.12.1)
fam.gam.av.num.serpulid
str(fam.gam.av.num.serpulid)
ilink.gam.av.num.serpulid<- fam.gam.av.num.serpulid$linkinv
ilink.gam.av.num.serpulid

mod.av.num.serpulid<-gam.av.nb.num.serpulid.12.1
ndata.av.num.serpulid <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                 length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.av.num.serpulid <- add_column(ndata.av.num.serpulid, fit = predict(mod.av.num.serpulid, newdata = ndata.av.num.serpulid, type = 'response'))

predict(mod.av.num.serpulid, newdata = ndata.av.num.serpulid, type = 'response')
ndata.av.num.serpulid <- bind_cols(ndata.av.num.serpulid, setNames(as_tibble(predict(mod.av.num.serpulid, ndata.av.num.serpulid, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.num.serpulid <- mutate(ndata.av.num.serpulid,
                         fit_resp  = ilink.gam.av.num.serpulid(fit_link),
                         right_upr = ilink.gam.av.num.serpulid(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.av.num.serpulid(fit_link - (2 * se_link)))


ndata.av.num.serpulid$av.pH.unscaled<-ndata.av.num.serpulid$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 
plt.av.num.serpulid <- ggplot(ndata.av.num.serpulid, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.serpulid, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop("Serpulid abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.num.serpulid,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.num.serpulid
ggsave("C:Graphs May 2019//num.serpulid_pred_av.png")


# GAM negbin orange sponge / gam.av.nb.orange_sponge.12.1 -------------------------------------------------------



poisson.12<-fitdistr(food.exp.data.12.2019_zscores$orange_sponge, "Poisson")
qqp(food.exp.data.12.2019_zscores$orange_sponge, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.orange_sponge <- fitdistr(food.exp.data.12.2019_zscores$orange_sponge, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$orange_sponge, "nbinom", size = nbinom12.orange_sponge$estimate[[1]], mu = nbinom12.orange_sponge$estimate[[2]])
#better, not great

glm.nb.orange_sponge.12<- glm.nb(orange_sponge~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

glm.poisson.orange_sponge.12.hydrogen<- glm(orange_sponge~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.orange_sponge.12.hydrogen)
#poisson is definitely overdispersed

#negative binomial first
gam.av.nb.orange_sponge.12<- gam(orange_sponge ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.orange_sponge$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.orange_sponge.12.1<- gam(orange_sponge ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.av.poisson.orange_sponge.12<- gam(orange_sponge ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.av.nb.orange_sponge.12, gam.av.nb.orange_sponge.12.1, glm.nb.orange_sponge.12, glm.poisson.orange_sponge.12.hydrogen, gam.av.poisson.orange_sponge.12)

###gam nb best
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.av.nb.orange_sponge.12.1)
qq_plot(gam.av.nb.orange_sponge.12.1, method = 'simulate')
plot(gam.av.nb.orange_sponge.12.1, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.av.nb.orange_sponge.12.1)
summary(gam.av.nb.orange_sponge.12.1)

#residuals a bit weird .... but they are the same in glm as in gam??? whic

gam.av.nb.orange_sponge.12.1.unordered<- gam(orange_sponge ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.av.orange_sponge <- family(gam.av.nb.orange_sponge.12.1)
fam.gam.av.orange_sponge
str(fam.gam.av.orange_sponge)
ilink.gam.av.orange_sponge<- fam.gam.av.orange_sponge$linkinv
ilink.gam.av.orange_sponge

mod.av.orange_sponge<-gam.av.nb.orange_sponge.12.1
ndata.av.orange_sponge <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                     length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.av.orange_sponge <- add_column(ndata.av.orange_sponge, fit = predict(mod.av.orange_sponge, newdata = ndata.av.orange_sponge, type = 'response'))

predict(mod.av.orange_sponge, newdata = ndata.av.orange_sponge, type = 'response')
ndata.av.orange_sponge <- bind_cols(ndata.av.orange_sponge, setNames(as_tibble(predict(mod.av.orange_sponge, ndata.av.orange_sponge, se.fit = TRUE)[1:2]),
                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.orange_sponge <- mutate(ndata.av.orange_sponge,
                             fit_resp  = ilink.gam.av.orange_sponge(fit_link),
                             right_upr = ilink.gam.av.orange_sponge(fit_link + (2 * se_link)),
                             right_lwr = ilink.gam.av.orange_sponge(fit_link - (2 * se_link)))


ndata.av.orange_sponge$av.pH.unscaled<-ndata.av.orange_sponge$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 
plt.av.orange_sponge <- ggplot(ndata.av.orange_sponge, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = orange_sponge, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop("Sponge abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.orange_sponge,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.orange_sponge
ggsave("C:Graphs May 2019//orange_sponge_pred_av.png")



# GAM negbin corella / gam.av.nb.num.corella.12 -------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$num.corella, "Poisson")
qqp(food.exp.data.12.2019_zscores$num.corella, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.num.corella <- fitdistr(food.exp.data.12.2019_zscores$num.corella, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$num.corella, "nbinom", size = nbinom12.num.corella$estimate[[1]], mu = nbinom12.num.corella$estimate[[2]])
#better, not great

glm.nb.num.corella.12<- glm.nb(num.corella~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

glm.poisson.num.corella.12.hydrogen<- glm(num.corella~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.num.corella.12.hydrogen)
#poisson is definitely overdispersed

#negative binomial first
gam.av.nb.num.corella.12<- gam(num.corella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.num.corella.12.1<- gam(num.corella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.av.nb.num.corella.12.2<- gam(num.corella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE)

gam.av.nb.num.corella.12.3<- gam(num.corella ~ s(av.pH, k=18)+ oFood.quality + s(av.pH, by=oFood.quality, k=20),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE)
gam.av.nb.num.corella.12.GCV<- gam(num.corella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="GCV.Cp")


gam.av.poisson.num.corella.12<- gam(num.corella ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.av.nb.num.corella.12.GCV, gam.av.nb.num.corella.12.3, gam.av.nb.num.corella.12.2, gam.av.nb.num.corella.12, gam.av.nb.num.corella.12.1, glm.nb.num.corella.12, glm.poisson.num.corella.12.hydrogen, gam.av.poisson.num.corella.12)

###gam negbinis best
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.av.nb.num.corella.12.2)
qq_plot(gam.av.nb.num.corella.12, method = 'simulate')
plot(gam.av.nb.num.corella.12.2, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.av.nb.num.corella.12.3)
summary(gam.av.nb.num.corella.12.3)

?gam

#residuals a bit patterny? 
###k significant

gam.av.nb.num.corella.12.unordered<- gam(num.corella ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.num.corella$estimate[[1]]), select=TRUE, method="REML")

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.av.num.corella <- family(gam.av.nb.num.corella.12)
fam.gam.av.num.corella
str(fam.gam.av.num.corella)
ilink.gam.av.num.corella<- fam.gam.av.num.corella$linkinv
ilink.gam.av.num.corella

mod.av.num.corella<-gam.av.nb.num.corella.12
ndata.av.num.corella <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                      length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.av.num.corella <- add_column(ndata.av.num.corella, fit = predict(mod.av.num.corella, newdata = ndata.av.num.corella, type = 'response'))

predict(mod.av.num.corella, newdata = ndata.av.num.corella, type = 'response')
ndata.av.num.corella <- bind_cols(ndata.av.num.corella, setNames(as_tibble(predict(mod.av.num.corella, ndata.av.num.corella, se.fit = TRUE)[1:2]),
                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.num.corella <- mutate(ndata.av.num.corella,
                              fit_resp  = ilink.gam.av.num.corella(fit_link),
                              right_upr = ilink.gam.av.num.corella(fit_link + (2 * se_link)),
                              right_lwr = ilink.gam.av.num.corella(fit_link - (2 * se_link)))


ndata.av.num.corella$av.pH.unscaled<-ndata.av.num.corella$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 
plt.av.num.corella <- ggplot(ndata.av.num.corella, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = num.corella, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(atop(italic("Corella")~ "abundance", paste("(# of individuals)"))))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.num.corella,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.num.corella
ggsave("C:Graphs May 2019//num.corella_pred_av.png")



# GAM poisson clam / gam.av.poisson.clam.12 ----------------------------------------------------------------

poisson.12<-fitdistr(food.exp.data.12.2019_zscores$clam, "Poisson")
qqp(food.exp.data.12.2019_zscores$clam, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.clam <- fitdistr(food.exp.data.12.2019_zscores$clam, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$clam, "nbinom", size = nbinom12.clam$estimate[[1]], mu = nbinom12.clam$estimate[[2]])
#better, not great

glm.nb.clam.12<- glm.nb(clam~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)

glm.poisson.clam.12.hydrogen<- glm(clam~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.poisson.clam.12.hydrogen)
#poisson is not overdispersed

#negative binomial first
gam.av.nb.clam.12<- gam(clam ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.clam$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.clam.12.1<- gam(clam ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")


gam.av.poisson.clam.12<- gam(clam ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")

AICtab(gam.av.nb.clam.12, gam.av.nb.clam.12.1, glm.nb.clam.12, glm.poisson.clam.12.hydrogen, gam.av.poisson.clam.12)

###gam poisson is best
#going with GAM across all is AIC <2
#easier for comparison with multiple models 
#easily interpretable outcome
#possibility for bayseian 

appraise(gam.av.poisson.clam.12)
qq_plot(gam.av.poisson.clam.12, method = 'simulate')
plot(gam.av.poisson.clam.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
k.check(gam.av.poisson.clam.12)
summary(gam.av.poisson.clam.12)

#residuals a bit patterny as well 

gam.av.poisson.clam.12.unordered<- gam(clam ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")


want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

fam.gam.av.clam <- family(gam.av.poisson.clam.12)
fam.gam.av.clam
str(fam.gam.av.clam)
ilink.gam.av.clam<- fam.gam.av.clam$linkinv
ilink.gam.av.clam

mod.av.clam<-gam.av.poisson.clam.12
ndata.av.clam <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                    length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))

## add the fitted values by predicting from the model for the new data
ndata.av.clam <- add_column(ndata.av.clam, fit = predict(mod.av.clam, newdata = ndata.av.clam, type = 'response'))

predict(mod.av.clam, newdata = ndata.av.clam, type = 'response')
ndata.av.clam <- bind_cols(ndata.av.clam, setNames(as_tibble(predict(mod.av.clam, ndata.av.clam, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.clam <- mutate(ndata.av.clam,
                            fit_resp  = ilink.gam.av.clam(fit_link),
                            right_upr = ilink.gam.av.clam(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.av.clam(fit_link - (2 * se_link)))


ndata.av.clam$av.pH.unscaled<-ndata.av.clam$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 
plt.av.clam <- ggplot(ndata.av.clam, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = clam, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Clam abundance\n(# of individuals)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.clam,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.clam
ggsave("C:Graphs May 2019//clam_pred_av.png")



# Fig 1 plot generation ---------------------------------------------------


library(cowplot)
fig.1.av<-plot_grid(plt.av.gam.av.hydroid,plt.av.alive.bot,plt.av.formicula,plt.caprellid.av.percent,plt.av.alive.mem,plt.av.didemnum,
          plt.av.mussel_complete,plt.av.num.barn.alive,plt.av.disporella,plt.av.schizo,plt.av.num.nudi,plt.av.num.serpulid,
          plt.av.orange_sponge,plt.av.num.corella,plt.av.clam,ncol=5, align='v', 
          labels=c('(a)', '(b)','(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)','(k)','(l)','(m)','(n)','(o)', 
          label_size=12))

fig.1.av

ggsave("C:For submission//Fig1.av.png", width=65, height=35, units="cm")



# Pulling model results to a table ----------------------------------------


#edf is low 0.0124
#You can often find spline terms with EDF <1 that are not linear. 
#In that case it seems that the null space (the linear, perfectly smooth bits)
#has been shrunk as well as the range space (the wiggly bits) but the smoothness parameter(s) 
#for the wiggly bit allow some small amount of wiggliness and the shrinkage to the null space drags 
#the EDF below 1. This is fine; the model is estimating a slightly wiggly function but uncertainty in 
#that estimate probably means that a linear function will reside in the 95% confidence interval for the smooth.
#https://stats.stackexchange.com/questions/179591/gam-with-low-e-d-f-estimated-degrees-of-freedom-value-in-main-effect-not-inte

hydroid.gam.av<-summary(gam.av.beta.hydroid.12.3)
botryllus.gam.av<-summary(gam.av.beta.alive.bot.12.3)
caprellid.gam.av<-summary(gam.beta.caprellid.av.percent.12)
formicula.gam.av<-summary(gam.av.beta.formicula.12)
alive.mem.gam.av<-summary(gam.av.beta.alive.mem.12)
didemnum.gam.av<-summary(gam.av.beta.didemnum.12)
mussel_complete.gam.av<-summary(gam.av.nb.mussel_complete.12)
alive.barn.gam.av<-summary(gam.av.nb.num.barn.alive.12)
disporella.gam.av<-summary(gam.av.nb.disporella.12)
schizo.gam.av<-summary(gam.av.nb.schizo.12)
num.nudi.gam.av<-summary(gam.av.nb.num.nudi.12)
num.serpulid.gam.av<-summary(gam.av.nb.num.serpulid.12.1)
clam.gam.av<-summary(gam.av.poisson.clam.12)
corella.gam.av<-summary(gam.av.nb.num.corella.12)
orange_sponge.gam.av<-summary(gam.av.nb.orange_sponge.12.1)

hydroid.gam.av.unordered<-summary(gam.av.beta.hydroid.12.3.unordered)
botryllus.gam.av.unordered<-summary(gam.av.beta.alive.bot.12.3.unordered)
caprellid.gam.av.unordered<-summary(gam.beta.caprellid.av.percent.12.unordered)
formicula.gam.av.unordered<-summary(gam.av.beta.formicula.12.unordered)
alive.mem.gam.av.unordered<-summary(gam.av.beta.alive.mem.12.unordered)
didemnum.gam.av.unordered<-summary(gam.av.beta.didemnum.12.unordered)
mussel_complete.gam.av.unordered<-summary(gam.av.nb.mussel_complete.12.unordered)
alive.barn.gam.av.unordered<-summary(gam.av.nb.num.barn.alive.12.unordered)
disporella.gam.av.unordered<-summary(gam.av.nb.disporella.12.unordered)
schizo.gam.av.unordered<-summary(gam.av.nb.schizo.12.unordered)
num.nudi.gam.av.unordered<-summary(gam.av.nb.num.nudi.12.unordered)
num.serpulid.gam.av.unordered<-summary(gam.av.nb.num.serpulid.12.1.unordered)
clam.gam.av.unordered<-summary(gam.av.poisson.clam.12.unordered)
corella.gam.av.unordered<-summary(gam.av.nb.num.corella.12.unordered)
orange_sponge.gam.av.unordered<-summary(gam.av.nb.orange_sponge.12.1.unordered)


hydroid.gam.av.p.table<-as.data.frame(hydroid.gam.av.unordered$p.table)
hydroid.gam.av.s.table<-as.data.frame(hydroid.gam.av$s.table)

botryllus.gam.av.p.table<-as.data.frame(botryllus.gam.av.unordered$p.table)
botryllus.gam.av.s.table<-as.data.frame(botryllus.gam.av$s.table)

caprellid.gam.av.p.table<-as.data.frame(caprellid.gam.av.unordered$p.table)
caprellid.gam.av.s.table<-as.data.frame(caprellid.gam.av$s.table)

formicula.gam.av.p.table<-as.data.frame(formicula.gam.av.unordered$p.table)
formicula.gam.av.s.table<-as.data.frame(formicula.gam.av$s.table)

alive.mem.gam.av.p.table<-as.data.frame(alive.mem.gam.av.unordered$p.table)
alive.mem.gam.av.s.table<-as.data.frame(alive.mem.gam.av$s.table)

didemnum.gam.av.p.table<-as.data.frame(didemnum.gam.av.unordered$p.table)
didemnum.gam.av.s.table<-as.data.frame(didemnum.gam.av$s.table)

mussel_complete.gam.av.p.table<-as.data.frame(mussel_complete.gam.av.unordered$p.table)
mussel_complete.gam.av.s.table<-as.data.frame(mussel_complete.gam.av$s.table)

alive.barn.gam.av.p.table<-as.data.frame(alive.barn.gam.av.unordered$p.table)
alive.barn.gam.av.s.table<-as.data.frame(alive.barn.gam.av$s.table)

disporella.gam.av.p.table<-as.data.frame(disporella.gam.av.unordered$p.table)
disporella.gam.av.s.table<-as.data.frame(disporella.gam.av$s.table)

schizo.gam.av.p.table<-as.data.frame(schizo.gam.av.unordered$p.table)
schizo.gam.av.s.table<-as.data.frame(schizo.gam.av$s.table)

num.nudi.gam.av.p.table<-as.data.frame(num.nudi.gam.av.unordered$p.table)
num.nudi.gam.av.s.table<-as.data.frame(num.nudi.gam.av$s.table)

num.serpulid.gam.av.p.table<-as.data.frame(num.serpulid.gam.av.unordered$p.table)
num.serpulid.gam.av.s.table<-as.data.frame(num.serpulid.gam.av$s.table)

orange_sponge.gam.av.p.table<-as.data.frame(orange_sponge.gam.av.unordered$p.table)
orange_sponge.gam.av.s.table<-as.data.frame(orange_sponge.gam.av$s.table)

corella.gam.av.p.table<-as.data.frame(corella.gam.av.unordered$p.table)
corella.gam.av.s.table<-as.data.frame(corella.gam.av$s.table)

clam.gam.av.p.table<-as.data.frame(clam.gam.av.unordered$p.table)
clam.gam.av.s.table<-as.data.frame(clam.gam.av$s.table)


#### Building the stats table
ptable.av<-rbind(hydroid.gam.av.p.table, 
               botryllus.gam.av.p.table, 
               caprellid.gam.av.p.table,
               formicula.gam.av.p.table,
               alive.mem.gam.av.p.table,
               didemnum.gam.av.p.table,
               mussel_complete.gam.av.p.table,
               alive.barn.gam.av.p.table,
               disporella.gam.av.p.table,
               schizo.gam.av.p.table,
               num.nudi.gam.av.p.table,
               num.serpulid.gam.av.p.table,
               orange_sponge.gam.av.p.table,
               corella.gam.av.p.table,
               clam.gam.av.p.table)


colnames(ptable.av) <- c("Estimate", "SE", "z", "p")
ptable.av$Factor<-rep(c("Intercept", "Low quality food", "High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.av %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, z, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia", 1, 3) %>%
  group_rows("Botryllus", 4, 6) %>% 
  group_rows("Formicula",7, 9) %>% 
  group_rows("Caprella", 10, 12) %>% 
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
  group_rows("Clams", 43, 45)%>% 
  save_kable(file = "C:For submission//ptable.av.html", self_contained = T)


### s table
stable.av<-rbind(hydroid.gam.av.s.table, 
              botryllus.gam.av.s.table, 
              caprellid.gam.av.s.table,
              formicula.gam.av.s.table,
              alive.mem.gam.av.s.table,
              didemnum.gam.av.s.table,
              mussel_complete.gam.av.s.table,
              alive.barn.gam.av.s.table,
              disporella.gam.av.s.table,
              schizo.gam.av.s.table,
              num.nudi.gam.av.s.table,
              num.serpulid.gam.av.s.table,
              orange_sponge.gam.av.s.table,
              corella.gam.av.s.table,
              clam.gam.av.s.table)


colnames(stable.av) <- c("Estimated_df", "Reference_df", "Chi_squared", "p_smooth")
stable.av$Smooth_terms<-rep(c("smooth av.pH", "smooth av.pH * Low quality food", "smooth av.pH * High quality food"))


stable.av %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, Chi_squared, p_smooth) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Obelia", 1, 3) %>%
  group_rows("Botryllus", 4, 6) %>% 
  group_rows("Formicula",7, 9) %>% 
  group_rows("Caprella", 10, 12) %>% 
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
  save_kable(file = "C:For submission//stable.av.html", self_contained = T)


pstable.av<-cbind(ptable.av, stable.av)

pstable.av %>% 
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
  group_rows("Caprella, negative binomial",7, 9) %>% 
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
  save_kable(file = "C:For submission//pstable.av.html", self_contained = T)











# Richness ----------------------------------------------------------------


poisson.12<-fitdistr(food.exp.data.12.2019_zscores$richness, "Poisson")
qqp(food.exp.data.12.2019_zscores$richness, "pois", lambda=poisson.12$estimate[[1]])


nbinom12.richness<- fitdistr(food.exp.data.12.2019_zscores$richness, "Negative Binomial")
qqp(food.exp.data.12.2019_zscores$richness, "nbinom", size = nbinom12.richness$estimate[[1]], mu = nbinom12.richness$estimate[[2]])
#worse than poiisson


glm.av.nb.richness.12.hydrogen<- glm.nb(richness~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores)
plot(glm.av.binomial.richness.12.hydrogen)

glm.av.poisson.richness.12.hydrogen<- glm(richness~ av.pH*oFood.quality, data = food.exp.data.12.2019_zscores, family="poisson")
overdisp_fun(glm.av.poisson.richness.12.hydrogen)
#not overdispersed

#negative binomial 
gam.av.nb.richness.12<- gam(richness ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = negbin(nbinom12.richness$estimate[[1]]), select=TRUE, method="REML")
gam.av.nb.richness.12.1<- gam(richness ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = nb(), select=TRUE, method="REML")
gam.av.poisson.richness.12<- gam(richness ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")



AICtab( gam.av.nb.richness.12.1, gam.av.poisson.richness.12)

plot(gam.av.poisson.richness.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.poisson.richness.12)
qq_plot(gam.av.poisson.richness.12, method = 'simulate')
k.check(gam.av.poisson.richness.12)
summary(gam.av.poisson.richness.12)


#a few outside the QQ plot on both ends... doesn't look great 


gam.av.poisson.richness.12.unordered<- gam(richness ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = poisson, select=TRUE, method="REML")


fam.gam.av.richness <- family(gam.av.poisson.richness.12)
fam.gam.av.richness
str(fam.gam.av.richness)
ilink.gam.av.richness<- fam.gam.av.richness$linkinv
ilink.gam.av.richness


mod.av.richness<-gam.av.poisson.richness.12
ndata.av.richness <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                 length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.richness <- add_column(ndata.av.richness, fit = predict(mod.av.richness, newdata = ndata.av.richness, type = 'response'))

predict(mod.av.richness, newdata = ndata.av.richness, type = 'response')
ndata.av.richness <- bind_cols(ndata.av.richness, setNames(as_tibble(predict(mod.av.richness, ndata.av.richness, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.richness <- mutate(ndata.av.richness,
                         fit_resp  = ilink.gam.av.richness(fit_link),
                         right_upr = ilink.gam.av.richness(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.av.richness(fit_link - (2 * se_link)))


ndata.av.richness$av.pH.unscaled<-ndata.av.richness$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.richness <- ggplot(ndata.av.richness, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = richness, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Species richness")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.richness,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.richness
ggsave("C:Graphs May 2019//richness_pred_av.png")







# Evenness ----------------------------------------------------------------
gamma.12.evenness<-fitdistr(food.exp.data.12.2019_zscores$evenness+0.01, "gamma")
qqp(food.exp.data.12.2019$evenness, "gamma", shape = gamma.12.evenness$estimate[[1]], rate = gamma.12.evenness$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$evenness)
#normal works

lm.av.evenness<-lm(evenness ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
glm.av.evenness<-glm(evenness ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.av.lm.evenness.12<- gam(evenness ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.evenness.12.1<- gam(evenness ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")



AICtab(lm.evenness, glm.av.evenness, gam.av.lm.evenness.12, gam.av.gamma.evenness.12.1)

plot(gam.av.lm.evenness.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.evenness.12)
qq_plot(gam.av.lm.evenness.12, method = 'simulate')
k.check(gam.av.lm.evenness.12)
summary(gam.av.lm.evenness.12)





gam.av.lm.evenness.12.unordered<- gam(evenness ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.av.evenness <- family(gam.av.lm.evenness.12)
fam.gam.av.evenness
str(fam.gam.av.evenness)
ilink.gam.av.evenness<- fam.gam.av.evenness$linkinv
ilink.gam.av.evenness


mod.av.evenness<-gam.av.lm.evenness.12
ndata.av.evenness <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                 length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.evenness <- add_column(ndata.av.evenness, fit = predict(mod.av.evenness, newdata = ndata.av.evenness, type = 'response'))

predict(mod.av.evenness, newdata = ndata.av.evenness, type = 'response')
ndata.av.evenness <- bind_cols(ndata.av.evenness, setNames(as_tibble(predict(mod.av.evenness, ndata.av.evenness, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.evenness <- mutate(ndata.av.evenness,
                         fit_resp  = ilink.gam.av.evenness(fit_link),
                         right_upr = ilink.gam.av.evenness(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.av.evenness(fit_link - (2 * se_link)))


ndata.av.evenness$av.pH.unscaled<-ndata.av.evenness$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.evenness <- ggplot(ndata.av.evenness, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = evenness, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Species evenness")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.evenness,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.evenness
ggsave("C:Graphs May 2019//evenness_pred_av.png")






# Occupied space ----------------------------------------------------------


beta.12<-fitdist(0.01*(food.exp.data.12.2019_zscores$occupied.space), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019_zscores$occupied.space), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])
#good! 

gamma.12.occupied.space<-fitdistr((food.exp.data.12.2019_zscores$occupied.space*0.01), "gamma")
qqp(food.exp.data.12.2019_zscores$occupied.space, "gamma", shape = gamma.12.occupied.space$estimate[[1]], rate = gamma.12.occupied.space$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$occupied.space)
#normal is okay not great 

lm.occupied.space<-lm(occupied.space ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
glm.av.occupied.space<-glm(occupied.space*0.01 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.av.lm.occupied.space.12<- gam(occupied.space ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.occupied.space.12.1<- gam(occupied.space*0.01 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

gam.av.binomial.occupied.space.12<- gam(formula = cbind(occupied.space, 100-occupied.space)~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = binomial, select=TRUE, method="REML")


food.exp.data.12.2019_zscores$occupied.space.001<-0.01*(food.exp.data.12.2019_zscores$occupied.space)

gam.av.beta.occupied.space.12<- gam(occupied.space.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.occupied.space.12.1<- gam(occupied.space.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.occupied.space.12.2<- gam(occupied.space.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
gam.av.beta.occupied.space.12.3<- gam(occupied.space.001~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cauchit"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(lm.occupied.space, glm.av.occupied.space, gam.av.lm.occupied.space.12, gam.av.gamma.occupied.space.12.1, gam.av.beta.occupied.space.12, gam.av.beta.occupied.space.12.1, gam.av.beta.occupied.space.12.2, gam.av.binomial.occupied.space.12)


#beta is best 

plot(gam.av.beta.occupied.space.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.occupied.space.12)
qq_plot(gam.av.beta.occupied.space.12, method = 'simulate')
k.check(gam.av.beta.occupied.space.12)
summary(gam.av.beta.occupied.space.12)





gam.av.beta.occupied.space.12.unordered<- gam(occupied.space.001~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
summary(gam.av.beta.occupied.space.12.unordered)

fam.gam.av.occupied.space <- family(gam.av.beta.occupied.space.12)
fam.gam.av.occupied.space
str(fam.gam.av.occupied.space)
ilink.gam.av.occupied.space<- fam.gam.av.occupied.space$linkinv
ilink.gam.av.occupied.space


mod.av.occupied.space<-gam.av.beta.occupied.space.12
ndata.av.occupied.space <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                       length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.occupied.space <- add_column(ndata.av.occupied.space, fit = predict(mod.av.occupied.space, newdata = ndata.av.occupied.space, type = 'response'))

predict(mod.av.occupied.space, newdata = ndata.av.occupied.space, type = 'response')
ndata.av.occupied.space <- bind_cols(ndata.av.occupied.space, setNames(as_tibble(predict(mod.av.occupied.space, ndata.av.occupied.space, se.fit = TRUE)[1:2]),
                                                                 c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.occupied.space <- mutate(ndata.av.occupied.space,
                               fit_resp  = ilink.gam.av.occupied.space(fit_link),
                               right_upr = ilink.gam.av.occupied.space(fit_link + (2 * se_link)),
                               right_lwr = ilink.gam.av.occupied.space(fit_link - (2 * se_link)))


ndata.av.occupied.space$av.pH.unscaled<-ndata.av.occupied.space$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.occupied.space <- ggplot(ndata.av.occupied.space, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = occupied.space.001, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Proportion of space on tile occupied")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.occupied.space,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.occupied.space
ggsave("C:Graphs May 2019//occupied.space_pred_av.png")



# Everything wet weight ---------------------------------------------------

gamma.12.everything.wet.weight<-fitdistr(food.exp.data.12.2019_zscores$everything.wet.weight+0.01, "gamma")
qqp(food.exp.data.12.2019$everything.wet.weight, "gamma", shape = gamma.12.everything.wet.weight$estimate[[1]], rate = gamma.12.everything.wet.weight$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$everything.wet.weight)

qqp(log(food.exp.data.12.2019_zscores$everything.wet.weight))
#log normal the best

qqp(food.exp.data.12.2019_zscores$everything.wet.weight, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.av.everything.wet.weight<-lm(everything.wet.weight ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.av.log.everything.wet.weight<-lm(log(everything.wet.weight) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.everything.wet.weight<-glm(everything.wet.weight ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.av.loglink.everything.wet.weight<-glm(everything.wet.weight ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.av.lm.everything.wet.weight.12<- gam(everything.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.everything.wet.weight.12<- gam(log(everything.wet.weight) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.everything.wet.weight.12.1<- gam(everything.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.everything.wet.weight.12.1, gam.av.lm.log.everything.wet.weight.12, gam.av.tweedie.everything.wet.weight.12.1, lm.log.everything.wet.weight, glm.av.loglink.everything.wet.weight, lm.everything.wet.weight, glm.av.everything.wet.weight, gam.av.lm.everything.wet.weight.12, gam.av.gamma.everything.wet.weight.12.1)

#gam.av.lm.log.everything.wet.weight.12 is best 


plot(gam.av.lm.log.everything.wet.weight.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.log.everything.wet.weight.12)
qq_plot(gam.av.lm.log.everything.wet.weight.12, method = 'simulate')
k.check(gam.av.lm.log.everything.wet.weight.12)
summary(gam.av.lm.log.everything.wet.weight.12)





gam.av.lm.log.everything.wet.weight.12.unordered<- gam(log(everything.wet.weight) ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.av.everything.wet.weight <- family(gam.av.lm.log.everything.wet.weight.12)
fam.gam.av.everything.wet.weight
str(fam.gam.av.everything.wet.weight)
ilink.gam.av.everything.wet.weight<- fam.gam.av.everything.wet.weight$linkinv
ilink.gam.av.everything.wet.weight


mod.av.everything.wet.weight<-gam.av.lm.log.everything.wet.weight.12
ndata.av.everything.wet.weight <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.everything.wet.weight <- add_column(ndata.av.everything.wet.weight, fit = predict(mod.av.everything.wet.weight, newdata = ndata.av.everything.wet.weight, type = 'response'))

predict(mod.av.everything.wet.weight, newdata = ndata.av.everything.wet.weight, type = 'response')
ndata.av.everything.wet.weight <- bind_cols(ndata.av.everything.wet.weight, setNames(as_tibble(predict(mod.av.everything.wet.weight, ndata.av.everything.wet.weight, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.everything.wet.weight <- mutate(ndata.av.everything.wet.weight,
                                      fit_resp  = ilink.gam.av.everything.wet.weight(fit_link),
                                      right_upr = ilink.gam.av.everything.wet.weight(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.av.everything.wet.weight(fit_link - (2 * se_link)))


ndata.av.everything.wet.weight$av.pH.unscaled<-ndata.av.everything.wet.weight$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.everything.wet.weight <- ggplot(ndata.av.everything.wet.weight, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(everything.wet.weight), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Total wet biomass per mesocosm \n(g, log scale)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.everything.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.everything.wet.weight
ggsave("C:Graphs May 2019//everything.wet.weight_pred_av.png")




# Everything wet weight per 1 ---------------------------------------------


gamma.12.everything.wet.weight.per.1<-fitdistr(food.exp.data.12.2019_zscores$everything.wet.weight.per.1+0.01, "gamma")
qqp(food.exp.data.12.2019$everything.wet.weight.per.1, "gamma", shape = gamma.12.everything.wet.weight.per.1$estimate[[1]], rate = gamma.12.everything.wet.weight.per.1$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$everything.wet.weight.per.1)

qqp(log(food.exp.data.12.2019_zscores$everything.wet.weight.per.1))
#log normal the best

qqp(food.exp.data.12.2019_zscores$everything.wet.weight.per.1, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.av.everything.wet.weight.per.1<-lm(everything.wet.weight.per.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.av.log.everything.wet.weight.per.1<-lm(log(everything.wet.weight.per.1) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.everything.wet.weight.per.1<-glm(everything.wet.weight.per.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.av.loglink.everything.wet.weight.per.1<-glm(everything.wet.weight.per.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.lm.av.everything.wet.weight.per.1.12<- gam(everything.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.everything.wet.weight.per.1.12<- gam(log(everything.wet.weight.per.1) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.everything.wet.weight.per.1.12.1<- gam(everything.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.everything.wet.weight.per.1.12.1, gam.av.lm.log.everything.wet.weight.per.1.12, gam.av.tweedie.everything.wet.weight.per.1.12.1, lm.av.log.everything.wet.weight.per.1, glm.av.loglink.everything.wet.weight.per.1, lm.av.everything.wet.weight.per.1, glm.av.everything.wet.weight.per.1, gam.lm.av.everything.wet.weight.per.1.12, gam.av.gamma.everything.wet.weight.per.1.12.1)

#gam.lm.log.everything.wet.weight.per.1.12 is best 


plot(gam.av.lm.log.everything.wet.weight.per.1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.log.everything.wet.weight.per.1.12)
qq_plot(gam.av.lm.log.everything.wet.weight.per.1.12, method = 'simulate')
k.check(gam.av.lm.log.everything.wet.weight.per.1.12)
summary(gam.av.lm.log.everything.wet.weight.per.1.12)





gam.av.lm.log.everything.wet.weight.per.1.12.unordered<- gam(log(everything.wet.weight.per.1) ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.av.gam.everything.wet.weight.per.1 <- family(gam.av.lm.log.everything.wet.weight.per.1.12)
fam.av.gam.everything.wet.weight.per.1
str(fam.av.gam.everything.wet.weight.per.1)
ilink.av.gam.everything.wet.weight.per.1<- fam.av.gam.everything.wet.weight.per.1$linkinv
ilink.av.gam.everything.wet.weight.per.1


mod.av.everything.wet.weight.per.1<-gam.av.lm.log.everything.wet.weight.per.1.12
ndata.av.everything.wet.weight.per.1 <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                                    length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.everything.wet.weight.per.1 <- add_column(ndata.av.everything.wet.weight.per.1, fit = predict(mod.av.everything.wet.weight.per.1, newdata = ndata.av.everything.wet.weight.per.1, type = 'response'))

predict(mod.av.everything.wet.weight.per.1, newdata = ndata.av.everything.wet.weight.per.1, type = 'response')
ndata.av.everything.wet.weight.per.1 <- bind_cols(ndata.av.everything.wet.weight.per.1, setNames(as_tibble(predict(mod.av.everything.wet.weight.per.1, ndata.av.everything.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.everything.wet.weight.per.1 <- mutate(ndata.av.everything.wet.weight.per.1,
                                            fit_resp  = ilink.av.gam.everything.wet.weight.per.1(fit_link),
                                            right_upr = ilink.av.gam.everything.wet.weight.per.1(fit_link + (2 * se_link)),
                                            right_lwr = ilink.av.gam.everything.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.av.everything.wet.weight.per.1$av.pH.unscaled<-ndata.av.everything.wet.weight.per.1$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.everything.wet.weight.per.1 <- ggplot(ndata.av.everything.wet.weight.per.1, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(everything.wet.weight.per.1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Log Biomass per 1 % cover \n(wet weight)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.everything.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.everything.wet.weight.per.1
ggsave("C:Graphs May 2019//everything.wet.weight.per.1_pred_av.png")



# Dry weight --------------------------------------------------------------


gamma.12.total_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$total_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$total_dry_biomass, "gamma", shape = gamma.12.total_dry_biomass$estimate[[1]], rate = gamma.12.total_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$total_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$total_dry_biomass))


qqp(food.exp.data.12.2019_zscores$total_dry_biomass, "lnorm")




lm.total_dry_biomass<-lm(total_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.total_dry_biomass<-lm(log(total_dry_biomass) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.total_dry_biomass<-glm(total_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
glm.av.loglink.total_dry_biomass<-glm(total_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.av.lm.total_dry_biomass.12<- gam(total_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.total_dry_biomass.12<- gam(log(total_dry_biomass) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.total_dry_biomass.12.1<- gam(total_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.total_dry_biomass.12.1, gam.av.lm.log.total_dry_biomass.12, gam.av.tweedie.total_dry_biomass.12.1, lm.log.total_dry_biomass, glm.av.loglink.total_dry_biomass, lm.total_dry_biomass, glm.av.total_dry_biomass, gam.av.lm.total_dry_biomass.12, gam.av.gamma.total_dry_biomass.12.1)

#gam.av.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.av.lm.total_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.total_dry_biomass.12 )
qq_plot(gam.av.lm.total_dry_biomass.12 , method = 'simulate')
k.check(gam.av.lm.total_dry_biomass.12 )
summary(gam.av.lm.total_dry_biomass.12 )





gam.av.lm.total_dry_biomass.12.unordered<- gam(total_dry_biomass ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
summary(gam.av.lm.total_dry_biomass.12.unordered)

fam.gam.av.total_dry_biomass <- family(gam.av.lm.total_dry_biomass.12)
fam.gam.av.total_dry_biomass
str(fam.gam.av.total_dry_biomass)
ilink.gam.av.total_dry_biomass<- fam.gam.av.total_dry_biomass$linkinv
ilink.gam.av.total_dry_biomass


mod.av.total_dry_biomass<-gam.av.lm.total_dry_biomass.12
ndata.av.total_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.total_dry_biomass <- add_column(ndata.av.total_dry_biomass, fit = predict(mod.av.total_dry_biomass, newdata = ndata.av.total_dry_biomass, type = 'response'))

predict(mod.av.total_dry_biomass, newdata = ndata.av.total_dry_biomass, type = 'response')
ndata.av.total_dry_biomass <- bind_cols(ndata.av.total_dry_biomass, setNames(as_tibble(predict(mod.av.total_dry_biomass, ndata.av.total_dry_biomass, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.total_dry_biomass <- mutate(ndata.av.total_dry_biomass,
                                  fit_resp  = ilink.gam.av.total_dry_biomass(fit_link),
                                  right_upr = ilink.gam.av.total_dry_biomass(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.av.total_dry_biomass(fit_link - (2 * se_link)))


ndata.av.total_dry_biomass$av.pH.unscaled<-ndata.av.total_dry_biomass$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.total_dry_biomass <- ggplot(ndata.av.total_dry_biomass, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (total_dry_biomass), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Total dry biomass per tile (g)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.total_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.total_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//total_dry_biomass_pred_av.png")

# hydroid biomass ---------------------------------------------------------



gamma.12.hydroid_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$hydroid_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$hydroid_dry_biomass, "gamma", shape = gamma.12.hydroid_dry_biomass$estimate[[1]], rate = gamma.12.hydroid_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$hydroid_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$hydroid_dry_biomass))


qqp(food.exp.data.12.2019_zscores$hydroid_dry_biomass, "lnorm")




lm.hydroid_dry_biomass<-lm(hydroid_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.hydroid_dry_biomass<-lm(log(hydroid_dry_biomass+0.1) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.hydroid_dry_biomass<-glm(hydroid_dry_biomass+0.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.av.lm.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.hydroid_dry_biomass.12<- gam(log(hydroid_dry_biomass+0.1) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.hydroid_dry_biomass.12<- gam(hydroid_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.hydroid_dry_biomass.12, gam.av.lm.log.hydroid_dry_biomass.12, gam.av.tweedie.hydroid_dry_biomass.12, lm.log.hydroid_dry_biomass, lm.hydroid_dry_biomass, glm.av.hydroid_dry_biomass, gam.av.lm.hydroid_dry_biomass.12, gam.av.gamma.hydroid_dry_biomass.12)



plot(gam.av.gamma.hydroid_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.gamma.hydroid_dry_biomass.12 )
qq_plot(gam.av.gamma.hydroid_dry_biomass.12 , method = 'simulate')
k.check(gam.av.gamma.hydroid_dry_biomass.12 )
summary(gam.av.gamma.hydroid_dry_biomass.12 )





gam.av.gamma.hydroid_dry_biomass.12.unordered<- gam(hydroid_dry_biomass+0.1 ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.av.hydroid_dry_biomass <- family(gam.av.gamma.hydroid_dry_biomass.12)
fam.gam.av.hydroid_dry_biomass
str(fam.gam.av.hydroid_dry_biomass)
ilink.gam.av.hydroid_dry_biomass<- fam.gam.av.hydroid_dry_biomass$linkinv
ilink.gam.av.hydroid_dry_biomass


mod.av.hydroid_dry_biomass<-gam.av.gamma.hydroid_dry_biomass.12
ndata.av.hydroid_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                            length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.hydroid_dry_biomass <- add_column(ndata.av.hydroid_dry_biomass, fit = predict(mod.av.hydroid_dry_biomass, newdata = ndata.av.hydroid_dry_biomass, type = 'response'))

predict(mod.av.hydroid_dry_biomass, newdata = ndata.av.hydroid_dry_biomass, type = 'response')
ndata.av.hydroid_dry_biomass <- bind_cols(ndata.av.hydroid_dry_biomass, setNames(as_tibble(predict(mod.av.hydroid_dry_biomass, ndata.av.hydroid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.hydroid_dry_biomass <- mutate(ndata.av.hydroid_dry_biomass,
                                    fit_resp  = ilink.gam.av.hydroid_dry_biomass(fit_link),
                                    right_upr = ilink.gam.av.hydroid_dry_biomass(fit_link + (2 * se_link)),
                                    right_lwr = ilink.gam.av.hydroid_dry_biomass(fit_link - (2 * se_link)))


ndata.av.hydroid_dry_biomass$av.pH.unscaled<-ndata.av.hydroid_dry_biomass$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.hydroid_dry_biomass <- ggplot(ndata.av.hydroid_dry_biomass, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (hydroid_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(italic("Obelia") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.hydroid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ylim(0,.75)
plt.av.hydroid_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//hydroid_dry_biomass_pred_av.png")




# botryllus biomass -------------------------------------------------------
food.exp.data.12.2019_zscores$tunicate_dry_biomass[food.exp.data.12.2019_zscores$tunicate_dry_biomass<0]<-0


gamma.12.tunicate_dry_biomass<-fitdistr((food.exp.data.12.2019_zscores$tunicate_dry_biomass+0.01), "gamma")
qqp(food.exp.data.12.2019_zscores$tunicate_dry_biomass, "gamma", shape = gamma.12.tunicate_dry_biomass$estimate[[1]], rate = gamma.12.tunicate_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$tunicate_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$tunicate_dry_biomass))


qqp(food.exp.data.12.2019_zscores$tunicate_dry_biomass, "lnorm")




lm.tunicate_dry_biomass<-lm(tunicate_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.tunicate_dry_biomass<-lm(log(tunicate_dry_biomass+0.1) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.tunicate_dry_biomass<-glm(tunicate_dry_biomass+0.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.av.lm.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.tunicate_dry_biomass.12<- gam(log(tunicate_dry_biomass+0.1) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.tunicate_dry_biomass.12<- gam(tunicate_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.tunicate_dry_biomass.12, gam.av.lm.log.tunicate_dry_biomass.12, gam.av.tweedie.tunicate_dry_biomass.12, lm.log.tunicate_dry_biomass, lm.tunicate_dry_biomass, glm.av.tunicate_dry_biomass, gam.av.lm.tunicate_dry_biomass.12, gam.av.gamma.tunicate_dry_biomass.12)

#gam.av.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.av.gamma.tunicate_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.gamma.tunicate_dry_biomass.12 )
qq_plot(gam.av.gamma.tunicate_dry_biomass.12 , method = 'simulate')
k.check(gam.av.gamma.tunicate_dry_biomass.12 )
summary(gam.av.gamma.tunicate_dry_biomass.12 )





gam.av.gamma.tunicate_dry_biomass.12.unordered<- gam(tunicate_dry_biomass+0.1 ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.av.tunicate_dry_biomass <- family(gam.av.gamma.tunicate_dry_biomass.12)
fam.gam.av.tunicate_dry_biomass
str(fam.gam.av.tunicate_dry_biomass)
ilink.gam.av.tunicate_dry_biomass<- fam.gam.av.tunicate_dry_biomass$linkinv
ilink.gam.av.tunicate_dry_biomass


mod.av.tunicate_dry_biomass<-gam.av.gamma.tunicate_dry_biomass.12
ndata.av.tunicate_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                             length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.tunicate_dry_biomass <- add_column(ndata.av.tunicate_dry_biomass, fit = predict(mod.av.tunicate_dry_biomass, newdata = ndata.av.tunicate_dry_biomass, type = 'response'))

predict(mod.av.tunicate_dry_biomass, newdata = ndata.av.tunicate_dry_biomass, type = 'response')
ndata.av.tunicate_dry_biomass <- bind_cols(ndata.av.tunicate_dry_biomass, setNames(as_tibble(predict(mod.av.tunicate_dry_biomass, ndata.av.tunicate_dry_biomass, se.fit = TRUE)[1:2]),
                                                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.tunicate_dry_biomass <- mutate(ndata.av.tunicate_dry_biomass,
                                     fit_resp  = ilink.gam.av.tunicate_dry_biomass(fit_link),
                                     right_upr = ilink.gam.av.tunicate_dry_biomass(fit_link + (2 * se_link)),
                                     right_lwr = ilink.gam.av.tunicate_dry_biomass(fit_link - (2 * se_link)))


ndata.av.tunicate_dry_biomass$av.pH.unscaled<-ndata.av.tunicate_dry_biomass$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.tunicate_dry_biomass <- ggplot(ndata.av.tunicate_dry_biomass, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (tunicate_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(italic("Botryllus") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.tunicate_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.tunicate_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//tunicate_dry_biomass_pred_av.png")



# caprellid biomass -------------------------------------------------------


gamma.12.caprellid_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$caprellid_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$caprellid_dry_biomass, "gamma", shape = gamma.12.caprellid_dry_biomass$estimate[[1]], rate = gamma.12.caprellid_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$caprellid_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$caprellid_dry_biomass))


qqp(food.exp.data.12.2019_zscores$caprellid_dry_biomass, "lnorm")




lm.caprellid_dry_biomass<-lm(caprellid_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.caprellid_dry_biomass<-lm(log(caprellid_dry_biomass+0.1) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.caprellid_dry_biomass<-glm(caprellid_dry_biomass+0.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.av.lm.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.caprellid_dry_biomass.12<- gam(log(caprellid_dry_biomass+0.1) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.caprellid_dry_biomass.12<- gam(caprellid_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.caprellid_dry_biomass.12, gam.av.lm.log.caprellid_dry_biomass.12, gam.av.tweedie.caprellid_dry_biomass.12, lm.log.caprellid_dry_biomass, lm.caprellid_dry_biomass, glm.av.caprellid_dry_biomass, gam.av.lm.caprellid_dry_biomass.12, gam.av.gamma.caprellid_dry_biomass.12)

#gam.av.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.av.gamma.caprellid_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.gamma.caprellid_dry_biomass.12 )
qq_plot(gam.av.gamma.caprellid_dry_biomass.12 , method = 'simulate')
k.check(gam.av.gamma.caprellid_dry_biomass.12 )
summary(gam.av.gamma.caprellid_dry_biomass.12 )





gam.av.gamma.caprellid_dry_biomass.12.unordered<- gam(caprellid_dry_biomass+0.1 ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")


fam.gam.av.caprellid_dry_biomass <- family(gam.av.gamma.caprellid_dry_biomass.12)
fam.gam.av.caprellid_dry_biomass
str(fam.gam.av.caprellid_dry_biomass)
ilink.gam.av.caprellid_dry_biomass<- fam.gam.av.caprellid_dry_biomass$linkinv
ilink.gam.av.caprellid_dry_biomass


mod.av.caprellid_dry_biomass<-gam.av.gamma.caprellid_dry_biomass.12
ndata.av.caprellid_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                              length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.caprellid_dry_biomass <- add_column(ndata.av.caprellid_dry_biomass, fit = predict(mod.av.caprellid_dry_biomass, newdata = ndata.av.caprellid_dry_biomass, type = 'response'))

predict(mod.av.caprellid_dry_biomass, newdata = ndata.av.caprellid_dry_biomass, type = 'response')
ndata.av.caprellid_dry_biomass <- bind_cols(ndata.av.caprellid_dry_biomass, setNames(as_tibble(predict(mod.av.caprellid_dry_biomass, ndata.av.caprellid_dry_biomass, se.fit = TRUE)[1:2]),
                                                                               c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.caprellid_dry_biomass <- mutate(ndata.av.caprellid_dry_biomass,
                                      fit_resp  = ilink.gam.av.caprellid_dry_biomass(fit_link),
                                      right_upr = ilink.gam.av.caprellid_dry_biomass(fit_link + (2 * se_link)),
                                      right_lwr = ilink.gam.av.caprellid_dry_biomass(fit_link - (2 * se_link)))


ndata.av.caprellid_dry_biomass$av.pH.unscaled<-ndata.av.caprellid_dry_biomass$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.caprellid_dry_biomass <- ggplot(ndata.av.caprellid_dry_biomass, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (caprellid_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(italic("Caprella") ~ "dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.caprellid_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.caprellid_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//caprellid_dry_biomass_pred_av.png")



# rest biomass ------------------------------------------------------------


gamma.12.rest_dry_biomass<-fitdistr(food.exp.data.12.2019_zscores$rest_dry_biomass+0.01, "gamma")
qqp(food.exp.data.12.2019$rest_dry_biomass, "gamma", shape = gamma.12.rest_dry_biomass$estimate[[1]], rate = gamma.12.rest_dry_biomass$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$rest_dry_biomass)
#normal good!!! 
qqp(log(food.exp.data.12.2019_zscores$rest_dry_biomass))


qqp(food.exp.data.12.2019_zscores$rest_dry_biomass, "lnorm")




lm.rest_dry_biomass<-lm(rest_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.log.rest_dry_biomass<-lm(log(rest_dry_biomass+0.1) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.rest_dry_biomass<-glm(rest_dry_biomass+0.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)

gam.av.lm.rest_dry_biomass.12<- gam(rest_dry_biomass ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.rest_dry_biomass.12<- gam(log(rest_dry_biomass+0.1) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.rest_dry_biomass.12<- gam(rest_dry_biomass+0.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.rest_dry_biomass.12, gam.av.lm.log.rest_dry_biomass.12, gam.av.tweedie.rest_dry_biomass.12, lm.log.rest_dry_biomass, lm.rest_dry_biomass, glm.av.rest_dry_biomass, gam.av.lm.rest_dry_biomass.12, gam.av.gamma.rest_dry_biomass.12)

#gam.av.lm.loglink is best but lm is only 0.1 so going with simpler.... 


plot(gam.av.gamma.rest_dry_biomass.12 , shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.gamma.rest_dry_biomass.12 )
qq_plot(gam.av.gamma.rest_dry_biomass.12 , method = 'simulate')
k.check(gam.av.gamma.rest_dry_biomass.12 )
summary(gam.av.gamma.rest_dry_biomass.12 )

###k check is significant




gam.av.gamma.rest_dry_biomass.12.unordered<- gam(rest_dry_biomass+0.1 ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
summary(gam.av.gamma.rest_dry_biomass.12.unordered)

fam.gam.av.rest_dry_biomass <- family(gam.av.gamma.rest_dry_biomass.12)
fam.gam.av.rest_dry_biomass
str(fam.gam.av.rest_dry_biomass)
ilink.gam.av.rest_dry_biomass<- fam.gam.av.rest_dry_biomass$linkinv
ilink.gam.av.rest_dry_biomass


mod.av.rest_dry_biomass<-gam.av.gamma.rest_dry_biomass.12
ndata.av.rest_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                         length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.rest_dry_biomass <- add_column(ndata.av.rest_dry_biomass, fit = predict(mod.av.rest_dry_biomass, newdata = ndata.av.rest_dry_biomass, type = 'response'))

predict(mod.av.rest_dry_biomass, newdata = ndata.av.rest_dry_biomass, type = 'response')
ndata.av.rest_dry_biomass <- bind_cols(ndata.av.rest_dry_biomass, setNames(as_tibble(predict(mod.av.rest_dry_biomass, ndata.av.rest_dry_biomass, se.fit = TRUE)[1:2]),
                                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.rest_dry_biomass <- mutate(ndata.av.rest_dry_biomass,
                                 fit_resp  = ilink.gam.av.rest_dry_biomass(fit_link),
                                 right_upr = ilink.gam.av.rest_dry_biomass(fit_link + (2 * se_link)),
                                 right_lwr = ilink.gam.av.rest_dry_biomass(fit_link - (2 * se_link)))


ndata.av.rest_dry_biomass$av.pH.unscaled<-ndata.av.rest_dry_biomass$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.rest_dry_biomass <- ggplot(ndata.av.rest_dry_biomass, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = (rest_dry_biomass+0.01), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression("Remaining dry weight per tile (g)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.rest_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.rest_dry_biomass
ggplot2::ggsave("C:Graphs May 2019//rest_dry_biomass_pred_av.png")




# Mussel wet weight -------------------------------------------------------


gamma.12.Mussel.wet.weight<-fitdistr(food.exp.data.12.2019_zscores$Mussel.wet.weight+0.01, "gamma")
qqp(food.exp.data.12.2019$Mussel.wet.weight, "gamma", shape = gamma.12.Mussel.wet.weight$estimate[[1]], rate = gamma.12.Mussel.wet.weight$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight)

qqp(log(food.exp.data.12.2019_zscores$Mussel.wet.weight+1))
#log normal the best

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight, "lnorm")


#none of these look good 
head(food.exp.data.12.2019_zscores)


lm.av.Mussel.wet.weight<-lm(Mussel.wet.weight ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.av.log.Mussel.wet.weight<-lm(log(Mussel.wet.weight+1) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

#glm.av.Mussel.wet.weight<-glm(Mussel.wet.weight+1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
#glm.av.loglink.Mussel.wet.weight<-glm(Mussel.wet.weight ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.av.lm.Mussel.wet.weight.12<- gam(Mussel.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight+1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.Mussel.wet.weight.12<- gam(log(Mussel.wet.weight+1) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.Mussel.wet.weight.12.1<- gam(Mussel.wet.weight ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.Mussel.wet.weight.12.1, gam.av.lm.log.Mussel.wet.weight.12, gam.av.tweedie.Mussel.wet.weight.12.1, lm.av.log.Mussel.wet.weight, lm.av.Mussel.wet.weight, gam.av.lm.Mussel.wet.weight.12, gam.av.gamma.Mussel.wet.weight.12.1)

#gam.av.lm.log.Mussel.wet.weight.12 is best 


plot(gam.av.lm.log.Mussel.wet.weight.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.log.Mussel.wet.weight.12)
qq_plot(gam.av.lm.log.Mussel.wet.weight.12, method = 'simulate')
k.check(gam.av.lm.log.Mussel.wet.weight.12)
summary(gam.av.lm.log.Mussel.wet.weight.12)

#residuals many in a straight line but otherwise good 



gam.av.lm.log.Mussel.wet.weight.12.unordered<- gam(log(Mussel.wet.weight+1) ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.av.Mussel.wet.weight <- family(gam.av.lm.log.Mussel.wet.weight.12)
fam.gam.av.Mussel.wet.weight
str(fam.gam.av.Mussel.wet.weight)
ilink.gam.av.Mussel.wet.weight<- fam.gam.av.Mussel.wet.weight$linkinv
ilink.gam.av.Mussel.wet.weight


mod.av.Mussel.wet.weight<-gam.av.lm.log.Mussel.wet.weight.12
ndata.av.Mussel.wet.weight <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                          length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.Mussel.wet.weight <- add_column(ndata.av.Mussel.wet.weight, fit = predict(mod.av.Mussel.wet.weight, newdata = ndata.av.Mussel.wet.weight, type = 'response'))

predict(mod.av.Mussel.wet.weight, newdata = ndata.av.Mussel.wet.weight, type = 'response')
ndata.av.Mussel.wet.weight <- bind_cols(ndata.av.Mussel.wet.weight, setNames(as_tibble(predict(mod.av.Mussel.wet.weight, ndata.av.Mussel.wet.weight, se.fit = TRUE)[1:2]),
                                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.Mussel.wet.weight <- mutate(ndata.av.Mussel.wet.weight,
                                  fit_resp  = ilink.gam.av.Mussel.wet.weight(fit_link),
                                  right_upr = ilink.gam.av.Mussel.wet.weight(fit_link + (2 * se_link)),
                                  right_lwr = ilink.gam.av.Mussel.wet.weight(fit_link - (2 * se_link)))


ndata.av.Mussel.wet.weight$av.pH.unscaled<-ndata.av.Mussel.wet.weight$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.Mussel.wet.weight <- ggplot(ndata.av.Mussel.wet.weight, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = log(Mussel.wet.weight+1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(italic("Mytilus") ~"wet weight (g, log scale)"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.Mussel.wet.weight,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.Mussel.wet.weight
ggsave("C:Graphs May 2019//Mussel.wet.weight_pred_av.png")






# Mussel wet weight per individual mussel ---------------------------------
#big outlier .... 5 mussels weigh 105.5 g?? so 17 g per mussel? 
View(food.exp.data.12.2019_zscores)

hist(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1)

food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1<-(food.exp.data.12.2019_zscores$Mussel.wet.weight)/(food.exp.data.12.2019_zscores$mussel_complete+1)

gamma.12.Mussel.wet.weight.per.1<-fitdistr(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1+0.01, "gamma")
qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1, "gamma", shape = gamma.12.Mussel.wet.weight.per.1$estimate[[1]], rate = gamma.12.Mussel.wet.weight.per.1$estimate[[2]])

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1)

qqp(log(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1+1))
#log normal the best

qqp(food.exp.data.12.2019_zscores$Mussel.wet.weight.per.1, "lnorm")


#none of these look good 
head(food.exp.data.12.2019_zscores)


lm.av.Mussel.wet.weight.per.1<-lm(Mussel.wet.weight.per.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.av.log.Mussel.wet.weight.per.1<-lm(log(Mussel.wet.weight.per.1+0.01) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

#glm.av.Mussel.wet.weight.per.1<-glm(Mussel.wet.weight.per.1+0.01 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=Gamma)
#glm.av.loglink.Mussel.wet.weight.per.1<-glm(Mussel.wet.weight.per.1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))

gam.av.lm.Mussel.wet.weight.per.1.12<- gam(Mussel.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.gamma.Mussel.wet.weight.per.1.12<- gam(Mussel.wet.weight.per.1+0.01 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")
gam.av.lm.log.Mussel.wet.weight.per.1.12<- gam(log(Mussel.wet.weight.per.1+0.01) ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.tweedie.Mussel.wet.weight.per.1.12.1<- gam(Mussel.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = tw, select=TRUE, method="REML")
gam.av.loglink.Mussel.wet.weight.per.1.12.1<- gam(Mussel.wet.weight.per.1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.Mussel.wet.weight.per.1.12.1, gam.av.lm.log.Mussel.wet.weight.per.1.12, gam.av.tweedie.Mussel.wet.weight.per.1.12.1, lm.av.log.Mussel.wet.weight.per.1, lm.av.Mussel.wet.weight.per.1, gam.av.lm.Mussel.wet.weight.per.1.12, gam.av.gamma.Mussel.wet.weight.per.1.12)

#gamma is best but produces very weird plkot so going with loglink


plot(gam.av.gamma.Mussel.wet.weight.per.1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.gamma.Mussel.wet.weight.per.1.12)
qq_plot(gam.av.gamma.Mussel.wet.weight.per.1.12, method = 'simulate')
k.check(gam.av.gamma.Mussel.wet.weight.per.1.12)
summary(gam.av.gamma.Mussel.wet.weight.per.1.12)

#don't look great... 


gam.av.gamma.Mussel.wet.weight.per.1.12.unordered<- gam(Mussel.wet.weight.per.1+0.01 ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = Gamma, select=TRUE, method="REML")

fam.gam.av.Mussel.wet.weight.per.1 <- family(gam.av.gamma.Mussel.wet.weight.per.1.12)
fam.gam.av.Mussel.wet.weight.per.1
str(fam.gam.av.Mussel.wet.weight.per.1)
ilink.gam.av.Mussel.wet.weight.per.1<- fam.gam.av.Mussel.wet.weight.per.1$linkinv
ilink.gam.av.Mussel.wet.weight.per.1


mod.av.Mussel.wet.weight.per.1<-gam.av.gamma.Mussel.wet.weight.per.1.12
ndata.av.Mussel.wet.weight.per.1 <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.Mussel.wet.weight.per.1 <- add_column(ndata.av.Mussel.wet.weight.per.1, fit = predict(mod.av.Mussel.wet.weight.per.1, newdata = ndata.av.Mussel.wet.weight.per.1, type = 'response'))

predict(mod.av.Mussel.wet.weight.per.1, newdata = ndata.av.Mussel.wet.weight.per.1, type = 'response')
ndata.av.Mussel.wet.weight.per.1 <- bind_cols(ndata.av.Mussel.wet.weight.per.1, setNames(as_tibble(predict(mod.av.Mussel.wet.weight.per.1, ndata.av.Mussel.wet.weight.per.1, se.fit = TRUE)[1:2]),
                                                                                   c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.Mussel.wet.weight.per.1 <- mutate(ndata.av.Mussel.wet.weight.per.1,
                                        fit_resp  = ilink.gam.av.Mussel.wet.weight.per.1(fit_link),
                                        right_upr = ilink.gam.av.Mussel.wet.weight.per.1(fit_link + (2 * se_link)),
                                        right_lwr = ilink.gam.av.Mussel.wet.weight.per.1(fit_link - (2 * se_link)))


ndata.av.Mussel.wet.weight.per.1$av.pH.unscaled<-ndata.av.Mussel.wet.weight.per.1$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.Mussel.wet.weight.per.1 <- ggplot(ndata.av.Mussel.wet.weight.per.1, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = Mussel.wet.weight.per.1, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Log biomass per mussel \n(wet weight)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.Mussel.wet.weight.per.1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.Mussel.wet.weight.per.1
ggsave("C:Graphs May 2019//Mussel.wet.weight.per.1_pred_av.png")




# Hydtobot ----------------------------------------------------------------
food.exp.data.12.2019_zscores$hydtobot[food.exp.data.12.2019_zscores$hydtobot==1]<-0.99
food.exp.data.12.2019_zscores$hydtobot[food.exp.data.12.2019_zscores$hydtobot==0]<-0.01

beta.12<-fitdist(0.01*(food.exp.data.12.2019$hydtobot), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$hydtobot), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

lm.av.hydtobot<-lm(hydtobot ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.av.log.hydtobot<-lm(log(hydtobot+0.01) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.loglink.hydtobot<-glm(hydtobot ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))


gam.av.lm.hydtobot.12<- gam(log(hydtobot+0.01)~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.log.lm.hydtobot.12<- gam(hydtobot~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


gam.av.tweedie.hydtobot.12<- gam(hydtobot~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family =tw(), select=TRUE, method="REML")

gam.av.beta.hydtobot.12<- gam(hydtobot~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.hydtobot.12.1<- gam(hydtobot~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.hydtobot.12.2<- gam(hydtobot~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(gam.av.tweedie.hydtobot.12 ,gam.av.beta.hydtobot.12, gam.av.lm.hydtobot.12, gam.av.log.lm.hydtobot.12, gam.av.beta.hydtobot.12.1, gam.av.beta.hydtobot.12.2, lm.av.hydtobot, lm.av.log.hydtobot)
#logit REML and cloglog same ... just logit.... 

### logit is the best 

plot(gam.av.beta.hydtobot.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.hydtobot.12)
qq_plot(gam.av.beta.hydtobot.12, method = 'simulate')
k.check(gam.av.beta.hydtobot.12)
gam.check(gam.av.beta.hydtobot.12)
summary(gam.av.beta.hydtobot.12)


#not the best qq plot
#resids are in rows

gam.av.beta.hydtobot.12.unordered<- gam(hydtobot~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")



fam.gam.av.hydtobot <- family(gam.av.beta.hydtobot.12)
fam.gam.av.hydtobot 
str(fam.gam.av.hydtobot )
ilink.gam.av.hydtobot <- fam.gam.av.hydtobot$linkinv
ilink.gam.av.hydtobot


food.exp.data.12.2019_zscores$av.pH.unscaled <-food.exp.data.12.2019_zscores$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')
head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

mod.av.hydtobot<-gam.av.beta.hydtobot.12
ndata.av.hydtobot <- with(food.exp.data.12.2019_zscores, data_frame(av.pH= seq(min(av.pH), max(av.pH),
                                                                                length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.hydtobot <- add_column(ndata.av.hydtobot, fit = predict(mod.av.hydtobot, newdata = ndata.av.hydtobot, type = 'response'))


ndata.av.hydtobot <- bind_cols(ndata.av.hydtobot, setNames(as_tibble(predict(mod.av.hydtobot, ndata.av.hydtobot, se.fit = TRUE)[1:2]),
                                                     c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.hydtobot <- mutate(ndata.av.hydtobot,
                         fit_resp  = ilink.gam.av.hydtobot(fit_link),
                         right_upr = ilink.gam.av.hydtobot(fit_link + (2 * se_link)),
                         right_lwr = ilink.gam.av.hydtobot(fit_link - (2 * se_link)))




ndata.av.hydtobot$av.pH.unscaled<-ndata.av.hydtobot$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.gam.av.hydtobot <- ggplot(ndata.av.hydtobot, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = hydtobot, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "cover ratio"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.hydtobot,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)+coord_cartesian(ylim = c(0, 1)) 
plt.av.gam.av.hydtobot 
ggsave("C:For submission//hydtobot_pred_av.png")







# Hydtobot_dry_biomass ----------------------------------------------------
food.exp.data.12.2019_zscores$hydtobot_dry_biomass[food.exp.data.12.2019_zscores$hydtobot_dry_biomass==1]<-0.99
food.exp.data.12.2019_zscores$hydtobot_dry_biomass[food.exp.data.12.2019_zscores$hydtobot_dry_biomass==0]<-0.01

beta.12<-fitdist(0.01*(food.exp.data.12.2019$hydtobot_dry_biomass), "beta", start=NULL)
qqp(0.01*(food.exp.data.12.2019$hydtobot_dry_biomass), "beta", shape1 = beta.12$estimate[[1]], shape2 = beta.12$estimate[[2]])

lm.av.hydtobot_dry_biomass<-lm(hydtobot_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)
lm.av.log.hydtobot_dry_biomass<-lm(log(hydtobot_dry_biomass+0.01) ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

glm.av.loglink.hydtobot_dry_biomass<-glm(hydtobot_dry_biomass ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores, family=gaussian(link="log"))


gam.av.lm.hydtobot_dry_biomass.12<- gam(log(hydtobot_dry_biomass+0.01)~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.log.lm.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


gam.av.tweedie.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family =tw(), select=TRUE, method="REML")

gam.av.beta.hydtobot_dry_biomass.12<- gam(hydtobot_dry_biomass~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")
gam.av.beta.hydtobot_dry_biomass.12.1<- gam(hydtobot_dry_biomass~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="probit"), select=TRUE, method="REML")
gam.av.beta.hydtobot_dry_biomass.12.2<- gam(hydtobot_dry_biomass~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="cloglog"), select=TRUE, method="REML")
#cauchit doesn't run with REML

AICtab(gam.av.tweedie.hydtobot_dry_biomass.12 ,gam.av.beta.hydtobot_dry_biomass.12, gam.av.lm.hydtobot_dry_biomass.12, gam.av.log.lm.hydtobot_dry_biomass.12, gam.av.beta.hydtobot_dry_biomass.12.1, gam.av.beta.hydtobot_dry_biomass.12.2, lm.av.hydtobot_dry_biomass, lm.av.log.hydtobot_dry_biomass)
#logit REML and cloglog same ... just logit.... 

### logit is the best 

plot(gam.av.beta.hydtobot_dry_biomass.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.beta.hydtobot_dry_biomass.12)
qq_plot(gam.av.beta.hydtobot_dry_biomass.12, method = 'simulate')
k.check(gam.av.beta.hydtobot_dry_biomass.12)
gam.check(gam.av.beta.hydtobot_dry_biomass.12)
summary(gam.av.beta.hydtobot_dry_biomass.12)


#not the best qq plot
#resids are in rows

gam.av.beta.hydtobot_dry_biomass.12.unordered<- gam(hydtobot_dry_biomass~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality), data = food.exp.data.12.2019_zscores, family = betar(link="logit"), select=TRUE, method="REML")



fam.gam.av.hydtobot_dry_biomass <- family(gam.av.beta.hydtobot_dry_biomass.12)
fam.gam.av.hydtobot_dry_biomass 
str(fam.gam.av.hydtobot_dry_biomass )
ilink.gam.av.hydtobot_dry_biomass <- fam.gam.av.hydtobot_dry_biomass$linkinv
ilink.gam.av.hydtobot_dry_biomass


food.exp.data.12.2019_zscores$av.pH.unscaled <-food.exp.data.12.2019_zscores$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')
head(food.exp.data.12.2019_zscores)

want <- seq(1, nrow(food.exp.data.12.2019_zscores), length.out = 100)

mod.av.hydtobot_dry_biomass<-gam.av.beta.hydtobot_dry_biomass.12
ndata.av.hydtobot_dry_biomass <- with(food.exp.data.12.2019_zscores, data_frame(av.pH= seq(min(av.pH), max(av.pH),
                                                                               length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.hydtobot_dry_biomass <- add_column(ndata.av.hydtobot_dry_biomass, fit = predict(mod.av.hydtobot_dry_biomass, newdata = ndata.av.hydtobot_dry_biomass, type = 'response'))


ndata.av.hydtobot_dry_biomass <- bind_cols(ndata.av.hydtobot_dry_biomass, setNames(as_tibble(predict(mod.av.hydtobot_dry_biomass, ndata.av.hydtobot_dry_biomass, se.fit = TRUE)[1:2]),
                                                           c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.hydtobot_dry_biomass <- mutate(ndata.av.hydtobot_dry_biomass,
                            fit_resp  = ilink.gam.av.hydtobot_dry_biomass(fit_link),
                            right_upr = ilink.gam.av.hydtobot_dry_biomass(fit_link + (2 * se_link)),
                            right_lwr = ilink.gam.av.hydtobot_dry_biomass(fit_link - (2 * se_link)))




ndata.av.hydtobot_dry_biomass$av.pH.unscaled<-ndata.av.hydtobot_dry_biomass$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')

# plot 

plt.av.gam.av.hydtobot_dry_biomass <- ggplot(ndata.av.hydtobot_dry_biomass, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y = hydtobot_dry_biomass, shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab(expression(italic("Botryllus")~ "to" ~ italic("Obelia") ~ "biomass ratio"))+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.hydtobot_dry_biomass,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")+ geom_hline(yintercept=0.5, linetype="dashed", color="black", size=1)+coord_cartesian(ylim = c(0, 1)) 
plt.av.gam.av.hydtobot_dry_biomass 
ggsave("C:For submission//hydtobot_dry_biomass_pred_av.png")



# CAP1 --------------------------------------------------------------------

qqp(food.exp.data.12.2019_zscores$CAP1)


qqp(food.exp.data.12.2019_zscores$CAP1, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.av.CAP1<-lm(CAP1 ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

gam.av.lm.CAP1.12<- gam(CAP1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.loglink.CAP1.12.1<- gam(CAP1 ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab( lm.av.CAP1,  gam.av.lm.CAP1.12)

#gam.av.lm.CAP1


plot(gam.av.lm.CAP1.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.CAP1.12)
qq_plot(gam.av.lm.CAP1.12, method = 'simulate')
k.check(gam.av.lm.CAP1.12)
summary(gam.av.lm.CAP1.12)





gam.av.lm.CAP1.12.unordered<- gam(CAP1 ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")


fam.gam.av.CAP1 <- family(gam.av.lm.CAP1.12)
fam.gam.av.CAP1
str(fam.gam.av.CAP1)
ilink.gam.av.CAP1<- fam.gam.av.CAP1$linkinv
ilink.gam.av.CAP1


mod.av.CAP1<-gam.av.lm.CAP1.12
ndata.av.CAP1 <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                             length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.CAP1 <- add_column(ndata.av.CAP1, fit = predict(mod.av.CAP1, newdata = ndata.av.CAP1, type = 'response'))

predict(mod.av.CAP1, newdata = ndata.av.CAP1, type = 'response')
ndata.av.CAP1 <- bind_cols(ndata.av.CAP1, setNames(as_tibble(predict(mod.av.CAP1, ndata.av.CAP1, se.fit = TRUE)[1:2]),
                                             c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.CAP1 <- mutate(ndata.av.CAP1,
                     fit_resp  = ilink.gam.av.CAP1(fit_link),
                     right_upr = ilink.gam.av.CAP1(fit_link + (2 * se_link)),
                     right_lwr = ilink.gam.av.CAP1(fit_link - (2 * se_link)))


ndata.av.CAP1$av.pH.unscaled<-ndata.av.CAP1$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.CAP1 <- ggplot(ndata.av.CAP1, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y =(CAP1), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Partial-dbRDA axis 1\n(36% of constrained variation)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.CAP1,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.CAP1
ggsave("C:Graphs May 2019//CAP1_pred_av.png")




# Distances ---------------------------------------------------------------

qqp(food.exp.data.12.2019_zscores$distances)


qqp(food.exp.data.12.2019_zscores$distances, "lnorm")

head(food.exp.data.12.2019_zscores)


lm.av.distances<-lm(distances ~ av.pH*oFood.quality, data=food.exp.data.12.2019_zscores)

gam.av.lm.distances.12<- gam(distances ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
gam.av.loglink.distances.12.1<- gam(distances ~ s(av.pH)+ oFood.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, family = gaussian(link="log"), select=TRUE, method="REML")


AICtab(gam.av.loglink.distances.12.1,lm.av.distances,  gam.av.lm.distances.12)

#gam.av.lm.distances


plot(gam.av.lm.distances.12, shade = TRUE, pages = 1, scale = 0, seWithMean = TRUE)
appraise(gam.av.lm.distances.12)
qq_plot(gam.av.lm.distances.12, method = 'simulate')
k.check(gam.av.lm.distances.12)
summary(gam.av.lm.distances.12)





gam.av.lm.distances.12.unordered<- gam(distances ~ s(av.pH)+ Food.quality + s(av.pH, by=oFood.quality),data = food.exp.data.12.2019_zscores, select=TRUE, method="REML")
summary(gam.av.lm.distances.12.unordered)

fam.gam.av.distances <- family(gam.av.lm.distances.12)
fam.gam.av.distances
str(fam.gam.av.distances)
ilink.gam.av.distances<- fam.gam.av.distances$linkinv
ilink.gam.av.distances


mod.av.distances<-gam.av.lm.distances.12
ndata.av.distances <- with(food.exp.data.12.2019_zscores, data_frame(av.pH = seq(min(av.pH), max(av.pH),
                                                                                  length = 100),  oFood.quality = oFood.quality[want],  CO2= CO2[want]))


## add the fitted values by predicting from the model for the new data
ndata.av.distances <- add_column(ndata.av.distances, fit = predict(mod.av.distances, newdata = ndata.av.distances, type = 'response'))

predict(mod.av.distances, newdata = ndata.av.distances, type = 'response')
ndata.av.distances <- bind_cols(ndata.av.distances, setNames(as_tibble(predict(mod.av.distances, ndata.av.distances, se.fit = TRUE)[1:2]),
                                                       c('fit_link','se_link')))

## create the interval and backtransform

ndata.av.distances <- mutate(ndata.av.distances,
                          fit_resp  = ilink.gam.av.distances(fit_link),
                          right_upr = ilink.gam.av.distances(fit_link + (2 * se_link)),
                          right_lwr = ilink.gam.av.distances(fit_link - (2 * se_link)))


ndata.av.distances$av.pH.unscaled<-ndata.av.distances$av.pH * attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:scale') + attr(food.exp.data.12.2019_zscores$av.pH, 'scaled:center')


# plot 
plt.av.distances <- ggplot(ndata.av.distances, aes(x = av.pH.unscaled, y = fit)) + 
  theme_classic()+
  geom_line(size=1.5, aes(colour=oFood.quality)) +
  geom_point(aes(y =(distances), shape=CO2, colour=oFood.quality), size=3, data = food.exp.data.12.2019_zscores)+
  xlab(expression("Average pH")) + ylab("Homogeneity of multivariate dispersions\n(distance to multivariate centroid)")+  
  scale_color_manual(values=colorset2)+
  scale_fill_manual(values=colorset2)+
  scale_shape_manual(values=c(19,17))+
  geom_ribbon(data = ndata.av.distances,aes(ymin = right_lwr, ymax = right_upr, fill=oFood.quality), alpha = 0.10)+
  theme(legend.position="none")
plt.av.distances
ggsave("C:Graphs May 2019//distances_pred_av.png")




# Community plotting ------------------------------------------------------
library(cowplot)
fig.av.biomass<-plot_grid(plt.av.total_dry_biomass,plt.av.everything.wet.weight ,plt.av.hydroid_dry_biomass,
                           plt.av.tunicate_dry_biomass,plt.av.Mussel.wet.weight,ncol=5, align='v', 
                           labels=c('(a)', '(b)','(c)', '(d)', '(e)', 
                                    label_size=12))

fig.av.biomass

ggplot2::ggsave("C:For submission//Supplemental//Fig.av.biomass.png", width=65, height=10, units="cm")


fig.av.community<-plot_grid(plt.av.richness, plt.av.evenness, plt.av.occupied.space,plt.av.CAP1, plt.av.distances,ncol=5, align='v', 
                         labels=c('(a)', '(b)','(c)', '(d)', '(e)', 
                                  label_size=12))

fig.av.community

ggplot2::ggsave("C:For submission//Fig.av.community.png", width=65, height=10, units="cm")


# Community level tables --------------------------------------------------


richness.gam.av<- summary(gam.av.poisson.richness.12)
richness.gam.av.unordered<- summary(gam.av.poisson.richness.12.unordered)

evenness.gam.av<-summary(gam.av.lm.evenness.12)
evenness.gam.av.unordered<-summary(gam.av.lm.evenness.12.unordered)


occupied.space.gam.av<- summary(gam.av.beta.occupied.space.12)
occupied.space.gam.av.unordered<- summary(gam.av.beta.occupied.space.12.unordered)

distances.gam.av<- summary(gam.av.lm.distances.12)
distances.gam.av.unordered<- summary(gam.av.lm.distances.12.unordered)


CAP1.gam.av <- summary(gam.av.lm.CAP1.12)
CAP1.gam.av.unordered <- summary(gam.av.lm.CAP1.12.unordered)

#dry biomass
total_dry_biomass.gam.av <- summary(gam.av.lm.total_dry_biomass.12)
total_dry_biomass.gam.av.unordered <- summary(gam.av.lm.total_dry_biomass.12.unordered)

hydroid_dry_biomass.gam.av.unordered <- summary(gam.av.gamma.hydroid_dry_biomass.12.unordered) 
hydroid_dry_biomass.gam.av<- summary(gam.av.gamma.hydroid_dry_biomass.12) 

caprellid_dry_biomass.gam.av.unordered <- summary(gam.av.gamma.caprellid_dry_biomass.12.unordered)
caprellid_dry_biomass.gam.av <- summary(gam.av.gamma.caprellid_dry_biomass.12)

tunicate_dry_biomass.gam.av <- summary(gam.av.gamma.tunicate_dry_biomass.12)
tunicate_dry_biomass.gam.av.unordered <- summary(gam.av.gamma.tunicate_dry_biomass.12.unordered)

rest_dry_biomass.gam.av.unordered <- summary(gam.av.gamma.rest_dry_biomass.12.unordered)
rest_dry_biomass.gam.av <- summary(gam.av.gamma.rest_dry_biomass.12)


#wet biomass
everything.wet.weight.gam.av <-summary(gam.av.lm.log.everything.wet.weight.12)
everything.wet.weight.gam.av.unordered <- summary(gam.av.lm.log.everything.wet.weight.12.unordered)

everything.wet.weight.per.1.gam.av <-summary(gam.av.lm.log.everything.wet.weight.per.1.12)
everything.wet.weight.per.1.gam.av.unordered <- summary(gam.av.lm.log.everything.wet.weight.per.1.12.unordered)

Mussel.wet.weight.gam.av <- summary(gam.av.lm.log.Mussel.wet.weight.12)
Mussel.wet.weight.gam.av.unordered <-summary(gam.av.lm.log.Mussel.wet.weight.12.unordered)

Mussel.wet.weight.per.1.gam.av <- summary(gam.av.gamma.Mussel.wet.weight.per.1.12)
Mussel.wet.weight.per.1.gam.av.unordered <-summary(gam.av.gamma.Mussel.wet.weight.per.1.12.unordered)


#competition metric 
hydtobot.gam.av <- summary(gam.av.beta.hydtobot.12)
hydtobot.gam.av.unordered <- summary(gam.av.beta.hydtobot.12.unordered)


hydtobot_dry_biomass.gam.av <- summary(gam.av.beta.hydtobot_dry_biomass.12)
hydtobot_dry_biomass.gam.av.unordered <- summary(gam.av.beta.hydtobot_dry_biomass.12.unordered)

#ptable building
richness.gam.av.p.table<-as.data.frame(richness.gam.av.unordered$p.table)
richness.gam.av.s.table<-as.data.frame(richness.gam.av$s.table)

evenness.gam.av.p.table<-as.data.frame(evenness.gam.av.unordered$p.table)
evenness.gam.av.s.table<-as.data.frame(evenness.gam.av$s.table)

occupied.space.gam.av.p.table<-as.data.frame(occupied.space.gam.av.unordered$p.table)
occupied.space.gam.av.s.table<-as.data.frame(occupied.space.gam.av$s.table)

distances.gam.av.p.table<-as.data.frame(distances.gam.av.unordered$p.table)
distances.gam.av.s.table<-as.data.frame(distances.gam.av$s.table)

CAP1.gam.av.p.table<-as.data.frame(CAP1.gam.av.unordered$p.table)
CAP1.gam.av.s.table<-as.data.frame(CAP1.gam.av$s.table)

total_dry_biomass.gam.av.p.table<-as.data.frame(total_dry_biomass.gam.av.unordered$p.table)
total_dry_biomass.gam.av.s.table<-as.data.frame(total_dry_biomass.gam.av$s.table)

hydroid_dry_biomass.gam.av.p.table<-as.data.frame(hydroid_dry_biomass.gam.av.unordered$p.table)
hydroid_dry_biomass.gam.av.s.table<-as.data.frame(hydroid_dry_biomass.gam.av$s.table)

caprellid_dry_biomass.gam.av.p.table<-as.data.frame(caprellid_dry_biomass.gam.av.unordered$p.table)
caprellid_dry_biomass.gam.av.s.table<-as.data.frame(caprellid_dry_biomass.gam.av$s.table)

tunicate_dry_biomass.gam.av.p.table<-as.data.frame(tunicate_dry_biomass.gam.av.unordered$p.table)
tunicate_dry_biomass.gam.av.s.table<-as.data.frame(tunicate_dry_biomass.gam.av$s.table)

rest_dry_biomass.gam.av.p.table<-as.data.frame(rest_dry_biomass.gam.av.unordered$p.table)
rest_dry_biomass.gam.av.s.table<-as.data.frame(rest_dry_biomass.gam.av$s.table)

everything.wet.weight.gam.av.p.table<-as.data.frame(everything.wet.weight.gam.av.unordered$p.table)
everything.wet.weight.gam.av.s.table<-as.data.frame(everything.wet.weight.gam.av$s.table)

everything.wet.weight.per.1.gam.av.p.table<-as.data.frame(everything.wet.weight.per.1.gam.av.unordered$p.table)
everything.wet.weight.per.1.gam.av.s.table<-as.data.frame(everything.wet.weight.per.1.gam.av$s.table)

Mussel.wet.weight.gam.av.p.table<-as.data.frame(Mussel.wet.weight.gam.av.unordered$p.table)
Mussel.wet.weight.gam.av.s.table<-as.data.frame(Mussel.wet.weight.gam.av$s.table)

Mussel.wet.weight.per.1.gam.av.p.table<-as.data.frame(Mussel.wet.weight.per.1.gam.av.unordered$p.table)
Mussel.wet.weight.per.1.gam.av.s.table<-as.data.frame(Mussel.wet.weight.per.1.gam.av$s.table)

hydtobot.gam.av.p.table<-as.data.frame(hydtobot.gam.av.unordered$p.table)
hydtobot.gam.av.s.table<-as.data.frame(hydtobot.gam.av$s.table)

hydtobot_dry_biomass.gam.av.p.table<-as.data.frame(hydtobot_dry_biomass.gam.av.unordered$p.table)
hydtobot_dry_biomass.gam.av.s.table<-as.data.frame(hydtobot_dry_biomass.gam.av$s.table)


#richness.gam.av.p.table and  hydtobot.gam.av.p.table, is with z value 
colnames(richness.gam.av.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot.gam.av.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
colnames(hydtobot_dry_biomass.gam.av.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")

colnames(occupied.space.gam.av.p.table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")


#### Building the stats table
ptable.av.community.t<-rbind(richness.gam.av.p.table,
                          evenness.gam.av.p.table,
                          occupied.space.gam.av.p.table,
                          CAP1.gam.av.p.table,
                          distances.gam.av.p.table,
                          total_dry_biomass.gam.av.p.table,
                          everything.wet.weight.gam.av.p.table,
                          hydroid_dry_biomass.gam.av.p.table,
                          tunicate_dry_biomass.gam.av.p.table,
                          Mussel.wet.weight.gam.av.p.table,
                          hydtobot.gam.av.p.table, 
                          hydtobot_dry_biomass.gam.av.p.table
)


colnames(ptable.av.community.t) <- c("Estimate", "SE", "t", "p")
ptable.av.community.t$Factor<-rep(c("Intercept", "Low quality food", "High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

ptable.av.community.t %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (z)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (z)", 7,9) %>% 
  group_rows("Partial dbRDA (1st axis), normal",10,12) %>% 
  group_rows("Homogeneity of dispersions, normal", 13,15) %>% 
  group_rows("Total dry biomass, normal", 16,18) %>% 
  group_rows("Total wet biomass, normal (log)", 19,21) %>% 
  group_rows("Obelia dry biomass, gamma", 22,24) %>% 
  group_rows("Botryllus dry biomass, gamma", 25,27) %>% 
  group_rows("Mussel wet biomass, normal (log)", 28,30) %>% 
  group_rows("Botryllus to Obelia dominance ratio, beta (z)", 31,33) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (z)", 34,36) %>% 
  save_kable(file = "ptable.av.community.t.html", self_contained = T, pagetitle="ptable.av.community.t")

#again hydtobot and richness
#richness.gam.av.p.table and  hydtobot.gam.av.p.table, is with Chisq
colnames(richness.gam.av.s.table) <- c("edf", "Ref.df", "F", "p-value")
colnames(hydtobot.gam.av.s.table) <- c("edf", "Ref.df",  "F", "p-value")
colnames(hydtobot_dry_biomass.gam.av.s.table) <- c("edf", "Ref.df",  "F", "p-value")

colnames(occupied.space.gam.av.s.table) <- c("edf", "Ref.df",  "F", "p-value")


### s table
stable.av.community.f<-rbind(richness.gam.av.s.table,
                          evenness.gam.av.s.table,
                          occupied.space.gam.av.s.table,
                          CAP1.gam.av.s.table,
                          distances.gam.av.s.table,
                          total_dry_biomass.gam.av.s.table,
                          everything.wet.weight.gam.av.s.table,
                          hydroid_dry_biomass.gam.av.s.table,
                          tunicate_dry_biomass.gam.av.s.table,
                          Mussel.wet.weight.gam.av.s.table,
                          hydtobot.gam.av.s.table, 
                          hydtobot_dry_biomass.gam.av.s.table
)


colnames(stable.av.community.f) <- c("Estimated_df", "Reference_df", "F", "p_smooth")
stable.av.community.f$Smooth_terms<-rep(c("smooth av.pH", "smooth av.pH * Low quality food", "smooth av.pH * High quality food"))


#development of kable will make it so that modified cells can apply to round - i.e. after "cel_spec"

stable.av.community.f %>% 
  mutate_if(is.numeric, round, 4) %>% 
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p) %>% 
  kable(escape=F, digits=2) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (Chi square)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (Chi square)", 7,9) %>% 
  group_rows("Partial dbRDA (1st axis), normal",10,12) %>% 
  group_rows("Homogeneity of dispersions, normal", 13,15) %>% 
  group_rows("Total dry biomass, normal", 16,18) %>% 
  group_rows("Total wet biomass, normal (log)", 19,21) %>% 
  group_rows("Obelia dry biomass, gamma", 22,24) %>% 
  group_rows("Botryllus dry biomass, gamma", 25,27) %>% 
  group_rows("Mussel wet biomass, normal (log)", 28,30) %>% 
  group_rows("Botryllus to Obelia dominance ratio, beta (Chi square)", 31,33) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi square)", 34,36) %>% 
  save_kable(file = "stable.av.community.f.html", self_contained = T)




pstable.av.community<-cbind(ptable.av.community.t, stable.av.community.f)


pstable.av.community %>% 
  mutate_if(is.numeric, round, 3) %>% 
  mutate(p = ifelse(p<0.001, "<0.001",p)) %>%
  mutate(p_smooth = ifelse(p_smooth<0.001, "<0.001",p_smooth)) %>%
  mutate(p_smooth = cell_spec(p_smooth, bold = ifelse(p_smooth < 0.05, "TRUE", "FALSE"))) %>% 
  mutate(p = cell_spec(p, bold = ifelse(p < 0.051, "TRUE", "FALSE"))) %>% 
  dplyr::select(Smooth_terms, Estimated_df, Reference_df, F, p_smooth, Factor, Estimate, SE, t, p) %>% 
  kable(escape=F, digits=2, row.names = FALSE) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>%
  group_rows("Richness, poisson (Chi-square, z)", 1,3) %>% 
  group_rows("Evenness, normal", 4,6) %>%
  group_rows("Occupied space, beta (Chi-square, z)", 7,9) %>% 
  group_rows("Partial dbRDA (1st axis), normal",10,12) %>% 
  group_rows("Homogeneity of dispersions, normal", 13,15) %>% 
  group_rows("Total dry biomass, normal", 16,18) %>% 
  group_rows("Total wet biomass, normal (log)", 19,21) %>% 
  group_rows("Obelia dry biomass, gamma", 22,24) %>% 
  group_rows("Botryllus dry biomass, gamma", 25,27) %>% 
  group_rows("Mussel wet biomass, normal (log)", 28,30) %>% 
  group_rows("Botryllus to Obelia dominance ratio, beta (Chi-square, z)", 31,33) %>% 
  group_rows("Botryllus to Obelia dominance ratio by biomass, beta (Chi square, z)", 34,36) %>% 
  
  save_kable(file = "C:For submission//pstable.av.community.html", self_contained = T)







# old code ----------------------------------------------------------------















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

beta.hydroid.12.1<- betareg(y.transf.betareg(hydroid/100) ~ av.pH*Food.quality, data =  food.exp.data.12.2019_zscores)
summary(beta.hydroid.12.1)
plot(beta.hydroid.12.1)
coeftest(beta.hydroid.12.1)





plot(beta.hydroid.12.1)
plot(beta.hydroid.12.1, which = 6, type = "deviance", sub.caption = "")
plot(beta.hydroid.12.1, which = 4, type = "deviance", sub.caption = "")
plot(beta.hydroid.12.1, which = 2, type = "deviance", sub.caption = "")
plot(beta.hydroid.12.1, which = 5, type = "deviance", sub.caption = "")


beta.hydroid.12.2<- betareg(y.transf.betareg(hydroid/100) ~ av.pH*Food.quality, data =  food.exp.data.12.2019_zscores, link="loglog")
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

beta.hydroid.12.2<- betareg(y.transf.betareg(hydroid/100) ~ av.pH*Food.quality , data =  food.exp.data.12.2019_zscores)

sapply(c("logit", "probit", "cloglog", "cauchit", "loglog"), function(x) logLik(update(beta.hydroid.12.2, link = x)))

beta.hydroid.12.3<- betareg(y.transf.betareg(hydroid/100) ~ av.pH*Food.quality , data =  food.exp.data.12.2019_zscores, link="cauchit")
beta.hydroid.12.4<- betareg(y.transf.betareg(hydroid/100) ~ hydrogen.concentration*Food.quality , data =  food.exp.data.12.2019_zscores,  link="cauchit")


AIC(beta.hydroid.12.2, beta.hydroid.12.1, beta.hydroid.12.3, beta.hydroid.12.4)

#best one: 
beta.hydroid.12.3<- betareg(y.transf.betareg(hydroid/100) ~ av.pH*Food.quality , data =  food.exp.data.12.2019_zscores, link="cauchit")



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
lm.dry_total<-lm(dry_total~hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.dry_total)
summary(lm.dry_total)
Anova(lm.dry_total, type="3")

lm.dry_total.av<-lm(dry_total~hydrogen.concentration.av*Food.quality, data=food.exp.data.12.2019)
plot(lm.dry_total.av)
summary(lm.dry_total.av)
Anova(lm.dry_total.av, type="3")

lm.dry_total_per1<-lm(dry_total_per1~hydrogen.concentration*Food.quality, data=food.exp.data.12.2019)
plot(lm.dry_total_per1)
summary(lm.dry_total_per1)
Anova(lm.dry_total_per1, type="3")

glm.gamma.dry_total_per1.12<-glm(formula = (dry_total_per1) ~ hydrogen.concentration*Food.quality, data = food.exp.data.12.2019, family = "Gamma")
summary(glm.gamma.dry_total_per1.12)
plot(glm.gamma.dry_total_per1.12)
Anova(glm.gamma.dry_total_per1.12, type="3")

#hydroid binomial

#### by hydrogen concentration

contrasts(food.exp.data.12.2019$Food.quality)=contr.poly(3) 

options(contrasts = c("contr.treatment", "contr.poly"))
options(contrasts = c("contr.sum", "contr.poly"))


fitBinom=fitdist(data=food.exp.data.12.2019$hydroid, dist="binom", fix.arg=list(size=100), start=list(prob=0.01))
qqp(food.exp.data.12.2019$hydroid, "binom", size = 100, prob = fitBinom$estimate[[1]])


food.exp.data.12.2019$hydrogen.concentration<- 10^(-food.exp.data.12.2019$av.pH)

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


beta.formicula.12.1<- betareg(y.transf.betareg(food.exp.data.12.2019$formicula/food.exp.data.12.2019$total)  ~ av.pH*Food.quality, data =  food.exp.data.12.2019)
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

glm.binomial.bowerbankia.12<- glm(formula = cbind(bowerbankia, 100-bowerbankia)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = binomial(link="logit"))
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

glm.neg.binomial.hydroid.12<- glm.nb(hydroid~ av.pH*Food.quality, data = food.exp.data.12.2019)
plot(glm.neg.binomial.hydroid.12)
summary(glm.neg.binomial.hydroid.12)
Anova(glm.neg.binomial.hydroid.12, type=3)

library(mgcv)
library(car)
influencePlot(glm.binomial.hydroid.12.1.hydrogen)


qqp(food.exp.data.12.2019$hydroid, "norm")
hydroid.lm<-lm(hydroid~av.pH*Food.quality, data = food.exp.data.12.2019)
summary(hydroid.lm)
anova(hydroid.lm, test = "Chi")

####### try a gam
gam.av.hydroid.12<- gam(0.01*hydroid ~ s(av.pH, fx=FALSE, k=-1, bs="cr")*Food.quality, data = food.exp.data.12.2019_zscores, family = binomial)
plot(gam.av.hydroid.12)
summary(gam.av.hydroid.12)


anova(gam.av.hydroid.12, glm.binomial.hydroid.12.hydrogen)

#### quasibinomial
glm.quasi.hydroid.12<- glm(0.01*hydroid ~ av.pH*Food.quality, data = food.exp.data.12.2019_zscores, family = quasibinomial(link="probit"))
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

glm.overdispersed.hydroid<- glmer(formula = cbind(hydroid, 100-hydroid)~ av.pH*Food.quality + (1|Mesocosm), data = food.exp.data.12.2019_zscores, family =binomial("logit"))
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

beta.hydroid.12.1<- betareg(hydroid.2 ~ av.pH*Food.quality, data =  food.exp.data.12.2019)
summary(beta.hydroid.12.1)
plot(beta.hydroid.12)
plot(beta.hydroid.12, which = 6, type = "deviance", sub.caption = "")
plot(beta.hydroid.12, which = 4, type = "deviance", sub.caption = "")
plot(beta.hydroid.12, which = 1, type = "deviance", sub.caption = "")

gy_logit4 <- update(beta.hydroid.12.1, subset = -3)

coef(beta.hydroid.12, model = "precision")
coef(gy_logit4, model = "precision")


### beta regression with av.pH as an additional regressor for the precision parameterto account for heteroskedasticit
beta.hydroid.12.2<- betareg(hydroid.2 ~ av.pH*Food.quality|av.pH, data =  food.exp.data.12.2019)
beta.hydroid.12.3<- betareg(hydroid.2 ~ av.pH*Food.quality|Food.quality, data =  food.exp.data.12.2019)
beta.hydroid.12.4<- betareg(hydroid.2 ~ av.pH*Food.quality|Food.quality, data =  food.exp.data.12.2019, link="probit")
beta.hydroid.12.5<- betareg(hydroid.2 ~ av.pH*Food.quality|Food.quality, data =  food.exp.data.12.2019, link="loglog")


lrtest(beta.hydroid.12.1, beta.hydroid.12.2)
lrtest(beta.hydroid.12.1, beta.hydroid.12.5)


beta.hydroid.12.space<- betareg(hydroid.3 ~ av.pH*Food.quality, data =  food.exp.data.12.2019)
plot(beta.hydroid.12.space, which = 1:5, type = "deviance", sub.caption = "")
summary(beta.hydroid.12.space)
cooks.distance(beta.hydroid.12.space)
AIC(beta.hydroid.12.1)
AIC(beta.hydroid.12.6)


beta.hydroid.12.6<- betareg(hydroid.2 ~ av.pH*Food.quality, data =  food.exp.data.12.2019)

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
cp = betareg(av.pH~hydroid, food.exp.data.12.2019)
vp = betareg(x~f|f, df)
lrtest(cp, vp)

glm.quasi.hydroid.12<- glm(hydroid.2 ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = quasibinomial(link="probit"))
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
  -sum(dbetabinom(food.exp.data.12.2019$av.pH,prob,size,theta,log=TRUE))
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
glm.binomial.disporella.12<- glm(formula = cbind(disporella, 1.0001-disporella)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.disporella.12)
summary(glm.binomial.disporella.12)
anova(glm.binomial.disporella.12, test = "Chi")


### Number Disporella
nbinom12 <- fitdistr(food.exp.data.12.2019$num.disporella, "Negative Binomial")
qqp(food.exp.data.12.2019$num.disporella, "nbinom", size = nbinom12$estimate[[1]], mu = nbinom12$estimate[[2]])


poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella, "Poisson")
qqp(food.exp.data.12.2019$num.disporella, "pois", poisson.12$estimate)

glm.poisson.num.disporella.12<-glm(formula = (num.disporella) ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.12)
plot(glm.poisson.num.disporella.12)
anova(glm.poisson.num.disporella.12, test = "Chi")


glm.nb.num.disporella.12<-glm.nb(formula = (num.disporella) ~ av.pH*Food.quality, data = food.exp.data.12.2019)
summary(glm.nb.num.disporella.12)
anova(glm.nb.num.disporella.12, test = "Chi")


visreg(glm.poisson.num.disporella.12, xvar = "av.pH", by= "Food.quality", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))


### Proportion erect. all
food.exp.data.12.2019$num.prop.disporella.erect<-(food.exp.data.12.2019$num.disporella.erect.all)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.erect.12<- glm(formula = cbind(num.prop.disporella.erect, 1-num.prop.disporella.erect)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.12)

summary(glm.binomial.num.prop.disporella.erect.12)
anova(glm.binomial.num.prop.disporella.erect.12, test = "Chi")



### Proportion erect.knobs
food.exp.data.12.2019$num.prop.disporella.erect.knobs<-(food.exp.data.12.2019$num.disporella.erect.knobs)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.erect.knobs.12<- glm(formula = cbind(num.prop.disporella.erect.knobs, 1-num.prop.disporella.erect.knobs)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.erect.knobs.12)

summary(glm.binomial.num.prop.disporella.erect.knobs.12)
anova(glm.binomial.num.prop.disporella.erect.knobs.12, test = "Chi")

### Proportion flat. all
food.exp.data.12.2019$num.prop.disporella.flat<-(food.exp.data.12.2019$num.disporella.flat)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.flat.12<- glm(formula = cbind(num.prop.disporella.flat, 1-num.prop.disporella.flat)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.flat.12)

summary(glm.binomial.num.prop.disporella.flat.12)
anova(glm.binomial.num.prop.disporella.flat.12, test = "Chi")

### Proportion ruffles. all
food.exp.data.12.2019$num.prop.disporella.ruffles<-(food.exp.data.12.2019$num.disporella.erect.ruffles)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.ruffles.12<- glm(formula = cbind(num.prop.disporella.ruffles, 1-num.prop.disporella.ruffles)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.ruffles.12)

summary(glm.binomial.num.prop.disporella.ruffles.12)
anova(glm.binomial.num.prop.disporella.ruffles.12, test = "Chi")


### Proportion fan. all
food.exp.data.12.2019$num.prop.disporella.fan<-(food.exp.data.12.2019$num.disporella.erect.fan)/(food.exp.data.12.2019$num.disporella)

glm.binomial.num.prop.disporella.fan.12<- glm(formula = cbind(num.prop.disporella.fan, 1-num.prop.disporella.fan)~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "binomial")
plot(glm.binomial.num.prop.disporella.fan.12)

summary(glm.binomial.num.prop.disporella.fan.12)
anova(glm.binomial.num.prop.disporella.fan.12, test = "Chi")

#### Num. erect fan

poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella.erect.fan, "Poisson")
qqp(food.exp.data.12.2019$num.disporella.erect.fan, "pois", poisson.12$estimate)

glm.poisson.num.disporella.erect.fan.12<-glm(formula = (num.disporella.erect.fan) ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.erect.fan.12)
plot(glm.poisson.num.disporella.erect.fan.12)
anova(glm.poisson.num.disporella.erect.fan.12, test = "Chi")

#### Num ruffled
poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella.erect.ruffles, "Poisson")
qqp(food.exp.data.12.2019$num.disporella.erect.ruffles, "pois", poisson.12$estimate)

glm.poisson.num.disporella.erect.ruffles.12<-glm(formula = (num.disporella.erect.ruffles) ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.erect.ruffles.12)
plot(glm.poisson.num.disporella.erect.ruffles.12)
anova(glm.poisson.num.disporella.erect.ruffles.12, test = "Chi")


###Num knobs
poisson.12<-fitdistr(food.exp.data.12.2019$num.disporella.erect.knobs, "Poisson")
qqp(food.exp.data.12.2019$num.disporella.erect.knobs, "pois", poisson.12$estimate)

glm.poisson.num.disporella.erect.knobs.12<-glm(formula = (num.disporella.erect.knobs) ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.disporella.erect.knobs.12)
plot(glm.poisson.num.disporella.erect.knobs.12)
anova(glm.poisson.num.disporella.erect.knobs.12, test = "Chi")

### Num didemnum
poisson.12<-fitdistr(food.exp.data.12.2019$num.didemnum, "Poisson")
qqp(food.exp.data.12.2019$num.didemnum, "pois", poisson.12$estimate)

glm.poisson.num.didemnum.12<-glm(formula = (num.didemnum) ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "poisson")
summary(glm.poisson.num.didemnum.12)
plot(glm.poisson.num.didemnum.12)
anova(glm.poisson.num.didemnum.12, test = "Chi")



####################### NORMAL
library(lme4)


qqp(log(food.exp.data.12.2019$richness+1), "norm")

qqp(food.exp.data.12.2019$richness, "lnorm")

qqp(sqrt(food.exp.data.12.2019$richness), "norm")

lm.hydroid<-lm(formula = log(hydroid+1) ~ av.pH*Food.quality, data = food.exp.data.12.2019)
summary(lm.hydroid)


summary(lm.richness)
plot(lm.richness)
anova(lm.richness, test = "F")


lognormal.12<-fitdistr(food.exp.data.12.2019$richness+0.01, "lognormal")
qqp(food.exp.data.12.2019$richness, "lnorm", lognormal.12$estimate)


glm.gamma.richness.12<- glm(formula = richness ~ av.pH*Food.quality, data = food.exp.data.12.2019, family = "lognormal")


qqp(food.exp.data.12.2019$richness, "lnorm")
?glm
?qqp

visreg(lm.richness, xvar = "av.pH", by= "Food.quality", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))




###might need gamlss pacage to do exponential functions
