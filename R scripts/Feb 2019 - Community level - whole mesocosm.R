setwd("C:/Users/Norah/Dropbox/Projects/Food-availability-experiment-2015/Data")

library(vegan)
library(ggplot2)
library(betapart)
library(bipartite)
library(car)
library(fitdistrplus)

## For sp div and richness use nudi.combined but for community level use nudi.separated
food.exp.data.mesocosm<-read.csv("c:Mesocosm inventory data//mesocosm_inventory_counts.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )
head(food.exp.data.mesocosm)

mesocosm.key<-read.csv("c:Percent cover data//Environmental files for community analysis//mesocosm_key.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )
food.exp.data.tile.all<-read.csv("C:Percent cover data//Percent cover files over time//tile_week_12_all.csv" ,stringsAsFactors = FALSE, na.strings = c("NA","") )
head(food.exp.data.tile.all)
food.exp.data.mesocosm<-read.csv("c:Mesocosm inventory data//mesocosm_inventory_counts.csv",stringsAsFactors = FALSE, na.strings = c("NA","") )

food.exp.data.mesocosm.12<-read.csv("c:Mesocosm inventory data//food.exp.data.mesocosm.12.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )
head(food.exp.data.mesocosm.12)
#combination of food.exp.data.mesocosm and meso.key and then combo tiles... 

#cover = formicula, hydroid, bortyllus, membranipora, caprellid, bowerbankia, brown sponge
which( colnames(food.exp.data.tile.all)=="formicula" )
which( colnames(food.exp.data.tile.all)=="hydroid" )
which( colnames(food.exp.data.tile.all)=="alive.bot" )
which( colnames(food.exp.data.tile.all)=="alive.mem" )
which( colnames(food.exp.data.tile.all)=="caprellid" )
which( colnames(food.exp.data.tile.all)=="bowerbankia" )


#counts right order
paste(
  which( colnames(food.exp.data.mesocosm.12)=="num.nudi.eggs" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.nudi" ),
  which( colnames(food.exp.data.mesocosm.12)=="mussel_complete" ),
  which( colnames(food.exp.data.mesocosm.12)=="brown_sponge" ),
  which( colnames(food.exp.data.mesocosm.12)=="didemnum" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.corella" ),
  which( colnames(food.exp.data.mesocosm.12)=="schizo" ),
  which( colnames(food.exp.data.mesocosm.12)=="disporella" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.spirorbid" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.serpulid" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.barn.alive" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.cribrilina" ),
  which( colnames(food.exp.data.mesocosm.12)=="num.seastar" ),
  which( colnames(food.exp.data.mesocosm.12)=="orange_sponge" ),
  which( colnames(food.exp.data.mesocosm.12)=="TRB" ),
  which( colnames(food.exp.data.mesocosm.12)=="clam" ),
  which( colnames(food.exp.data.mesocosm.12)=="scallop" ),
  which( colnames(food.exp.data.mesocosm.12)=="pres.limpet" ),
  which( colnames(food.exp.data.mesocosm.12)=="pres.bubble.snail" ),
  which( colnames(food.exp.data.mesocosm.12)=="num_anemone" ),
  sep=",")




#don't include these b/c 0 at week 12
#which( colnames(food.exp.data.tile.all)=="corambe.nudis" )
#which( colnames(food.exp.data.tile.all)=="num.flatworm" )
#which( colnames(food.exp.data.tile.all)=="num.isopod" )
#which( colnames(food.exp.data.tile.all)=="white.worm.1" )

#changed these two:
#which( colnames(food.exp.data.mesocosm.12)=="pres.brown.sponge" )
#which( colnames(food.exp.data.mesocosm.12)=="seastar.eva" )

#53,52,51,46,50,42,39,40,49,38,48,47,55,45,17,16,19,23,26

species.rec_cover <- cbind(food.exp.data.tile.all[,c(1,8,9,11,12,21,38)], food.exp.data.mesocosm.12[,c(53,52,51,46,50,42,39,40,49,43,38,48,55,45,17,16,19,23,26,47)])
head(species.rec_cover)
row.names(species.rec_cover)<-species.rec_cover$Mesocosm
species.rec_cover<-species.rec_cover[,-1]

just.species.rec_cover<- species.rec_cover
head(just.species.rec_cover)


# New dataframes for other analyses ---------------------------------------



##########making a new one for use in PRC analysis - need sto be right ordeR:

species.rec_cover_2 <- cbind(food.exp.data.tile.all[,c(1,8,9,11,12,21)], 
                             food.exp.data.mesocosm.12[, c(51,50)], 
                             food.exp.data.tile.all[,c(25)],
                             food.exp.data.mesocosm.12[,c(49,54,48)],
                             food.exp.data.tile.all[,c(38)],
                             food.exp.data.mesocosm.12[,c(41,38,39,47,42,37,46)],
                             food.exp.data.tile.all[,c(66)],
                             food.exp.data.mesocosm.12[,c(45)],
                             food.exp.data.tile.all[,c(68,69)],
                             food.exp.data.mesocosm.12[,c(53,44,16,15,18,22,25)])

head(species.rec_cover_2)  
names(species.rec_cover_2)[9]<-"corambe.nudis"
names(species.rec_cover_2)[7]<-"nudi.eggs"
names(species.rec_cover_2)[8]<-"nudi"
names(species.rec_cover_2)[10]<-"mussel"
names(species.rec_cover_2)[11]<-"sponge.brown"
names(species.rec_cover_2)[13]<-"bowerbankia"
names(species.rec_cover_2)[15]<-"num.schizo"
names(species.rec_cover_2)[16]<-"num.disporella"
names(species.rec_cover_2)[21]<-"num.flatworm"
names(species.rec_cover_2)[22]<-"num.anemone"
names(species.rec_cover_2)[24]<-"white.worm"
names(species.rec_cover_2)[26]<-"num.orange.sponge"



write.csv(species.rec_cover_2,"C:Mesocosm inventory data/species.rec_cover_mesocosm_12.csv")



##### making a newdataframe for % cover only - to be used for evenness and shannon diversity ... but has to have 0.5 for the right ones. 
head(food.exp.data.tile.all)
paste(
  which( colnames(food.exp.data.tile.all)=="formicula" ),
  which( colnames(food.exp.data.tile.all)=="hydroid" ),
  which( colnames(food.exp.data.tile.all)=="alive.bot" ),
  which( colnames(food.exp.data.tile.all)=="alive.mem" ),
  which( colnames(food.exp.data.tile.all)=="caprellid" ),
  which( colnames(food.exp.data.tile.all)=="nudi.eggs" ),
  which( colnames(food.exp.data.tile.all)=="nudi" ),
  which( colnames(food.exp.data.tile.all)=="corambe.nudis" ),
  which( colnames(food.exp.data.tile.all)=="mussel" ),
  which( colnames(food.exp.data.tile.all)=="sponge.brown" ),
  which( colnames(food.exp.data.tile.all)=="didemnum" ),
  which( colnames(food.exp.data.tile.all)=="bowerbankia" ),
  which( colnames(food.exp.data.tile.all)=="corella" ),
  which( colnames(food.exp.data.tile.all)=="schizo" ),
  which( colnames(food.exp.data.tile.all)=="disporella" ),
  which( colnames(food.exp.data.tile.all)=="serpulid" ),
  which( colnames(food.exp.data.tile.all)=="alive.barn" ),
  which( colnames(food.exp.data.tile.all)=="cribrilina" ),
  which( colnames(food.exp.data.tile.all)=="flatworm" ),
  which( colnames(food.exp.data.tile.all)=="anemone" ),
  which( colnames(food.exp.data.tile.all)=="isopod" ),
  which( colnames(food.exp.data.tile.all)=="white.worm" ),
  which( colnames(food.exp.data.tile.all)=="seastar" ),
  sep=","
)

species.cover <- food.exp.data.tile.all[,c(8,9,11,12,21,23,24,25,28,29,31,38,10,13,14,20,26,30,33,34,35,36,37)]
head(species.cover)
just.species.cover<-species.cover

species.cover$richness<-specnumber(just.species.cover)
species.cover$shannon.diversity<-diversity(just.species.cover, index="shannon")
species.cover$evenness<-species.cover$shannon.diversity/(log(species.cover$richness))
food.exp.data.mesocosm.12$evenness<-species.cover$evenness

food.exp.data.mesocosm.12$richness <- specnumber(just.species.rec_cover)

write.csv(food.exp.data.mesocosm.12,"C:Mesocosm inventory data/food.exp.data.mesocosm.12.csv")


# MDS ---------------------------------------------------------------------


##read in environment dataset##

compiled.data <- mesocosm.key
head(compiled.data)

row.names(compiled.data)<-compiled.data$Mesocosm

head(compiled.data)
compiled.data$Combined.Treatment<-as.factor(compiled.data$Combined.Treatment) 



all.data.rec_cover<-cbind(just.species.rec_cover,compiled.data)



################################


cbbPalette.all.2<- c( "#F8766D", "#F8766D", "#00BA38" , "#00BA38", "#619CFF", "#619CFF")
cbbPalette.all.3<- c( "#F8766D","#00BA38" , "#619CFF")

cbbPalette.all.2.control<- c( "#F8766D", "#FFFFFF", "#00BA38" , "#FFFFFF", "#619CFF", "#FFFFFF")

cbbPalette.all<- c( "#F8766D", "#00BA38", "#619CFF", "#F8766D", "#00BA38", "#619CFF")

###CONSTRAINED Ordination

capscale_plot<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- cbbPalette.all.2#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordisurf(m ~ min.10.pH, data=compiled.data_zscores, method = "REML", select = TRUE)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Food  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=1)
  comm=species.rec_cover_jacc

}

#can I addd arrows? 
#scores(m,display="sp")
#arrows(m, "species") # not this one exactly


##October 11, 2018
# need to rescazle variables.... 
# need to have zscores at least for pH ... otherwise evaluating at 0 but not meaningful ... need to do something to resp. variables... 
compiled.data_zscores<-compiled.data
compiled.data_zscores$min.10.pH<-scale(compiled.data$min.10.pH, center=TRUE, scale=TRUE)
head(compiled.data_zscores)

#shoulod do for spp too ... but not for multivariate, just for univariate stuff. 
#head(species.rec_cover)
#species.rec_cover_zscores<-species.rec_cover
#species.rec_cover_zscores[1:26]<-lapply(species.rec_cover_zscores, function(x) {
#  y<-scale(x, center=TRUE, scale=TRUE)
#}
#)
#head(species.rec_cover_zscores)


# Jaccard -----------------------------------------------------------------



# Jaccard - answers the question - of community composition (regardless of abundance)
species.rec_cover_jacc<-decostand(species.rec_cover, method="pa")

?capscale


model.meso.jac<-capscale(species.rec_cover_jacc ~ min.10.pH*Food.quality,compiled.data_zscores , distance="jaccard", binary = TRUE)
capscale_plot(model.meso.jac, colorby=compiled.data$Combined.Treatment)
adonis(species.rec_cover_jacc ~ min.10.pH*Food.quality, method="jaccard", permutations = 9999, data=compiled.data_zscores, binary=TRUE)

model.meso.jac.scores<- as.data.frame(scores(model.meso.jac)$sites)
head(model.meso.jac.scores)
model.meso.jac.scores$Mesocosm<-row.names(model.meso.jac.scores)
model.meso.jac.scores.CAP<-merge(model.meso.jac.scores, compiled.data, by="Mesocosm")
head(model.meso.jac.scores.CAP)

plot.total.CAP1.12.hydrogen<- ggplot(model.meso.jac.scores.CAP, aes(x=min.10.pH, y=CAP1, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.total.CAP1.12.hydrogen<- plot.total.CAP1.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("CAP1")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.total.CAP1.12.hydrogen<- plot.total.CAP1.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.total.CAP1.12.hydrogen<- plot.total.CAP1.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.total.CAP1.12.hydrogen <-plot.total.CAP1.12.hydrogen +  scale_x_reverse( breaks=c(10^-7.0,10^-7.2,10^-7.4,10^-7.6, 10^-7.8), labels=c(7.0, 7.2, 7.4, 7.6, 7.8))
plot.total.CAP1.12.hydrogen 
library(ggplot2)

plot.total.CAP2.12.hydrogen<- ggplot(model.meso.jac.scores.CAP, aes(x=min.10.pH, y=CAP2, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.total.CAP2.12.hydrogen<- plot.total.CAP2.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("CAP2")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.total.CAP2.12.hydrogen<- plot.total.CAP2.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.total.CAP2.12.hydrogen<- plot.total.CAP2.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.total.CAP2.12.hydrogen <-plot.total.CAP2.12.hydrogen +  scale_x_reverse( breaks=c(10^-7.0,10^-7.2,10^-7.4,10^-7.6, 10^-7.8), labels=c(7.0, 7.2, 7.4, 7.6, 7.8))
plot.total.CAP2.12.hydrogen 



# Betadiversity -----------------------------------------------------------

#using betadiver
data(sipoo)
m <- betadiver(species.rec_cover_jacc)
plot(m)


#betadispersion
dist_jac <- vegdist(species.rec_cover_jacc, method = "jaccard", binary=TRUE)

mod.meso.jac.Food.quality<-betadisper(dist_jac, compiled.data$Food.quality, type="centroid")
anova(mod.meso.jac.Food.quality)


mod.meso.jac.CO2<-betadisper(dist_jac, compiled.data$CO2, type="centroid")
anova(mod.meso.jac.CO2)



mod.meso.jac.Combined.Treatment<-betadisper(dist_jac, compiled.data$Combined.Treatment, type="centroid", bias.adjust=TRUE)
anova(mod.meso.jac.Combined.Treatment)

#plot_net?
#web methods?
install.packages("bipartite")
library(bipartite)
nestedrank(web, method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=FALSE)
#web is a matrix with elements of a set (e.g., plants) as rows, elements of a second set (e.g., pollinators) as columns and number of interactions as entries.

head(species.rec_cover)
nested.matrix<-nestedrank(species.rec_cover, method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=FALSE)
nested.matrixnon.normal<-nestedrank(species.rec_cover, method = "NODF", weighted=TRUE, normalise=FALSE, return.matrix=FALSE)

?nestedrank

nested.matrix$'lower level'


model.meso.jac.scores.CAP$nestedness<-nested.matrix$'lower level'
head(model.meso.jac.scores.CAP)
model.meso.jac.scores.CAP$nestedness.nonnormal<-nested.matrixnon.normal$'lower level'

plot.total.nestedness.12.hydrogen<- ggplot(model.meso.jac.scores.CAP, aes(x=min.10.pH, y=nestedness, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.total.nestedness.12.hydrogen<- plot.total.nestedness.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("nestedness")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.total.nestedness.12.hydrogen<- plot.total.nestedness.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.total.nestedness.12.hydrogen<- plot.total.nestedness.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.total.nestedness.12.hydrogen <-plot.total.nestedness.12.hydrogen +  scale_x_reverse( breaks=c(10^-7.0,10^-7.2,10^-7.4,10^-7.6, 10^-7.8), labels=c(7.0, 7.2, 7.4, 7.6, 7.8))
plot.total.nestedness.12.hydrogen 




multipart.jac<-beta.multi(species.rec_cover_jacc, index.family = "jaccard")
head(multipart.jac)

betasample.jac<-beta.sample(species.rec_cover_jacc, index.family="jaccard", sites=10, samples = 100)

## Permutation test for F
permutest(mod.meso.jac.Combined.Treatment, pairwise = TRUE, permutations = 999)

## Tukey's Honest Significant Differences
(mod.HSD <- TukeyHSD(mod.meso.jac.Combined.Treatment))
plot(mod.HSD)

#can't do by hydrogen b/c it's by groups.... 

plot(dist.part.jac$beta.jac~compiled.data$min.10.pH)



#betadispersion partitioned - do this by Food and CO2... 
library(betapart)
dist.part.jac<-beta.pair(species.rec_cover_jacc, index.family = "jaccard")

#returns a distance matrix, pairwise between site values of each component of beta diversitity 
bd.jac<-betadisper(dist.part.jac[[3]],compiled.data$Combined.Treatment)
bd.nestedness.jac<-betadisper(dist.part.jac[[2]],compiled.data$Combined.Treatment)
bd.turnover.jac<-betadisper(dist.part.jac[[1]],compiled.data$Combined.Treatment)

dist.part.jac.core<-betapart.core(species.rec_cover_jacc)
str(dist.part.jac.core)


plot(bd.jac)
anova(bd.jac)
boxplot(bd.jac)

plot(bd.nestedness.jac)
anova(bd.nestedness.jac)
boxplot(bd.nestedness.jac)

plot(bd.turnover.jac)
anova(bd.turnover.jac)
boxplot(bd.turnover.jac)

str(bd.turnover.jac)

bd.turnover.jac$distances


bd.turnover.jac.distances<- as.data.frame(bd.turnover.jac$distances)
head(bd.turnover.jac.distances)
bd.turnover.jac.distances$distcentroid<-bd.turnover.jac.distances$`bd.turnover.jac$distances`
bd.turnover.jac.distances$Mesocosm<-row.names(bd.turnover.jac.distances)
bd.turnover.jac.distances.2<-merge(bd.turnover.jac.distances, compiled.data, by="Mesocosm")
head(bd.turnover.jac.distances.2)

plot.jac.turnover.distcentroid.12.hydrogen<- ggplot(bd.turnover.jac.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.jac.turnover.distcentroid.12.hydrogen<- plot.jac.turnover.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("jaccard turnover")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.jac.turnover.distcentroid.12.hydrogen<- plot.jac.turnover.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.jac.turnover.distcentroid.12.hydrogen<- plot.jac.turnover.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.jac.turnover.distcentroid.12.hydrogen <-plot.jac.turnover.distcentroid.12.hydrogen +  scale_x_reverse( breaks=c(10^-7.0,10^-7.2,10^-7.4,10^-7.6, 10^-7.8), labels=c(7.0, 7.2, 7.4, 7.6, 7.8))
plot.jac.turnover.distcentroid.12.hydrogen 

bd.nestedness.jac.distances<- as.data.frame(bd.nestedness.jac$distances)
head(bd.nestedness.jac.distances)
bd.nestedness.jac.distances$distcentroid<-bd.nestedness.jac.distances$`bd.nestedness.jac$distances`
bd.nestedness.jac.distances$Mesocosm<-row.names(bd.nestedness.jac.distances)
bd.nestedness.jac.distances.2<-merge(bd.nestedness.jac.distances, compiled.data, by="Mesocosm")
head(bd.nestedness.jac.distances.2)

plot.jac.nestedness.distcentroid.12.hydrogen<- ggplot(bd.nestedness.jac.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.jac.nestedness.distcentroid.12.hydrogen<- plot.jac.nestedness.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Jaccard Nestedness")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.jac.nestedness.distcentroid.12.hydrogen<- plot.jac.nestedness.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.jac.nestedness.distcentroid.12.hydrogen<- plot.jac.nestedness.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.jac.nestedness.distcentroid.12.hydrogen <-plot.jac.nestedness.distcentroid.12.hydrogen +  scale_x_reverse( breaks=c(10^-7.0,10^-7.2,10^-7.4,10^-7.6, 10^-7.8), labels=c(7.0, 7.2, 7.4, 7.6, 7.8))
plot.jac.nestedness.distcentroid.12.hydrogen 


bd.overall.jac.distances<- as.data.frame(bd.jac$distances)
head(bd.overall.jac.distances)
bd.overall.jac.distances$distcentroid<-bd.overall.jac.distances$`bd.jac$distances`
bd.overall.jac.distances$Mesocosm<-row.names(bd.overall.jac.distances)
bd.overall.jac.distances.2<-merge(bd.overall.jac.distances, compiled.data, by="Mesocosm")
head(bd.overall.jac.distances.2)

plot.jac.overall.distcentroid.12.hydrogen<- ggplot(bd.overall.jac.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.jac.overall.distcentroid.12.hydrogen<- plot.jac.overall.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Jaccard distcentroid")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.jac.overall.distcentroid.12.hydrogen<- plot.jac.overall.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.jac.overall.distcentroid.12.hydrogen<- plot.jac.overall.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.jac.overall.distcentroid.12.hydrogen <-plot.jac.overall.distcentroid.12.hydrogen +  scale_x_reverse( breaks=c(10^-7.0,10^-7.2,10^-7.4,10^-7.6, 10^-7.8), labels=c(7.0, 7.2, 7.4, 7.6, 7.8))
plot.jac.overall.distcentroid.12.hydrogen 


### Overall bd result ... 



#the overall bd result is more due to turnover than nestedness....#Baselga paper - species replacement not loss 
#elevated CO2 has more turnover. 






# Bray-Curtis -------------------------------------------------------------



#Standardizing by total of the species ... makes each one on their own scale... 
head(standardized.species.rec_cover)

standardized.species.rec_cover<-decostand(species.rec_cover, method="total", MARGIN=2)

model.meso.bray<-capscale(standardized.species.rec_cover ~ min.10.pH*Food.quality,compiled.data_zscores , distance="bray")
capscale_plot(model.meso.bray, colorby=compiled.data$Combined.Treatment)

model.meso.bray.sf<-ordisurf(model.meso.bray ~ min.10.pH, data=compiled.data_zscores, method = "REML", select = TRUE)
summary(model.meso.bray.sf)

adonis(standardized.species.rec_cover ~ min.10.pH*Food.quality, method="bray", permutations = 9999, data=compiled.data_zscores)

summary(model.meso.bray)

model.meso.bray.scores<- as.data.frame(scores(model.meso.bray)$sites)
head(model.meso.bray.scores)
model.meso.bray.scores$Mesocosm<-row.names(model.meso.bray.scores)
model.meso.bray.scores.CAP<-merge(model.meso.bray.scores, compiled.data, by="Mesocosm")
head(model.meso.bray.scores.CAP)

plot.CAP1.12.hydrogen<- ggplot(model.meso.bray.scores.CAP, aes(x=min.10.pH, y=CAP1, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.CAP1.12.hydrogen<- plot.CAP1.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("CAP1 (36% of constrained variation)")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.CAP1.12.hydrogen<- plot.CAP1.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.CAP1.12.hydrogen<- plot.CAP1.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.CAP1.12.hydrogen 


head(standardized.species.rec_cover)



str(compiled.data)
str(compiled.data_zscores)

### betadispersion 
dist_bray <- vegdist(standardized.species.rec_cover, method = "bray")
mod.meso.bray.Combined.Treatment<-betadisper(dist_bray, compiled.data_zscores$Combined.Treatment, type="centroid")
permutest(mod.meso.bray.Combined.Treatment)

boxplot(mod.meso.bray.Combined.Treatment$distances ~ compiled.data_zscores$Combined.Treatment)

head(mod.meso.bray.Combined.Treatment)

mod.meso.bray.CO2<-betadisper(dist_bray, compiled.data_zscores$CO2, type="centroid")
permutest(mod.meso.bray.CO2)

boxplot(mod.meso.bray.CO2$distances ~ compiled.data_zscores$CO2)

mod.meso.bray.Food.quality<-betadisper(dist_bray, compiled.data_zscores$Food.quality, type="centroid")
permutest(mod.meso.bray.Food.quality)


#betadispersion partitioned

dist.part.bray<-bray.part(standardized.species.rec_cover)
#returns a distance matrix, pairwise between site values of each component of beta diversitity 
#this is the right order .. but tunrover needs to be called "balanced vairation in species abundances" and nestedness "abundance gradient" 
bd.bray<-betadisper(dist.part.bray[[3]],compiled.data_zscores$Combined.Treatment, type="centroid" )
bd.nestedness.bray<-betadisper(dist.part.bray[[2]],compiled.data_zscores$Combined.Treatment, type="centroid")
bd.turnover.bray<-betadisper(dist.part.bray[[1]],compiled.data_zscores$Combined.Treatment, type="centroid")

head(bd.bray)
plot(bd.bray, hull=FALSE, ellipse = TRUE)
anova(bd.bray)
boxplot(bd.bray)
permutest(bd.bray)

plot(bd.nestedness.bray)
anova(bd.nestedness.bray)
boxplot(bd.nestedness.bray)

plot(bd.turnover.bray)
anova(bd.turnover.bray)
boxplot(bd.turnover.bray)

str(bd.turnover.bray)

bd.turnover.bray$distances


bd.turnover.bray.distances<- as.data.frame(bd.turnover.bray$distances)
head(bd.turnover.bray.distances)
bd.turnover.bray.distances$distcentroid<-bd.turnover.bray.distances$`bd.turnover.bray$distances`
bd.turnover.bray.distances$Mesocosm<-row.names(bd.turnover.bray.distances)
bd.turnover.bray.distances.2<-merge(bd.turnover.bray.distances, compiled.data, by="Mesocosm")
head(bd.turnover.bray.distances.2)

plot.turnover.distcentroid.12.hydrogen<- ggplot(bd.turnover.bray.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.turnover.distcentroid.12.hydrogen<- plot.turnover.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Turnover")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.turnover.distcentroid.12.hydrogen<- plot.turnover.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.turnover.distcentroid.12.hydrogen<- plot.turnover.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.turnover.distcentroid.12.hydrogen

bd.nestedness.bray.distances<- as.data.frame(bd.nestedness.bray$distances)
head(bd.nestedness.bray.distances)
bd.nestedness.bray.distances$distcentroid<-bd.nestedness.bray.distances$`bd.nestedness.bray$distances`
bd.nestedness.bray.distances$Mesocosm<-row.names(bd.nestedness.bray.distances)
bd.nestedness.bray.distances.2<-merge(bd.nestedness.bray.distances, compiled.data, by="Mesocosm")
head(bd.nestedness.bray.distances.2)

plot.nestedness.distcentroid.12.hydrogen<- ggplot(bd.nestedness.bray.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.nestedness.distcentroid.12.hydrogen<- plot.nestedness.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Nestedness")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.nestedness.distcentroid.12.hydrogen<- plot.nestedness.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.nestedness.distcentroid.12.hydrogen<- plot.nestedness.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.nestedness.distcentroid.12.hydrogen


bd.overall.bray.distances<- as.data.frame(bd.bray$distances)
head(bd.overall.bray.distances)
bd.overall.bray.distances$distcentroid<-bd.overall.bray.distances$`bd.bray$distances`
bd.overall.bray.distances$Mesocosm<-row.names(bd.overall.bray.distances)
bd.overall.bray.distances.2<-merge(bd.overall.bray.distances, compiled.data, by="Mesocosm")
head(bd.overall.bray.distances.2)

plot.overall.distcentroid.12.hydrogen<- ggplot(bd.overall.bray.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Distance to centroid")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.overall.distcentroid.12.hydrogen 


# mvabund -----------------------------------------------------------------
library(mvabund)
par(mar=c(2,10,2,2)) # adjusts the margins
head(standardized.species.rec_cover)
Norah_spp <- mvabund(standardized.species.rec_cover)
plot(Norah_spp~model.meso.bray.scores.CAP$Combined.Treatment, cex.axis=0.8, cex=0.8)

meanvar.plot(Norah_spp)

mod1 <- manylm(Norah_spp ~ compiled.data_zscores$min.10.pH*compiled.data_zscores$Food.quality)

plot(mod1)
anova(mod1)
#problem wit hthe residuals vs fitter - but alot the plot doesn't quite fit.... 





# Correlation -------------------------------------------------------------


library(ggcorrplot)
library(corrr)
library(PerformanceAnalytics)
library(cooccur)




corr_species<-round(cor(species.rec_cover, use="pairwise.complete.obs"), 1)
head(corr_species[, 1:6])
p.mat_species <- cor_pmat(corr_species)
head(p.mat_species[, 1:4])

ggcorrplot(corr_species)
ggcorrplot(corr_species, hc.order = TRUE, type = "lower",
           outline.col = "white")
ggcorrplot(corr_species, hc.order = TRUE,
           type = "lower", p.mat = p.mat_species)

ggcorrplot(corr_species, hc.order = TRUE, type = "lower",
           lab = TRUE)


#cooccur package
head(species.rec_cover_jacc)
species.rec_cover_jacc_transpose<-t(species.rec_cover_jacc)

cooccur.species <- cooccur(mat = species.rec_cover_jacc_transpose, type = "spp_site", thresh = TRUE, spp_names = TRUE)
summary(cooccur.species)
prob.table(cooccur.species)
plot(cooccur.species)

pair(mod=cooccur.species, "3")

pair.attributes(mod=cooccur.species)
pair.profile(mod=cooccur.species)

obs.v.exp(cooccur.species)


#netassoc package
library(netassoc)
head(species.rec_cover)
str(species.rec_cover)
species.rec_cover_int <- species.rec_cover %>% transform(hydroid = as.integer(hydroid), 
                                alive.bot= as.integer(alive.bot),
                                alive.mem= as.integer(alive.mem),
                                caprellid= as.integer(caprellid))
str(species.rec_cover_int)

obs<-t(species.rec_cover_int)
str(obs)

n<-make_netassoc_network(obs, nul=vegan::permatfull(obs)$perm[[1]],
                      method="partial_correlation", args=list(method="shrinkage",verbose=FALSE),
                      p.method="fdr", alpha=0.05, numnulls=1000,
                      plot=TRUE,plot.legend=TRUE, plot.title=TRUE, verbose=TRUE)



plot_netassoc_network (n$network_all)

plot(n$matrix_spsp_obs)
n$network_all

spxsp<-pairwise_association(obs, method = "condentropy")
image(spxsp)

#############################



#Simper
#The simper functions performs pairwise comparisons of groups of sampling units 
#and finds the average contributions of each species to the average overall Bray-Curtis dissimilarity between sites

#The results of simper can be very difficult to interpret. 
#The method very badly confounds the mean between group differences and within group variation, 
#and seems to single out variable species instead of distinctive species (Warton et al. 2012). 
#Even if you make groups that are copies of each other, the method will single out species with high contribution, 
#but these are not contributions to non-existing between-group differences 
#but to within-group variation in species abundance.


sim <- with(compiled.data, simper(species.rec_cover,CO2))
summary(sim)
simper(species.rec_cover,compiled.data$CO2, permutations = 0, trace = FALSE)

?mvabund
#mvabund - not sure how useful
install.packages("mvabund")
library(mvabund)
species.rec_cover.mvabund<-mvabund(species.rec_cover)
plot(species.rec_cover.mvabund ~ compiled.data$Combined.Treatment)


##### Indval
###"The indval approach looks for species.rec_cover that are
#both necessary and sufficient, i.e. if you find that species.rec_cover you should be in that type,
#and if you are in that type you should find that species.rec_cover
install.packages("labdsv")
library(labdsv)
iva<-indval(species.rec_cover,compiled.data$Combined.Treatment)

#### indv val is indicative of fidelity of that species.rec_cover to a particular site
#### 1 - mv abunda - 2- indval? 
 

iva$relfrq
iva$relabu
iva$indval

gr<- iva$maxcls[iva$pval<=0.05]
iv<- iva$indcls[iva$pval<=0.05]
pv<- iva$pval[iva$pval<=0.05]
fr<-apply(species.rec_cover>0, 2, sum)[iva$pval<=0.05]
fidg<-data.frame(group=gr, indval=iv, pvalue=pv, freq=fr)
fidg<-fidg[order(fidg$group,-fidg$indval), ]
fidg




?isamic
#Indicator Species Analysis Minimizing Intermediate Occurrences
#Calculates the degree to which species are either always present or always absent within clusters or types.
###do to break down into each part of combined. treatment ... Or just co2 vs just air
isamic(species.rec_cover,compiled.data$Combined.Treatment, sort=TRUE)

 
#october 7, 2018
# Gower - is good for dimensionally heterogeneous variables
model.meso.gower<-capscale(species.rec_cover ~ min.10.pH*Food.quality,compiled.data , distance="gower")
capscale_plot(model.meso.gower, colorby=compiled.data$Combined.Treatment)
adonis(species.rec_cover ~ min.10.pH*Food.quality, method="gower", permutations = 9999, data=compiled.data)
