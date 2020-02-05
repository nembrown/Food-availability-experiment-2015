#read in packages
library(vegan)
library(ggplot2)
library(betapart)
library(bipartite)
library(car)
library(fitdistrplus)

## read in data from both mesocosms and tiles
food.exp.data.mesocosm<-read.csv("C:Data//Mesocosm inventory data//mesocosm_inventory_counts.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )
head(food.exp.data.mesocosm)

mesocosm.key<-read.csv("C:Data//Percent cover data//Environmental files for community analysis//mesocosm_key.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )
food.exp.data.tile.all<-read.csv("C:Data//Percent cover data//Percent cover files over time//tile_week_12_all.csv" ,stringsAsFactors = FALSE, na.strings = c("NA","") )
head(food.exp.data.tile.all)

food.exp.data.mesocosm.12<-read.csv("c:Data//Mesocosm inventory data//food.exp.data.mesocosm.12.original.csv", stringsAsFactors = FALSE, na.strings = c("NA","") )
head(food.exp.data.mesocosm.12)


#Select data to be used: 
# some cover and some counts: 


#cover 
names_food_exp_tile<-c("Mesocosm", "formicula", "hydroid", "alive.bot", "alive.mem" , "caprellid" ,"bowerbankia" )

#counts
names_counts_mesocosm<-c("Mesocosm","num.nudi.eggs" ,
                          "num.nudi" ,
                          "mussel_complete" ,
                          "brown_sponge" ,
                          "didemnum" ,
                          "num.corella" ,
                          "schizo" ,
                          "disporella" ,
                          "num.spirorbid",
                          "num.serpulid" ,
                          "num.barn.alive" ,
                          "num.cribrilina" ,
                          "num.seastar" ,
                          "orange_sponge" ,
                          "TRB" ,
                          "clam" ,
                          "scallop" ,
                          "pres.limpet" ,
                          "pres.bubble.snail" ,
                          "num_anemone")
                          


#For sp div and richness use nudi.combined but for community level use nudi.separated
#don't include these b/c 0 at week 12
#which( colnames(food.exp.data.tile.all)=="corambe.nudis" )
#which( colnames(food.exp.data.tile.all)=="num.flatworm" )
#which( colnames(food.exp.data.tile.all)=="num.isopod" )
#which( colnames(food.exp.data.tile.all)=="white.worm.1" )
#changed these two:
#which( colnames(food.exp.data.mesocosm.12)=="pres.brown.sponge" )
#which( colnames(food.exp.data.mesocosm.12)=="seastar.eva" )

food.exp.data.tile.selected<-food.exp.data.tile.all[,colnames(food.exp.data.tile.all) %in% names_food_exp_tile]
head(food.exp.data.tile.selected)
food.exp.data.mesocosm.selected<-food.exp.data.mesocosm.12[,colnames(food.exp.data.mesocosm.12) %in% names_counts_mesocosm]
head(food.exp.data.mesocosm.selected)

#community level data combined counts and percent cover
species.rec_cover <- merge(food.exp.data.tile.selected, food.exp.data.mesocosm.selected)
head(species.rec_cover)

just.species.rec_cover<- species.rec_cover[,-1]
head(just.species.rec_cover)



##### making a newdataframe for % cover only - to be used for evenness and shannon diversity
head(food.exp.data.tile.all)

names_cover_food_exp_tile<-c("Mesocosm","formicula" ,
                            "hydroid" ,
                            "alive.bot" ,
                            "alive.mem" ,
                            "caprellid" ,
                            "nudi.eggs" ,
                            "nudi" ,
                            "corambe.nudis",
                            "mussel" ,
                            "sponge.brown" ,
                            "didemnum" ,
                            "bowerbankia" ,
                            "corella" ,
                            "schizo" ,
                            "disporella" ,
                            "serpulid" ,
                            "alive.barn" ,
                            "cribrilina" ,
                            "flatworm" ,
                            "anemone" ,
                            "isopod" ,
                            "white.worm" ,
                            "seastar" )



species.cover <- food.exp.data.tile.all[,colnames(food.exp.data.tile.all) %in% names_cover_food_exp_tile]
head(species.cover)
just.species.cover<-species.cover[,-1]

species.cover$richness<-specnumber(just.species.cover)
species.cover$shannon.diversity<-diversity(just.species.cover, index="shannon")
species.cover$evenness<-species.cover$shannon.diversity/(log(species.cover$richness))

#evenness has to be created from just cover data
food.exp.data.mesocosm.12$evenness<-species.cover$evenness

#richness can be from count data - it's just pres/abs of given species
food.exp.data.mesocosm.12$richness <- specnumber(just.species.rec_cover)


# MDS ---------------------------------------------------------------------

compiled.data <- mesocosm.key
head(compiled.data)

row.names(compiled.data)<-compiled.data$Mesocosm

head(compiled.data)
compiled.data$Combined.Treatment<-as.factor(compiled.data$Combined.Treatment) 

#Combinging species and environment
all.data.rec_cover<-merge(species.rec_cover,compiled.data)
head(all.data.rec_cover)

cbbPalette.all.2<- c( "#F8766D", "#F8766D", "#00BA38" , "#00BA38", "#619CFF", "#619CFF")

###CONSTRAINED Ordination

capscale_plot<- function(m, colorby){
  colorby<-factor(colorby) #Convert to factor (just in case it isn't already)
  cols <- colorset6#vector of colors needed
  shapesies<-c( 16,2,16,2,16,2)
  ordiplot(m, display = c("sites"), type = "n")
  #ordisurf(m ~ min.10.pH, data=compiled.data_zscores, method = "REML", select = TRUE)
  points(m, col = cols[colorby], pch = shapesies[colorby], cex=1.5)
  legend("topright", title ="Food  CO2", legend=levels(colorby), col=cols, pch = shapesies, cex=.5)
}


# need to have zscores for pH ... otherwise evaluating at ph=0 which is not meaningful 
compiled.data_zscores<-compiled.data
compiled.data_zscores$min.10.pH<-scale(compiled.data$min.10.pH, center=TRUE, scale=TRUE)
head(compiled.data_zscores)



# Bray-Curtis Capscale / constrained ordination -------------------------------------------------------------
colorset2 = c("High"="#F8A02E" ,"Low"="#439E5F","None"= "#666666")
colorset6 = c("HighAmbient"="#F8A02E", "HighElevated"="#F8A02E" ,"LowAmbient"="#439E5F", "LowElevated"="#439E5F","NoneAmbient"= "#666666", "NoneElevated"= "#666666")



#Standardizing by total of the species either percent or count
# This makes each species on their own scale, so the mesocosm got x % of the total mussels for eg.
standardized.species.rec_cover<-decostand(just.species.rec_cover, method="total", MARGIN=2)
head(standardized.species.rec_cover)

model.meso.bray<-capscale(standardized.species.rec_cover ~ min.10.pH*Food.quality,compiled.data_zscores , distance="bray")
capscale_plot(model.meso.bray, colorby=compiled.data$Combined.Treatment)

model.meso.bray.sf<-ordisurf(model.meso.bray ~ min.10.pH, data=compiled.data_zscores, method = "REML", select = TRUE)
summary(model.meso.bray.sf)

adonis(standardized.species.rec_cover ~ min.10.pH*Food.quality, method="bray", permutations = 9999, data=compiled.data_zscores)
summary(model.meso.bray)

model.meso.bray.scores<- as.data.frame(scores(model.meso.bray)$sites)
head(model.meso.bray.scores)
model.meso.bray.scores$Mesocosm<-species.rec_cover$Mesocosm
#Can't do it by row number because missing mesocosm #46

model.meso.bray.scores.CAP<-merge(model.meso.bray.scores, compiled.data, by="Mesocosm")
head(model.meso.bray.scores.CAP)

write.csv(model.meso.bray.scores,"C:Data//Mesocosm inventory data//model.meso.bray.scores.csv", row.names=FALSE)

food.exp.data.mesocosm.12<-merge(model.meso.bray.scores,food.exp.data.mesocosm.12)



# betadispersion partitioned ----------------------------------------------


dist.part.bray<-bray.part(standardized.species.rec_cover)
#returns a distance matrix, pairwise between site values of each component of beta diversitity 
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


bd.overall.bray.distances<- as.data.frame(bd.bray$distances)
head(bd.overall.bray.distances)
bd.overall.bray.distances$distcentroid<-bd.overall.bray.distances$`bd.bray$distances`
bd.overall.bray.distances$Mesocosm<-species.rec_cover$Mesocosm
bd.overall.bray.distances.2<-merge(bd.overall.bray.distances, compiled.data, by="Mesocosm")
head(bd.overall.bray.distances.2)

plot.overall.distcentroid.12.hydrogen<- ggplot(bd.overall.bray.distances.2, aes(x=min.10.pH, y=distcentroid, colour=Food.quality)) + geom_point(size=5,aes(colour=factor(Food.quality), shape=CO2)) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(aes(col=Food.quality, fill=Food.quality), alpha=0.15,size=1.5, method="lm") +scale_shape_manual(values=c(19,17))
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme_bw() +  xlab(expression("Minimum" ~"10"^"th"~"percentile pH")) + ylab("Distance to centroid")  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))#+ylim(0,0.75)
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm") )
plot.overall.distcentroid.12.hydrogen<- plot.overall.distcentroid.12.hydrogen+ theme(legend.position="none")+ scale_colour_discrete(name = "Food.quality")+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))+ theme(axis.text.x = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")), axis.text.y = element_text(margin=margin(0.5, 0.5, 0.5, 0.5, "cm")))
plot.overall.distcentroid.12.hydrogen 

write.csv(bd.overall.bray.distances.2,"C:Data//Mesocosm inventory data/bd.overall.bray.distances.2.csv",row.names=FALSE )

food.exp.data.mesocosm.12<-merge(bd.overall.bray.distances,food.exp.data.mesocosm.12)
head(food.exp.data.mesocosm.12)

write.csv(food.exp.data.mesocosm.12,"C:Data//Mesocosm inventory data/food.exp.data.mesocosm.12.csv", row.names=FALSE)


