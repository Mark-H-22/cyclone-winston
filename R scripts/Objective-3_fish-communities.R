# Objective 3 - Trajectories of reef fish communities

library(tidyverse); theme_set(theme_classic())
library(sciplot)   #for calculating standard errors
library(ggpubr)
library(lattice)
library(performance)
library(lme4)
library(sjPlot)

MyTheme <- theme(axis.title=element_text(size=13), axis.text=element_text(size=10),
                 legend.title=element_text(size=11), legend.text=element_text(size=10),
                 plot.title = element_text(size=13),
                 panel.background=element_rect(fill="transparent"), legend.box="vertical")



### Fish data cleaning & biomass plots

fish <- read.csv("Data/fish_full_new.csv")

# see observers in each survey year:
unique(fish[which(fish$Year==2013),]$Observers)
unique(fish[which(fish$Year==2014),]$Observers)
unique(fish[which(fish$Year==2016),]$Observers)
unique(fish[which(fish$Year==2018),]$Observers)
unique(fish[which(fish$Year==2020),]$Observers)


table(fish$Year)
#2013  2014  2016  2018  2020 
#2194  3182  1723  2238  3863


table(fish$Site)
# C3    C5  KB04  KB05  KB06  KB07   N19   N20  VIR1 VIR10 VIR11   
#614   666   902   742   807   693   706   841   705   379   585   

#VIR2  VIR3  VIR4  VIR5 VIR6  VIR7  VIR8  VIR9 
# 801   726   440   483 698   779  1035   598 


table(fish$Management)
#Namena_LMMA   Namena_MPA Vatuira_LMMA  Vatuira_MPA 
#       2780         3191         1856         5373

# Namena_MPA  = Namena Marine Reserve       = NMR     (this is a tabu, example of MCA)
# Vatuira_MPA = Vatu-i-Ra Conservation Park = ViRCP   (this is a tabu, example of MCA)
#                                             ViRCP_old = sites VIR 1 - 5, VIR11
#                                             ViRCP_new = sites VIR 8 & 9
# Namena_LMMA = Kubulau Fishing Ground      = KFG     (NMR is in this fishing ground, but unfished)
# Vatuira_LMMA= Nakorotubu Fishing Ground   = NFG     (ViRCP is in this fishing ground, but unfished)


# Split ViRCP into "old" and "new"
VIR.old <- c("VIR1","VIR2","VIR3","VIR4","VIR5","VIR11")
VIR.new <- c("VIR8","VIR9")
fish$Management <- ifelse(fish$Site %in% VIR.old, "ViRCP_old", fish$Management)
fish$Management <- ifelse(fish$Site %in% VIR.new, "ViRCP_new", fish$Management)

fish$Management <- recode_factor(fish$Management, 
                                 Namena_MPA="NMR",
                                 Namena_LMMA="KFG",
                                 Vatuira_LMMA="NFG")

table(fish$Management)
# NMR    KFG    NFG   ViRCP_old    ViRCP_new
#3191   2780   1856        3740         1633


fish$Site             <- as.factor(fish$Site)
fish$Reef.slope       <- as.factor(fish$Reef.slope)
fish$Management.rules <- as.factor(fish$Management.rules) 
fish$Management       <- as.factor(fish$Management)


# Merge some trophic groups together:
unique(fish$Functional.group)

# Excavators & scrapers
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "excavator/scraper" = c("excavator", "scraper"))

# Micro- & Macro-invertivores
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "invertivore" = c("micro_invertivore", "macro_invertivore",
                                                        "spongivore"))

# Browsers & croppers/grazers
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "browser/cropper/grazer" = c("browser", "cropper/grazer"))

# piscivores & pisci-invertivores
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "piscivore" = c("piscivore", "pisci-invertivore"))


# Only fish in following families in standard WCS protocol: 
# Acanthuridae, Balistidae, Caesionidae, Carangidae,
# Chaetodontidae, Ephippidae, Haemulidae, Kyphosidae, Labridae, Lethrinidae,
# Lutjanidae, Mullidae, Nemipteridae, Pomacanthidae, Priacanthidae, Scaridae,
# Serranidae (minus Anthias/Pseudanthias), Siganidae, Scombridae, Sphyraenidae,
# all sharks.

# find additional families
table(fish[which(fish$Year==2013),]$Fish.family) # Tetraodontidae
table(fish[which(fish$Year==2014),]$Fish.family) # Aulostomidae, Myliobatidae
table(fish[which(fish$Year==2016),]$Fish.family) # Diodontidae, Echeneidae, Pinguipedidae
table(fish[which(fish$Year==2018),]$Fish.family) # Aulostomidae, Clupeidae, Ginglymostomatidae, Holocentridae, 
#Monocanthidae, Myliobatidae, Pinguipedidae, Pomacentridae, Tetraodontidae
table(fish[which(fish$Year==2020),]$Fish.family) # Apogonidae, Aulostomidae, Holocentridae, 
#Monocanthidae, Muraenidae, Scorpaenidae, Sparidae, Tetraodontidae

# Remove families not in WCS survey list. 
# List to keep:
wcs.fams <- c("Acanthuridae", "Balistidae", "Caesionidae", "Carangidae",
              "Chaetodontidae", "Ephippidae", "Haemulidae", "Kyphosidae", "Labridae", 
              "Lethrinidae", "Lutjanidae", "Mullidae", "Nemipteridae", "Pomacanthidae", 
              "Priacanthidae", "Scaridae", "Serranidae", "Siganidae", "Scombridae", 
              "Sphyraenidae", "Zanclidae")
# don't include "Carcharhinidae" - UVC not good and also observer differences

fish1 <- fish[which(fish$Fish.family %in% wcs.fams),]
head(fish1)


# Address the difference in length estimates between observers in 2014:

fish14 <- data.frame(aggregate(fish1[,c("Count", "Biomass_kgha")], 
                               by=list(fish1$Year, fish1$Fish.family, 
                                       fish1$Observers, fish1$Site), 
                               sum, na.rm=T))
new.names2 <- c("Year","Fish.family","Observers","Site","Count","Biomass_kgha")
names(fish14)[1:6] <- new.names2
head(fish14)

# average over sites
fish14ii <- data.frame(aggregate(fish14[,c("Count","Biomass_kgha")], 
                                 by=list(fish14$Year, fish14$Fish.family, 
                                         fish14$Observers), 
                                 mean, na.rm=T))
new.names3 <- c("Year","Fish.family","Observers","Mean.count","Mean.biomass")
names(fish14ii)[1:5] <- new.names3
head(fish14ii)

# plot by fish count
fams14 <- ggplot(data=fish14ii[which(fish14ii$Year==2014),], 
                 aes(x=Fish.family, y=Mean.count, fill=as.factor(Observers))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_y_discrete(expand = c(0,0)) + ylim(0,80) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position=c(0.9,0.95)) + labs(fill=NULL) + MyTheme +
  xlab(NULL) + ylab("Mean count") + ggtitle("Reef fish counts, 2014")
fams14

# plot by fish biomass
fams14.bio <- ggplot(data=fish14ii[which(fish14ii$Year==2014),], 
                     aes(x=Fish.family, y=Mean.biomass, fill=as.factor(Observers))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_y_discrete(expand = c(0,0)) + ylim(0,5000) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position=c(0.9,0.95)) + labs(fill=NULL) + MyTheme +
  xlab(NULL) + ylab("Mean biomass (kg/ha)") + ggtitle("Reef fish biomass estimates, 2014")
fams14.bio


# Remove Observer 2 from 2014:
fish1 <- fish1[-which(fish1$Observers=="Observer 2"),]



### Biomass Plots


## whole fish assemblage

#sum abundance and biomass per transect:
f2 <- data.frame(aggregate(fish1[,c("Count","Biomass_kgha")], 
                           by=list(fish1$Year, fish1$Management, fish1$Site,
                                   fish1$Transect), 
                           sum, na.rm=T))
new.names <- c("Year","Management","Site","Transect")
names(f2)[1:4] <- new.names
head(f2)

# mean per site:
f3 <- data.frame(aggregate(f2[,c("Count","Biomass_kgha")], 
                           by=list(f2$Management, 
                                   f2$Site, f2$Year), 
                           mean, na.rm=T))
new.names2 <- c("Management","Site","Year")
names(f3)[1:3] <- new.names2
head(f3)

# means per mgmt area:
f4 <- data.frame(aggregate(f3[,c("Count","Biomass_kgha")], 
                           by=list(f3$Management, 
                                   f3$Year), 
                           mean, na.rm=T))
new.names3 <- c("Management","Year")
names(f4)[1:2] <- new.names3
head(f4)

#Add SEM across sites:
# (variables must be in same order as aggregating code above)
f4$bioSE <- aggregate(Biomass_kgha ~ Management + Year, 
                      f3, se)$Biomass_kgha

f4$abuSE <- aggregate(Count ~ Management + Year, 
                      f3, se)$Count


# relevel so Mgmt order in legend makes sense with colours
f4$Management <- fct_relevel(f4$Management, c("NMR","KFG","ViRCP_old",
                                              "ViRCP_new","NFG"))
# colours:
# NMR =           #AA2020
# KFG =           #F26A60
# ViRCP_old =     #095A62
# ViRCP_new =     #28A7B5
# NFG =           #9AE3EC


# Total biomass trends:
tot.bio <- ggplot(f4) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=1.25, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  ggtitle("Total biomass") +
  MyTheme + theme(legend.position="none",
                  plot.title=element_text(face="bold")) +
  ylab(expression(Biomass~(kg~ha^-1))) + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) 

tot.bio2 <- tot.bio + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                                    width=0.3, size=0.5, position=position_dodge(0.1)) 
tot.bio2



## Biomass per feeding group

#sum abundance and biomass per Functional group in each transect:
fish2 <- data.frame(aggregate(fish1[,c("Count","Biomass_kgha")], 
                              by=list(fish1$Year, fish1$Management, fish1$Site,
                                      fish1$Functional.group, fish1$Transect), 
                              sum, na.rm=T))
new.names <- c("Year","Management","Site","Func.group","Transect")
names(fish2)[1:5] <- new.names
head(fish2)


# mean per site:
fish3 <- data.frame(aggregate(fish2[,c("Count","Biomass_kgha")], 
                              by=list(fish2$Func.group, fish2$Management, 
                                      fish2$Site, fish2$Year), 
                              mean, na.rm=T))
new.names2 <- c("Func.group","Management","Site","Year")
names(fish3)[1:4] <- new.names2
head(fish3)


# means per mgmt area:
fish4 <- data.frame(aggregate(fish3[,c("Count","Biomass_kgha")], 
                              by=list(fish3$Func.group, fish3$Management, 
                                      fish3$Year), 
                              mean, na.rm=T))
new.names3 <- c("Func.group","Management","Year")
names(fish4)[1:3] <- new.names3
head(fish4)


#Add SEM across sites:
# (variables must be in same order as aggregating code above)
fish4$bioSE <- aggregate(Biomass_kgha ~ Func.group + Management + Year, 
                         fish3, se)$Biomass_kgha

fish4$abuSE <- aggregate(Count ~ Func.group + Management + Year, 
                         fish3, se)$Count


# relevel so Mgmt order in legend makes sense with colours
fish4$Management <- fct_relevel(fish4$Management, c("NMR","KFG","ViRCP_old",
                                                    "ViRCP_new","NFG"))


# Browser/cropper/grazer:
bro <- ggplot(fish4[which(fish4$Func.group=="browser/cropper/grazer"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Browser/cropper/grazer")

bro2 <- bro + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
bro2


# Corallivore:
cor <- ggplot(fish4[which(fish4$Func.group=="corallivore"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) +
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Corallivore")

cor2 <- cor + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
cor2


# Detritivore:
det <- ggplot(fish4[which(fish4$Func.group=="detritivore"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) +
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Detritivore")

det2 <- det + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
det2


# Excavator/scraper:
exc <- ggplot(fish4[which(fish4$Func.group=="excavator/scraper"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(expression(Biomass~(kg~ha^-1))) + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Excavator/scraper")

exc2 <- exc + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
exc2


# Invertivore:
inv <- ggplot(fish4[which(fish4$Func.group=="invertivore"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) +
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Invertivore")

inv2 <- inv + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
inv2


# Piscivores:
pis <- ggplot(fish4[which(fish4$Func.group=="piscivore"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) +
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("  ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Piscivore")

pis2 <- pis + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
pis2


# Planktivore:
pla <- ggplot(fish4[which(fish4$Func.group=="planktivore"),]) + 
  aes(x=Year, y=Biomass_kgha, colour=Management) + 
  geom_line(linewidth=0.75, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) +
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("  ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Planktivore")

pla2 <- pla + geom_errorbar(aes(ymin=Biomass_kgha-bioSE, ymax=Biomass_kgha+bioSE), 
                            width=0.3, size=0.5, position=position_dodge(0.1)) 
pla2


# Arrange biomass panels
fish.bio <- ggarrange(tot.bio2, bro2, cor2, det2, exc2, inv2, pis2, pla2, 
                      nrow=2, ncol=4, align="v", labels=c("a","b","c","d","e","f","g","h"),
                      common.legend=T, legend="top")
fish.bio

#pdf("plots/FG-biomass2.pdf", width=13, height=6)
#fish.bio
#dev.off()



### Plots of mean length and abundance trends (for supplementary material).

## Mean abundance of fish per FG
head(fish1)

# Sum abundance per FG per transect
FG_abu <- data.frame(aggregate(fish1[,c("Count")], 
                               by=list(fish1$Functional.group, fish1$Transect, 
                                       fish1$Site,
                                       fish1$Year, fish1$Management), 
                               sum, na.rm=T))
new.names <- c("Func.group","Transect","Site","Year","Management","Count")
names(FG_abu)[1:6] <- new.names
head(FG_abu)

# convert counts to abundance per hectare
# fish transects were 250 m^2
FG_abu$count_ha <- (FG_abu$Count/250)*10000


# Calculate mean abundance per site:
FG_abu2 <- data.frame(aggregate(FG_abu[,c("Count", "count_ha")], 
                                by=list(FG_abu$Func.group, FG_abu$Site,
                                        FG_abu$Year, FG_abu$Management), 
                                mean, na.rm=T))
new.names2 <- c("Func.group","Site","Year","Management","MeanCount","count_ha")
names(FG_abu2)[1:6] <- new.names2
head(FG_abu2)


# Calculate mean length per management level (across sites):
FG_abu3 <- data.frame(aggregate(FG_abu2[,c("MeanCount","count_ha")], 
                                by=list(FG_abu2$Func.group, 
                                        FG_abu2$Year, FG_abu2$Management), 
                                mean, na.rm=T))
new.names3 <- c("Func.group","Year","Management","MeanCount","count_ha")
names(FG_abu3)[1:5] <- new.names3
head(FG_abu3)

# Add in SEM
FG_abu3$Abu_SE <- aggregate(count_ha ~ Func.group + Year + Management, 
                            FG_abu2, se)$count_ha
head(FG_abu3)


# Plot mean length over time per management level, for each FG
# relevel so Mgmt order in legend makes sense with colours
FG_abu3$Management <- fct_relevel(FG_abu3$Management, c("NMR","KFG","ViRCP_old",
                                                        "ViRCP_new","NFG"))

# Browser/cropper/grazer:
bcg <- FG_abu3[which(FG_abu3$Func.group=="browser/cropper/grazer"),]
head(bcg)

bcg_abu <- ggplot(bcg) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Browser/cropper/grazer")

bcg_abu2 <- bcg_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                        ymax=count_ha + Abu_SE), 
                                    width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
bcg_abu2


# Corallivore:
crl <- FG_abu3[which(FG_abu3$Func.group=="corallivore"),]
head(crl)


crl_abu <- ggplot(crl) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Corallivore")

crl_abu2 <- crl_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                        ymax=count_ha + Abu_SE), 
                                    width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
crl_abu2


# Detritivore:
dtr <- FG_abu3[which(FG_abu3$Func.group=="detritivore"),]
head(dtr)


dtr_abu <- ggplot(dtr) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(expression(Abundance~(ha^-1))) + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Detritivore")

dtr_abu2 <- dtr_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                        ymax=count_ha + Abu_SE), 
                                    width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
dtr_abu2


# Excavator/scraper:
exsc <- FG_abu3[which(FG_abu3$Func.group=="excavator/scraper"),]
head(exsc)


exsc_abu <- ggplot(exsc) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Excavator/scraper")

exsc_abu2 <- exsc_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                          ymax=count_ha + Abu_SE), 
                                      width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
exsc_abu2


# Mobile invertivore
mobi <- FG_abu3[which(FG_abu3$Func.group=="invertivore"),]
head(mobi)


mobi_abu <- ggplot(mobi) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Invertivore")

mobi_abu2 <- mobi_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                          ymax=count_ha + Abu_SE), 
                                      width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
mobi_abu2


# Piscivore
pis <- FG_abu3[which(FG_abu3$Func.group=="piscivore"),]
head(pis)


pis_abu <- ggplot(pis) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Piscivore")

pis_abu2 <- pis_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                        ymax=count_ha + Abu_SE), 
                                    width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
pis_abu2


# Planktivore
plk <- FG_abu3[which(FG_abu3$Func.group=="planktivore"),]
head(plk)

plk_abu <- ggplot(plk) + 
  aes(x=Year, y=count_ha, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab(" ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle("Planktivore")

plk_abu2 <- plk_abu + geom_errorbar(aes(ymin=count_ha - Abu_SE, 
                                        ymax=count_ha + Abu_SE), 
                                    width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
plk_abu2



## Mean size of fish per FG

# MEAN LENGTH = sum(fish length * biomass of that length) / Total FG biomass

head(fish1)

# First, multiply length and biomass of length
fish1$lngth_bio <- fish1$Size * fish1$Biomass_kgha

# Second, sum lngth_bio per FG per transect
FG_bio <- data.frame(aggregate(fish1[,c("Biomass_kgha","lngth_bio","Count")], 
                               by=list(fish1$Functional.group, fish1$Transect, 
                                       fish1$Site,
                                       fish1$Year, fish1$Management), 
                               sum, na.rm=T))
new.names <- c("Func.group","Transect","Site","Year","Management","FG_Biomass",
               "FG_lngthbio","Count")
names(FG_bio)[1:8] <- new.names
head(FG_bio)


# Third, divide by total biomass of each FG per transect to get weighted mean
# length
FG_bio$wtMeanLength <- FG_bio$FG_lngthbio / FG_bio$FG_Biomass

head(FG_bio)


# Calculate (regular) mean length per site, from weighted means:
FG_bio2 <- data.frame(aggregate(FG_bio[,c("wtMeanLength")], 
                                by=list(FG_bio$Func.group, FG_bio$Site,
                                        FG_bio$Year, FG_bio$Management), 
                                mean, na.rm=T))
new.names2 <- c("Func.group","Site","Year","Management","MeanLength")
names(FG_bio2)[1:5] <- new.names2
head(FG_bio2)

# Calculate mean length per management level (across sites):
FG_bio3 <- data.frame(aggregate(FG_bio2[,c("MeanLength")], 
                                by=list(FG_bio2$Func.group, 
                                        FG_bio2$Year, FG_bio2$Management), 
                                mean, na.rm=T))
new.names3 <- c("Func.group","Year","Management","MeanLength")
names(FG_bio3)[1:4] <- new.names3
head(FG_bio3)

# Add in SEM
FG_bio3$ML_SE <- aggregate(MeanLength ~ Func.group + Year + Management, 
                           FG_bio2, se)$MeanLength
head(FG_bio3)


# Plot mean length over time per management level, for each FG
# relevel so Mgmt order in legend makes sense with colours
FG_bio3$Management <- fct_relevel(FG_bio3$Management, c("NMR","KFG","ViRCP_old",
                                                        "ViRCP_new","NFG"))


# Browser/cropper/grazer:
bcg <- FG_bio3[which(FG_bio3$Func.group=="browser/cropper/grazer"),]
head(bcg)

bcg_lngth <- ggplot(bcg) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\n ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

bcg_lngth2 <- bcg_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                            ymax=MeanLength + ML_SE), 
                                        width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
bcg_lngth2


# Corallivore:
crl <- FG_bio3[which(FG_bio3$Func.group=="corallivore"),]
head(crl)

crl_lngth <- ggplot(crl) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\n ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

crl_lngth2 <- crl_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                            ymax=MeanLength + ML_SE), 
                                        width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
crl_lngth2


# Detritivore:
dtr <- FG_bio3[which(FG_bio3$Func.group=="detritivore"),]
head(dtr)

dtr_lngth <- ggplot(dtr) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\nMean total length (cm)") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

dtr_lngth2 <- dtr_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                            ymax=MeanLength + ML_SE), 
                                        width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
dtr_lngth2


# Excavator/scraper:
exsc <- FG_bio3[which(FG_bio3$Func.group=="excavator/scraper"),]
head(exsc)

exsc_lngth <- ggplot(exsc) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\n ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

exsc_lngth2 <- exsc_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                              ymax=MeanLength + ML_SE), 
                                          width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
exsc_lngth2


# Invertivore
mobi <- FG_bio3[which(FG_bio3$Func.group=="invertivore"),]
head(mobi)

mobi_lngth <- ggplot(mobi) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\n ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

mobi_lngth2 <- mobi_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                              ymax=MeanLength + ML_SE), 
                                          width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
mobi_lngth2


# Piscivore
pis <- FG_bio3[which(FG_bio3$Func.group=="piscivore"),]
head(pis)

pis_lngth <- ggplot(pis) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\n ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

pis_lngth2 <- pis_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                            ymax=MeanLength + ML_SE), 
                                        width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
pis_lngth2


# Planktivore
plk <- FG_bio3[which(FG_bio3$Func.group=="planktivore"),]
head(plk)

plk_lngth <- ggplot(plk) + 
  aes(x=Year, y=MeanLength, colour=Management) + 
  geom_line(linewidth=1, aes(group=Management)) + 
  scale_linetype_manual(values=c("solid")) +
  geom_point(size=3.5, position=position_dodge(0.1)) + 
  scale_color_manual(values=c("NMR"="#AA2020", "KFG"="#F26A60", 
                              "ViRCP_old"="#095A62", "ViRCP_new"="#28A7B5", 
                              "NFG"="#9AE3EC")) +
  MyTheme + theme(legend.position="none") +
  ylab("\n ") + xlab(NULL) +
  geom_vline(xintercept=2015.8, col="gray20",lty=3,lwd=0.5) +
  ggtitle(" ")

plk_lngth2 <- plk_lngth + geom_errorbar(aes(ymin=MeanLength - ML_SE, 
                                            ymax=MeanLength + ML_SE), 
                                        width=0.3, linewidth=0.5, position=position_dodge(0.1)) 
plk_lngth2



# Arrange abundance & length panels

left <- ggarrange(bcg_abu2, crl_abu2, dtr_abu2, exsc_abu2, mobi_abu2, pis_abu2, plk_abu2, 
                  nrow=7, align="v", labels=c("a","c","e","g","i","k","m"),
                  common.legend=T, legend="top")

right <- ggarrange(bcg_lngth2, crl_lngth2, dtr_lngth2, exsc_lngth2, mobi_lngth2, pis_lngth2, plk_lngth2, 
                   nrow=7, align="v", labels=c("b","d","f","h","j","l","n"),
                   common.legend=T, legend="top")

all <- ggarrange(left, right, align="h")
all

#pdf("plots/FG-Abundance&MeanLength.pdf", width=11, height=14.5)
#all
#dev.off()


#####################################################################

### GLMMs for reef fish associations with benthic state




# Handy functions from Highland Statistics:
panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
    cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex * r)
}


Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) points(x, y, 
                                             pch = 16, cex = 0.8, 
                                             col = gray(0.1)))
}




### DATA MANIPULATION ###

fish <- read.csv("data/fish_full_new.csv") 
head(fish)


# Regroup management levels
VIR.old <- c("VIR1","VIR2","VIR3","VIR4","VIR5","VIR11")
VIR.new <- c("VIR8","VIR9")
fish$Management <- ifelse(fish$Site %in% VIR.old, "ViRCP_old", fish$Management)
fish$Management <- ifelse(fish$Site %in% VIR.new, "ViRCP_new", fish$Management)

fish$Management <- recode_factor(fish$Management, 
                                 Namena_MPA="NMR",
                                 Namena_LMMA="KFG",
                                 Vatuira_LMMA="NFG")

table(fish$Management)
# NMR    KFG    NFG   ViRCP_old    ViRCP_new
#3191   2780   1856        3740         1633

fish$Site             <- as.factor(fish$Site)
fish$Management       <- as.factor(fish$Management)


# Merge some trophic groups together:
unique(fish$Functional.group)

# Excavators & scrapers
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "excavator/scraper" = c("excavator", "scraper"))

# Micro- & Macro-invertivores & Spongivores
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "invertivore" = c("micro_invertivore", 
                                                        "macro_invertivore", 
                                                        "spongivore"))

# Browsers & croppers/grazers
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "browser/cropper/grazer" = c("browser", "cropper/grazer"))

# Piscivores & pisci-invertivores
fish$Functional.group <- fct_collapse(fish$Functional.group, 
                                      "piscivore" = c("piscivore", "pisci-invertivore"))



# SUM count and biomass per transect for each FG:
fish2 <- data.frame(aggregate(fish1[,c("Count", "Biomass_kgha")], 
                              by=list(fish1$Transect, fish1$Observers, 
                                      fish1$Year, fish1$Site, fish1$Management, 
                                      fish1$Functional.group), sum, na.rm=T)) 
new.names <- c("Transect","Observer","Year","Site","Management","Func.group")
names(fish2)[1:6] <- new.names


# Expand data to include zeros when a FG was not recorded in a transect:
vals <- expand.grid(Transect=unique(fish2$Transect),
                    Site=unique(fish2$Site),
                    Management=unique(fish2$Management),
                    Year=unique(fish2$Year),
                    Func.group=unique(fish2$Func.group))
vals

# Merge into dataset:
fish3 <- merge(fish2, vals, all=T)

head(fish3)

# Merging "vals" into the data adds a lot of unwanted NAs (transects that did 
# not exist).
# Need to remove sites that weren't in certain management areas:
unique(fish2$Management)
table(fish2[which(fish2$Management=="NMR"),]$Site) # KB04, KB05, N19, N20
table(fish2[which(fish2$Management=="KFG"),]$Site) # C3, C5, KB06, KB07
table(fish2[which(fish2$Management=="NFG"),]$Site) # VIR6, VIR7, VIR10
table(fish2[which(fish2$Management=="ViRCP_old"),]$Site) # VIR1, VIR2, VIR3, VIR4, VIR5, VIR11
table(fish2[which(fish2$Management=="ViRCP_new"),]$Site) # VIR8, VIR9

nmr <- c("KB04", "KB05", "N19", "N20")
kfg <- c("C3", "C5", "KB06", "KB07")
nfg <- c("VIR6", "VIR7", "VIR10")
vircp_old <- c("VIR1", "VIR2", "VIR3", "VIR4", "VIR5", "VIR11")
vircp_new <- c("VIR8", "VIR9")

fish3 <- fish3[-which(fish3$Management=="NMR" & !(fish3$Site %in% nmr)),]
fish3 <- fish3[-which(fish3$Management=="KFG" & !(fish3$Site %in% kfg)),]
fish3 <- fish3[-which(fish3$Management=="NFG" & !(fish3$Site %in% nfg)),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & !(fish3$Site %in% vircp_old)),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_new" & !(fish3$Site %in% vircp_new)),]

# No KFG or NMR sites were surveyed in 2018, so remove:
fish3 <- fish3[-which(fish3$Management=="KFG" & fish3$Year==2018),]
fish3 <- fish3[-which(fish3$Management=="NMR" & fish3$Year==2018),]

# Still some sites with all NA values:

# 2013:
# VIR10, NFG
# VIR11, ViRCP_old
table(fish2[which(fish2$Management=="NFG" & fish2$Year==2013),]$Site) 
table(fish2[which(fish2$Management=="ViRCP_old" & fish2$Year==2013),]$Site) 

# 2014:
# VIR1, VIR5, VIR11; ViRCP_old
# VIR10, NFG
# VIR9; ViRCP_new
table(fish2[which(fish2$Management=="ViRCP_old" & fish2$Year==2014),]$Site) 
table(fish2[which(fish2$Management=="NFG" & fish2$Year==2014),]$Site) 
table(fish2[which(fish2$Management=="ViRCP_new" & fish2$Year==2014),]$Site) 

# 2016:
# VIR4; ViRCP_old
table(fish2[which(fish2$Management=="ViRCP_old" & fish2$Year==2016),]$Site) 

# 2018:
# VIR4; ViRCP_old
table(fish2[which(fish2$Management=="ViRCP_old" & fish2$Year==2018),]$Site) 

# 2020:
# KB04; NMR
# VIR5; ViRCP_old
table(fish2[which(fish2$Management=="NMR" & fish2$Year==2020),]$Site) 
table(fish2[which(fish2$Management=="ViRCP_old" & fish2$Year==2020),]$Site) 


# None of the sites were surveyed in those years, so remove all:
#2013
fish3 <- fish3[-which(fish3$Management=="NFG" & fish3$Year==2013 & fish3$Site=="VIR10"),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2013 & fish3$Site=="VIR11"),]
#2014
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2014 & fish3$Site=="VIR1"),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2014 & fish3$Site=="VIR5"),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2014 & fish3$Site=="VIR11"),]
fish3 <- fish3[-which(fish3$Management=="NFG" & fish3$Year==2014 & fish3$Site=="VIR10"),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_new" & fish3$Year==2014 & fish3$Site=="VIR9"),]
#2016
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2016 & fish3$Site=="VIR4"),]
#2018
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2018 & fish3$Site=="VIR4"),]
#2020
fish3 <- fish3[-which(fish3$Management=="NMR" & fish3$Year==2020 & fish3$Site=="KB04"),]
fish3 <- fish3[-which(fish3$Management=="ViRCP_old" & fish3$Year==2020 & fish3$Site=="VIR5"),]


# Replace any biomass NAs with 0 (i.e. a feeding group was not observed):
fish3$Count <- ifelse(is.na(fish3$Count), 0, fish3$Count)
fish3$Biomass_kgha <- ifelse(is.na(fish3$Biomass_kgha), 0, fish3$Biomass_kgha)


benthic <- read.csv("data/benthic_category.csv", sep = ";")
head(benthic)
#100 points per transect (along 50m)
# bare = bare substrate
# CCA = crustose coralline algae
# BCM = cyanobacteria
# HC = hard coral
# MA = macroalgae
# other = other invertebrates
# rubble
# sand
# SC = soft coral
# TA = Turf algae


#Recode management names so fish and benthic datasets match:

table(benthic$Management)

benthic$Management <- recode_factor(benthic$Management, 
                                    KubulauFG="KFG",
                                    NakorotubuFG="NFG")
table(benthic$Management)
table(fish3$Management)


# Remove "Time","Latitude","Longitude","Exposure" and "Reef.type" from benthic data:
benthic2 <- benthic %>% select(-one_of(c("Latitude","Longitude",
                                         "Exposure","Reef.type")))
head(benthic2)


df <- left_join(fish3, benthic2, by=c("Year","Transect","Management","Site"))
head(df)

# Add in fished vs protected variable:
fishing <- c("KFG","NFG")
df$Fishery <- ifelse(df$Management %in% fishing, "Fished", "Protected")


#write.csv(df, "data/Fish-Benthos-analysis.csv", row.names=F)

# Observer name gaps filled in manually in Excel. And removed NA rows.

################################################################################


### DATA EXPLORATION ###

df <- read.csv("data/Fish-Benthos-analysis_edit_all-yrs.csv")
head(df)

unique(df$Year)
# 2013   2014   2016   2018   2020     


# GLMM for each FG:  fish biomass ~ HC + MA + TA + BCM + rubble + Management + 
#                                   (1|Year) + (1|Site)

# We are comparing between the effect sizes of different variables - best
# to scale numeric covariates ( x - mean(x) / sd(x) )
df2$HC.s <- scale(df2$HC)
df2$MA.s <- scale(df2$MA)
df2$TA.s <- scale(df2$TA)
df2$BCM.s <- scale(df2$BCM)
df2$rubble.s <- scale(df2$rubble)

# Create separate datasets for each FG:
bcg <- df2[which(df2$Func.group=="browser/cropper/grazer"),]
crl <- df2[which(df2$Func.group=="corallivore"),]
det <- df2[which(df2$Func.group=="detritivore"),]
exc <- df2[which(df2$Func.group=="excavator/scraper"),]
inv <- df2[which(df2$Func.group=="invertivore"),]
pis <- df2[which(df2$Func.group=="piscivore"),]
plk <- df2[which(df2$Func.group=="planktivore"),]

# Check for outliers in response:
dotchart(bcg$Biomass_kgha) # ok, a couple of high values >2000
dotchart(crl$Biomass_kgha) # ok
dotchart(det$Biomass_kgha) # ok, 1 high value ~500
dotchart(exc$Biomass_kgha) # ok, 3 high values >2000
dotchart(inv$Biomass_kgha) # ok, 3 high values >2500
dotchart(pis$Biomass_kgha) # ok
dotchart(plk$Biomass_kgha) # ok, 2 high values >4000

# outliers in explanatory variables (all FG datasets will be the same):
dotchart(bcg$HC.s)
dotchart(bcg$MA.s)
dotchart(bcg$TA.s)
dotchart(bcg$BCM.s)
dotchart(bcg$rubble.s)
# All good.


table(bcg$Year)
#2013 2014 2016 2018 2020 
#  75   39   54   60  101 

table(bcg$Site)

table(bcg$Fishery)
#Fished  Protected 
#   129        200


b.bio <- ggplot(data=bcg) + 
  geom_point(aes(y=Biomass_kgha, x=Year, group=Site, col=Site)) + 
  facet_grid(~Fishery)
b.bio # some higher values in protected areas

c.bio <- ggplot(data=crl) + 
  geom_point(aes(y=Biomass_kgha, x=Year, group=Site, col=Site)) + 
  facet_grid(~Fishery)
c.bio 

d.bio <- ggplot(data=det) + 
  geom_point(aes(y=Biomass_kgha, x=Year, group=Site, col=Site)) + 
  facet_grid(~Fishery)
d.bio  # 1 high value in 2020

e.bio <- ggplot(data=exc) + 
  geom_point(aes(y=Biomass_kgha, x=Year, group=Site, col=Site)) + 
  facet_grid(~Fishery)
e.bio 

i.bio <- ggplot(data=inv) + 
  geom_point(aes(y=Biomass_kgha, x=Year, group=Management, col=Management)) + 
  facet_grid(~Fishery)
i.bio 


# Check homogeneity of variance

#Year
boxplot(Biomass_kgha ~ Year, data=bcg) 
boxplot(Biomass_kgha ~ Year, data=crl) 
boxplot(Biomass_kgha ~ Year, data=det) 
boxplot(Biomass_kgha ~ Year, data=exc) 
boxplot(Biomass_kgha ~ Year, data=inv) 
boxplot(Biomass_kgha ~ Year, data=pis) 
boxplot(Biomass_kgha ~ Year, data=plk) 

#Management
boxplot(Biomass_kgha ~ Fishery, data=bcg) 
boxplot(Biomass_kgha ~ Fishery, data=crl) 
boxplot(Biomass_kgha ~ Fishery, data=det) 
boxplot(Biomass_kgha ~ Fishery, data=exc) 
boxplot(Biomass_kgha ~ Fishery, data=inv) 
boxplot(Biomass_kgha ~ Fishery, data=pis) 
boxplot(Biomass_kgha ~ Fishery, data=plk) 

# Normally distributed?
ggplot(data=bcg, aes(x=Biomass_kgha)) +
  geom_histogram(alpha=0.6, position = 'identity', binwidth=20) +
  labs(fill="") + facet_wrap(~Year)

ggplot(data=crl, aes(x=Biomass_kgha)) +
  geom_histogram(alpha=0.6, position = 'identity', binwidth=20) +
  labs(fill="") + facet_wrap(~Year)

ggplot(data=det, aes(x=Biomass_kgha)) +
  geom_histogram(alpha=0.6, position = 'identity', binwidth=20) +
  labs(fill="") + facet_wrap(~Year)

ggplot(data=exc, aes(x=Biomass_kgha)) +
  geom_histogram(alpha=0.6, position = 'identity', binwidth=20) +
  labs(fill="") + facet_wrap(~Year)

ggplot(data=inv, aes(x=Biomass_kgha)) +
  geom_histogram(alpha=0.6, position = 'identity', binwidth=20) +
  labs(fill="") + facet_wrap(~Year)

ggplot(data=plk, aes(x=Biomass_kgha)) +
  geom_histogram(alpha=0.6, position = 'identity', binwidth=20) +
  labs(fill="") + facet_wrap(~Year)

#all right-skewed (lots of 0s) - better to log biomass values


# Zeros in data:

length(bcg[which(bcg$Biomass_kgha==0),]$Biomass_kgha)
#  13
length(bcg$Biomass_kgha)
# 329

13/329
# 4.0% are 0s


length(crl[which(crl$Biomass_kgha==0),]$Biomass_kgha)
# 42
length(crl$Biomass_kgha)
# 329

42/329
# 12.8% are 0s


length(det[which(det$Biomass_kgha==0),]$Biomass_kgha)
# 87
length(det$Biomass_kgha)
# 329

87/329
# 26.4% are 0s


length(exc[which(exc$Biomass_kgha==0),]$Biomass_kgha)
# 10
length(exc$Biomass_kgha)
# 329

10/329
# 3.0% are 0s


length(inv[which(inv$Biomass_kgha==0),]$Biomass_kgha)
# 1


length(pis[which(pis$Biomass_kgha==0),]$Biomass_kgha)
# 8
length(pis$Biomass_kgha)
# 329

8/329
# 2.4% are 0s


length(plk[which(plk$Biomass_kgha==0),]$Biomass_kgha)
# 61
length(plk$Biomass_kgha)
# 329

61/329
# 18.5% are 0s



# Check collinearity between variables:

Mypairs(bcg[,c("HC.s","MA.s","TA.s","BCM.s","rubble.s","Year")]) 
# the only strong correlation is between hard coral and rubble.
# could remove rubble if this causes problems (check VIF).


# Relationships between X and Y:

# Hard coral
xyplot(Biomass_kgha ~ HC.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=bcg,
       panel=function(x,y){
         panel.points(x, y, col = 1)  #Add points
         panel.abline(lm(y~x))        #Add regression line
       })

xyplot(Biomass_kgha ~ HC.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=crl,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
# increase in abundance with coral cover

xyplot(Biomass_kgha ~ HC.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=det,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ HC.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=exc,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
# slight increase on fished and protected sites with higher coral cover

xyplot(Biomass_kgha ~ HC.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=inv,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ HC.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=plk,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })



# Macroalgae
xyplot(Biomass_kgha ~ MA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=bcg,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ MA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=crl,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ MA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=det,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
#decrease with increasing macroalgae - especially fished reefs

xyplot(Biomass_kgha ~ MA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=exc,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
#decrease with higher macroalgae 

xyplot(Biomass_kgha ~ MA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=inv,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ MA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=plk,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })



# Turf algae
xyplot(Biomass_kgha ~ TA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=bcg,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ TA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=crl,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ TA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=det,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
#increase with higher turf on fished reefs, slight decrease on protected.

xyplot(Biomass_kgha ~ TA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=exc,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ TA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=inv,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ TA.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=plk,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
# decrease with increasing TA on fished reefs, slight increase on protected.



# Cyanobacteria
xyplot(Biomass_kgha ~ BCM.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=bcg,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ BCM.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=crl,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ BCM.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=det,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
#increase with higher BCM, especially on fished reefs.

xyplot(Biomass_kgha ~ BCM.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=exc,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ BCM.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=inv,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ BCM.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=plk,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })



# Rubble
xyplot(Biomass_kgha ~ rubble.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=bcg,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ rubble.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=crl,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
# biomass decrease with increasing rubble

xyplot(Biomass_kgha ~ rubble.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=det,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
# biomass decrease with higher rubble on fished reefs

xyplot(Biomass_kgha ~ rubble.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=exc,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })
# biomass decrease with increasing rubble

xyplot(Biomass_kgha ~ rubble.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=inv,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })

xyplot(Biomass_kgha ~ rubble.s | Fishery,
       scales=list(y=list(relation="free")), auto.key=T, data=plk,
       panel=function(x,y){
         panel.points(x, y, col = 1)  
         panel.abline(lm(y~x))        
       })



### Running the GLMMs


# benthic effects on fish biomass

# BROWSER/CROPPER/GRAZER

hist(bcg$Biomass_kgha, breaks=50) # right-skewed
hist(log(bcg$Biomass_kgha+1), breaks=50) #better

bcg1 <- lmer(Biomass_kgha ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + Fishery + 
               (1|Year) + (1|Site), data=bcg)
summary(bcg1)
plot(bcg1)  # unequal variances - try log of response

# also include Observer as a FIXED effect (random not appropriate as some
# samplers were only involved in one or two years)
bcg2 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + Observer + (1|Year) + (1|Site), data=bcg)
summary(bcg2)
( bcg.res <- plot(bcg2) )  # good

Eb2 <- resid(bcg2, type="pearson")
hist(Eb2, breaks=50) # good

par(mar = rep(3, 4))
qqnorm(resid(bcg2))
qqline(resid(bcg2)) 
( bcg.qq <- recordPlot() )
dev.off()

check_collinearity(bcg2)
# all VIF <3

( bcg.f <- plot_model(bcg2, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects:
bcg.r <- plot_model(bcg2, type="re", line.size=1, dot.size=3,
                    vline.color="gray70") 
bcg.r[1]
bcg.r[2]


# CORALLIVORE

hist(crl$Biomass_kgha, breaks=50) #right-skewed
hist(log(crl$Biomass_kgha+1), breaks=50) # lots of 0s

crl1 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + (1|Year) + (1|Site) + Observer, data=crl)
summary(crl1) 
( crl.res <- plot(crl1) ) # ok (line due to zeros)

Ec1 <- resid(crl1, type="pearson")
hist(Ec1, breaks=50) # good

qqnorm(resid(crl1))
qqline(resid(crl1))
( crl.qq <- recordPlot() )
dev.off()

check_collinearity(crl1)
#all VIF <3

( crl.f <- plot_model(crl1, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects
crl.r <- plot_model(crl1, type="re", line.size=1, dot.size=3,
                    vline.color="gray70") 
crl.r[1]
crl.r[2]


# DETRITIVORE

hist(det$Biomass_kgha, breaks=50) #right-skewed
hist(log(det$Biomass_kgha+1), breaks=50)  # lots of 0s

det1 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + (1|Year) + (1|Site) + Observer, data=det)
summary(det1)
( det.res <- plot(det1) )  # ok - line from the zeros

Ed1 <- resid(det1, type="pearson")
hist(Ed1, breaks=50) # ok 

qqnorm(resid(det1))
qqline(resid(det1)) 
( det.qq <- recordPlot() )
dev.off()

check_collinearity(det1)
# all VIF <3

( det.f <- plot_model(det1, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects
det.r <- plot_model(det1, type="re", line.size=1, dot.size=3,
                    vline.color="gray70")
det.r[1]
det.r[2]


# EXCAVATOR/SCRAPER

hist(exc$Biomass_kgha, breaks=50) #right-skewed
hist(log(exc$Biomass_kgha+1), breaks=50) #ok

exc1 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + (1|Year) + (1|Site) + Observer, data=exc)
summary(exc1)
( exc.res <- plot(exc1) )  # good

Ee1 <- resid(exc1, type="pearson")
hist(Ee1, breaks=50) # ok

qqnorm(resid(exc1))
qqline(resid(exc1)) 
( exc.qq <- recordPlot() )
dev.off()

check_collinearity(exc1)
# all VIF <3

( exc.f <- plot_model(exc1, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects
exc.r <- plot_model(exc1, type="re", line.size=1, dot.size=3,
                    vline.color="gray70")
exc.r[1]
exc.r[2]


# INVERTIVORE

hist(inv$Biomass_kgha, breaks=50) #right-skewed
hist(log(inv$Biomass_kgha+1), breaks=50) #ok

inv1 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + (1|Year) + (1|Site) + Observer, data=inv)
summary(inv1)
( inv.res <- plot(inv1) )  # good

Ei1 <- resid(inv1, type="pearson")
hist(Ei1, breaks=50) # good

qqnorm(resid(inv1))
qqline(resid(inv1)) 
( inv.qq <- recordPlot() )
dev.off()

check_collinearity(inv1)
# all VIF <3

( inv.f <- plot_model(inv1, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects
inv.r <- plot_model(inv1, type="re", line.size=1, dot.size=3,
                    vline.color="gray70")
inv.r[1]
inv.r[2]


# PISCIVORE

hist(pis$Biomass_kgha, breaks=50) #right-skewed
hist(log(pis$Biomass_kgha+1), breaks=50) #ok

pis1 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + (1|Year) + (1|Site) + Observer, data=pis)
summary(pis1)
( pis.res <- plot(pis1) )  # good

Ei1 <- resid(pis1, type="pearson")
hist(Ei1, breaks=50) # good

qqnorm(resid(pis1))
qqline(resid(pis1))
( pis.qq <- recordPlot() )
dev.off()

check_collinearity(pis1)
# all VIF <3

( pis.f <- plot_model(pis1, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects
pis.r <- plot_model(pis1, type="re", line.size=1, dot.size=3,
                    vline.color="gray70")
pis.r[1]
pis.r[2]


# PLANKTIVORE

hist(plk$Biomass_kgha, breaks=50) #right-skewed
hist(log(plk$Biomass_kgha+1), breaks=50) # lots of zeros

plk1 <- lmer(log(Biomass_kgha+1) ~ HC.s + MA.s + TA.s + BCM.s + rubble.s + 
               Fishery + (1|Year) + (1|Site) + Observer, data=plk)

summary(plk1)
( plk.res <- plot(plk1) )  # ok - line from zeros

Ep1 <- resid(plk1, type="pearson")
hist(Ep1, breaks=50) # left-skewed

qqnorm(resid(plk1))
qqline(resid(plk1))
( plk.qq <- recordPlot() )
dev.off()

check_collinearity(plk1)
# all VIF <3

( plk.f <- plot_model(plk1, show.values=F,
                      value.offset = .2, vline.color = "gray") )

# random effects
plk.r <- plot_model(plk1, type="re", line.size=1, dot.size=3,
                    vline.color="gray70")
plk.r[1]
plk.r[2] 




### Plot fixed effects for Figure 7:

unique(bcg$Fishery) # "Fished"    "Protected"

bcg.fixed <- plot_model(bcg2, line.size=1, dot.size=3, 
                        value.offset=0.2, title="Browser/cropper/grazer", 
                        vline.color="gray70", 
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected"),  # need factor LEVEL too
                        axis.labels=c("MCA","Rubble","Cyanobacteria",
                                      "Turf algae","Macroalgae","Hard coral")) +
  font_size(labels.x=10, labels.y=13) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
bcg.fixed 


crl.fixed <- plot_model(crl1, line.size=1, 
                        dot.size=3, value.offset=0.2, title="Corallivore", 
                        vline.color="gray70",
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected")) +
  font_size(labels.x=10) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
crl.fixed


det.fixed <- plot_model(det1, line.size=1, 
                        dot.size=3, value.offset=0.2, title="Detritivore", 
                        vline.color="gray70",
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected")) +
  font_size(labels.x=10) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
det.fixed


exc.fixed <- plot_model(exc1, line.size=1, 
                        dot.size=3, value.offset=0.2, title="Excavator/scraper", 
                        vline.color="gray70",
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected")) +
  font_size(labels.x=10) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
exc.fixed


inv.fixed <- plot_model(inv1, line.size=1, 
                        dot.size=3, value.offset=0.2, title="Invertivore", 
                        vline.color="gray70",
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected"),
                        axis.labels=c("MCA","Rubble","Cyanobacteria",
                                      "Turf algae","Macroalgae","Hard coral")) +
  font_size(labels.x=10,labels.y=13) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
inv.fixed


pis.fixed <- plot_model(pis1, line.size=1, 
                        dot.size=3, value.offset=0.2, title="Piscivore", 
                        vline.color="gray70",
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected")) +
  font_size(labels.x=10) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
pis.fixed


plk.fixed <- plot_model(plk1, line.size=1, 
                        dot.size=3, value.offset=0.2, title="Planktivore", 
                        vline.color="gray70",
                        terms=c("HC.s","MA.s","TA.s","BCM.s","rubble.s",
                                "FisheryProtected")) +
  font_size(labels.x=10) + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank()) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
plk.fixed


# Arrange plot panels:

all.FGs <- ggarrange(bcg.fixed, crl.fixed, det.fixed, exc.fixed, 
                     inv.fixed, pis.fixed, plk.fixed,
                     nrow=2, ncol=4,  widths = c(1.4,1,1,1,
                                                 1.4,1,1))
all.FGs

#pdf("plots/Fish-benthos-GLMM.pdf", width=10, height=6)
#all.FGs
#dev.off()


# Random effects - supplementary material only.
# Years
bcg.r[2]
crl.r[2]
det.r[2]
exc.r[2]
inv.r[2]
pis.r[2]
plk.r[2] 



##### END OF SCRIPT #####