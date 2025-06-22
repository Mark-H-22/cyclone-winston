#Benthic models to test for significant changes over time in all benthic groups

library(multcomp)
require(ggeffects)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(Matrix)
library(lmerTest)
library(lme4)
library(glmmTMB)
library(dplyr)
library(ggplot2)
library(nlme)
library(car)
library(emmeans)
library(DHARMa)
library(tidyr)             
library(reshape2) 

setwd("C:/Users/akfor/OneDrive/IN DEVELOPMENT/TC Winston Recovery/R-files/Github")
dat <- read.delim("benthic_full.txt")

head(dat)
str(dat)

#Make factors where necessary
dat <- dat %>%
  mutate_at(vars(Country, Site, Month, Day, Year, Exposure, Reef.slope, Reef.type, Reef.zone, Relative.depth, Transect.number, Area, Benthic.category), as.factor)

dat$Benthic.category   
length(unique(dat$Benthic.category))
list(unique(dat$Benthic.category))

#Break up ViRCP into old and new
dat <- dat %>%
  mutate(Area = as.factor(case_when(
    grepl("^VIR[1-4]|VIR10|VIR11", Site) ~ "ViRCP_old",
    grepl("^VIR[89]", Site) ~ "ViRCP_new",
    TRUE ~ as.character(Area)
  )))

#Convert dataframe to wide format so that benthic categories have columns with counts
library(data.table)
library(stringr)

dat2 <- dcast(setDT(dat), Year+Site+Area+Transect.number ~ dat$Benthic.category, length) # replace colon with plus
#dat2$dat <- as.character(dat2$dat)
#dat3 <- separate(data = dat2, col = dat, into = c("Year", "Site", "Management", "Transect"), sep = "\\:")
dat4 <- rename(dat2, bare = "Bare substrate", CCA = "Crustose coralline algae", BCM = "Cyanobacteria", HC = "Hard coral", MA = "Macroalgae", other = "Other invertebrates", rubble = "Rubble", sand = "Sand", SC = "Soft coral", TA = "Turf algae")  

dat_check <- mutate(dat4, total.points = bare + CCA + BCM + HC + MA + other + rubble + sand + SC + TA)
#All rows add up to 100 so values represent percentage!

#Unsure how to get rid of the 0 data so export and reimport
write.table(dat4, file = "benthic_edited.csv", sep = ";", dec = ".")

dat <- read.csv("benthic_edited.csv", sep = ";")

head(dat)
str(dat)
dat$Site <- as.factor(dat$Site)
dat$Area <- as.factor(dat$Area)

dat <- dat %>%
  group_by(Site) %>%
  mutate(Period = case_when(
    # If both 2013 and 2014 data are present for a site, assign pre-cyclone to 2014 data
    Year == "2014" & any(Year == "2013") ~ "pre-cyclone",
    # Assign pre-cyclone to 2013 data if there is no 2014 data
    Year == "2013" & !any(Year == "2014") ~ "pre-cyclone",
    # Assign other periods based on the year
    Year == "2016" ~ "post-cyclone",
    Year == "2018" ~ "recovery-2018",
    Year == "2020" ~ "recovery-2020"
  )) %>%
  ungroup()


dat$Year <- as.factor(dat$Year)
dat$Period <- as.factor(dat$Period)

dat$Period <- factor(dat$Period, levels = c('pre-cyclone', 'post-cyclone', 'recovery-2018', 'recovery-2020'))
dat <- dat[complete.cases(dat$Period), ]

dat <- subset(dat, dat$Year!="2018")

#Models for hard coral (HC)
hist(dat$HC) #not skewed

plot_coral <- ggplot(dat, aes(x = Period, y = HC, colour = Area)) + 
  ylab("Hard Coral Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_coral

# Fit a linear mixed-effects model with all periods as factors
coral_model <- lmer(HC ~ Period * Area + (1 | Site), data = dat)
summary(coral_model)

plot(coral_model) #looks good
qqnorm(resid(coral_model)) #looks good
qqline(resid(coral_model)) #looks good

plot_model(coral_model, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic()

Anova(coral_model, type = "III")

marg_means <- emmeans(coral_model, ~ Period:Area) # Compute marginal means
pairwise_comp <- pairs(marg_means, by="Area") # Perform pairwise comparisons by Management
print(pairwise_comp)

############################################################################################################################################

#Models for rubble

hist(dat$rubble) #right skewed

plot_rubble <- ggplot(dat, aes(x = Period, y = rubble, colour = Area)) + 
  ylab("Rubble Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_rubble

# Fit a linear mixed-effects model with all periods as factors
rubble_model <- lmer(rubble ~ Period * Area + (1 | Site), data = dat)
summary(rubble_model) #issues

plot(rubble_model) #definitely some patterns in the residuals
qqnorm(resid(rubble_model)) 
qqline(resid(rubble_model)) 

#Should fit a model with a different distribution - negative binomial?
library(glmmTMB)
library(DHARMa)

rubble_model <- glmmTMB(rubble ~ Period*Area + (1|Site), family = "nbinom1", data= dat) #can change to nbinom1 but AIC lower for nbinom1
summary(rubble_model)
res = simulateResiduals(rubble_model) 
plot(res, asFactor = T) #issues with levene's test, tried multiple different distributions, add dispformula


rubble_model1 <- glmmTMB(rubble ~ Period*Area + (1|Site), family = "nbinom1", data= dat, dispformula = ~ Area) #can change to nbinom1 but AIC lower for nbinom1
summary(rubble_model1)
res = simulateResiduals(rubble_model1) 
plot(res, asFactor = T) #looks good

rubble_model2 <- glmmTMB(rubble ~ Period*Area + (1|Site), family = "nbinom2", data= dat, dispformula = ~ Area) #can change to nbinom1 but AIC lower for nbinom1
AIC(rubble_model1, rubble_model2) #1 has the lowest AIC so go with nbinom1

plot_model(rubble_model1, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic() #can't figure out the error here

Anova(rubble_model1, type = "III")

marg_means <- emmeans(rubble_model1, ~ Period:Area) # Compute marginal means
pairwise_comp <- pairs(marg_means, by="Area") # Perform pairwise comparisons by Management
print(pairwise_comp)

############################################################################################################################################

#Models for turf algae
hist(dat$TA) #right skewed

plot_turf <- ggplot(dat, aes(x = Period, y = TA, colour = Area)) + 
  ylab("Turf Algal Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_turf

# Fit a linear mixed-effects model with all periods as factors
turf_model <- lmer(TA ~ Period * Area + (1 | Site), data = dat)
summary(turf_model)

plot(turf_model) #definite pattern
qqnorm(resid(turf_model))
qqline(resid(turf_model))

#Should fit a model with a different distribution - negative binomial?
turf_model2 <- glmmTMB(TA ~ Period*Area + (1|Site), dispformula = ~Period*Area, family = tweedie(), data= dat) #tried nbinom1, nbinom2, tweedie and different dispformulas - this is the only one that works
summary(turf_model2)
res = simulateResiduals(turf_model2) #no issues
plot(res, asFactor = T)  #no issues

turf_model3 <- glmmTMB(TA ~ Period*Area + (1|Site), family = tweedie(), data= dat) #tried nbinom1, nbinom2, tweedie and different dispformulas - this is the only one that works
summary(turf_model3)

plot_model(turf_model3, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic() #also plots dispersion



Anova(turf_model2, type = "III")

marg_means <- emmeans(turf_model2, ~ Period:Area) # Compute marginal means
pairwise_comp <- pairs(marg_means, by="Area") # Perform pairwise comparisons by Management
print(pairwise_comp)

############################################################################################################################################

#Models for fleshy algae

hist(dat$MA) #right skewed

plot_macroalgae <- ggplot(dat, aes(x = Period, y = MA, colour = Area)) + 
  ylab("Macroalgal Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_macroalgae

# Fit a linear mixed-effects model with all periods as factors
macroalgal_model <- lme(MA ~ Period * Area, random = ~1 | Site, data = dat)
summary(macroalgal_model)

plot(macroalgal_model) #definitely some patterns in the residuals
qqnorm(resid(macroalgal_model)) 
qqline(resid(macroalgal_model)) 

#Should fit a model with a different distribution - negative binomial?
macroalgal_model1 <- glmmTMB(MA ~ Period*Area + (1|Site), family = "nbinom1", data= dat) 
summary(macroalgal_model1)
res = simulateResiduals(macroalgal_model1) 
plot(res, asFactor = T)#no issues

macroalgal_model2 <- glmmTMB(MA ~ Period*Area + (1|Site), family = "nbinom2", data= dat) 
summary(macroalgal_model2)
res = simulateResiduals(macroalgal_model2) 
plot(res, asFactor = T)#no issues

AIC(macroalgal_model1, macroalgal_model2) #nbinom2 is better

plot_model(macroalgal_model2, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic()

Anova(macroalgal_model2, type = "III")

#Pairwise test to see where 
marg_means <- emmeans(macroalgal_model2, ~ Period:Area) 
pairwise_comp <- pairs(marg_means, by="Area") 
print(pairwise_comp)

############################################################################################################################################

#Models for CCA

hist(dat$CCA) #right skewed

plot_CCA <- ggplot(dat, aes(x = Period, y = CCA, colour = Area)) + 
  ylab("CCA Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_CCA

# Fit a linear mixed-effects model with all periods as factors
CCA_model <- lme(CCA ~ Period * Area, random = ~1 | Site, data = dat)
summary(CCA_model)

plot(CCA_model) #definitely some patterns in the residuals
qqnorm(resid(CCA_model)) 
qqline(resid(CCA_model)) 

#Should fit a model with a different distribution - negative binomial?
CCA_model1 <- glmmTMB(CCA ~ Period*Area + (1|Site), family = "nbinom1", data= dat) #can change to nbinom1 but AIC lower for nbinom1
summary(CCA_model1)
res = simulateResiduals(CCA_model1) #no issues
plot(res, asFactor = T)

CCA_model2 <- glmmTMB(CCA ~ Period*Area + (1|Site), family = "nbinom2", data= dat) #can change to nbinom1 but AIC lower for nbinom1
summary(CCA_model2)
res = simulateResiduals(CCA_model2) #no issues
plot(res, asFactor = T)

AIC(CCA_model1, CCA_model2) #1 is lower

plot_model(CCA_model1, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic()

Anova(CCA_model1, type = "III")

#Pairwise test to see where 
marg_means <- emmeans(CCA_model1, ~ Period:Area) # Compute marginal means
pairwise_comp <- pairs(marg_means, by="Area") # Perform pairwise comparisons by Management
print(pairwise_comp)

############################################################################################################################################

#Models for soft coral

hist(dat$SC) #right skewed

plot_SC <- ggplot(dat, aes(x = Period, y = SC, colour = Area)) + 
  ylab("Soft Coral Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_SC

# Fit a linear mixed-effects model with all periods as factors
SC_model <- lme(SC ~ Period * Area, random = ~1 | Site, data = dat)
summary(SC_model)

plot(SC_model) #definitely some patterns in the residuals

#Should fit a model with a different distribution - negative binomial?
SC_model1 <- glmmTMB(SC ~ Period*Area + (1|Site), family = "nbinom1", data= dat) #can change to nbinom1 but AIC lower for nbinom1
summary(SC_model1)
res = simulateResiduals(SC_model1) #no issues
plot(res, asFactor = T)

SC_model2 <- glmmTMB(SC ~ Period*Area + (1|Site), family = "nbinom2", data= dat) #can change to nbinom1 but AIC lower for nbinom1
summary(SC_model2)
res = simulateResiduals(SC_model2) #no issues
plot(res, asFactor = T)

AIC(SC_model1, SC_model2) #1 is lower but it violates levene's test

plot_model(SC_model2, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic()

Anova(SC_model2, type = "III")

marg_means <- emmeans(SC_model2, ~ Period:Area) 
pairwise_comp <- pairs(marg_means, by="Area") 
print(pairwise_comp)

############################################################################################################################################

#Models for BCM

hist(dat$BCM) #right skewed

plot_BCM <- ggplot(dat, aes(x = Period, y = BCM, colour = Area)) + 
  ylab("BCM Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_BCM

# Fit a linear mixed-effects model with all periods as factors
BCM_model <- lme(BCM ~ Period * Area, random = ~1 | Site, data = dat)
summary(BCM_model)

plot(BCM_model) #definitely some patterns in the residuals

#Should fit a model with a different distribution - negative binomial?
BCM_model1 <- glmmTMB(BCM ~ Period*Area + (1|Site), family = "nbinom1", data= dat) #can change to nbinom1 but AIC lower for nbinom1
summary(BCM_model1)
res = simulateResiduals(BCM_model1) #no issues
plot(res, asFactor = T)

BCM_model2 <- glmmTMB(BCM ~ Period*Area + (1|Site), family = "nbinom2", data= dat) 
summary(BCM_model2)
res = simulateResiduals(BCM_model2) #no issues
plot(res, asFactor = T)

AIC(BCM_model1, BCM_model2) #1 is lower

plot_model(BCM_model1, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic() 

Anova(BCM_model1, type = "III")

#Pairwise test to see where 
marg_means <- emmeans(BCM_model1, ~ Period:Area) 
pairwise_comp <- pairs(marg_means, by="Area") 
print(pairwise_comp)


############################################################################################################################################

#Models for bare substrate

hist(dat$bare) #right skewed

plot_bare <- ggplot(dat, aes(x = Period, y = bare, colour = Area)) + 
  ylab("Bare Substrate Cover (%)") +
  facet_wrap(~Area) +
  geom_boxplot(alpha = 0.5) +
  #geom_line(data = cbind(dat, pred = predict(coral)), aes(y = pred), size = 1) +
  theme_classic()  
plot_bare

# Fit a linear mixed-effects model with all periods as factors
bare_model <- lme(bare ~ Period * Area, random = ~1 | Site, data = dat)
summary(bare_model)

plot(bare_model) #definitely some patterns in the residuals

#Should fit a model with a different distribution - negative binomial?
bare_model2 <- glmmTMB(bare ~ Period*Area + (1|Site), family = "nbinom2", data= dat) #nbinom1 violates Levene's test
summary(bare_model2)
res = simulateResiduals(bare_model2) 
plot(res, asFactor = T) #no issues

plot_model(bare_model2, show.values = TRUE,
           value.offset = .2, vline.color = "gray") + theme_classic() #can't figure out the error here

#Use type III ANOVA to test for general significance of interactive terms
Anova(bare_model2, type = "III")
#The significant p-value for the interaction term (Period:Management) suggests that the effect of Period on HC varies significantly across 
#different Management levels. In other words, the change in MA over time is not consistent across all Management levels; there are variations 
#in how Management influences the change in bare substrate.

#Pairwise test to see where 
marg_means <- emmeans(bare_model2, ~ Period:Area) # Compute marginal means
pairwise_comp <- pairs(marg_means, by="Area") # Perform pairwise comparisons by Management
print(pairwise_comp)
#no apparent impact of the cyclone apart from at NMR and to a lesser extend KubulauFG

#################################################################################

# Script for analysing changes in composition of genera following TC Winston
# Genus level analysis

library(dplyr)     
library(ggplot2)   
library(tidyr)             
library(reshape2) 
library(ggpubr)
library(gtsummary)

#setwd("C:/Users/akfor/OneDrive/IN DEVELOPMENT/TC Winston Recovery/R-files/Github")
dat <- read.delim("benthic_full.txt")

#Make factors where necessary
dat <- dat %>%
  mutate_at(vars(Country, Site, Month, Day, Year, Exposure, Reef.slope, Reef.type, Reef.zone, Relative.depth, Transect.number, Area, Benthic.category), as.factor)

dat$Benthic.attribute   
length(unique(dat$Benthic.attribute))
list(unique(dat$Benthic.attribute))

#Convert dataframe to wide format so that benthic categories have columns with counts
library(data.table)
library(stringr)

dat2 <- dcast(setDT(dat), Year+Site+Area+Transect.number ~ dat$Benthic.category, length) 
#dat2$dat <- as.character(dat2$dat)
#dat3 <- separate(data = dat2, col = dat, into = c("Year", "Site", "Management", "Transect"), sep = "\\:")
dat4 <- rename(dat2, bare = "Bare substrate", CCA = "Crustose coralline algae", BCM = "Cyanobacteria", HC = "Hard coral", MA = "Macroalgae", other = "Other invertebrates", rubble = "Rubble", sand = "Sand", SC = "Soft coral", TA = "Turf algae")  

dat_check <- mutate(dat4, total.points = bare + CCA + BCM + HC + MA + other + rubble + sand + SC + TA)
#All rows add up to 100 so values represent percentage!

write.table(dat4, file = "benthic_edited.csv", sep = ";", dec = ".")

dat <- read.csv("benthic_edited.csv", sep = ";")

head(dat)
str(dat)
dat$Site <- as.factor(dat$Site)
dat$Area <- as.factor(dat$Area)

# Filter rows where 'Benthic.category' is 'Hard coral'
hard_coral_df <- dat %>% filter(Benthic.category == 'Hard coral')

dat2 <- dcast(setDT(hard_coral_df), Year+Site+Area+Transect.number ~ hard_coral_df$Benthic.attribute, length) # replace colon with plus

# We will set up a new column called 'Period' where we will use the 2014 data for 'pre-cyclone' UNLESS there is on 2013 data
#in which case we will use this for the 'pre-cyclone' data for that site. This avoids duplication of sites in the pre-cyclone
#category. 2016 will become 'post-cyclone', 2018 'recovery-2018' and 2020 'recovery-2020'.
dat4 <- dat2 %>%
  group_by(Area, Site) %>%
  filter(!(Year == "2013" & "2014" %in% Year)) %>%  # Remove 2013 data if 2014 exists
  mutate(Period = case_when(
    Year == "2013" ~ "pre-cyclone",
    Year == "2014" ~ "pre-cyclone",
    Year == "2016" ~ "post-cyclone",
    Year == "2018" ~ "recovery-2018",
    Year == "2020" ~ "recovery-2020",
    TRUE ~ NA_character_
  )) %>%
  ungroup()

#Break up ViRCP into old and new
dat4 <- dat4 %>%
  mutate(Area = as.factor(case_when(
    grepl("^VIR[1-4]|VIR10|VIR11", Site) ~ "ViRCP_old",
    grepl("^VIR[89]", Site) ~ "ViRCP_new",
    TRUE ~ as.character(Area)
  )))

#Take a mean of each site/period
mean_genus <- dat4 %>%
  filter(!is.na(Period)) %>%
  group_by(Site, Area, Period) %>%
  summarize(across(`Acanthastrea`:`Zoopilus`, mean, na.rm = TRUE), .groups = 'drop')

mean_genus

long <- melt(setDT(mean_genus), id.vars = c("Period","Area", "Site"), variable.name = "Genus")

long <- long %>%
  mutate(Genus = as.character(Genus),
         Genus = ifelse(Genus == "Hard coral", "Other", Genus))

long <- long %>%
  group_by(Area, Genus, Period) %>%
  summarize(overall_mean = mean(value, na.rm = TRUE))

long <- long %>%
  group_by(Area, Genus, Period) %>%
  filter(any(overall_mean != 0)) %>%
  ungroup()

long_1 <- long %>%
  group_by(Genus) %>%
  mutate(Genus = ifelse(all(overall_mean < 1), 'Other', as.character(Genus))) %>%
  ungroup()


library(RColorBrewer)
library(viridis)
library(scales)

long_1 <- long_1 %>%
  mutate(Period = factor(Period, levels = c('pre-cyclone', 'post-cyclone', 'recovery-2018', 'recovery-2020')))

#install.packages("forcats") 
library(forcats)

#remove recovery-2018 to simplify graph
long_1 <- long_1 %>%
  mutate(Period = fct_recode(Period, 'pre' = 'pre-cyclone', 'post' = 'post-cyclone', 'recovery' = 'recovery-2020')) %>%
  filter(Period != 'recovery-2018')

# Now create the plot
ggplot(long_1, aes(x = Period, y = overall_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Overall Means of Each Genus Over the Years (*genera with at least 1% cover at a site in any year)",
       x = "Period",
       y = "Overall Mean Hard Coral Cover (%)") +
  facet_wrap(~Area, scales = "free_y") +
  theme_light()

# Define the desired order for Area
management_order <- c("NakorotubuFG", "ViRCP_new", "ViRCP_old", "KubulauFG", "NMR")

# Convert Area to a factor with custom order
long_1$Area <- factor(long_1$Area, levels = management_order)

# Define custom colors
custom_colors <- c(
  "Acropora" = "#33cccc",  # Blue
  "Astreopora" = "5F9EA0",
  "Diploastrea" = "#F0E68C",
  "Echinopora" = "violetred",
  "Favia" = "#9999ff",
  "Favites" = "#6633cc",
  "Goniastrea" = "#3399ff",
  "Isopora" = "#33ff99",
  "Merulina" = "#33a02c",
  'Millepora' = "#b2df8a",
  "Montipora" = "#1f78b4",
  "Other" = "plum1",  
  "Pavona" = "ivory4",
  "Pocillopora" = "#bcbddc",
  "Porites" = "salmon2"
)

ggplot(long_1, aes(x = Period, y = overall_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(
    #title = "Overall Means of Each Genus Over the Years (*genera with at least 1% cover at a site in any year)",
    x = "Time",
    y = "Overall Mean Hard Coral Cover (%)"
  ) +
  facet_wrap(~Area, ncol=5) +
  scale_fill_manual(values = custom_colors) +  
  theme_bw()

sum_overall_mean <- sum(long_1$overall_mean, na.rm = TRUE)
# Display the result
sum_overall_mean


#same but for those with >5% cover only
long_5 <- long_1 %>%
  group_by(Genus) %>%
  mutate(Genus = ifelse(all(overall_mean < 5), 'Other', as.character(Genus))) %>% #everything with <5% becomes 'other'
  ungroup()

sum_overall_mean <- sum(long_5$overall_mean, na.rm = TRUE)
# Display the result
sum_overall_mean

# Now create the plot
ggplot(long_5, aes(x = Period, y = overall_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(title = "Overall Means of Each Genus Over the Years (*genera with at least 5% cover at an area in any period)",
       x = "Period",
       y = "Overall Mean Hard Coral Cover (%)") +
  facet_wrap(~Area, scales = "free_y") +
  theme_light()

# Convert Area to a factor with custom order
long_5$Area <- factor(long_5$Area, levels = management_order)

# Define custom colors
custom_colors <- c(
  "Acropora" = "#33cccc",  
  "Montipora" = "#1f78b4",
  "Other" = "plum1",  
  "Pocillopora" = "#bcbddc",
  "Pavona" = "ivory4",
  "Porites" = "salmon2"
  # Add more colors as needed
)

# Plot the data
ggplot(long_5, aes(x = Period, y = overall_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Overall Means of Each Genus Over the Years (*genera with at least 5% cover at an area in any period)",
    x = "Period",
    y = "Overall Mean Hard Coral Cover (%)"
  ) +
  facet_wrap(~Area, scales = "free_y") +
  scale_fill_manual(values = custom_colors) +  
  theme_light()

#Add on the average wave exposure during winston as a panel
#the following file includes the routine/winston wave exposure
library(patchwork)

exposure <- read.csv("storm_impact.csv", sep = ",") # open to merge with management categories again
dat_exposure <- merge(dat4, exposure, by = "Site")

mean_by_period <- dat_exposure %>%
  group_by(Area) %>%
  summarize(
    mean_FN_durHs = mean(FN_durHs, na.rm = TRUE),
    mean_Win_v_Routine = mean(Win_v_Routine, na.rm = TRUE)
  )

mean_by_period$Area <- factor(mean_by_period$Area, levels = management_order)

font_size <- 12  

# Plot the data (Plot 5)
Plot_5 <- ggplot(long_5, aes(x = Period, y = overall_mean, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Period",
    y = "Overall Mean Hard Coral Cover (%)"
  ) +
  facet_wrap(~Area, nrow = 1) +  
  scale_fill_manual(values = custom_colors) +  
  theme_classic() +
  theme(
    title = element_text(size = font_size), 
    axis.text = element_text(size = 10),  
    axis.title = element_text(size = font_size),  
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    strip.text = element_text(size = font_size, face = "bold"),  
    strip.background = element_blank(),  
    plot.margin = margin(10, 10, 10, 10)  
  )

# Add the average wave exposure during Winston (Plot 2)
Plot_2 <- ggplot(mean_by_period, aes(x = "", y = mean_FN_durHs)) +  
  geom_bar(stat = "identity", fill = "white", color = "black", size = 0.5, width = 0.3) +  
  labs(
    y = "WaveDUR"
  ) +
  facet_wrap(~Area, nrow = 1) +  
  theme_classic() +
  theme(
    title = element_text(size = font_size),  
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),  
    axis.text = element_text(size = 10),  
    strip.text = element_blank(),  
    strip.background = element_blank(),  
    axis.title.x = element_blank(),  
    plot.margin = margin(10, 10, 10, 10)  
  )

# Combine the two plots into one figure using patchwork
library(patchwork)

combined_plot <- Plot_5 / Plot_2 +  
  plot_annotation(
    tag_levels = "A"  
  )

# Display the combined plot
combined_plot

#################################################################

#Add plot of PCA to growth form
dat_morphology <- read.csv("benthic_morphology.csv", sep = ";") 

#remove 2018 data
dat_morphology <- subset(dat_morphology, dat_morphology$Year!="2018")

dat_morphology <- dat_morphology %>%
  group_by(Site) %>%
  mutate(Period = case_when(
    # If both 2013 and 2014 data are present for a site, assign pre-cyclone to 2014 data
    Year == "2014" & any(Year == "2013") ~ "pre-cyclone",
    # Assign pre-cyclone to 2013 data if there is no 2014 data
    Year == "2013" & !any(Year == "2014") ~ "pre-cyclone",
    # Assign other periods based on the year
    Year == "2016" ~ "post-cyclone",
    Year == "2020" ~ "recovery-2020"
  )) %>%
  ungroup()

#Remove VIR5 as it a much shallower different reef system than all the others used in the study
dat_morphology$Site <- as.factor(dat_morphology$Site)
dat_morphology <- dat_morphology[dat_morphology$Site != "VIR5", ]
dat_morphology$Site <- droplevels(dat_morphology$Site)
levels(dat_morphology$Site)

# Drop rows with NA in Year
dat_morphology <- dat_morphology %>%
  filter(!is.na(Period))

dat_morphology$Period <- as.factor(dat_morphology$Period)

#load relevant packages
#library(vegan)
#library(ape)
library(dplyr)

#Can have issues with plyr here so if not working try this:
detach(package:ggbiplot)
detach(package:plyr)

#Take a mean of each site/year
str(dat_morphology)

means <- dat_morphology %>% 
  group_by(Site,Period, Management) %>% 
  summarize(mean(Arborescent), mean(Branching), mean(Columnar),
            mean(Corymbose), mean(Digitate), mean(Encrusting), mean(Foliose), mean(Massive), mean(Mushroom),
            mean(Plates), mean(Submassive))
means

means2 <- means %>%
  rename(Arborescent = "mean(Arborescent)", 
         Branching = "mean(Branching)", 
         Columnar = "mean(Columnar)", 
         Corymbose = "mean(Corymbose)", 
         Digitate = "mean(Digitate)", 
         Encrusting = "mean(Encrusting)", 
         Foliose = "mean(Foliose)", 
         Massive = "mean(Massive)", 
         Mushroom = "mean(Mushroom)", 
         Plates = "mean(Plates)",
         Submassive = "mean(Submassive)"
  )

pca.dat <- means2[,c(4:14)]
str(pca.dat)

pca<-prcomp(pca.dat,center=T,scale=T, na.rm=FALSE) # run the PCA : code is standardising data
summary(pca) # proportion of variance explained

screeplot(pca,main="Scree Plot", xlab="Component")
screeplot(pca,type="line",main="Scree Plot")

eigen<-pca$sdev^2 
eigen # retain the ones > 1

load<-pca$rotation
load # loadings of each component

## PCA plot
library(devtools) 
install_github("vqv/ggbiplot")
library(ggbiplot)

fPeriod<-factor(means$Period)
plot1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, ellipse = T, varname.size = 5, groups=fPeriod) +
  geom_point(aes(colour=fPeriod), size=4,  pch=16)
plot1
plot2 <- plot1 + theme_bw() + theme_classic() + 
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_blank()) + scale_colour_manual(values=c("#CCCCCC", "#FF3333", "#666666"))
plot2

mean(means2$Massive)

###################################################################

#and now summarising to management level

#Take a mean of each site/year
#Can have issues with plyr here so if not working try this:
detach(package:ggbiplot)
detach(package:plyr)

means <- dat_morphology %>% 
  group_by(Period, Management) %>% 
  summarize(mean(Arborescent), mean(Branching), mean(Columnar),
            mean(Corymbose), mean(Digitate), mean(Encrusting), mean(Foliose), mean(Massive), mean(Mushroom),
            mean(Plates), mean(Submassive))

means

means2 <- means %>%
  rename(Arborescent = "mean(Arborescent)", 
         Branching = "mean(Branching)", 
         Columnar = "mean(Columnar)", 
         Corymbose = "mean(Corymbose)", 
         Digitate = "mean(Digitate)", 
         Encrusting = "mean(Encrusting)", 
         Foliose = "mean(Foliose)", 
         Massive = "mean(Massive)", 
         Mushroom = "mean(Mushroom)", 
         Plates = "mean(Plates)",
         Submassive = "mean(Submassive)"
  )


pca.dat <- means2[,c(3:13)]
str(pca.dat)

pca<-prcomp(pca.dat,center=T,scale=T, na.rm=FALSE) # run the PCA : code is standardising data
summary(pca) # proportion of variance explained

screeplot(pca,main="Scree Plot", xlab="Component")
screeplot(pca,type="line",main="Scree Plot")

eigen<-pca$sdev^2 
eigen 

load<-pca$rotation
load 

## PCA plot
library(devtools) # don't forget to install Rtools first
install_github("vqv/ggbiplot")
library(ggbiplot)

fPeriod<-factor(means2$Period)
fManagement <-factor(means2$Management)

plot1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, ellipse = T, varname.size = 4, var.axes = TRUE, groups=fPeriod, linetype = "dashed") +
  geom_point(aes(colour=fManagement, pch=fPeriod), size=5)
plot1
plot2 <- plot1 + theme_bw() + theme_classic() + 
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.position = "top") + scale_colour_manual(values=c("#333333", "#666666", "#CCCCCC", "TOMATO", "TURQUOISE", "TOMATO4", "TURQUOISE3", "TURQUOISE4"))
plot2

#LOOK FOR COLLINEAR TERMS
pca.dat <- means2[,c(3:13)]
library(GGally)
ggpairs(pca.dat)

#Remove arborescent and columnar as almost all zero values

pca.dat <- means2[,c(4,6:13)]
str(pca.dat)

pca<-prcomp(pca.dat,center=T,scale=T, na.rm=FALSE) # run the PCA : code is standardising data
summary(pca) # proportion of variance explained

screeplot(pca,main="Scree Plot", xlab="Component")
screeplot(pca,type="line",main="Scree Plot")

eigen<-pca$sdev^2 
eigen 

load<-pca$rotation
load 

## PCA plot

fPeriod<-factor(means2$Period)
fManagement <-factor(means2$Management)

plot1 <- ggbiplot(pca, obs.scale = 1, var.scale = 1, ellipse = T, varname.size = 4, var.axes = TRUE, groups=fPeriod, linetype = "dashed") +
  geom_point(aes(colour=fManagement, pch=fPeriod), size=5)
plot1
plot2 <- plot1 + theme_bw() + theme_classic() + 
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 14), legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(), legend.position = "top") + scale_colour_manual(values=c("#333333", "#666666", "#CCCCCC", "TOMATO", "TURQUOISE", "TOMATO4", "TURQUOISE3", "TURQUOISE4"))
plot2

#Change so that recovery is with the triangle shape, make the arrows and text gray: 
plot_pca <- ggbiplot(
  pca,
  obs.scale = 1,
  var.scale = 1,
  ellipse = TRUE,
  varname.size = 4,
  var.axes = TRUE,
  groups = fPeriod,
  linetype = "dashed"
) +
  geom_point(aes(colour = fManagement, pch = fPeriod), size = 5) +
  geom_line(aes(group = interaction(fManagement, fPeriod), colour = fManagement), size = 1) +  # connector lines
  theme_bw() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  scale_colour_manual(values = c("#333333", "#666666", "#CCCCCC", "TOMATO", "TURQUOISE", "TOMATO4", "TURQUOISE3", "TURQUOISE4")) +
  scale_shape_manual(values = c("recovery-2020" = 17, "post-cyclone" = 16, "pre-cyclone" = 15)) +
  scale_fill_manual(values = c("#33333380", "#66666680", "#CCCCCC80")) 

# Modify the arrows and text colour for variables
plot_pca$layers[[1]]$aes_params$colour <- "grey70" # Change arrow colour
plot_pca$layers[[1]]$aes_params$alpha <- 0.8       # Adjust arrow transparency
plot_pca$layers[[3]]$aes_params$colour <- "grey70" # Change text colour on arrows
# Display the updated plot
plot_pca


plot_pca <- ggbiplot(
  pca,
  obs.scale = 1,
  var.scale = 1,
  ellipse = TRUE,
  varname.size = 4,
  var.axes = TRUE,
  groups = fPeriod,
  linetype = "dashed"
) +
  geom_point(aes(colour = fManagement, pch = fPeriod), size = 5) +
  geom_line(aes(group = interaction(fManagement, fPeriod), colour = fManagement), size = 1) +
  theme_bw() +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  ) +
  scale_colour_manual(values = c("#333333", "#666666", "#CCCCCC", "TOMATO", "TURQUOISE", "TOMATO4", "TURQUOISE3", "TURQUOISE4")) +
  scale_shape_manual(values = c("recovery-2020" = 17, "post-cyclone" = 16, "pre-cyclone" = 15)) +
  scale_fill_manual(values = c("#33333380", "#66666680", "#CCCCCC80"))  # 80 = 50% opacity in hex

# Modify the arrows and text colours
plot_pca$layers[[1]]$aes_params$colour <- "grey70"
plot_pca$layers[[1]]$aes_params$alpha <- 0.8
plot_pca$layers[[3]]$aes_params$colour <- "grey70"

# Display
plot_pca
