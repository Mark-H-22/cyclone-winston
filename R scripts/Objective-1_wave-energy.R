# Script for analysing objective 1 in Ford et al.
# Assessing the link between wave exposure and hard coral cover

library(dplyr)     
library(ggplot2)   
library(tidyr)             
library(reshape2) 
library(ggpubr)
library(gtsummary)

setwd("C:/Users/akfor/OneDrive/IN DEVELOPMENT/TC Winston Recovery/R-files/Github")
dat <- read.delim("benthic_full.txt")

head(dat)
str(dat)

#Make factors where necessary
dat <- dat %>%
  mutate_at(vars(Country, Site, Month, Day, Year, Exposure, Reef.slope, Reef.type, Reef.zone, Relative.depth, Transect.number, Area, Benthic.category), as.factor)

# Filter rows where 'Benthic.category' is 'Hard coral'
hard_coral_df <- dat %>% filter(Benthic.category == 'Hard coral')

#Convert dataframe to wide format so that growth forms have columns with counts
library(data.table)
library(stringr)

# Transform data to wide format by benthic category 
dat2 <- dcast(setDT(hard_coral_df), Year + Site + Area + Transect.number ~ hard_coral_df$Benthic.category, length)

# Set up data by Period
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

#Break up ViRCP into old and new (management)
dat4 <- dat4 %>%
  mutate(Area = as.factor(case_when(
    grepl("^VIR[1-4]|VIR10|VIR11", Site) ~ "ViRCP_old",
    grepl("^VIR[89]", Site) ~ "ViRCP_new",
    TRUE ~ as.character(Area)
  )))

dat4 %>%
  select(Period, Area, Site) %>%
  tbl_summary(by = Period) %>%
  bold_labels()

# Calculate the average HC per site and period
HC_summary <- dat4 %>%
  group_by(Site, Period, Area) %>%  # Group by Site, Period, and Area
  summarise(mean_hard_coral = mean(`Hard coral`, na.rm = TRUE), .groups = "drop") %>%
  arrange(Site, Period, Area)  
print(HC_summary)

HC_summary2 <- HC_summary %>%
  pivot_wider(names_from = Period, values_from = mean_hard_coral) %>%
  mutate(
    impact = `post-cyclone` - `pre-cyclone`, 
    recovery_2018 = `recovery-2018` - `post-cyclone`,
    recovery_2020 = `recovery-2020` - `post-cyclone`
  )
print(HC_summary2)

# Calculate relative changes
HC_summary2$impact_rel <- 100/HC_summary2$`pre-cyclone`*HC_summary2$impact
HC_summary2$recovery_2018_rel <- 100/HC_summary2$`recovery-2018`*HC_summary2$recovery_2018
HC_summary2$recovery_2020_rel <- 100/HC_summary2$`recovery-2020`*HC_summary2$recovery_2020

write.table(HC_summary2, file = "benthic_for_model.csv", sep = ";", dec = ".")


##################################################################################
# Script to test wave exposure on benthos

library(lmerTest)
library(ggrepel)
library(nlme)
library(lmtest)

dat <- read.csv("benthic_for_model.csv", sep = ";")

#the following file includes the routine/winston wave exposure
exposure <- read.csv("storm_impact.csv", sep = ",") # open to merge with management categories again
dat_exposure <- merge(dat, exposure, by = "Site")

str(dat_exposure)

hist(dat_exposure$impact) #normal
hist(dat_exposure$recovery_2018) #kind of normal
hist(dat_exposure$recovery_2020) #normal

plot(impact ~ FN_durHs, data = dat_exposure)

# Set theme and color palette
theme_custom <- theme_minimal() +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.text.y = element_text(angle = 90, vjust = 0.5),  # Rotate y-axis text to be horizontal
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    #plot.title = element_text(size = 16, hjust = 0.5, face = "bold", color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    panel.grid = element_line(color = "gray", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()  # Remove minor vertical lines
  )

dat_exposure$Area <- factor(dat_exposure$Area,
                                       levels = c("KubulauFG", "NMR", "NakorotubuFG", "ViRCP_new", "ViRCP_old"))

color_palette <- c("tomato", "tomato4", "#6defe2",  "#059398", "#014a4d")

# Create the base plot
impact <- ggplot(dat_exposure, aes(x = FN_durHs, y = impact, color = Area)) +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = color_palette) +  
  geom_hline(yintercept = c(-60, -40, -20, 0, 20), linetype = "dashed", color = "gray") +
  labs(
    x = "",
    y = "Absolute Hard Coral Loss (%) \n",
    color = "Area"  
  ) +
  theme_custom

plot(impact)

# Display the plot

# Add repelling labels without showing in the legend
impact_labeled <- impact +
  geom_text_repel(aes(label = Site), size = 4, box.padding = 0.5, show.legend = FALSE) +  # Customize x-axis labels
  #annotate("text", x = 0.08, y = -62, label = "Low Impact Expected", size = 4, color = "grey", vjust = 1) +
  #annotate("text", x = 0.22, y = -62, label = "High Impact Expected", size = 4, color = "grey", vjust = 1) +
  labs(
    x = "",
    y = "Absolute Hard Coral Loss (%) \n"
  )

# Display the plot
print(impact_labeled)

#####################################
#And for relative loss of coral cover - same script to plot

# Create the base plot
impact_rel <- ggplot(dat_exposure, aes(x = FN_durHs, y = impact_rel, color = Area)) +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = color_palette, breaks = c("KubulauFG", "NMR", "NakorotubuFG", "ViRCP_new", "ViRCP_old")) +  # Specify the order in the legend
  geom_hline(yintercept = c(-100, -75, -50, -25, 0, 25), linetype = "dashed", color = "gray") +
  labs(
    x = "",
    y = "Absolute Hard Coral Loss (%) \n",
    color = "Area"  
  ) +
  theme_custom

# Add repelling labels without showing in the legend
impact_rel_labeled <- impact_rel +
  geom_text_repel(aes(label = Site), size = 4, box.padding = 0.5, show.legend = FALSE) +
  
  # Customize x-axis labels
  #annotate("text", x = 0.08, y = -103, label = "Low Impact Expected", size = 4, color = "darkgrey", vjust = 0.5) +
  #annotate("text", x = 0.22, y = -103, label = "High Impact Expected", size = 4, color = "darkgrey", vjust = 0.5) +
  
  # Align y-axis labels with tick marks and lines
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25),
                     labels = c("-100", "-75", "-50", "-25", "0", "25")) +
  
  labs(
    x = "WaveDUR",
    y = "Relative Hard Coral Loss (%) \n"
  )

# Display the plot
print(impact_rel_labeled)


# Combine the two plots into one figure using patchwork
library(patchwork)

combined_plot <- impact_labeled / impact_rel_labeled + 
  plot_annotation(
    tag_levels = "a"  
  )

# Display the combined plot
combined_plot

#Export manually

# Run LMER to check for relationship between wave exposure and hard coral loss
library(lme4)

model <- lmer(impact ~ FN_durHs + (1 | Area), na.action = na.exclude, data = dat_exposure)
summary(model) 
plot(model)

plot(fitted(model), resid(model), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
qqnorm(resid(model))
qqline(resid(model), col = "red")

model_null <- lmer(impact ~ (1 | Area), na.action = na.exclude, data = dat_exposure)
summary(model_null) 

# Perform the likelihood ratio test
lr_test <- anova(model_null, model)
print(lr_test)

# Calculate R-squared values
library(MuMIn)
r_squared <- r.squaredGLMM(model)
print(r_squared)

hist(residuals(model))
qqnorm(resid(model)) 
qqline(resid(model)) 

model <- lmer(impact_rel ~ FN_durHs + (1 | Area), na.action = na.exclude, data = dat_exposure)
summary(model) 

plot(fitted(model), resid(model), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red") #pattern here - fanning out

hist(dat_exposure$impact_rel)

#invert data so it isn't negative
dat_exposure2 <- dat_exposure %>%
  filter(!is.na(impact_rel)) %>%
  mutate(impact_rel_sqrt = sqrt(impact_rel + abs(min(impact_rel, na.rm = TRUE)) + 1))

dat_exposure2$impact_rel_sqrt <- sqrt(dat_exposure2$impact_rel + abs(min(dat_exposure2$impact_rel)) + 1)

model <- lmer(impact_rel_sqrt ~ FN_durHs + (1 | Area), na.action = na.exclude, data = dat_exposure2)
summary(model)

plot(fitted(model), resid(model), 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
qqnorm(resid(model))
qqline(resid(model), col = "red")

model_null <- lmer(impact_rel_sqrt ~ (1 | Area), na.action = na.exclude, data = dat_exposure2)
summary(model_null) 

# Perform the likelihood ratio test
lr_test <- anova(model_null, model)
print(lr_test)

# Calculate R-squared values
library(MuMIn)
r_squared <- r.squaredGLMM(model)
print(r_squared)

#########################################################
