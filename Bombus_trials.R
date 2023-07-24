# Analysis of Bombus appositus and B. rufocinctus behavior in response to thermopsine.
# The number of visits to an artificial flower and the time spent on the flower were recorded...
# during 5 minutes in flight cages for trials with or without thermopsine added to pollen.
# Sébastien Rivest, 16/02/2023

# Load packages
Packages <- c("readxl", "lme4", "lmerTest", "blmeco", "DHARMa",
              "ggplot2", "ggthemes", "glmmTMB", "dplyr")
lapply(Packages, library, character.only = TRUE)

# Import dataset
Bombus_trials <- read_excel("C:/Users/seb69/OneDrive/Documents/Doc/lupinus_figs/Bombus_trials.xlsx")
Bombus_trials <- mutate(Bombus_trials, Colony = as.factor(Colony))

# Model for the Number of visits to the artificial flower
hist(Bombus_trials$Number_visits)
visit_mod1 <- glmer(Number_visits ~ Species + Treatment + (1|Colony),
                           data = Bombus_trials, poisson(link = "log")) 
# Test for overdispersion
dispersion_glmer(visit_mod1) 
r_visits_mod1 <- simulateResiduals(visit_mod1)
plot(r_visits_mod1)
summary(visit_mod1)

# Model for the time (in seconds) spent on the artificial flower
hist(Bombus_trials$Time_on_flower)
visit_mod2 <- glmmTMB(Time_on_flower ~ Species + Treatment + (1|Colony),
                     data = Bombus_trials, Gamma(link = "log")) 
r_visits_mod2 <- simulateResiduals(visit_mod2)
plot(r_visits_mod2)
summary(visit_mod2)

# Plot for the number of visits to the artificial flower
p1 <- ggplot(Bombus_trials, aes(x=Treatment, y=Number_visits, fill=Treatment)) + 
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  geom_jitter(shape=16, position=position_jitter(0.15), alpha = 0.4) +
  scale_fill_manual(values = c("#ca2e3090", "#f69f5690")) +
  theme_classic() +
  labs(x="Treatment", y = "Number of visits") +
  theme(legend.position = "none")
p1

# Plot for the time spent on the artificial flower
p2 <- ggplot(Bombus_trials, aes(x=Treatment, y=Time_on_flower, fill=Treatment)) + 
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  geom_jitter(shape=16, position=position_jitter(0.15), alpha = 0.4) +
  scale_fill_manual(values = c("#ca2e3090", "#f69f5690")) +
  theme_classic() +
  labs(x="Treatment", y = "Time on flower") +
  theme(legend.position = "none")
p2
