# Analysis of Lupinus argenteus pollen and petal alkaloids.
# This data contains the values of alkaloid concentrations for the quinolizidine alkaloids in petals and pollen.
# Lupinus argenteus individuals were sampled at and near the RMBL in 2021 and 2022 in three sites: 401, GT, and EL. 
# Sébastien Rivest, 16/02/2023

# Load packages
library(readxl)
library(vegan)
library(ggplot2)
library(tidyverse)

# Loading dataset
poll_pet_alk <- read_excel("C:/Users/seb69/OneDrive/Documents/Doc/lupinus_figs/Alkaloid_composition.xlsx")
poll_pet_alk <- subset(poll_pet_alk, is.na(Thermopsine)==FALSE) # remove missing rows with missing data
# Keep only columns with information on alkaloids for the alkaloid composition matrix.
data_lup <- subset(poll_pet_alk,rowSums(poll_pet_alk[,6:19])>0)

# The GT population belongs to a different chemotype than 401 and EL. 
# The two chemotypes are analysed separately.
# Start with the 401 and EL chemotype ("upvalley").
data_lup_upvalley <- data_lup %>% subset(., Site == "401" | Site == "EL") 
# Keep only alkaloids present in this chemotype
data_lup_upvalley <- data_lup_upvalley[,1:10]
data_lup_upvalley <- subset(data_lup_upvalley, rowSums(data_lup_upvalley[,6:10])>0)
# Nonmetric multidimensional scaling
ord1 <- metaMDS(data_lup_upvalley[,6:10], distance="bray", k=1)
ord2 <- metaMDS(data_lup_upvalley[,6:10], distance="bray", k=2)
ord3 <- metaMDS(data_lup_upvalley[,6:10], distance="bray", k=3)
ord4 <- metaMDS(data_lup_upvalley[,6:10], distance="bray", k=4)
ord5 <- metaMDS(data_lup_upvalley[,6:10], distance="bray", k=5)
# Stress plot
stress <- c(ord1$stress,ord2$stress,ord3$stress,ord4$stress,ord5$stress)
plot(stress~c(1:5)) 

layout(matrix(seq(1,2),nrow=2))
par(mar=c(5,5,1,1))
stressplot(ord3)
text(0.3,1,"k = 3 axes")
stressplot(ord2)
text(0.3,1,"k = 2 axes")

dev.off()
# Fit explanatory variables in the ordination
site.dat <- data_lup_upvalley[, 1:5] %>% mutate(Year = as.factor(Year))
envfit(ord2,site.dat[, c(1,2,5)], permu=999)

NMDS1 <- ord2$points[,1] 
NMDS2 <- ord2$points[,2]
data.scores <- cbind(NMDS1, NMDS2, site.dat)

# Plot
colors <- c(rep("#7442c8",length(which(data_lup_upvalley$Tissue == "Petal")) ), 
          rep("#ffcf48",length(which(data_lup_upvalley$Tissue == "Pollen"))))
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Tissue)) + 
  geom_point(data = data.scores, aes(fill = Tissue),shape = 21, col = "black", size = 2, alpha = 0.7) + 
  scale_fill_manual(values = c("steelblue", "#bd8200")) + 
  stat_ellipse(type='t',size =1, linetype = 2) +
  scale_colour_manual(values = c("steelblue", "#bd8200")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
        legend.key = element_blank()) +
  labs(colour = "Tissue")
gg

# Analyse of the GT site chemotype.
data_lup_GT <- data_lup %>% subset(., Site == "GT") 

ord1 <- metaMDS(data_lup_GT[,11:19], distance="bray", k=1)
ord2 <- metaMDS(data_lup_GT[,11:19], distance="bray", k=2)
ord3 <- metaMDS(data_lup_GT[,11:19], distance="bray", k=3)
ord4 <- metaMDS(data_lup_GT[,11:19], distance="bray", k=4)
ord5 <- metaMDS(data_lup_GT[,11:19], distance="bray", k=5)

stress <- c(ord1$stress,ord2$stress,ord3$stress,ord4$stress,ord5$stress)
plot(stress~c(1:5)) 

layout(matrix(seq(1,2),nrow=2))
par(mar=c(5,5,1,1))
stressplot(ord3)
text(0.3,1,"k = 3 axes")
stressplot(ord2)
text(0.3,1,"k = 2 axes")

dev.off()

site.dat<-data_lup_GT[, 1:5] %>% mutate(Year = as.factor(Year))
envfit(ord2,site.dat[, c(1, 3, 5)],permu=999)

NMDS1 <- ord2$points[,1] 
NMDS2 <- ord2$points[,2]
data.scores <- cbind(NMDS1, NMDS2, site.dat)

colors<-c(rep("#7442c8", length(which(data_lup_upvalley$Tissue == "Petal")) ), 
          rep("#ffcf48", length(which(data_lup_upvalley$Tissue == "Pollen"))))

gg2 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Tissue)) + 
  geom_point(data = data.scores, aes(fill = Tissue),shape = 21, col = "black", size = 2, alpha = 0.7) + 
  scale_fill_manual(values = c("steelblue", "#bd8200")) + 
  stat_ellipse(type='t', size =1, linetype = 2) +
  scale_colour_manual(values = c("steelblue", "#bd8200")) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
        legend.key = element_blank()) +
  labs(colour = "Tissue")
gg2
