## Compare Rob Cooke's data to mine
## except volancy: I have not collected the data

library(dplyr)
library(ggplot2)
library(ggpubr)
library(EnvStats)
source("Functions_for_comparison_with_RC_data.R")

## Load data

# Robert Cooke's data, imputed and collected (need to be standardised within class and scaled to 0-mean and unit variance)
RC_Collected <- read.csv("../../Data/RCooke_Bird_and_Mammal_imputed_data/Cooke_et_al_2018_GEB_collected_trait_data.csv") %>%
  mutate(activity=ifelse(activity==1, "Other", ifelse(activity==2, "Nocturnal", NA)))
colnames(RC_Collected)[ c(3:6, 8)] <- c("Body_mass_g", "Litter_size", "Diel_activity", "Habitat_breadth_IUCN", "Primary_diet")

RC_Imputed <- read.csv("../../Data/RCooke_Bird_and_Mammal_imputed_data/Cooke_et_al_2018_GEB_single_imputed.csv") %>%
  mutate(activity=ifelse(activity==1, "Other", ifelse(activity==2, "Nocturnal", NA)))
colnames(RC_Imputed)[ c(3:6, 8)] <- c("Body_mass_g", "Litter_size", "Diel_activity", "Habitat_breadth_IUCN", "Primary_diet")

# My data, collected and imputed (one dataframe selected randomly) - not transformed and standardised
# Collected
X <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv") %>% 
  select(Class, Best_guess_binomial, Diel_activity, Body_mass_g, Litter_size, Habitat_breadth_IUCN, Primary_diet)

Y <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv") %>%
  select(Class, Best_guess_binomial, Diel_activity, Body_mass_g, Litter_size, Habitat_breadth_IUCN, Primary_diet)

AE_Collected <- rbind(X,Y) %>%  mutate(Best_guess_binomial=as.character(Best_guess_binomial))
AE_Collected <- AE_Collected[order(AE_Collected$Best_guess_binomial),]
rm(X,Y)

# Imputed 
# select one randomly
i <- sample(1:8,1)
AE_Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_not_standardised/List_of_8_sets.rds")[[i]]

X <- AE_Imputed$M$Imputed.Dataset %>%
  select(Class, Best_guess_binomial, Diel_activity, Body_mass_g, Litter_size, Habitat_breadth_IUCN, Primary_diet,Phylo_info)
Y <- AE_Imputed$B$Imputed.Dataset %>%
  select(Class, Best_guess_binomial, Diel_activity, Body_mass_g, Litter_size, Habitat_breadth_IUCN, Primary_diet, Phylo_info) 

AE_Imputed <- rbind(X,Y) %>%  mutate(Best_guess_binomial=as.character(Best_guess_binomial))
AE_Imputed <- AE_Imputed[order(AE_Imputed$Best_guess_binomial),] 
rm(X,Y)


# # Reprocess primary diet to align to Rob's categories
RC_Collected$Primary_diet <- as.character(RC_Collected$Primary_diet)
AE_Collected <- AE_Collected %>% Reprocess_diet()
RC_Imputed$Primary_diet <- as.character(RC_Imputed$Primary_diet)
AE_Imputed <- AE_Imputed %>% Reprocess_diet()

# Body mass in kg
RC_Collected$log10_Body_mass_kg <- log10(RC_Collected$Body_mass_g / 1000)
RC_Imputed$log10_Body_mass_kg <- log10(RC_Imputed$Body_mass_g / 1000)
AE_Collected$log10_Body_mass_kg <-log10(AE_Collected$Body_mass_g / 1000)
AE_Imputed$log10_Body_mass_kg <- log10(AE_Imputed$Body_mass_g / 1000)



## 1.a. Compare initial trait coverage
pdf(file="../../Results/Plots/Comparison_with_RCooke/Coverage.pdf", width=7, height=4.5, family="Times", pointsize=16)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(2,7,2,2), oma=c(1,1,1,1))
par(mfrow=c(2,2))
Plot.Cov(RC_Collected[RC_Collected$class=="Aves",], c("Body_mass_g","Litter_size", "Habitat_breadth_IUCN","Diel_activity", "Primary_diet"),
         Main="A. RC data, birds")
Plot.Cov(RC_Collected[RC_Collected$class=="Mammalia",], c("Body_mass_g","Litter_size", "Habitat_breadth_IUCN","Diel_activity", "Primary_diet"),
         Main="B. RC data, mammals")
Plot.Cov(AE_Collected[AE_Collected$Class=="Aves",], c("Body_mass_g","Litter_size", "Habitat_breadth_IUCN","Diel_activity", "Primary_diet"),
         Main="C. AE data, birds")
Plot.Cov(AE_Collected[AE_Collected$Class=="Mammalia",], c("Body_mass_g","Litter_size", "Habitat_breadth_IUCN","Diel_activity", "Primary_diet"),
         Main="D. AE data, mammals")
mtext(at=50, line=-7.5, "% coverage", cex=0.8)
mtext(at=-210, line=-7.5, "% coverage", cex=0.8)
dev.off()


## 1.b. Plot all traits against each other
p1 <- Plot_all_values(RC_Imputed, AE_Imputed, "log10_Body_mass_kg", FALSE, AxisX = "Body mass (kg, log10, AE)", AxisY = "Body mass (kg, log10, RC)", Diet = FALSE) + 
  labs(tag = "A") + theme(plot.tag.position = c(0.3, 0.92))
p2 <- Plot_all_values(RC_Imputed, AE_Imputed, "Litter_size", FALSE, AxisX = "Litter size (AE)", AxisY = "Litter size (RC)", Diet = FALSE)+ 
  labs(tag = "B") + theme(plot.tag.position = c(0.3, 0.92))
p3 <- Plot_all_values(RC_Imputed, AE_Imputed, "Habitat_breadth_IUCN", FALSE, AxisX = "Habitat breadth (AE)", AxisY = "Habitat breadth (RC)", Diet = FALSE)+ 
  labs(tag = "C") + theme(plot.tag.position = c(0.3, 0.92))
p4 <- Plot_all_values(RC_Imputed, AE_Imputed, "Diel_activity", TRUE, AxisX = "Diel activity", AxisY = "Diel activity", Diet = FALSE)+ 
  labs(tag = "D") + theme(plot.tag.position = c(0.3, 0.92))
p5 <- Plot_all_values(RC_Imputed, AE_Imputed, "Primary_diet", TRUE, AxisX = "Primary diet", AxisY = "Pirmary diet", Diet = TRUE)
p5p <- p5$p + labs(tag = "E") + theme(plot.tag.position = c(0.3, 0.92))
p <- ggarrange(p1,p2,p3,p4,p5p,common.legend = TRUE)
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_all_values.pdf", width =8, height =6)


## 2. Comparison of collected data for mammals and birds
p1 <- Compare(RC_Collected, NULL, AE_Collected, NULL, FALSE,"log10_Body_mass_kg", FALSE, "BM ","collected", "collected", FALSE, FALSE, FALSE) + 
  labs(tag = "B") + theme(plot.tag.position = c(0.2, 0.94))
p2 <- Compare(RC_Collected, NULL, AE_Collected, NULL, FALSE, "Litter_size", FALSE,"LCS ", "collected", "collected", FALSE, FALSE, FALSE) +
  labs(tag = "A") + theme(plot.tag.position = c(0.2, 0.94))
p3 <- Compare(RC_Collected, NULL, AE_Collected, NULL, FALSE,"Habitat_breadth_IUCN", FALSE, "HB ", "collected", "collected", FALSE, FALSE, FALSE) +
  labs(tag = "C") + theme(plot.tag.position = c(0.2, 0.94))
p4 <- Compare(RC_Collected, NULL,  AE_Collected, NULL, FALSE, "Diel_activity", TRUE, "DA", "collected", "collected", FALSE, FALSE, FALSE) +
  labs(tag = "D") + theme(plot.tag.position = c(0.22, 0.94))
X <- Compare(RC_Collected, NULL,  AE_Collected, NULL, FALSE, "Primary_diet", TRUE, "PD", "collected", "collected", FALSE, FALSE, TRUE)
p5 <- X$p + labs(tag = "E") + theme(plot.tag.position = c(0.22, 0.94))
p <- ggarrange(p2,p1,p3,p4,p5,common.legend = TRUE, legend="top")
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_collected.pdf", width =10, height =6)


## 3. Comparison of imputed data for mammals and birds
# p1 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, TRUE, "Body_mass_g", FALSE,"Body mass","imputed","imputed", FALSE, FALSE, FALSE)
p2 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, TRUE, "Litter_size", FALSE, "LCS ", "imputed", "imputed", FALSE, FALSE, FALSE) +
  labs(tag = "A") + theme(plot.tag.position = c(0.17, 0.94))
p3 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, TRUE, "Habitat_breadth_IUCN", FALSE, "HB ", "imputed", "imputed", FALSE, FALSE, FALSE) +
  labs(tag = "C") + theme(plot.tag.position = c(0.17, 0.94))
p4 <- Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,TRUE, "Diel_activity", TRUE, "DA", "imputed","imputed", FALSE, FALSE, FALSE) +
  labs(tag = "B") + theme(plot.tag.position = c(0.19, 0.94))
X <- Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed, TRUE, "Primary_diet", TRUE, "PD", "imputed","imputed", FALSE, FALSE, TRUE)
p5 <- X$p + labs(tag = "D") + theme(plot.tag.position = c(0.19, 0.94))
p <- ggarrange(p2,p4,p3,p5,common.legend = TRUE)
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_imputed.pdf", width =8, height =6)


## 4.Imputed VS collected, collected VS imputed
p1 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, FALSE, "log10_Body_mass_kg", FALSE, "BM, log10 ", "imputed", "imputed", TRUE, TRUE, FALSE) +
  labs(tag = "E") + theme(plot.tag.position = c(0.20, 0.94))

p2 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,  FALSE, "Litter_size", FALSE, "LCS ", "imputed", "imputed", TRUE, TRUE, FALSE) +
  labs(tag = "C") + theme(plot.tag.position = c(0.20, 0.94))
p2b <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, FALSE, "Litter_size", FALSE, "LCS ", "imputed", "imputed", TRUE, FALSE, FALSE)+
  labs(tag = "D") + theme(plot.tag.position = c(0.20, 0.94))

p3 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, FALSE, "Habitat_breadth_IUCN", FALSE, "HB ","imputed", "imputed", TRUE, TRUE, FALSE)+ labs(tag = "B") +
  theme(plot.tag.position = c(0.20, 0.94))
p3b <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, FALSE, "Habitat_breadth_IUCN", FALSE, "HB ","imputed", "imputed", TRUE, FALSE, FALSE)+ labs(tag = "A") + 
  theme(plot.tag.position = c(0.20, 0.94))

p4 <- ggarrange(
  Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,FALSE, "Diel_activity", TRUE, "DA", "imputed","imputed", TRUE, TRUE, FALSE)+ labs(tag = "F") + 
    theme(plot.tag.position = c(0.25, 0.94)),
  Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed, FALSE, "Diel_activity", TRUE, "DA", "imputed","imputed", TRUE, FALSE, FALSE)+ labs(tag = "G")
  + theme(plot.tag.position = c(0.25, 0.94))
)

p5 <- ggarrange(
  Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed, FALSE, "Primary_diet", TRUE, "PD", "imputed","imputed", TRUE, TRUE, TRUE)$p
  + theme(plot.tag.position = c(0.25, 0.94))+ labs(tag = "H"),
  Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed, FALSE, "Primary_diet", TRUE, "PD", "imputed","imputed", TRUE, FALSE, TRUE)$p + 
  theme(plot.tag.position = c(0.25, 0.94))+ labs(tag = "I")
)

p <- ggarrange(p3b,p3,p4,p2,p2b,p5,p1, common.legend = TRUE, widths = c(1/4,1/4, 1/2))
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_imputed_VS_collected.pdf", width =13, height =9)

