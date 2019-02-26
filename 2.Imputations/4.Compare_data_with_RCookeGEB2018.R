## Compare Rob Cooke's data to mine
# diet?
# volancy: I have not collected the data

library(dplyr)
library(ggplot2)
library(ggpubr)
source("Functions_for_comparison_with_RC_data.R")

## Load data

# Robert Cooke's data, imputed and collected (need to be standardised within class and scaled to 0-mean and unit variance)
RC_Collected <- read.csv("../../Data/RCooke_Bird_and_Mammal_imputed_data/Cooke_et_al_2018_GEB_collected_trait_data.csv") %>%
  mutate(activity=ifelse(activity==1, "Other", ifelse(activity==2, "Nocturnal", NA)))

RC_Imputed <- read.csv("../../Data/RCooke_Bird_and_Mammal_imputed_data/Cooke_et_al_2018_GEB_single_imputed.csv") %>%
  mutate(activity=ifelse(activity==1, "Other", ifelse(activity==2, "Nocturnal", NA)))

# My data, collected and imputed (one dataframe selected randomly)

# Collected
X <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv") %>% 
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN, Primary_diet) %>%
  mutate(Class="Mammalia") 

Y <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv") %>%
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN, Primary_diet) %>%
  mutate(Class="Aves") 

AE_Collected <- rbind(X,Y) %>%  mutate(Best_guess_binomial=as.character(Best_guess_binomial))
AE_Collected <- AE_Collected[order(AE_Collected$Best_guess_binomial),]
rm(X,Y)

# Imputed 
# select one randomly
i <- sample(1:8,1)
AE_Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_datasets/List_8_imputed_sets_v2.rds")[[i]]

X <- AE_Imputed$M$Imputed.Dataset %>%
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN, Primary_diet,Phylo_info) %>%
  mutate(Class="Mammalia")

Y <- AE_Imputed$B$Imputed.Dataset %>%
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN, Primary_diet, Phylo_info) %>%
  mutate(Class="Aves")

AE_Imputed <- rbind(X,Y) %>%  mutate(Best_guess_binomial=as.character(Best_guess_binomial))
AE_Imputed <- AE_Imputed[order(AE_Imputed$Best_guess_binomial),] 
rm(X,Y)


# Reprocess primary diet to align to Rob's categories
AE_Imputed <- Reprocess_diet(AE_Imputed)
AE_Collected <- Reprocess_diet(AE_Collected)

RC_Collected$diet_5cat <- as.character(RC_Collected$diet_5cat)
AE_Collected$Reprocessed_PD <- as.character(AE_Collected$Reprocessed_PD)
RC_Imputed$diet_5cat <- as.character(RC_Imputed$diet_5cat)
AE_Imputed$Reprocessed_PD <- as.character(AE_Imputed$Reprocessed_PD)


## 1. Compare initial trait coverage
pdf(file="../../Results/Plots/Comparison_with_RCooke/Coverage.pdf", width=7, height=4.5, family="Times", pointsize=16)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(2,7,2,2), oma=c(1,1,1,1))
par(mfrow=c(2,2))

Plot.Cov(RC_Collected[RC_Collected$class=="Aves",], c("body_mass_median","litter_clutch_size", "hab_breadth","activity"),
         Main="R.C., birds")
Plot.Cov(RC_Collected[RC_Collected$class=="Mammalia",], c("body_mass_median","litter_clutch_size", "hab_breadth","activity"),
         Main="R.C., mammals")
Plot.Cov(AE_Collected[AE_Collected$Class=="Aves",], c("log10_Body_mass_g","log10_Litter_size", "sqrt_Habitat_breadth_IUCN","Diel_activity"),
         Main="A.E., birds")
Plot.Cov(AE_Collected[AE_Collected$Class=="Mammalia",], c("log10_Body_mass_g","log10_Litter_size", "sqrt_Habitat_breadth_IUCN","Diel_activity"),
         Main="A.E., mammals")
mtext(at=50, line=-7.5, "% coverage", cex=0.8)
mtext(at=-210, line=-7.5, "% coverage", cex=0.8)
dev.off()

## 2. Comparison of collected data for mammals and birds
p1 <- Compare(RC_Collected, NULL, AE_Collected, NULL, "body_mass_median","log10_Body_mass_g", FALSE,"BM","collected", "collected", FALSE, FALSE, FALSE, FALSE)  + labs(tag = "B") + theme(plot.tag.position = c(0.17, 0.94))
p2 <- Compare(RC_Collected, NULL, AE_Collected, NULL,  "litter_clutch_size","log10_Litter_size", FALSE,"LCS", "collected", "collected",FALSE,FALSE, FALSE,FALSE)+ labs(tag = "A") + theme(plot.tag.position = c(0.17, 0.94))
p3 <- Compare(RC_Collected, NULL, AE_Collected, NULL,"hab_breadth","sqrt_Habitat_breadth_IUCN", FALSE, NULL, "collected", "collected", FALSE,FALSE, FALSE,FALSE)+ labs(tag = "C") + theme(plot.tag.position = c(0.17, 0.94))
p4 <- Compare(RC_Collected, NULL,  AE_Collected, NULL,  "activity","Diel_activity", TRUE, "DA", "collected", "collected", FALSE,FALSE, FALSE,FALSE)+ labs(tag = "D") + theme(plot.tag.position = c(0.17, 0.94))
X <- Compare(RC_Collected, NULL,  AE_Collected, NULL,  "diet_5cat","Reprocessed_PD", TRUE, "PD", "collected", "collected", FALSE,FALSE, FALSE, TRUE)
p5 <- X$p + labs(tag = "E") + theme(plot.tag.position = c(0.17, 0.94))
p <- ggarrange(p2,p1,p3,p4,p5,common.legend = TRUE, legend="top")
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_collected.pdf", width =8, height =6)

## 3. Comparison of imputed data for mammals and birds
#p1 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, "body_mass_median","log10_Body_mass_g", FALSE,"Body mass","imputed","imputed", TRUE, FALSE, FALSE)
p2 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,  "litter_clutch_size","log10_Litter_size", FALSE,"LCS", "imputed", "imputed",TRUE,FALSE, FALSE, FALSE)+ labs(tag = "A") + theme(plot.tag.position = c(0.17, 0.94))
p3 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,"hab_breadth","sqrt_Habitat_breadth_IUCN", FALSE, NULL, "imputed", "imputed", TRUE, FALSE, FALSE, FALSE)+ labs(tag = "C") + theme(plot.tag.position = c(0.17, 0.94))
p4 <- Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,  "activity","Diel_activity", TRUE, "DA", "imputed","imputed", TRUE, FALSE, FALSE, FALSE)+ labs(tag = "B") + theme(plot.tag.position = c(0.17, 0.94))
X <- Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,  "diet_5cat","Reprocessed_PD", TRUE, "PD", "imputed","imputed", TRUE, FALSE, FALSE, TRUE)
p5 <- X$p + labs(tag = "D") + theme(plot.tag.position = c(0.95, 0.94))
p <- ggarrange(p2,p4,p3,p5,common.legend = TRUE)
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_imputed.pdf", width =8, height =6)


## 4.Imputed VS collected, collected VS imputed

# Collected RC vs Imputed AE
p1 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, "body_mass_median","log10_Body_mass_g", FALSE,"BM","imputed","imputed", FALSE, TRUE, TRUE) + labs(tag = "E") + theme(plot.tag.position = c(0.10, 0.94))
# p1 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, "body_mass_median","log10_Body_mass_g", FALSE,"Body mass","imputed","imputed", FALSE, TRUE, FALSE)

p2 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,  "litter_clutch_size","log10_Litter_size", FALSE,"LCS", "imputed", "imputed", FALSE, TRUE, FALSE) + labs(tag = "A") + theme(plot.tag.position = c(0.10, 0.94))
p2b <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,  "litter_clutch_size","log10_Litter_size", FALSE,"LCS", "imputed", "imputed", FALSE, TRUE, TRUE) + labs(tag = "C") + theme(plot.tag.position = c(0.10, 0.94))

p3 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,"hab_breadth","sqrt_Habitat_breadth_IUCN", FALSE, NULL, "imputed", "imputed", FALSE, TRUE, FALSE)+ labs(tag = "B") + theme(plot.tag.position = c(0.10, 0.94))
p3b <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,"hab_breadth","sqrt_Habitat_breadth_IUCN", FALSE, NULL, "imputed", "imputed", FALSE, TRUE, TRUE)+ labs(tag = "D") + theme(plot.tag.position = c(0.10, 0.94))

p4 <- ggarrange(
  Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,  "activity","Diel_activity", TRUE, "DA", "imputed","imputed", FALSE, TRUE, FALSE)+ labs(tag = "F") + theme(plot.tag.position = c(0.2, 0.94)),
  Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,  "activity","Diel_activity", TRUE, "DA", "imputed","imputed", FALSE, TRUE, TRUE)+ labs(tag = "G") + theme(plot.tag.position = c(0.2, 0.94)),
  common.legend = TRUE
)

X <- Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,  
              "diet_5cat","Reprocessed_PD", TRUE, "PD", "imputed","imputed", 
              FALSE, TRUE, FALSE, TRUE)
Xt <- X$outputs
p5 <- X$p + theme(plot.tag.position = c(0.95, 0.94))+ labs(tag = "H")

p <- ggarrange(p2,p3,p2b,p3b,p1,p4,p5,common.legend = TRUE, ncol=2, nrow=3)
ggsave(p, file="../../Results/Plots/Comparison_with_RCooke/Comparison_imputed_VS_collected.pdf", width =13, height =11)


