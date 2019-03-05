## Imputation errors across 8 imputed datasets for each class

library(reshape)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape)
source("Functions_for_results_congruence.R")

Cattraits <- c("Trophic_level", "Diel_activity", "Specialisation", "Primary_diet")

## Load imputed datasets (8); retrieve errors
Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_not_standardised/List_of_8_sets.rds")

## Amphibians
Amphibians <- Get_all_results(Imputed, "A")$Errors
Errors_amp <- data.table::rbindlist(Amphibians) %>% as.data.frame()

# Categorical
Y <- which(colnames(Errors_amp) %in% paste(Cattraits, "PFC", sep=" "))
PFC_amp <- Errors_amp[,Y]*100
PFC_amp <- melt(PFC_amp)
PFC_amp_means <- PFC_amp %>% group_by(variable) %>% summarise(Mean=mean(value))

ggplot(PFC_amp, aes(variable, value)) + geom_point() + ylim(0,100) +
  geom_point(aes(PFC_amp_means$variable, PFC_amp_means$Mean))


## Reptiles
Reptiles <- Get_all_results(Imputed, "R")$Errors
Errors_rep <- data.table::rbindlist(Reptiles)

## Mammals
Mammals <- Get_all_results(Imputed, "M")$Errors
Errors_mam <- data.table::rbindlist(Mammals)

## Birds
Birds <- Get_all_results(Imputed, "B")$Errors
Errors_bir <- data.table::rbindlist(Birds)


