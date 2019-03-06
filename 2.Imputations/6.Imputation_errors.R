## Imputation errors across 8 imputed datasets for each class

library(reshape)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape)
source("Functions_for_results_congruence.R")

Cattraits <- c("Trophic_level", "Diel_activity", "Specialisation", "Primary_diet")

## Function to plot imputation errors for categorical imputed traits: PFC
PlotPFC <- function(PFCErrors) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))

  PFC <- melt(PFCErrors)
  p <- ggplot(PFC, aes(variable, value, col=Class)) + GGPoptions +
    geom_point(alpha=0.5) + 
    ylim(0,25) +
    scale_y_continuous(trans="log10", breaks=c(0.01,0.10,1.00,10.00,20.00)) +
    scale_x_discrete(labels=c("DA", "TL", "PD", "Sp")) +
    xlab("") + ylab("% falsely classified")
  
  return(p)
}


GetPFC <- function(Errors, isReptile){

  if(!isReptile){Cattraits <- c("Trophic_level", "Diel_activity", "Specialisation", "Primary_diet")}
  else{Cattraits <- c("Trophic_level", "Diel_activity", "Specialisation")}
  
  Y <- which(colnames(Errors) %in% paste(Cattraits, "PFC", sep=" "))
  PFC <- Errors[,Y]*100
    return(PFC)
}

## Load imputed datasets (8); retrieve errors
Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_not_standardised/List_of_8_sets.rds")

## Amphibians
Amphibians <- Get_all_results(Imputed, "A")$Errors
Errors_amp <- data.table::rbindlist(Amphibians)%>% as.data.frame()

## Reptiles
Reptiles <- Get_all_results(Imputed, "R")$Errors
Errors_rep <- data.table::rbindlist(Reptiles)%>% as.data.frame()

## Mammals
Mammals <- Get_all_results(Imputed, "M")$Errors
Errors_mam <- data.table::rbindlist(Mammals)%>% as.data.frame()

## Birds
Birds <- Get_all_results(Imputed, "B")$Errors
Errors_bir <- data.table::rbindlist(Birds)%>% as.data.frame()


## PFC for all classes
PFC_amp <- GetPFC(Errors_amp, FALSE); PFC_amp$Class <- "Amphibians"
PFC_rep <- GetPFC(Errors_rep, FALSE); PFC_rep$`Primary_diet PFC` <- NA;  PFC_rep$Class <- "Reptiles"
PFC_mam <- GetPFC(Errors_mam, FALSE); PFC_mam$Class <- "Mammals"
PFC_bir <- GetPFC(Errors_bir, FALSE); PFC_bir$Class <- "Birds"
PFC_all <- rbind(PFC_amp, PFC_rep, PFC_mam, PFC_bir)

## Plot PFC for all classes
p <- PlotPFC(PFC_all)
ggsave(p, file="../../Results/Plots/Imputation_errors/PFC.pdf", width=4, height=2.5)
p$data


