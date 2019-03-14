## Imputation errors across 8 imputed datasets for each class

library(reshape)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape)
source("Functions_for_results_congruence.R")

Cattraits <- c("Trophic_level", "Diel_activity", "Specialisation", "Primary_diet")

## Functions to get and to plot imputation errors for categorical imputed traits: PFC
GetPFC <- function(Errors){
  
  Cattraits <- c("Trophic_level", "Diel_activity", "Specialisation", "Primary_diet")
  Y <- which(colnames(Errors) %in% paste(Cattraits, "PFC", sep=" "))
  PFC <- Errors[,Y]*100
  return(PFC)
}

PlotPFC <- function(PFCErrors) {
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  PFC <- melt(PFCErrors)
  
  PFC <- PFC %>% group_by(Class, variable) %>%
    summarise(Mean=mean(value, na.rm=TRUE), 
              Max=max(value, na.rm = TRUE),
              Min=min(value, na.rm = TRUE))
  PFC$Max[15] <- NA
  PFC$Min[15] <- NA
  
  p <- ggplot(PFC, aes(variable, Mean, col=Class)) + GGPoptions +
    geom_point() + 
    scale_x_discrete(labels=c("DA", "TL", "PD", "Sp")) +
    geom_segment(aes(x=variable, xend=variable, y=PFC$Min, yend=PFC$Max, col=Class)) +
    ylim(0,20) +
    xlab("") + ylab("OOB % falsely classified") + coord_flip()
  
  return(p)
}

## Functions to get and to plot imputation errors for continuous traits: MSE
GetrootMSE <- function(Errors){
  
  Contraits <- c("Body_mass_g", "Longevity_d", "Litter_size", "Diet_breadth","Habitat_breadth_IUCN","Range_size_m2")

  Y <- which(colnames(Errors) %in% paste(Contraits, "MSE", sep=" "))
  Y <- Errors[,Y]
  
  Y <- sqrt(Y)
  
  Y[, "Body_mass_g MSE"] <- Y[, "Body_mass_g MSE"] / 1000  # in kg
  Y[, "Longevity_d MSE"] <- Y[, "Longevity_d MSE"] / 365.25 # in years
  Y[, "Range_size_m2 MSE"] <- Y[, "Range_size_m2 MSE"] / 1000000000000 # in Mega-meters
  
  return(Y)
}

PlotMSE <- function(MSEErrors) {
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  MSE <- melt(MSEErrors)
  
  MSE <- MSE %>% group_by(Class, variable) %>%
    summarise(Mean=mean(value, na.rm=TRUE), 
              Max=max(value, na.rm = TRUE),
              Min=min(value, na.rm = TRUE))
  MSE$Max[23] <- NA
  MSE$Min[23] <- NA

  p <- ggplot(MSE, aes(variable, Mean, col=Class)) + GGPoptions +
    geom_point() + 
    scale_x_discrete(labels=c("BM (kg)", "LG (years)", "LCS", expression("RS (mega-m"^2*")"), "DB", "HB")) +
    scale_y_continuous(labels=c(0.1,1,10,100,1000), breaks=c(0.1,1,10,100,1000), trans="log10") +
    geom_segment(aes(x=variable, xend=variable, y=MSE$Min, yend=MSE$Max, col=Class)) +
    xlab("") + ylab(expression(sqrt("OOB mean-square error"))) + coord_flip()
  
  return(p)
}


## Load imputed datasets (8); retrieve errors
Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_corrected_trees/List_of_8_sets.rds")

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
PFC_amp <- GetPFC(Errors_amp); PFC_amp$Class <- "Amphibians"
PFC_rep <- GetPFC(Errors_rep); PFC_rep$`Primary_diet PFC` <- NA;  PFC_rep$Class <- "Reptiles"
PFC_mam <- GetPFC(Errors_mam); PFC_mam$Class <- "Mammals"
PFC_bir <- GetPFC(Errors_bir); PFC_bir$Class <- "Birds"
PFC_all <- rbind(PFC_amp, PFC_rep, PFC_mam, PFC_bir)

## MSE for all classes
MSE_amp <- GetrootMSE(Errors_amp); MSE_amp$Class <- "Amphibians"
MSE_rep <- GetrootMSE(Errors_rep);  MSE_rep$`Diet_breadth MSE` <- NA;  MSE_rep$Class <- "Reptiles"
MSE_mam <- GetrootMSE(Errors_mam);  MSE_mam$Class <- "Mammals"
MSE_bir <- GetrootMSE(Errors_bir); MSE_bir$Class <- "Birds"
MSE_all <- rbind(MSE_amp, MSE_rep, MSE_mam, MSE_bir)


## Plot PFC for all classes
pPFC <- PlotPFC(PFC_all)
## Plot MSE for all classes
pMSE <- PlotMSE(MSE_all)

p <- ggarrange(pMSE + labs(tag = "A") + theme(plot.tag.position = c(0.95, 0.94)),
               pPFC + labs(tag = "B") + theme(plot.tag.position = c(0.95, 0.94)),
               common.legend = TRUE, widths=c(0.555,0.445), legend="right")
ggsave(p, file="../../Results/Plots/Imputation_errors/MSE_PFC.pdf", width=8.5, height=3)



