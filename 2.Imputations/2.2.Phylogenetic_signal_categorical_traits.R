## Phylogenetic signal for categorical traits: testing the difference from null distribution of delta values
library(gmodels)
library(ggplot2)

Mammals <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalMammals.rds")

## Function to plot the distribution of null values, calculate mean and 95% CI for all categorcial traits
Comparison_null_dist <- function(List, Catraits) {
  
  # Dataset of observed signal
  Delta <- vector()
  for (i in 1:length(Catraits)) { 
    Delta <- c(Delta, List[[Catraits[i]]][1]) %>%
      as.data.frame() %>%
      t 
  }
  rownames(Delta) <- Catraits
  colnames(Delta) <- "delta_obs"
  Delta <- as.data.frame(Delta)
  
  # Add simulation means and 95% CI
  n <- List[[Catraits[1]]][2]$Delta0 %>% nrow
  Simulations <- as.data.frame(matrix(nrow=n, ncol=length(Catraits)))
  colnames(Simulations) <- Catraits
  
  for (i in 1:length(Catraits)) {
    Sim <- List[[Catraits[i]]][2]$Delta0 %>% as.vector()
    colnames(Sim) <- "Simulated"
    Delta$MeanSim[i] <- mean(Sim$Simulated, na.rm=TRUE)
    Delta$CI_up[i] <- ci(Sim$Simulated, confidence=0.95, na.rm = TRUE)[3]
    Delta$CI_low[i]<- ci(Sim$Simulated, confidence=0.95, na.rm = TRUE)[2]
    Simulations[,i] <- Sim
  }
  
  return(list(Delta, Simulations))
  
}

## Function to plot results
Plot_Delta <- function(Delta) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=13), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=13),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  Delta$Trait <- rownames(Delta)
  
  p <- ggplot(Delta, aes(Trait, delta_obs)) + 
    geom_point() + ylim(min(Delta$CI_low), max(Delta$delta_obs)+1) +
    geom_point(aes(Trait, MeanSim), col="blue") +
    geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=0.3, col="blue") +
    scale_y_continuous(trans="log10") +
    GGPoptions
  
  return(p)
}

# Mammals
Results_Mammals <- Comparison_null_dist(Mammals, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaMammal <- Results_Mammals[[1]]
SimMammal <- Results_Mammals[[2]]
pMammals <- Plot_Delta(DeltaMammal) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
ggsave(pMammals, file="../../Results/Plots/Phylosignal_categorical/Mammals.pdf", height=2, width=3)
plot(density(SimMammal$Specialisation, na.rm=TRUE))
plot(density(SimMammal$Trophic_level, na.rm=TRUE))
plot(density(SimMammal$Diel_activity, na.rm=TRUE))





