## Phylogenetic signal for categorical traits: testing the difference from null distribution of delta values
library(gmodels)
library(ggplot2)
library(dplyr)
library(boot)
library(ggpubr)

## Function to plot the distribution of null values, calculate mean and 95% CI for all categorical traits
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
    Delta$MedianSim[i] <- median(Sim$Simulated, na.rm=TRUE)
    
    # 95% ci for the median (bootstrapping)
    bootmed <- apply(matrix(sample(Sim$Simulated, rep=TRUE, 10^5*length(Sim$Simulated)), nrow=10^5), 1, median, na.rm=TRUE)
    X <- quantile(bootmed, c(.025, 0.975), na.rm=TRUE)
    Delta$CI_up[i] <- X[2]
    Delta$CI_low[i]<- X[1]
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
    geom_point(aes(Trait, MedianSim), col="blue") +
    geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=0.3, col="blue") +
    scale_y_continuous(trans="log10") +
    GGPoptions
  
  return(p)
}

## Function to run Wilcoxon Rank Sum and Signed rank tests
## Testing whether the observed median is significantly different from the null distribution
Wilcox_Test <- function(Simulations, Delta) {

  Delta <- Delta %>% 
    select(delta_obs)
  
  for (i in 1:nrow(Delta)) {
    
    Delta$W_pvalue[i] <- wilcox.test(Simulations[,i], mu=Delta$delta_obs[i], alternative="less")["p.value"] %>% as.numeric()
    Delta$W_statistic[i] <- wilcox.test(Simulations[,i], mu=Delta$delta_obs[i], alternative="less")["statistic"]%>% as.numeric()
    
  }
  
  return(Delta)
}



## Load data
Mammals <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalMammals.rds")
Birds <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalBirds.rds")
Reptiles <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalReptiles.rds")
Amphibians <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalAmphibians.rds")


# Mammals
Results_Mammals <- Comparison_null_dist(Mammals, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaMammal <- Results_Mammals[[1]]
SimMammal <- Results_Mammals[[2]]
pMammals <- Plot_Delta(DeltaMammal) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimMammal$Specialisation, na.rm=TRUE))
plot(density(SimMammal$Trophic_level, na.rm=TRUE))
plot(density(SimMammal$Diel_activity, na.rm=TRUE))

# Birds
Results_Birds <- Comparison_null_dist(Birds, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaBirds <- Results_Birds[[1]]
SimBirds <- Results_Birds[[2]]
pBirds <- Plot_Delta(DeltaBirds) + scale_x_discrete(labels=c("DA", "Sp", "TL", "PD")) + xlab("") + ylab(expression(delta))
plot(density(SimBirds$Specialisation, na.rm=TRUE))
plot(density(SimBirds$Trophic_level, na.rm=TRUE))
plot(density(SimBirds$Diel_activity, na.rm=TRUE))
plot(density(SimBirds$Primary_diet, na.rm=TRUE))

# Reptiles
Results_Reptiles <- Comparison_null_dist(Reptiles, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaReptiles <- Results_Reptiles[[1]]
SimReptiles <- Results_Reptiles[[2]]
pReptiles <- Plot_Delta(DeltaReptiles) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimReptiles$Specialisation, na.rm=TRUE))
plot(density(SimReptiles$Trophic_level, na.rm=TRUE))
plot(density(SimReptiles$Diel_activity, na.rm=TRUE))

# Amphibians
Results_Amphibians <- Comparison_null_dist(Amphibians, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaAmphibians <- Results_Amphibians[[1]]
SimAmphibians <- Results_Amphibians[[2]]
pAmphibians <- Plot_Delta(DeltaAmphibians) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimAmphibians$Specialisation, na.rm=TRUE))
plot(density(SimAmphibians$Trophic_level, na.rm=TRUE))
plot(density(SimAmphibians$Diel_activity, na.rm=TRUE))
plot(density(SimAmphibians$Primary_diet, na.rm=TRUE))

p <- ggarrange(pMammals + labs(tag = "A") + theme(plot.tag.position = c(0.95, 0.94)), 
               pBirds + labs(tag = "B") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""),
               pReptiles + labs(tag = "C") + theme(plot.tag.position = c(0.95, 0.94)), 
               pAmphibians + labs(tag = "D") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""))
p
ggsave(p, file="../../Results/Plots/Phylosignal_categorical/Allclasses.pdf", height=4, width=5)



# Wilcoxon Tests
WBirds <- Wilcox_Test(SimBirds, DeltaBirds); WBirds$Class <- "Birds"
WReptiles <- Wilcox_Test(SimReptiles, DeltaReptiles); WReptiles$Class <- "Reptiles"
WMammals <- Wilcox_Test(SimMammal, DeltaMammal); WMammals$Class <- "Mammals"
WAmphibians <- Wilcox_Test(SimAmphibians, DeltaAmphibians); WAmphibians$Class <- "Amphibians"

WTests <- rbind(WBirds, WReptiles, WMammals, WAmphibians)
write.csv(WTests,"../../Results/1.Traits_before_imputations/Phylogenetic_signal/Categorical_traits_significance_Wilcoxon.csv", row.names = FALSE)


