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
    geom_point() +
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
# with corrected tree
Mammals_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalMammals.rds")
Birds_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalBirds.rds")
Reptiles_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalReptiles.rds")
Amphibians_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalAmphibians.rds")

# with original tree
Mammals_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalMammals.rds")
Birds_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalBirds.rds")
Reptiles_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalReptiles.rds")
Amphibians_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalAmphibians.rds")


## Phylogenetic signal in categorical traits with corrected phylogenies
# Mammals_c
Results_Mammals_c <- Comparison_null_dist(Mammals_c, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaMammals_c <- Results_Mammals_c[[1]]
SimMammals_c <- Results_Mammals_c[[2]]
pMammals_c <- Plot_Delta(DeltaMammals_c) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimMammals_c$Specialisation, na.rm=TRUE))
plot(density(SimMammals_c$Trophic_level, na.rm=TRUE))
plot(density(SimMammals_c$Diel_activity, na.rm=TRUE))

# Birds_c
Results_Birds_c <- Comparison_null_dist(Birds_c, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaBirds_c <- Results_Birds_c[[1]]
SimBirds_c <- Results_Birds_c[[2]]
pBirds_c <- Plot_Delta(DeltaBirds_c) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimBirds_c$Specialisation, na.rm=TRUE))
plot(density(SimBirds_c$Trophic_level, na.rm=TRUE))
plot(density(SimBirds_c$Diel_activity, na.rm=TRUE))
plot(density(SimBirds_c$Primary_diet, na.rm=TRUE))

# Reptiles_c
Results_Reptiles_c <- Comparison_null_dist(Reptiles_c, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaReptiles_c <- Results_Reptiles_c[[1]]
SimReptiles_c <- Results_Reptiles_c[[2]]
pReptiles_c <- Plot_Delta(DeltaReptiles_c) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimReptiles_c$Specialisation, na.rm=TRUE))
plot(density(SimReptiles_c$Trophic_level, na.rm=TRUE))
plot(density(SimReptiles_c$Diel_activity, na.rm=TRUE))

# Amphibians_c
Results_Amphibians_c <- Comparison_null_dist(Amphibians_c, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaAmphibians_c <- Results_Amphibians_c[[1]]
SimAmphibians_c <- Results_Amphibians_c[[2]]
pAmphibians_c <- Plot_Delta(DeltaAmphibians_c) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimAmphibians_c$Specialisation, na.rm=TRUE))
plot(density(SimAmphibians_c$Trophic_level, na.rm=TRUE))
plot(density(SimAmphibians_c$Diel_activity, na.rm=TRUE))
plot(density(SimAmphibians_c$Primary_diet, na.rm=TRUE))


p <- ggarrange(pMammals_c + labs(tag = "A") + theme(plot.tag.position = c(0.95, 0.94)), 
               pBirds_c + labs(tag = "B") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""),
               pReptiles_c + labs(tag = "C") + theme(plot.tag.position = c(0.95, 0.94)), 
               pAmphibians_c + labs(tag = "D") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""))
p
ggsave(p, file="../../Results/Plots/Phylosignal_categorical/Allclasses_correctedtree.pdf", height=4, width=5)


## Phylogenetic signal in categorical traits with original phylogenies
# Mammals_o
Results_Mammals_o <- Comparison_null_dist(Mammals_o, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaMammals_o <- Results_Mammals_o[[1]]
SimMammals_o <- Results_Mammals_o[[2]]
pMammals_o <- Plot_Delta(DeltaMammals_o) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimMammals_o$Specialisation, na.rm=TRUE))
plot(density(SimMammals_o$Trophic_level, na.rm=TRUE))
plot(density(SimMammals_o$Diel_activity, na.rm=TRUE))

# Birds_o
Results_Birds_o <- Comparison_null_dist(Birds_o, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaBirds_o <- Results_Birds_o[[1]]
SimBirds_o <- Results_Birds_o[[2]]
pBirds_o <- Plot_Delta(DeltaBirds_o) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimBirds_o$Specialisation, na.rm=TRUE))
plot(density(SimBirds_o$Trophic_level, na.rm=TRUE))
plot(density(SimBirds_o$Diel_activity, na.rm=TRUE))
plot(density(SimBirds_o$Primary_diet, na.rm=TRUE))

# Reptiles_o
Results_Reptiles_o <- Comparison_null_dist(Reptiles_o, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaReptiles_o <- Results_Reptiles_o[[1]]
SimReptiles_o <- Results_Reptiles_o[[2]]
pReptiles_o <- Plot_Delta(DeltaReptiles_o) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimReptiles_o$Specialisation, na.rm=TRUE))
plot(density(SimReptiles_o$Trophic_level, na.rm=TRUE))
plot(density(SimReptiles_o$Diel_activity, na.rm=TRUE))

# Amphibians_o
Results_Amphibians_o <- Comparison_null_dist(Amphibians_o, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaAmphibians_o <- Results_Amphibians_o[[1]]
SimAmphibians_o <- Results_Amphibians_o[[2]]
pAmphibians_o <- Plot_Delta(DeltaAmphibians_o) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
plot(density(SimAmphibians_o$Specialisation, na.rm=TRUE))
plot(density(SimAmphibians_o$Trophic_level, na.rm=TRUE))
plot(density(SimAmphibians_o$Diel_activity, na.rm=TRUE))
plot(density(SimAmphibians_o$Primary_diet, na.rm=TRUE))


p <- ggarrange(pMammals_o + labs(tag = "A") + theme(plot.tag.position = c(0.95, 0.94)), 
               pBirds_o + labs(tag = "B") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""),
               pReptiles_o + labs(tag = "C") + theme(plot.tag.position = c(0.95, 0.94)), 
               Amphibians_o + labs(tag = "D") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""))


# Wilcoxon Tests
WBirds <- Wilcox_Test(SimBirds_c, DeltaBirds_c); WBirds$Class <- "Birds"
WReptiles <- Wilcox_Test(SimReptiles_c, DeltaReptiles_c); WReptiles$Class <- "Reptiles"
WMammals <- Wilcox_Test(SimMammals_c, DeltaMammals_c); WMammals$Class <- "Mammals"
WAmphibians <- Wilcox_Test(SimAmphibians_c, DeltaAmphibians_c); WAmphibians$Class <- "Amphibians"

WTests <- rbind(WBirds, WReptiles, WMammals, WAmphibians)
min(WTests$W_pvalue)
write.csv(WTests,"../../Results/1.Traits_before_imputations/Phylogenetic_signal/Categorical_traits_significance_Wilcoxon.csv", row.names = FALSE)


