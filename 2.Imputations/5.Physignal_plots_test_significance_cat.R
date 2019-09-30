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

## Function to plot results for categorical traits
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

## Function to plot results for categorical traits (delta) for both original and corrected phylogenies
Plot_Delta_All <- function(Delta) {
  
  Delta$Phylogeny <- factor(Delta$Phylogeny, levels=c("Original", "Modified"))
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=13), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=13),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  p <- ggplot(Delta, aes(Trait, delta_obs, group=Phylogeny, col=Phylogeny)) + 
    geom_point(position = position_dodge(width=0.2), shape=24, aes(fill=Phylogeny)) +
    geom_point(aes(Trait, MedianSim), position = position_dodge(width=0.2)) +
    geom_errorbar(aes(ymin=CI_low, ymax=CI_up), width=0.1, position = position_dodge(width=0.2)) +
    scale_y_continuous(trans="log10") +
    ylab(expression(delta)) +
    GGPoptions
  
  return(p)
}

## Function to plot results for continuous traits (lambda) for both original and corrected phylogenies
Plot_Lambda_All <- function(Lambda_c, Lambda_o) {
  
  Lambda <- rbind(Lambda_c, Lambda_o)[c(1,5),]
  Lambda<- reshape::melt(Lambda)
  
  Lambda$Phylogeny <- factor(Lambda$Phylogeny, levels=c("Original", "Modified"))
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=13), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=13),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  p <- ggplot(Lambda, aes(variable, value, group=Phylogeny, col=Phylogeny)) + 
    geom_point(shape=24, aes(fill=Phylogeny), alpha=0.8) + #position = position_dodge(width=0.2), 
    geom_hline(yintercept= 0.9, linetype="dashed") +
   #scale_y_continuous(trans="log10", breaks = c(0.7,0.8,0.9,1)) +
    ylab(expression(lambda)) +
    GGPoptions
  
  return(p)
}


Plot_Lambda <- function(Lambda_c) {
  
  Lambda <- Lambda_c[1,]
  Lambda<- reshape::melt(Lambda)
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=13), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=13),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  p <- ggplot(Lambda, aes(variable, value)) + 
    geom_point(shape=24, alpha=0.8, col="blue", fill="blue") + #position = position_dodge(width=0.2), 
    geom_hline(yintercept= 0.9, linetype="dashed") +
    #scale_y_continuous(trans="log10", breaks = c(0.7,0.8,0.9,1)) +
    ylab(expression(lambda)) +
    GGPoptions
  
  return(p)
}


## Function to run Wilcoxon Rank Sum and Signed rank tests
## Testing whether the observed median is significantly different from the null distribution for phylogenetic signal in categorcial traits
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
# with corrected tree -- categorical
Mammals_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalMammals.rds")
Birds_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalBirds.rds")
Reptiles_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalReptiles.rds")
Amphibians_c <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/CategoricalAmphibians.rds")

# with original tree -- categorical
Mammals_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalMammals.rds")
Birds_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalBirds.rds")
Reptiles_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalReptiles.rds")
Amphibians_o <- readRDS("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalAmphibians.rds")


## Continuous traits (for plotting)
Lambda_M_c <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/ContinuousMammals_log.csv")
Lambda_B_c <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/ContinuousBirds_log.csv")
Lambda_R_c <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/ContinuousReptiles_log.csv")
Lambda_A_c <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_corrected_phylogenies/ContinuousAmphibians_log.csv")

Lambda_M_o <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousMammals_log10.csv")
Lambda_B_o <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousBirds_log10.csv")
Lambda_R_o <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousReptiles_log10.csv")
Lambda_A_o <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousAmphibians_log10.csv")



## Phylogenetic signal in categorical traits with corrected phylogenies
# Mammals_c
Results_Mammals_c <- Comparison_null_dist(Mammals_c, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaMammals_c <- Results_Mammals_c[[1]]
SimMammals_c <- Results_Mammals_c[[2]]
pMammals_c <- Plot_Delta(DeltaMammals_c) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimMammals_c$Specialisation, na.rm=TRUE))
# plot(density(SimMammals_c$Trophic_level, na.rm=TRUE))
# plot(density(SimMammals_c$Diel_activity, na.rm=TRUE))

# Birds_c
Results_Birds_c <- Comparison_null_dist(Birds_c, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaBirds_c <- Results_Birds_c[[1]]
SimBirds_c <- Results_Birds_c[[2]]
pBirds_c <- Plot_Delta(DeltaBirds_c) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimBirds_c$Specialisation, na.rm=TRUE))
# plot(density(SimBirds_c$Trophic_level, na.rm=TRUE))
# plot(density(SimBirds_c$Diel_activity, na.rm=TRUE))
# plot(density(SimBirds_c$Primary_diet, na.rm=TRUE))

# Reptiles_c
Results_Reptiles_c <- Comparison_null_dist(Reptiles_c, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaReptiles_c <- Results_Reptiles_c[[1]]
SimReptiles_c <- Results_Reptiles_c[[2]]
pReptiles_c <- Plot_Delta(DeltaReptiles_c) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimReptiles_c$Specialisation, na.rm=TRUE))
# plot(density(SimReptiles_c$Trophic_level, na.rm=TRUE))
# plot(density(SimReptiles_c$Diel_activity, na.rm=TRUE))

# Amphibians_c
Results_Amphibians_c <- Comparison_null_dist(Amphibians_c, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaAmphibians_c <- Results_Amphibians_c[[1]]
SimAmphibians_c <- Results_Amphibians_c[[2]]
pAmphibians_c <- Plot_Delta(DeltaAmphibians_c) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimAmphibians_c$Specialisation, na.rm=TRUE))
# plot(density(SimAmphibians_c$Trophic_level, na.rm=TRUE))
# plot(density(SimAmphibians_c$Diel_activity, na.rm=TRUE))
# plot(density(SimAmphibians_c$Primary_diet, na.rm=TRUE))


pcatc <- ggarrange(pMammals_c + theme(plot.tag.position = c(0.95, 0.94)) + ggtitle("Mammals"),
               pBirds_c  + ggtitle("Birds") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""),
               pReptiles_c  + ggtitle("Reptiles") + theme(plot.tag.position = c(0.95, 0.94)),
               pAmphibians_c  + ggtitle("Amphibians") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""))
# ggsave(pcatc, file="../../Results/Plots/Phylosignal_categorical/Allclasses_correctedtree.pdf", height=4, width=5)
ggsave(pcatc, file="../../Results/Plots_CBER_talk_20119/Pysicat.png", height=4, width=5)

## Phylogenetic signal in categorical traits with original phylogenies
# Mammals_o
Results_Mammals_o <- Comparison_null_dist(Mammals_o, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaMammals_o <- Results_Mammals_o[[1]]
SimMammals_o <- Results_Mammals_o[[2]]
pMammals_o <- Plot_Delta(DeltaMammals_o) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimMammals_o$Specialisation, na.rm=TRUE))
# plot(density(SimMammals_o$Trophic_level, na.rm=TRUE))
# plot(density(SimMammals_o$Diel_activity, na.rm=TRUE))

# Birds_o
Results_Birds_o <- Comparison_null_dist(Birds_o, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaBirds_o <- Results_Birds_o[[1]]
SimBirds_o <- Results_Birds_o[[2]]
pBirds_o <- Plot_Delta(DeltaBirds_o) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimBirds_o$Specialisation, na.rm=TRUE))
# plot(density(SimBirds_o$Trophic_level, na.rm=TRUE))
# plot(density(SimBirds_o$Diel_activity, na.rm=TRUE))
# plot(density(SimBirds_o$Primary_diet, na.rm=TRUE))

# Reptiles_o
Results_Reptiles_o <- Comparison_null_dist(Reptiles_o, c("Specialisation", "Trophic_level", "Diel_activity")) 
DeltaReptiles_o <- Results_Reptiles_o[[1]]
SimReptiles_o <- Results_Reptiles_o[[2]]
pReptiles_o <- Plot_Delta(DeltaReptiles_o) + scale_x_discrete(labels=c("DA", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimReptiles_o$Specialisation, na.rm=TRUE))
# plot(density(SimReptiles_o$Trophic_level, na.rm=TRUE))
# plot(density(SimReptiles_o$Diel_activity, na.rm=TRUE))

# Amphibians_o
Results_Amphibians_o <- Comparison_null_dist(Amphibians_o, c("Specialisation", "Trophic_level", "Diel_activity", "Primary_diet")) 
DeltaAmphibians_o <- Results_Amphibians_o[[1]]
SimAmphibians_o <- Results_Amphibians_o[[2]]
pAmphibians_o <- Plot_Delta(DeltaAmphibians_o) + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) + xlab("") + ylab(expression(delta))
# plot(density(SimAmphibians_o$Specialisation, na.rm=TRUE))
# plot(density(SimAmphibians_o$Trophic_level, na.rm=TRUE))
# plot(density(SimAmphibians_o$Diel_activity, na.rm=TRUE))
# plot(density(SimAmphibians_o$Primary_diet, na.rm=TRUE))


# pcato <- ggarrange(pMammals_o + labs(tag = "A") + theme(plot.tag.position = c(0.95, 0.94)), 
#                pBirds_o + labs(tag = "B") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""),
#                pReptiles_o + labs(tag = "C") + theme(plot.tag.position = c(0.95, 0.94)), 
#                pAmphibians_o + labs(tag = "D") + theme(plot.tag.position = c(0.95, 0.94)) + ylab(""))


## Arrange plots together (original + corrected) for categorical traits
pMammals_c$data$Class <- "Mammals"
pBirds_c$data$Class <- "Birds"
pReptiles_c$data$Class <- "Reptiles"
pAmphibians_c$data$Class <- "Amphibians"
Cat_corrected <- rbind(pMammals_c$data, pBirds_c$data, pReptiles_c$data, pAmphibians_c$data) %>%
  mutate(Phylogeny="Modified")

pMammals_o$data$Class <- "Mammals"
pBirds_o$data$Class <- "Birds"
pReptiles_o$data$Class <- "Reptiles"
pAmphibians_o$data$Class <- "Amphibians"
Cat_original <- rbind(pMammals_o$data, pBirds_o$data, pReptiles_o$data, pAmphibians_o$data) %>%
  mutate(Phylogeny="Original")


Cat_signal_all <- rbind(Cat_corrected, Cat_original)
Cat_signal_all[31, ] <- c(NA, NA, NA, NA, "Primary_diet", "Reptiles", "Original")
Cat_signal_all[32, ] <- c(NA, NA, NA, NA, "Primary_diet", "Reptiles", "Modified")
Cat_signal_all$delta_obs <- as.numeric(Cat_signal_all$delta_obs)
Cat_signal_all$MedianSim <- as.numeric(Cat_signal_all$MedianSim)
Cat_signal_all$CI_low <- as.numeric(Cat_signal_all$CI_low)
Cat_signal_all$CI_up <- as.numeric(Cat_signal_all$CI_up)


pCat <- ggarrange(
Plot_Delta_All(Cat_signal_all[Cat_signal_all$Class=="Mammals",]) + xlab("") + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL")) +
  labs(tag = "A") + theme(plot.tag.position = "topleft"),
Plot_Delta_All(Cat_signal_all[Cat_signal_all$Class=="Birds",])+ xlab("") + ylab("") + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL"))+
  labs(tag = "B") + theme(plot.tag.position = "topleft") + scale_y_continuous(trans="log10", labels=scales::number_format(decimal.mark = '.', accuracy = 0.01)),
Plot_Delta_All(Cat_signal_all[Cat_signal_all$Class=="Reptiles",])+ xlab("") + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL"))+
  labs(tag = "C") + theme(plot.tag.position = "topleft"),
Plot_Delta_All(Cat_signal_all[Cat_signal_all$Class=="Amphibians",])+ xlab("")+ylab("") + scale_x_discrete(labels=c("DA", "PD", "Sp", "TL"))+
  labs(tag = "D") + theme(plot.tag.position = "topleft"),
common.legend = TRUE, legend="right"
)

ggsave(pCat, filename="../../Results/Plots/Phylosignal/Categorical.pdf", width=8, height=6)

## Plotting continuous signals
Lambda_B_c$Class <- "Birds"
Lambda_B_c$Phylogeny <- "Modified"
Lambda_M_c$Class <- "Mammals"
Lambda_M_c$Phylogeny <- "Modified"
Lambda_R_c$Class <- "Reptiles"
Lambda_R_c$Phylogeny <- "Modified"
Lambda_A_c$Class <- "Amphibians"
Lambda_A_c$Phylogeny <- "Modified"

Lambda_B_o$Class <- "Birds"
Lambda_B_o$Phylogeny <- "Original"
Lambda_M_o$Class <- "Mammals"
Lambda_M_o$Phylogeny <- "Original"
Lambda_R_o$Class <- "Reptiles"
Lambda_R_o$Phylogeny <- "Original"
Lambda_A_o$Class <- "Amphibians"
Lambda_A_o$Phylogeny <- "Original"


pCont <- ggarrange(
Plot_Lambda_All(Lambda_M_c, Lambda_M_o) + scale_x_discrete(labels=c("BM", "L", "LCS", "DB", "RS", "HB", "GL", "BL")) + xlab("")+
  labs(tag = "A") + theme(plot.tag.position ="topleft"),
Plot_Lambda_All(Lambda_B_c, Lambda_B_o) + scale_x_discrete(labels=c("BM", "L", "LCS", "DB", "RS", "HB", "GL")) + xlab("") + ylab("")+
  labs(tag = "B") + theme(plot.tag.position ="topleft"),
Plot_Lambda_All(Lambda_R_c, Lambda_R_o) + scale_x_discrete(labels=c("BM", "L", "LCS", "RS", "HB", "BL", "SM"))+ xlab("")+
  labs(tag = "C") + theme(plot.tag.position = "topleft"),
Plot_Lambda_All(Lambda_A_c, Lambda_A_o) + scale_x_discrete(labels=c("BM", "L", "LCS", "DB", "RS", "HB", "BL"))+ xlab("") + ylab("")+
  labs(tag = "D") + theme(plot.tag.position = "topleft"),
common.legend=TRUE, legend="right"
)
ggsave(pCont, filename="../../Results/Plots/Phylosignal/Continuous.pdf", width=8, height=6)



# Wilcoxon Tests of significance for categorical traits

# corrected phylogenies
WBirds <- Wilcox_Test(SimBirds_c, DeltaBirds_c); WBirds$Class <- "Birds"
WReptiles <- Wilcox_Test(SimReptiles_c, DeltaReptiles_c); WReptiles$Class <- "Reptiles"
WMammals <- Wilcox_Test(SimMammals_c, DeltaMammals_c); WMammals$Class <- "Mammals"
WAmphibians <- Wilcox_Test(SimAmphibians_c, DeltaAmphibians_c); WAmphibians$Class <- "Amphibians"

WTests <- rbind(WBirds, WReptiles, WMammals, WAmphibians)
max(WTests$W_pvalue)
write.csv(WTests,"../../Results/1.Traits_before_imputations/Phylogenetic_signal/Categorical_traits_significance_Wilcoxon_corrected.csv", row.names = TRUE)

# uncorrected phylogenies
WBirds <- Wilcox_Test(SimBirds_o, DeltaBirds_o); WBirds$Class <- "Birds"
WReptiles <- Wilcox_Test(SimReptiles_o, DeltaReptiles_o); WReptiles$Class <- "Reptiles"
WMammals <- Wilcox_Test(SimMammals_o, DeltaMammals_o); WMammals$Class <- "Mammals"
WAmphibians <- Wilcox_Test(SimAmphibians_o, DeltaAmphibians_o); WAmphibians$Class <- "Amphibians"

WTestsO <- rbind(WBirds, WReptiles, WMammals, WAmphibians)
max(WTestsO$W_pvalue)
write.csv(WTestsO,"../../Results/1.Traits_before_imputations/Phylogenetic_signal/Categorical_traits_significance_Wilcoxon_original.csv", row.names = TRUE)

pcomtc<-ggarrange(
Plot_Lambda(Lambda_M_c) + scale_x_discrete(labels=c("BM", "L", "LCS", "DB", "RS", "HB", "GL", "BL")) + xlab("") + ggtitle("Mammals"),
Plot_Lambda(Lambda_B_c) + scale_x_discrete(labels=c("BM", "L", "LCS", "DB", "RS", "HB", "GL")) + xlab("") + ylab("") + ggtitle("Birds"),
Plot_Lambda(Lambda_R_c) + scale_x_discrete(labels=c("BM", "L", "LCS", "RS", "HB", "BL", "SM"))+ xlab("")+ ggtitle("Reptiles"),
Plot_Lambda(Lambda_A_c) + scale_x_discrete(labels=c("BM", "L", "LCS", "DB", "RS", "HB", "BL"))+ xlab("")+ ylab("") + ggtitle("Amphibians")
)
ggsave(pcomtc, file="../../Results/Plots_CBER_talk_20119/Pysicont.png", height=4, width=6)

