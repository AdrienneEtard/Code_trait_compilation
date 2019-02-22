## Comparison of basic trait coverage with VS without taxonomic correction
## and of number of species for which traits were compiled


# # #   P  R  E  A  M  B  L  E

X <- c("dplyr", "ggplot2", "ggpubr", "grid", "cowplot")
lapply(X, library, character.only=TRUE); rm(X)
source("Comparison_with_without_taxonomic_corrections_functions.R")

## Load trait files

# No taxonomic correction
UN_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv")
UN_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv")
UN_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv")
UN_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv")

# With taxonomic correction
C_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv")
C_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv")
C_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv")
C_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv")

## Load PREDICTS
UN_Predicts <- readRDS("../../Data/PREDICTS_database.rds") %>% 
  filter(Class %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"))

C_Predicts <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")


## 1. Differences in species number

# Amphibians
Delta_Number(UN_Amphibians, C_Amphibians, "Amphibia", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Amphibians, C_Amphibians, "Amphibia", "Predicts", C_Predicts, UN_Predicts)

# Mammals
Delta_Number(UN_Mammals, C_Mammals, "Mammalia", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Mammals, C_Mammals, "Mammalia", "Predicts", C_Predicts, UN_Predicts)

# Birds
Delta_Number(UN_Birds, C_Birds, "Aves", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Birds, C_Birds, "Aves", "Predicts", C_Predicts, UN_Predicts)

# Reptiles
Delta_Number(UN_Reptiles, C_Reptiles, "Reptilia", "all species", C_Predicts, UN_Predicts)
Delta_Number(UN_Reptiles, C_Reptiles, "Reptilia", "Predicts", C_Predicts, UN_Predicts)


## 2. Plotting coverage

# 2.1. % Species representation in phylogenies
CovPhyloAll <- Phylo_Delta(C_Mammals, C_Birds, C_Reptiles, C_Amphibians, C_Predicts, FALSE,
                           UN_Mammals, UN_Birds, UN_Reptiles, UN_Amphibians, UN_Predicts)

CovPhyloPredicts <- Phylo_Delta(C_Mammals, C_Birds, C_Reptiles, C_Amphibians, C_Predicts, TRUE,
                                UN_Mammals, UN_Birds, UN_Reptiles, UN_Amphibians, UN_Predicts)

p <- PlotPhyloCov(CovPhyloAll,CovPhyloPredicts, 15)
ggsave(p, file="../../Results/Plots/Coverage/Phylogenies/Species_representation.pdf", width=8, height=4.2)


# 2.2. % Trait coverage

## Plotting coverage, basic or superposed (with and without taxonomic correction)

TraitsBirds <-  c("log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Maturity_d", "log10_Longevity_d",
                  "log10_Litter_size", "Range_size_m2", "Diel_activity", "Primary_diet","Trophic_level","Specialisation",
                  "sqrt_Habitat_breadth_IUCN")

TraitsMammals <-  c("log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Forearm_length_mm","log10_Head_length_mm",
                    "log10_Generation_length_d", "log10_Maturity_d", "log10_Longevity_d", "log10_AFR_d",
                    "log10_Litter_size", "Range_size_m2", "Diel_activity", "Primary_diet","Trophic_level","Specialisation",
                    "sqrt_Habitat_breadth_IUCN")

TraitsAmphibians <-  c("log10_Body_mass_g", "log10_Body_length_mm","log10_Svl_length_mm", "log10_Maturity_d", "log10_Longevity_d",
                       "log10_Litter_size", "Range_size_m2", "Diel_activity","Primary_diet", "Trophic_level","Specialisation",
                       "sqrt_Habitat_breadth_IUCN")

TraitsReptiles <-  c("log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Maturity_d", "log10_Longevity_d",
                     "log10_Litter_size", "Range_size_m2", "Diel_activity", "Primary_diet","Trophic_level", "Specialisation",
                     "sqrt_Habitat_breadth_IUCN")



# # For all species -- corrected
# #pdf(file="../../Results/Plots/Trait_coverage/For_trait_selection/All_species.pdf", width=8, height=5, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot.Cov(C_Birds, TraitsBirds, "Birds", FALSE, C_Predicts)
# Plot.Cov(C_Mammals, TraitsMammals, "Mammals", FALSE, C_Predicts)
# Plot.Cov(C_Amphibians, TraitsAmphibians, "Amphibians", FALSE, C_Predicts)
# Plot.Cov(C_Reptiles, TraitsReptiles, "Reptiles", FALSE, C_Predicts)
# #dev.off()
# 
# # For PREDICTS species -- corrected
# #pdf(file="../../Results/Plots/Trait_coverage/For_trait_selection/Predicts_species.pdf", width=8, height=5, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot.Cov(C_Birds, TraitsBirds, "Birds", TRUE, C_Predicts)
# Plot.Cov(C_Mammals, TraitsMammals, "Mammals", TRUE, C_Predicts)
# Plot.Cov(C_Amphibians, TraitsAmphibians, "Amphibians", TRUE, C_Predicts)
# Plot.Cov(C_Reptiles, TraitsReptiles, "Reptiles", TRUE, C_Predicts)
# #dev.off()


## Delta in trait coverage before and after taxonomic correction -- FOR ALL TRAITS

# For all species

# pdf(file="../../Results/Plots/Trait_coverage/All_Delta_taxonomy.pdf", width=8, height=5, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TraitsMammals, FALSE, UN_Predicts, C_Predicts, "Mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TraitsBirds, FALSE, UN_Predicts, C_Predicts, "Birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TraitsReptiles, FALSE, UN_Predicts, C_Predicts, "Reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TraitsAmphibians, FALSE, UN_Predicts, C_Predicts, "Amphibians")
box("outer")
# dev.off()

## For PREDICTS species
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TraitsMammals, TRUE, UN_Predicts, C_Predicts, "Mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TraitsBirds, TRUE, UN_Predicts, C_Predicts, "Birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TraitsReptiles, TRUE, UN_Predicts, C_Predicts, "Reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TraitsAmphibians, TRUE, UN_Predicts, C_Predicts, "Amphibians")
box("outer")


## Delta in trait coverage before and after taxonomic correction -- FOR TARGET TRAITS

TargetTraits <- c("log10_Body_mass_g",
                  "log10_Longevity_d",
                  "log10_Litter_size", 
                  "Range_size_m2", 
                  "Diel_activity",
                  "Trophic_level",
                  "Specialisation",
                  "sqrt_Habitat_breadth_IUCN",
                  "Primary_diet")

# For all species

pdf(file="../../Results/Plots/Coverage/Target_traits/All_species.pdf", width=7, height=5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TargetTraits, FALSE, UN_Predicts, C_Predicts, "A. Mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TargetTraits, FALSE, UN_Predicts, C_Predicts, "B. Birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TargetTraits, FALSE, UN_Predicts, C_Predicts, "C. Reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TargetTraits, FALSE, UN_Predicts, C_Predicts, "D. Amphibians")
box("outer")

mtext(at=50, line=-10, "% coverage", cex=0.8)
mtext(at=-110, line=-10, "% coverage", cex=0.8)

legend(x=-110, y=-6, xpd="NA", title="Trait coverage across all species",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("lightgrey", "deepskyblue3"),  box.lty=0)

dev.off()


# For PREDICTS species
pdf(file="../../Results/Plots/Coverage/Target_traits/Predicts.pdf", width=7, height=5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TargetTraits, TRUE, UN_Predicts, C_Predicts, "A. PREDICTS mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TargetTraits, TRUE, UN_Predicts, C_Predicts, "B. PREDICTS birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TargetTraits, TRUE, UN_Predicts, C_Predicts, "C. PREDICTS reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TargetTraits, TRUE, UN_Predicts, C_Predicts, "D. PREDICTS amphibians")
box("outer")

mtext(at=50, line=-10, "% coverage", cex=0.8)
mtext(at=-110, line=-10, "% coverage", cex=0.8)

legend(x=-110, y=-6, xpd="NA", title="Trait coverage across all species",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("lightgrey", "deepskyblue3"),  box.lty=0)

dev.off()



## Plotting percent information across species

Percent_info <- function(TraitDF, Traits) {
  
  TraitDF <- TraitDF[, c("Best_guess_binomial", Traits)]
  Results <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y))) %>% as.data.frame() %>%
    setNames(., "Percent")
  Results$Percent <- Results$Percent/length(Traits)*100
  
  return(Results)
  
}

Test <- Percent_info(Mammals, TraitsMammals)
Test <- Percent_info(Birds, TraitsBirds)
Test <- Percent_info(Amphibians, TraitsAmphibians)
Test <- Percent_info(Reptiles, TraitsReptiles)


plot(density(Test$Percent))
hist(Test$Percent, breaks = 10)



## On the basis of these plots, I exclude: 
#    -> head and forearm length; sexual maturity age and age at first reproduction (for mammals); 
#    -> svl length and sexual maturity age (for birds); 
#    -> svl length and sexual maturity age (for birds); 
#    -> svl length and sexual maturity age (for birds); 



# ## Plotting trait distributions
# 
# DistCont <- function(DF, Trait, Breaks, Main) {
#   par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0))
#   hist(DF[,Trait], breaks = Breaks, main=Main, xlab="")    
# }
# 
# 
# par(mfrow=c(2,2))
# DistCont(Reptiles, "log10_Body_mass_g", 100, "Reptiles")
# DistCont(Birds, "log10_Body_mass_g", 100, "Birds")
# DistCont(Amphibians, "log10_Body_mass_g", 100, "Amphibians")
# DistCont(Mammals, "log10_Body_mass_g", 100, "Mammals")


# hist(log10(Reptiles$Adult_svl_cm), breaks=100)
# 
# 
# ## log10 generation length, log10 longevity
# par(mfrow=c(2,2))
# hist(log10(Mammals$Generation_length_d), breaks=100)
# hist(log10(Birds$Generation_length_d), breaks=100)
# hist(log10(Amphibians$Generation_length_d), breaks=100)
# hist(log10(Reptiles$Generation_length_d), breaks=100)
# 
# hist(log10(Mammals$Longevity_d), breaks=100)
# hist(log10(Birds$Longevity_d), breaks=100)
# hist(log10(Amphibians$Longevity_d), breaks=100)
# hist(log10(Reptiles$Longevity_d), breaks=100)
# 
# ## log10 litter and clutch size
# par(mfrow=c(2,2))
# hist(log10(Mammals$Litter_size), breaks=100)
# hist(log10(Birds$Litter_size), breaks=100)
# hist(log10(Amphibians$Litter_size), breaks=100)
# hist(log10(Reptiles$Litter_size), breaks=100)
# 
# par(mfrow=c(2,2))
# hist(Mammals$Litter_size, breaks=100)
# hist(Birds$Litter_size, breaks=100)
# hist(Amphibians$Litter_size, breaks=100)
# hist(Reptiles$Litter_size, breaks=100)
# 
# ## SQRT habitat breadth
# par(mfrow=c(2,2))
# hist(sqrt(Mammals$Habitat_breadth_IUCN), breaks=100)
# hist(sqrt(Birds$Habitat_breadth_IUCN), breaks=100)
# hist(sqrt(Amphibians$Habitat_breadth_IUCN), breaks=100)
# hist(sqrt(Reptiles$Habitat_breadth_IUCN), breaks=100)
# 
# par(mfrow=c(2,2))
# hist(Mammals$Habitat_breadth_IUCN, breaks=100)
# hist(Birds$Habitat_breadth_IUCN, breaks=100)
# hist(Amphibians$Habitat_breadth_IUCN, breaks=100)
# hist(Reptiles$Habitat_breadth_IUCN, breaks=100)


# ## SQRT diet breadth
# par(mfrow=c(2,2))
# hist(sqrt(Mammals$Diet_breadth), breaks=100)
# hist(sqrt(Birds$Diet_breadth), breaks=100)
# hist(sqrt(Amphibians$Diet_breadth), breaks=100)
# 
# par(mfrow=c(2,2))
# hist(Mammals$Diet_breadth, breaks=100)
# hist(Birds$Diet_breadth, breaks=100)
# hist(Amphibians$Diet_breadth, breaks=100)







