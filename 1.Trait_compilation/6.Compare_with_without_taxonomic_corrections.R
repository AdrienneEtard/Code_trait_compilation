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



# ## Delta in trait coverage before and after taxonomic correction -- FOR ALL TRAITS
# 
# # For all species
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot.Delta.Cov(UN_Mammals, C_Mammals, TraitsMammals, FALSE, UN_Predicts, C_Predicts, "Mammals")
# Plot.Delta.Cov(UN_Birds, C_Birds, TraitsBirds, FALSE, UN_Predicts, C_Predicts, "Birds")
# Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TraitsReptiles, FALSE, UN_Predicts, C_Predicts, "Reptiles")
# Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TraitsAmphibians, FALSE, UN_Predicts, C_Predicts, "Amphibians")
# box("outer")
# 
# ## For PREDICTS species
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(1,2,1,1))
# par(mfrow=c(2,2))
# Plot.Delta.Cov(UN_Mammals, C_Mammals, TraitsMammals, TRUE, UN_Predicts, C_Predicts, "Mammals")
# Plot.Delta.Cov(UN_Birds, C_Birds, TraitsBirds, TRUE, UN_Predicts, C_Predicts, "Birds")
# Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TraitsReptiles, TRUE, UN_Predicts, C_Predicts, "Reptiles")
# Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TraitsAmphibians, TRUE, UN_Predicts, C_Predicts, "Amphibians")
# box("outer")


## Delta in trait coverage before and after taxonomic correction -- FOR TARGET TRAITS

TargetTraits <- c("log10_Body_mass_g",
                  "log10_Longevity_d",
                  "log10_Litter_size", 
                  "Range_size_m2", 
                  "Diel_activity",
                  "Trophic_level",
                  "sqrt_Diet_breadth",
                  "Specialisation",
                  "sqrt_Habitat_breadth_IUCN",
                  "Primary_diet")

# For all species

pdf(file="../../Results/Plots/Coverage/Target_traits/All_species.pdf", width=7, height=5.5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TargetTraits, FALSE, UN_Predicts, C_Predicts, "A. Mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TargetTraits, FALSE, UN_Predicts, C_Predicts, "B. Birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TargetTraits, FALSE, UN_Predicts, C_Predicts, "C. Reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TargetTraits, FALSE, UN_Predicts, C_Predicts, "D. Amphibians")
box("outer")

mtext(at=50, line=-11, "% coverage", cex=0.8)
mtext(at=-130, line=-11, "% coverage", cex=0.8)

legend(x=-110, y=-6, xpd="NA", title="Trait coverage across all species",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("lightgrey", "deepskyblue3"),  box.lty=0)

dev.off()


# For PREDICTS species
pdf(file="../../Results/Plots/Coverage/Target_traits/Predicts.pdf", width=7, height=5.5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TargetTraits, TRUE, UN_Predicts, C_Predicts, "A. PREDICTS mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TargetTraits, TRUE, UN_Predicts, C_Predicts, "B. PREDICTS birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TargetTraits, TRUE, UN_Predicts, C_Predicts, "C. PREDICTS reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TargetTraits, TRUE, UN_Predicts, C_Predicts, "D. PREDICTS amphibians")
box("outer")

mtext(at=50, line=-11, "% coverage", cex=0.8)
mtext(at=-130, line=-11, "% coverage", cex=0.8)

legend(x=-110, y=-6, xpd="NA", title="Trait coverage across all species",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("lightgrey", "deepskyblue3"),  box.lty=0)

dev.off()



## Delta in trait coverage before and after taxonomic correction -- FOR TRAITS USED IN MISSFOREST IMPUTATIONS only

Traits_cont <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", "Range_size_m2", "sqrt_Habitat_breadth_IUCN", "sqrt_Diet_breadth")
Traits_cat <- c("Specialisation", "Diel_activity","Trophic_level")

TMammalsI <- c(Traits_cont, Traits_cat, "log10_Generation_length_d", "log10_Adult_svl_cm")
TBirdsI <- c(Traits_cont, Traits_cat)
TReptilesI <- c(Traits_cont, Traits_cat, "log10_Adult_svl_cm", "log10_Maturity_d")
TAmphibiansI <- c(Traits_cont, Traits_cat, "log10_Body_length_mm")

# For all species

pdf(file="../../Results/Plots/Coverage/Predictor_traits/All_species.pdf", width=7, height=5.5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TMammalsI, FALSE, UN_Predicts, C_Predicts, "A. Mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TBirdsI, FALSE, UN_Predicts, C_Predicts, "B. Birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TReptilesI, FALSE, UN_Predicts, C_Predicts, "C. Reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TAmphibiansI, FALSE, UN_Predicts, C_Predicts, "D. Amphibians")
box("outer")

mtext(at=50, line=-11, "% coverage", cex=0.8)
mtext(at=-130, line=-11, "% coverage", cex=0.8)

legend(x=-110, y=-6, xpd="NA", title="Trait coverage across all species",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("lightgrey", "deepskyblue3"),  box.lty=0)
dev.off()


# For PREDICTS species
pdf(file="../../Results/Plots/Coverage/Predictor_traits/Predicts.pdf", width=7, height=5.5, family="Times", pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(5,2,1,1))
par(mfrow=c(2,2))
Plot.Delta.Cov(UN_Mammals, C_Mammals, TMammalsI, TRUE, UN_Predicts, C_Predicts, "A. PREDICTS mammals")
Plot.Delta.Cov(UN_Birds, C_Birds, TBirdsI, TRUE, UN_Predicts, C_Predicts, "B. PREDICTS birds")
Plot.Delta.Cov(UN_Reptiles, C_Reptiles, TReptilesI, TRUE, UN_Predicts, C_Predicts, "C. PREDICTS reptiles")
Plot.Delta.Cov(UN_Amphibians, C_Amphibians, TAmphibiansI, TRUE, UN_Predicts, C_Predicts, "D. PREDICTS amphibians")
box("outer")

mtext(at=50, line=-11, "% coverage", cex=0.8)
mtext(at=-130, line=-11, "% coverage", cex=0.8)

legend(x=-110, y=-6, xpd="NA", title="Trait coverage across all species",
       legend = c("Without taxonomic correction", "With taxonomic correction"),
       fill = c("lightgrey", "deepskyblue3"),  box.lty=0)
dev.off()



