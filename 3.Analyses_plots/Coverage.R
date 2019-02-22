## Plotting coverage

library(dplyr)

# # # # # # #
#           # 
# Load data # 
#           #
# # # # # # # 

## Mammals
TraitsMammal <- read.csv("../Results/1.Traits_before_imputations/TraitsMammals.csv")
TraitsMammal.Predicts <- read.csv("../Results/1.Traits_before_imputations/TraitsMammalsPredicts.csv")
TraitsMammal.Predicts_Imputed <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsMammal.Predicts_ContImputed.csv")
TraitsMammal.Predicts_Imputed_Randomadd <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsMammal.Predicts_ContImputed_Randomadd.csv")
TraitsMammals.Predicts_imputed_correlated <- read.csv("../Results/3.TraitsContinuous_after_correlations/TraitsMammals.Predicts_imputed_correlated.csv")


## Amphibians -- all
TraitsAmphibian_basic <- read.csv("../Results/1.Traits_before_imputations/Amphibians/TraitsAmphibians_basic_compilation.csv")
TraitsAmphibian_correlated <- read.csv("../Results/1.Traits_before_imputations/Amphibians/TraitsAmphibians_correlated.csv")

#             -- PREDICTS
TraitsAmphibian.Predicts_basic <- read.csv("../Results/1.Traits_before_imputations/Amphibians/TraitsAmphibiansPredicts_basic_compilation.csv")
TraitsAmphibian.Predicts_correlated <- read.csv("../Results/1.Traits_before_imputations/Amphibians/TraitsAmphibiansPredicts_correlated.csv")

# #             -- PREDICTS imputed
TraitsAmphibian.Predicts_Imputed <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsAmphibians.Predicts_ContImputed.csv")
TraitsAmphibian.Predicts_Imputed_Randomadd <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsAmphibians.Predicts_ContImputed_Randomadd.csv")
# TraitsAmphibian.Predicts_imputed_correlated <- read.csv("../Results/3.TraitsContinuous_after_correlations/TraitsAmphibians.Predicts_imputed_correlated.csv")



## Birds
TraitsBird <- read.csv("../Results/1.Traits_before_imputations/TraitsBird.csv")
TraitsBird.Predicts <- read.csv("../Results/1.Traits_before_imputations/TraitsBirdPredicts.csv")


## Reptiles -- all
TraitsReptile_basic <- read.csv("../Results/1.Traits_before_imputations/Reptiles/TraitsReptiles_basic_compilation.csv")
TraitsReptile_correlated <- read.csv("../Results/1.Traits_before_imputations/Reptiles/TraitsReptiles_correlated.csv")

#           -- Predicts 
TraitsReptile.Predicts_basic <- read.csv("../Results/1.Traits_before_imputations/Reptiles/TraitsReptilesPredicts_basic_compilation.csv")
TraitsReptile.Predicts_correlated <- read.csv("../Results/1.Traits_before_imputations/Reptiles/TraitsReptilePredicts_correlated.csv")


#           -- after imputations
TraitsReptile.Predicts_Imputed <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsReptile.Predicts_ContImputed.csv")




#  #  #  #  #  #                          F  u  n  c  t  i  o  n  s                 #  #  #  #  #  # 

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Plot coverage per trait --------------------------------------------------------------------------

Plot.Cov.Basic <- function(TraitData, Traits_name, Main, Xlab) {
  
  Completeness <- apply(TraitData[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness <- as.data.frame(Completeness/nrow(TraitData)*100)
  colnames(Completeness) <- "Completeness"
  Completeness <- Completeness[order(Completeness, decreasing=FALSE), , drop=FALSE]
  
  barplot(Completeness$Completeness, horiz = TRUE, 
          xlim = c(0,100), names.arg=rownames(Completeness), las=1, 
          main=paste(Main), xlab=paste(Xlab))
  
  abline(v=100, lty="dotted")
  
}


# Coverages on continuous traits before VS after correlations --------------

Plot.Cov <- function(TraitData1, TraitData2, Traits_name, Main) {
  
  Completeness1 <- apply(TraitData1[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness1 <- as.data.frame(Completeness1/nrow(TraitData1)*100)
  colnames(Completeness1) <- "Completeness"
  
  Completeness2 <- apply(TraitData2[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness2 <- as.data.frame(Completeness2/nrow(TraitData2)*100)
  colnames(Completeness2) <- "Completeness"
  
  Completeness2 <- Completeness2[order(Completeness2, decreasing=FALSE), , drop=FALSE]
  X <- row.names(Completeness2)
  
  Completeness1 <- Completeness1[X,] %>% as.data.frame()
  
  row.names(Completeness1) <- X
  colnames(Completeness1) <- "Completeness"
  
  barplot(Completeness2$Completeness, horiz = TRUE, 
          xlim = c(0,100), names.arg=rownames(Completeness2), las=1, col="blue", main=Main, xlab="% Coverage")
  
  barplot(Completeness1$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="grey")
  
  abline(v=100, lty="dotted")
  
}



# Coverages on continuous traits before VS after imputations --------------

Plot.Cov.Delta <- function(TraitData1, TraitData2, Traits_name, Main) {
  
  Completeness1 <- apply(TraitData1[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness1 <- as.data.frame(Completeness1/nrow(TraitData1)*100)
  colnames(Completeness1) <- "Completeness"
  
  Completeness2 <- apply(TraitData2[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness2 <- as.data.frame(Completeness2/nrow(TraitData2)*100)
  colnames(Completeness2) <- "Completeness"
  
  Completeness2 <- Completeness2[order(Completeness2, decreasing=FALSE), , drop=FALSE]
  X <- row.names(Completeness2)
  
  Completeness1 <- Completeness1[X,] %>% as.data.frame()

  row.names(Completeness1) <- X
  colnames(Completeness1) <- "Completeness"
  
  barplot(Completeness2$Completeness, horiz = TRUE, 
          xlim = c(0,100), names.arg=rownames(Completeness2), las=1, col="black", main=Main, xlab="% Coverage")
  
  barplot(Completeness1$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="grey")
  
  abline(v=100, lty="dotted")
  
}

# Coverages on continuous traits before VS after imputations with random addition of species at genus level in the phylogeny  --------------

Plot.Cov.Delta2 <- function(TraitData1, TraitData2, TraitData3, Traits_name, Main, x_lim) {
  
  Completeness1 <- apply(TraitData1[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness1 <- as.data.frame(Completeness1/nrow(TraitData1)*100)
  colnames(Completeness1) <- "Completeness"
  
  Completeness2 <- apply(TraitData2[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness2 <- as.data.frame(Completeness2/nrow(TraitData2)*100)
  colnames(Completeness2) <- "Completeness"
  
  Completeness3 <- apply(TraitData3[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness3 <- as.data.frame(Completeness3/nrow(TraitData3)*100)
  colnames(Completeness3) <- "Completeness"
  
  Completeness3 <- Completeness3[order(Completeness3, decreasing=FALSE), , drop=FALSE]
  X <- row.names(Completeness3)
  
  Completeness2 <- Completeness2[X,] %>% as.data.frame()
  row.names(Completeness2) <- X
  colnames(Completeness2) <- "Completeness"
  
  Completeness1 <- Completeness1[X,] %>% as.data.frame()
  row.names(Completeness1) <- X
  colnames(Completeness1) <- "Completeness"
  
  barplot(Completeness3$Completeness, horiz = TRUE, 
          names.arg=rownames(Completeness2), las=1, col="red", main=Main, xlab="% Coverage", xlim=x_lim)
  
  barplot(Completeness2$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="black")
  
  barplot(Completeness1$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="grey")
  
  abline(v=100, lty="dotted")
  
}


# Coverage after correlations ---------------------------------------------
Plot.Cov.Delta3 <- function(TraitData1, TraitData2, TraitData3, TraitData4, Traits_name, Main, x_lim) {
  
  Completeness1 <- apply(TraitData1[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness1 <- as.data.frame(Completeness1/nrow(TraitData1)*100)
  colnames(Completeness1) <- "Completeness"
  
  Completeness2 <- apply(TraitData2[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness2 <- as.data.frame(Completeness2/nrow(TraitData2)*100)
  colnames(Completeness2) <- "Completeness"
  
  Completeness3 <- apply(TraitData3[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness3 <- as.data.frame(Completeness3/nrow(TraitData3)*100)
  colnames(Completeness3) <- "Completeness"
  
  Completeness4 <- apply(TraitData4[, Traits_name], 2,  function(y) sum(!is.na(y))) 
  Completeness4 <- as.data.frame(Completeness4/nrow(TraitData4)*100)
  colnames(Completeness4) <- "Completeness"
  
  Completeness4 <- Completeness4[order(Completeness4, decreasing=FALSE), , drop=FALSE]
  X <- row.names(Completeness4)
  
  Completeness3 <- Completeness3[X,] %>% as.data.frame()
  row.names(Completeness3) <- X
  colnames(Completeness3) <- "Completeness"
  
  Completeness2 <- Completeness2[X,] %>% as.data.frame()
  row.names(Completeness2) <- X
  colnames(Completeness2) <- "Completeness"
  
  Completeness1 <- Completeness1[X,] %>% as.data.frame()
  row.names(Completeness1) <- X
  colnames(Completeness1) <- "Completeness"
  
  barplot(Completeness4$Completeness, horiz = TRUE, 
          names.arg=rownames(Completeness1), las=1, col="red", main=Main, xlab="% Coverage", xlim=x_lim)
  
  barplot(Completeness3$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="black")
  
  barplot(Completeness2$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="blue")
  
  barplot(Completeness1$Completeness, horiz=TRUE, add=T, names.arg=rownames(Completeness1), las=1, col="grey")
  
  abline(v=100, lty="dotted")
  
}


#  #  #  #  #  #  #  #  #  #  #  #  #  #      R  u  n  s      #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 

# # # # # # # # 
#             #
#  Mammals    #
#             #
# # # # # # # # 

## Coverage before imputation
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation")
pdf(file="../Results/Plots/Coverage_mammal_compil2.pdf", width=6, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
par(mfrow=c(2,1))
Plot.Cov(TraitsMammal, Trait_names, "Trait coverage across all mammals", "% Coverage")
Plot.Cov(TraitsMammal.Predicts, Trait_names, "Trait coverage across Predicts mammals", "% Coverage")
dev.off()


## Coverage after imputations - continuous (Molina-Venegas)
pdf(file="../Results/Plots/Coverage_Predicts_mammal_continuous_imputed.pdf", width=6, height=3, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,0,2,2), oma=c(2,9,1,1))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN")
Plot.Cov.Delta2(TraitsMammal.Predicts, TraitsMammal.Predicts_Imputed, TraitsMammal.Predicts_Imputed_Randomadd, Trait_names, "PREDICTS Mammals", c(80,100))
dev.off()


## Coverage after imputations and correlations - continuous (Molina-Venegas)
pdf(file="../Results/Plots/Coverage_Predicts_mammal_continuous_imputed.pdf", width=6, height=3, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,0,2,2), oma=c(2,9,1,1))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN")
Plot.Cov.Delta3(TraitsMammal.Predicts, TraitsMammal.Predicts_Imputed, TraitsMammal.Predicts_Imputed_Randomadd, TraitsMammals.Predicts_imputed_correlated,
                Trait_names, "PREDICTS Mammals", c(80,100))
dev.off()


pdf(file="../Results/Plots/Coverage_Predicts_mammal_continuous_imputed.pdf", width=6, height=3, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,0,2,2), oma=c(2,9,1,1))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation", "Diet_breadth")
Plot.Cov.Delta3(TraitsMammal.Predicts, TraitsMammal.Predicts_Imputed, TraitsMammal.Predicts_Imputed_Randomadd, TraitsMammals.Predicts_imputed_correlated,
                Trait_names, "PREDICTS Mammals", c(80,100))
dev.off()



# # # # # # # # 
#             #
#  Amphibians #
#             #
# # # # # # # # 


## 1.Coverage before imputation

## Basic coverage
Trait_names <- c("Body_mass_g", "Litter_size", "Maturity_d", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation")
pdf(file="../Results/Plots/Amphibian_coverage/Coverage_amphibians_basic.pdf", width=6, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
par(mfrow=c(2,1))
Plot.Cov.Basic(TraitsAmphibian_basic, Trait_names, "Basic trait coverage across amphibians","% Coverage")
Plot.Cov.Basic(TraitsAmphibian.Predicts_basic, Trait_names, "Basic trait coverage across Predicts amphibians", "% Coverage")
dev.off()

## Coverage after correlations
Trait_names <- c("Body_mass_g", "Litter_size", "Maturity_d", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation")
pdf(file="../Results/Plots/Amphibian_coverage/Coverage_amphibians_correlated.pdf", width=6, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
par(mfrow=c(2,1))
Plot.Cov(TraitsAmphibian_basic, TraitsAmphibian_correlated, Trait_names, "Trait coverage across amphibians after correlations on continuous traits")
Plot.Cov(TraitsAmphibian.Predicts_basic, TraitsAmphibian.Predicts_correlated, Trait_names, "Trait coverage across Predicts amphibians after correlations on continuous traits")
dev.off()


# ## 2.Coverage after imputations - continuous (Molina-Venegas) without random addition
# Trait_names <- c("Body_mass_g", "Litter_size", "Maturity_d", "Range_size_m2", "Habitat_breadth_IUCN")
# pdf(file="../Results/Plots/Coverage_Predicts_amphibians_continuous_imputed.pdf", width=6, height=3, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
# Plot.Cov.Delta(TraitsAmphibian.Predicts_correlated, TraitsAmphibian.Predicts_Imputed,Trait_names, "PREDICTS Amphibians")
# dev.off()

## 2.Coverage after imputations - continuous (Molina-Venegas) after random addition
Trait_names <- c("Body_mass_g", "Litter_size", "Maturity_d", "Range_size_m2", "Habitat_breadth_IUCN")
pdf(file="../Results/Plots/Amphibian_coverage/Coverage_Predicts_amphibians_correlated_imputed.pdf", width=6, height=3, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
Plot.Cov.Delta3(TraitsAmphibian.Predicts_basic,
                TraitsAmphibian.Predicts_correlated,
                TraitsAmphibian.Predicts_Imputed,
                TraitsAmphibian.Predicts_Imputed_Randomadd,
                Trait_names, "PREDICTS Amphibians", c(0,100))
dev.off()


# # # # # # # # 
#             #
#    Birds    #
#             #
# # # # # # # # 


Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation")
pdf(file="../Results/Plots/Coverage_birds_compil2.pdf", width=6, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
par(mfrow=c(2,1))
Plot.Cov(TraitsBird, Trait_names, "Trait coverage across birds \n before imputations", "% Coverage")
Plot.Cov(TraitsBird.Predicts, Trait_names, "Trait coverage across Predicts birds", "% Coverage")
dev.off()


# # # # # # # # 
#             #
#  Reptiles   #
#             #
# # # # # # # # 


## Basic compil coverage
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation")
pdf(file="../Results/Plots/Coverage_reptiles_compil.pdf", width=6, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,2,2,2), oma=c(2,7,1,1))
par(mfrow=c(2,1))
Plot.Cov(TraitsReptile, Trait_names, "Trait coverage across reptiles \n before imputations", "% Coverage")
Plot.Cov(TraitsReptile.Predicts, Trait_names, "Trait coverage across Predicts reptiles", "% Coverage")
dev.off()

## Coverage after imputations on continuous traits
Plot.Cov.Delta(TraitsReptile.Predicts, TraitsReptile.Predicts_Imputed, Trait_names, "Bla")



## All four classes -- basic compilation
pdf(file="../Results/Plots/Basic_coverage_PREDICTS.pdf", width=10, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(2,3,1,1))
par(mfrow=c(2,2))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Trophic_level", "Diet_breadth")
Plot.Cov.Basic(TraitsMammal.Predicts, Trait_names, "Mammals", "% Coverage")
Plot.Cov.Basic(TraitsBird.Predicts, Trait_names, "Birds", "% Coverage")
Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d","Trophic_level")
Plot.Cov.Basic(TraitsReptile.Predicts_basic, Trait_names, "Reptiles", "% Coverage")
Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d", "Trophic_level", "Diet_breadth")
Plot.Cov.Basic(TraitsAmphibian.Predicts_basic, Trait_names, "Amphibians", "% Coverage")
dev.off()


## All four classes -- after correlation
pdf(file="../Results/Plots/Basic_coverage_correlated_PREDICTS.pdf", width=10, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(2,3,1,1))
par(mfrow=c(2,2))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Trophic_level", "Diet_breadth")
Plot.Cov.Basic(TraitsMammal.Predicts, Trait_names, "Mammals", "% Coverage")
Plot.Cov.Basic(TraitsBird.Predicts, Trait_names, "Birds", "% Coverage")
Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d","Trophic_level")
Plot.Cov(TraitsReptile.Predicts_basic, TraitsReptile.Predicts_correlated, Trait_names, "Reptiles")
Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d", "Trophic_level", "Diet_breadth")
Plot.Cov(TraitsAmphibian.Predicts_basic, TraitsAmphibian.Predicts_correlated, Trait_names, "Amphibians")
dev.off()



## All four classes -- after correlation
pdf(file="../Results/Plots/Basic_coverage_correlated_PREDICTS.pdf", width=10, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(2,3,1,1))
par(mfrow=c(2,2))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Trophic_level", "Diet_breadth")
Plot.Cov.Basic(TraitsMammal.Predicts, Trait_names, "Mammals", "% Coverage")
Plot.Cov.Basic(TraitsBird.Predicts, Trait_names, "Birds", "% Coverage")
Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d","Trophic_level")
Plot.Cov(TraitsReptile.Predicts_basic, TraitsReptile.Predicts_correlated, Trait_names, "Reptiles")
Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d", "Trophic_level", "Diet_breadth")
Plot.Cov(TraitsAmphibian.Predicts_basic, TraitsAmphibian.Predicts_correlated, Trait_names, "Amphibians")
dev.off()



## All four classes -- after correlation and imputations
pdf(file="../Results/Plots/Coverage_PREDICTS_imputed_correlated.pdf", width=10, height=7, family="Times", pointsize=11)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,7,2,2), oma=c(2,3,1,1))
par(mfrow=c(2,2))
Trait_names <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Trophic_level", "Diet_breadth")
Plot.Cov.Delta2(TraitsMammal.Predicts,
                TraitsMammal.Predicts_Imputed,
                TraitsMammal.Predicts_Imputed_Randomadd,
                Trait_names, "PREDICTS Mammals", c(0,100))

Plot.Cov.Basic(TraitsBird.Predicts, Trait_names, "PREDICTS Birds", "% Coverage")

Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d","Trophic_level")
Plot.Cov.Delta3(TraitsReptile.Predicts_basic,
                TraitsReptile.Predicts_correlated,
                TraitsReptile.Predicts_Imputed,
                TraitsReptile.Predicts_Imputed,
                Trait_names, "PREDICTS Reptiles", c(0,100))

Trait_names <- c("Body_mass_g", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Specialisation", "Maturity_d", "Max_longevity_d", "Trophic_level", "Diet_breadth")
Plot.Cov.Delta3(TraitsAmphibian.Predicts_basic,
                TraitsAmphibian.Predicts_correlated,
                TraitsAmphibian.Predicts_Imputed,
                TraitsAmphibian.Predicts_Imputed_Randomadd,
                Trait_names, "PREDICTS Amphibians", c(0,100))
dev.off()


