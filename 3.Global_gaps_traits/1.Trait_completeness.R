## Trait completeness

## Preamble
X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ggplot2", "ggpubr", "plyr", "reshape", "reshape2", "mi", "cowplot")
lapply(X, library, character.only=TRUE); rm(X)

# # # # # # # 

## Functions

## Count target traits (same traits across all classes, with the GL-L equivalence as well as the BM-BL equivalence)
Count_traits <- function(TraitDF, Traits) {
  TraitDF$CountTraits <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y)))
  TraitDF$PropTraits <- TraitDF$CountTraits / length(Traits) * 100
  return(TraitDF[, c("Best_guess_binomial", "CountTraits", "PropTraits")])
}

## Plotting percent information across species ("trait filling" or histogram or number of traits sampled across species)
Completeness_plot <- function(TraitDFC, TraitDFU, FontSize, Hist_or_Density=c("histogram", "density"), BW) {
  
  GGPoptions <- theme_classic()+ theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=FontSize, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=FontSize), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=FontSize),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=FontSize)) 
  
  TraitDFC$Taxonomy <- "Corrected"
  TraitDFU$Taxonomy <- "Uncorrected"
  
  Results <- rbind(TraitDFC, TraitDFU)
  Results$Taxonomy <- factor(Results$Taxonomy, levels = c("Uncorrected", "Corrected"))

  medN <- ddply(Results, "Taxonomy", summarise, grp.median=median(CountTraits))
  medP <- ddply(Results, "Taxonomy", summarise, grp.median=median(PropTraits))
  print(paste("median count: ", medN))
  print(paste("median proportion: ", medP))
  
  if(Hist_or_Density=="histogram") {
    
    p <-ggplot(Results, aes(CountTraits, fill=Taxonomy)) +
       geom_histogram(alpha=0.8,position="dodge2", stat="count") + GGPoptions +
      geom_vline(data=medN[medN$Taxonomy=="Uncorrected", ], aes(xintercept=grp.median), linetype="dashed", col="red", alpha=0.7) +
      geom_vline(data=medN[medN$Taxonomy=="Corrected", ], aes(xintercept=grp.median), linetype="dashed", col="blue", alpha=0.7) +
      xlab("Number of sampled traits") + ylab("Count")
    
  }
  
  if(Hist_or_Density=="density") {
    p <- ggplot(Results, aes(PropTraits, fill=Taxonomy)) +
      geom_density(alpha=0.5, adjust=BW) + GGPoptions +
      geom_vline(data=medP[medP$Taxonomy=="Uncorrected", ], aes(xintercept=grp.median), linetype="dashed", col="red", alpha=0.7) +
      geom_vline(data=medP[medP$Taxonomy=="Corrected", ], aes(xintercept=grp.median), linetype="dashed", col="blue", alpha=0.7) +
      xlab("Completeness (%)") + ylab("Density of species")
  }
  
  return(p)
  
}

# For PREDICTS species only
GetPredictsSpecies <- function(Predicts, DF, Class){
  Psp <- unique(Predicts$Best_guess_binomial[Predicts$Class==Class])
  DF <- DF %>% filter(Best_guess_binomial %in% Psp)
  return(DF)}

# To test whether median completeness differs across classes: Kruskall-Wallis test 
KTest <- function(M, B, R, A, Pairwise){
  
  ResultsM <- M %>%
    mutate(Class="Mammalia")
  
  ResultsB <- B %>%
    mutate(Class="Birds")
  
  ResultsR <- R %>%
    mutate(Class="Reptiles")
  
  ResultsA <- A %>%
    mutate(Class="Amphibians")
  
  Results <- rbind(ResultsB,ResultsA, ResultsR, ResultsM)
  
  Results$Class <- as.factor(Results$Class)
  
  if(!Pairwise) {
    return(kruskal.test(Results$CountTraits, Results$Class))
  }
  
  else{
    return(pairwise.wilcox.test(Results$CountTraits, Results$Class,
                                p.adjust.method = "BH"))
  }
  
}


# # # # # # # #

## Load trait data, with and without taxonomic correction
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
UMammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
UBirds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
UAmphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")
UReptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

Predicts_C <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")
Predicts_U <- readRDS("../../Data/Vertebrates.rds")

## Target traits:

# equivalence body mass/body length and longevity/generation length
Common_traits <- c("Litter_size", 
                   "Habitat_breadth_IUCN",
                   "Specialisation",
                   "Trophic_level",
                   "Diel_activity")

TraitsReptiles <- c(Common_traits, "Adult_svl_cm", "Longevity_d")
TraitsAmphibians <- c(Common_traits, "Body_length_mm", "Longevity_d")
TraitsMammals <- c(Common_traits, "Body_mass_g", "Generation_length_d")
TraitsBirds <- c(Common_traits, "Body_mass_g", "Generation_length_d")

# For all species
T_Mammals_C <- Count_traits(Mammals, TraitsMammals)
T_Mammals_U <- Count_traits(UMammals, TraitsMammals)

T_Birds_C <- Count_traits(Birds, TraitsBirds)
T_Birds_U <- Count_traits(UBirds, TraitsBirds)

T_Reptiles_C <- Count_traits(Reptiles, TraitsReptiles)
T_Reptiles_U <- Count_traits(UReptiles, TraitsReptiles)

T_Amphibians_C <- Count_traits(Amphibians, TraitsAmphibians)
T_Amphibians_U <- Count_traits(UAmphibians, TraitsAmphibians)

# plot as histogram
p <- ggarrange(
Completeness_plot(T_Mammals_C, T_Mammals_U, 12, "histogram", 2) + xlab(""),
Completeness_plot(T_Birds_C, T_Birds_U, 12, "histogram", 2)+ xlab("") + ylab(""),
Completeness_plot(T_Reptiles_C, T_Reptiles_U, 12, "histogram", 2),
Completeness_plot(T_Amphibians_C, T_Amphibians_U, 12, "histogram", 2) + ylab(""),
common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/All_species/TARGET_traits_count.pdf", width=6.5, height=5)

# plot as density
p <- ggarrange(
  Completeness_plot(T_Mammals_C, T_Mammals_U, 12, "density", 4) + xlab(""),
  Completeness_plot(T_Birds_C, T_Birds_U, 12, "density", 4.5)+ xlab("") + ylab(""),
  Completeness_plot(T_Reptiles_C, T_Reptiles_U, 12, "density", 4),
  Completeness_plot(T_Amphibians_C, T_Amphibians_U, 12, "density", 4.5) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/All_species/TARGET_traits_density.pdf", width=6.5, height=5)



# For PREDICTS species only
Predicts_Mammals_TC <- GetPredictsSpecies(Predicts_C, T_Mammals_C, "Mammalia")
Predicts_Mammals_TU <- GetPredictsSpecies(Predicts_U, T_Mammals_U, "Mammalia")
Predicts_Birds_TC <- GetPredictsSpecies(Predicts_C, T_Birds_C, "Aves")
Predicts_Birds_TU <- GetPredictsSpecies(Predicts_U, T_Birds_U, "Aves")
Predicts_Reptiles_TC <- GetPredictsSpecies(Predicts_C, T_Reptiles_C, "Reptilia")
Predicts_Reptiles_TU <- GetPredictsSpecies(Predicts_U, T_Reptiles_U, "Reptilia")
Predicts_Amphibians_TC <- GetPredictsSpecies(Predicts_C, T_Amphibians_C, "Amphibia")
Predicts_Amphibians_TU <- GetPredictsSpecies(Predicts_U, T_Amphibians_U, "Amphibia")

# plot as histogram
p <- ggarrange(
  Completeness_plot(Predicts_Mammals_TC, Predicts_Mammals_TU, 12, "histogram", 2) + xlab(""),
  Completeness_plot(Predicts_Birds_TC, Predicts_Birds_TU, 12, "histogram", 2)+ xlab("") + ylab(""),
  Completeness_plot(Predicts_Reptiles_TC, Predicts_Reptiles_TU, 12, "histogram", 2),
  Completeness_plot(Predicts_Amphibians_TC, Predicts_Amphibians_TU, 12, "histogram", 2) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/PREDICTS_species/TARGET_traits_count.pdf", width=6.5, height=5)

# plot as density
p <- ggarrange(
  Completeness_plot(Predicts_Mammals_TC, Predicts_Mammals_TU, 12, "density", 2) + xlab(""),
  Completeness_plot(Predicts_Birds_TC, Predicts_Birds_TU, 12, "density", 4)+ xlab("") + ylab(""),
  Completeness_plot(Predicts_Reptiles_TC, Predicts_Reptiles_TU, 12, "density", 2),
  Completeness_plot(Predicts_Amphibians_TC, Predicts_Amphibians_TU, 12, "density", 2) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/PREDICTS_species/TARGET_traits_density.pdf", width=6.5, height=5)


## Test whether taxonomic corrections increase median trait completeness for target traits
wilcox.test(x=T_Mammals_U$CountTraits, y=T_Mammals_C$CountTraits, alternative = "less") # yes
wilcox.test(x=T_Birds_U$CountTraits, y=T_Birds_C$CountTraits, alternative = "less") # yes
wilcox.test(x=T_Reptiles_U$CountTraits, y=T_Reptiles_C$CountTraits, alternative = "less") # no
wilcox.test(x=T_Amphibians_U$CountTraits, y=T_Amphibians_C$CountTraits, alternative = "less") # yes

## Test whether class affects trait completeness
KTest(M=T_Mammals_C, 
      R=T_Reptiles_C, 
      B=T_Birds_C, 
      A=T_Amphibians_C, 
      Pairwise = TRUE)

KTest(M=T_Mammals_C, 
      R=T_Reptiles_C, 
      B=T_Birds_C, 
      A=T_Amphibians_C, 
      Pairwise = FALSE)


## PREDICTORS traits - same as above, but with all traits included as predictors in the RF algorithm.
Traits <- c("Body_mass_g",
            "Longevity_d",
            "Litter_size", 
            "Range_size_m2", 
            "Habitat_breadth_IUCN",
            "Specialisation",
            "Trophic_level",
            "Diel_activity",
            "Primary_diet",
            "Diet_breadth")

TraitsReptiles <- c(Traits[1:8], "Adult_svl_cm", "Maturity_d")
TraitsAmphibians <- c(Traits, "Body_length_mm")
TraitsMammals <- c(Traits, "Adult_svl_cm", "Generation_length_d")
TraitsBirds <- c(Traits, "Generation_length_d")

# For all species
P_Mammals_C <- Count_traits(Mammals, TraitsMammals)
P_Mammals_U <- Count_traits(UMammals, TraitsMammals)

P_Birds_C <- Count_traits(Birds, TraitsBirds)
P_Birds_U <- Count_traits(UBirds, TraitsBirds)

P_Reptiles_C <- Count_traits(Reptiles, TraitsReptiles)
P_Reptiles_U <- Count_traits(UReptiles, TraitsReptiles)

P_Amphibians_C <- Count_traits(Amphibians, TraitsAmphibians)
P_Amphibians_U <- Count_traits(UAmphibians, TraitsAmphibians)

# plot as histogram
p <- ggarrange(
  Completeness_plot(P_Mammals_C, P_Mammals_U, 12, "histogram", 2) + xlab(""),
  Completeness_plot(P_Birds_C, P_Birds_U, 12, "histogram", 2)+ xlab("") + ylab(""),
  Completeness_plot(P_Reptiles_C, P_Reptiles_U, 12, "histogram", 2),
  Completeness_plot(P_Amphibians_C, P_Amphibians_U, 12, "histogram", 2) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/All_species/PREDICTOR_traits_count.pdf", width=6.5, height=5)

# plot as density
p <- ggarrange(
  Completeness_plot(P_Mammals_C, P_Mammals_U, 12, "density", 4) + xlab(""),
  Completeness_plot(P_Birds_C, P_Birds_U, 12, "density", 4.5)+ xlab("") + ylab(""),
  Completeness_plot(P_Reptiles_C, P_Reptiles_U, 12, "density", 4),
  Completeness_plot(P_Amphibians_C, P_Amphibians_U, 12, "density", 4.5) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/All_species/PREDICTOR_traits_density.pdf", width=6.5, height=5)



# For PREDICTS species only
Predicts_Mammals_PC <- GetPredictsSpecies(Predicts_C, P_Mammals_C, "Mammalia")
Predicts_Mammals_PU <- GetPredictsSpecies(Predicts_U, P_Mammals_U, "Mammalia")
Predicts_Birds_PC <- GetPredictsSpecies(Predicts_C, P_Birds_C, "Aves")
Predicts_Birds_PU <- GetPredictsSpecies(Predicts_U, P_Birds_U, "Aves")
Predicts_Reptiles_PC <- GetPredictsSpecies(Predicts_C, P_Reptiles_C, "Reptilia")
Predicts_Reptiles_PU <- GetPredictsSpecies(Predicts_U, P_Reptiles_U, "Reptilia")
Predicts_Amphibians_PC <- GetPredictsSpecies(Predicts_C, P_Amphibians_C, "Amphibia")
Predicts_Amphibians_PU <- GetPredictsSpecies(Predicts_U, P_Amphibians_U, "Amphibia")

# plot as histogram
p <- ggarrange(
  Completeness_plot(Predicts_Mammals_PC, Predicts_Mammals_PU, 12, "histogram", 2) + xlab(""),
  Completeness_plot(Predicts_Birds_PC, Predicts_Birds_PU, 12, "histogram", 2)+ xlab("") + ylab(""),
  Completeness_plot(Predicts_Reptiles_PC, Predicts_Reptiles_PU, 12, "histogram", 2),
  Completeness_plot(Predicts_Amphibians_PC, Predicts_Amphibians_PU, 12, "histogram", 2) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/PREDICTS_species/PREDICTOR_traits_count.pdf", width=6.5, height=5)

# plot as density
p <- ggarrange(
  Completeness_plot(Predicts_Mammals_PC, Predicts_Mammals_PU, 12, "density", 2) + xlab(""),
  Completeness_plot(Predicts_Birds_PC, Predicts_Birds_PU, 12, "density", 4)+ xlab("") + ylab(""),
  Completeness_plot(Predicts_Reptiles_PC, Predicts_Reptiles_PU, 12, "density", 2),
  Completeness_plot(Predicts_Amphibians_PC, Predicts_Amphibians_PU, 12, "density", 2) + ylab(""),
  common.legend = TRUE, labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-21,-21), vjust=2, font.label = list(family="serif")
)
ggsave(p, file="../../Results/Plots/Trait_missing_values/PREDICTS_species/PREDICTOR_traits_density.pdf", width=6.5, height=5)


## Test whether taxonomic corrections increase median trait completeness for target traits
wilcox.test(x=P_Mammals_U$CountTraits, y=P_Mammals_C$CountTraits, alternative = "less") # yes
wilcox.test(x=P_Birds_U$CountTraits, y=P_Birds_C$CountTraits, alternative = "less") # yes
wilcox.test(x=P_Reptiles_U$CountTraits, y=P_Reptiles_C$CountTraits, alternative = "less") # yes
wilcox.test(x=P_Amphibians_U$CountTraits, y=P_Amphibians_C$CountTraits, alternative = "less") # yes

## Test whether class affects trait completeness
KTest(M=P_Mammals_C, 
      R=P_Reptiles_C, 
      B=P_Birds_C, 
      A=P_Amphibians_C, 
      Pairwise = TRUE)

KTest(M=P_Mammals_C, 
      R=P_Reptiles_C, 
      B=P_Birds_C, 
      A=P_Amphibians_C, 
      Pairwise = FALSE)
## class significantly affects % completeness


## Save completeness datasets (for corrected taxonomy)
CompletenessMammals_corrected <- cbind(T_Mammals_C, P_Mammals_C)[-4]
colnames(CompletenessMammals_corrected) <- c("Best_guess_binomial", "N_target", "prop_target", "N_predictors", "prop_predictors")

CompletenessBirds_corrected <- cbind(T_Birds_C, P_Birds_C)[-4]
colnames(CompletenessBirds_corrected) <- c("Best_guess_binomial", "N_target", "prop_target", "N_predictors", "prop_predictors")

CompletenessReptiles_corrected <- cbind(T_Reptiles_C, P_Reptiles_C)[-4]
colnames(CompletenessReptiles_corrected) <- c("Best_guess_binomial", "N_target", "prop_target", "N_predictors", "prop_predictors")

CompletenessAmphibians_corrected <- cbind(T_Amphibians_C, P_Amphibians_C)[-4]
colnames(CompletenessAmphibians_corrected) <- c("Best_guess_binomial", "N_target", "prop_target", "N_predictors", "prop_predictors")

write.csv(CompletenessMammals_corrected, "../../Results/3.Global_gaps_traits/Completeness_mammals.csv", row.names = FALSE)
write.csv(CompletenessBirds_corrected, "../../Results/3.Global_gaps_traits/Completeness_birds.csv", row.names = FALSE)
write.csv(CompletenessReptiles_corrected, "../../Results/3.Global_gaps_traits/Completeness_reptiles.csv", row.names = FALSE)
write.csv(CompletenessAmphibians_corrected, "../../Results/3.Global_gaps_traits/Completeness_amphibians.csv", row.names = FALSE)





