############################################################
# AIM: PHYLOGENETIC SIGNAL IN CONTINUOUS TRAITS, and IMPUTATIONS IF > 0.9
# given a trait dataset and a phylogeny.


## Preamble
X <- c("dplyr", "phytools", "picante", "adephylo", "geiger", "phylobase")
invisible(lapply(X, library, character.only=TRUE)); rm(X)

source("../Code/1.Sourced_v.2/2.Phylogenetic_imputations_functions_Molina-Venegas.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Load and prepare phylogenies; data obtained from www.biodiversitycenter.org/ttol


## Mammals; downloaded 06 July 2018
Mammal_phylo <- read.newick("../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk")
Mammal_phylo$tip.label[Mammal_phylo$tip.label=="Micoureus_regina"] <- "Marmosa_regina"
Mammal_phylo$tip.label[Mammal_phylo$tip.label=="Thylamys_macrura"] <- "Thylamys_macrurus"
Mammal_phylo$tip.label[Mammal_phylo$tip.label=="Herpailurus_yaguarondi"] <- "Herpailurus_yagouaroundi"
Mammal_phylo$tip.label[Mammal_phylo$tip.label=="Uncia_uncia"] <- "Panthera_uncia"
Mammal_phylo$tip.label[Mammal_phylo$tip.label=="Galerella_sanguinea"] <- "Herpestes_sanguineus"

## Adding the species that do not match at random in the phylogeny
# For mammals, the species are: 
# "Cryptonanus agricolai"      "Falsistrellus petersi"    "Falsistrellus tasmaniensis"   "Haeromys margarettae"    "Smutsia gigantea"      
# Placing into synonymous genuses:
Mammal_phylo_Randomadd <- Mammal_phylo
Mammal_phylo_Randomadd <- add.species.to.genus(Mammal_phylo_Randomadd, "Smutsia gigantea", genus="Manis", where="random")
Mammal_phylo_Randomadd <- add.species.to.genus(Mammal_phylo_Randomadd, "Falsistrellus tasmaniensis", genus="Pipistrellus", where="random")
Mammal_phylo_Randomadd <- add.species.to.genus(Mammal_phylo_Randomadd, "Falsistrellus petersi", genus="Pipistrellus", where="random")
Mammal_phylo_Randomadd <- add.species.to.genus(Mammal_phylo_Randomadd, "Cryptonanus agricolai", genus="Gracilinanus", where="random")
# Mammal_phylo_Randomadd <- add.species.to.genus(Mammal_phylo, "Haeromys margarettae", where="random")

Mammal_phylo <- .Format_tiplabels(Mammal_phylo)
Mammal_phylo_Randomadd <- .Format_tiplabels(Mammal_phylo_Randomadd)


## Amphibians
Amphibian_phylo <- read.newick("../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk")
Amphibian_phylo$tip.label[Amphibian_phylo$tip.label=="Rhamphophryne_festae"] <- "Rhinella_festae"

## Adding genuses that don't match the phylogeny to it -- for Body_mass_g imputations
Amphibian_phylo_Randomadd <- Amphibian_phylo
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Sclerophrys togoensis", genus="Amietophrynus", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Pristimantis uisae", genus="Pristimantis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Pristimantis kelephas", genus="Pristimantis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Incilius cristatus", genus="Incilius", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Hypodactylus mantipus", genus="Hypodactylus", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Boophis burgeri", genus="Boophis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Boophis guibei", genus="Boophis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Boophis lichenoides", genus="Boophis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Boophis reticulatus", genus="Boophis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Boophis rufioculis", genus="Boophis", where="random")
Amphibian_phylo_Randomadd <- add.species.to.genus(Amphibian_phylo_Randomadd, "Duellmanohyla eutisanota", genus="Duellmanohyla", where="random")

## Adding genuses that don't match the phylogeny to it -- for Habitat breadth imputations (only)
Amphibian_phylo_Randomadd_HB <- Amphibian_phylo
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Duellmanohyla eutisanota", genus="Duellmanohyla", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Hylarana occidentalis", genus="Hylarana", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Hypopachus pictiventris", genus="Hypopachus", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Hypopachus ustus", genus="Hypopachus", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Lithobates brownorum", genus="Lithobates", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Megophrys brachykolos", genus="Megophrys", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Phrynoidis aspera", genus="Phrynoidis", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Pristimantis kelephas", genus="Pristimantis", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Sclerophrys togoensis", genus="Amietophrynus", where="random")
Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Raorchestes viridis", genus="Raorchestes", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Pseudophilautus alto", genus="Pseudophilautus", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Pseudophilautus sarasinorum", genus="Pseudophilautus", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Pseudophilautus silus", genus="Pseudophilautus", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Pseudophilautus sordidus", genus="Pseudophilautus", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Phrynobatrachus intermedius", genus="Phrynobatrachus", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Agalychnis buckleyi", genus="Agalychnis", where="random")
# Amphibian_phylo_Randomadd_HB <- add.species.to.genus(Amphibian_phylo_Randomadd_HB, "Boophis calcaratus", genus="Boophis", where="random")

Amphibian_phylo <- .Format_tiplabels(Amphibian_phylo)
Amphibian_phylo_Randomadd <- .Format_tiplabels(Amphibian_phylo_Randomadd)
Amphibian_phylo_Randomadd_HB <- .Format_tiplabels(Amphibian_phylo_Randomadd_HB)


## Birds
Bird_phylo <- read.newick("../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk")
Bird_phylo <- .Format_tiplabels(Bird_phylo)


## Reptiles
Reptile_phylo <- read.newick("../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk")
Reptile_phylo <- .Format_tiplabels(Reptile_phylo)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Load trait data

TraitsMammal.Predicts <- read.csv("../Results/1.Traits_before_imputations/TraitsMammalsPredicts.csv")
TraitsAmphibian.Predicts <- read.csv("../Results/1.Traits_before_imputations/Amphibians/TraitsAmphibiansPredicts_correlated.csv")
TraitsBird.Predicts <- read.csv("../Results/1.Traits_before_imputations/TraitsBirdPredicts.csv")
TraitsReptile.Predicts <- read.csv("../Results/1.Traits_before_imputations/TraitsReptilePredicts.csv")

# # # # Check NAs and species names that don't match with the phylogeny
# TraitNA <- TraitsAmphibian.Predicts$Best_guess_binomial[is.na(TraitsAmphibian.Predicts$Habitat_breadth_IUCN)]
# X <- setdiff(TraitNA, Amphibian_phylo$tip.label)
# Bird_phylo$tip.label[grepl(" curvirostris", Bird_phylo$tip.label)]

#
# setdiff(TraitNA, Amphibian_phylo_Randomadd$tip.label)
# TraitsAmphibian.Predicts$Best_guess_binomial[grepl("Blommersia", TraitsAmphibian.Predicts$Best_guess_binomial)]
# TraitsAmphibian$Body_length_mm[grepl("Hypodactylus mantipus", TraitsAmphibian$Binomial_name)]
# setdiff(TraitsAmphibian.Predicts$Best_guess_binomial, Amphibian_phylo$tip.label)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Run the function on continuous traits -- imputations 

## Mammals
TraitsMammalContinuous <- c("Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN")

## Mammals -- Predicts subset
Imputed_mammals_phylopars <- Predictions_Molina_venegas(TraitsMammal.Predicts, Mammal_phylo, TraitsMammalContinuous)
Imputed_mammals_phylopars_Randomadd <- Predictions_Molina_venegas(TraitsMammal.Predicts, Mammal_phylo_Randomadd, TraitsMammalContinuous)

TraitsMammal.Predicts_Imputed <- Add_Results(Imputed_mammals_phylopars, TraitsMammal.Predicts)
TraitsMammal.Predicts_Imputed_Randomadd <- Add_Results(Imputed_mammals_phylopars_Randomadd, TraitsMammal.Predicts)


## Amphibians
TraitsAmphibianContinuous <- c("Body_mass_g", "Maturity_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN")

Imputed_amphibians_phylopars <- Predictions_Molina_venegas(TraitsAmphibian.Predicts, Amphibian_phylo, TraitsAmphibianContinuous)
# Body_mass
Imputed_amphibians_phylopars_Randomadd <- Predictions_Molina_venegas(TraitsAmphibian.Predicts, Amphibian_phylo_Randomadd, TraitsAmphibianContinuous)
# Habitat_breadth
Imputed_amphibians_phylopars_Randomadd_HB <- Predictions_Molina_venegas(TraitsAmphibian.Predicts, Amphibian_phylo_Randomadd_HB, TraitsAmphibianContinuous)

## Adding results
TraitsAmphibian.Predicts_Imputed <- Add_Results(Imputed_amphibians_phylopars, TraitsAmphibian.Predicts)
TraitsAmphibian.Predicts_Imputed_Randomadd <- Add_Results(Imputed_amphibians_phylopars_Randomadd, TraitsAmphibian.Predicts)
TraitsAmphibian.Predicts_Imputed_Randomadd_HB <- Add_Results(Imputed_amphibians_phylopars_Randomadd_HB, TraitsAmphibian.Predicts)
TraitsAmphibian.Predicts_Imputed_Randomadd$Habitat_breadth_IUCN <- TraitsAmphibian.Predicts_Imputed_Randomadd_HB$Habitat_breadth_IUCN


## Birds 
TraitsBirdContinuous <- c("Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Body_mass_g")
Imputed_birds_phylopars <- Predictions_Molina_venegas(TraitsBird.Predicts, Bird_phylo, TraitsBirdContinuous)
TraitsBird.Predicts_Imputed <- Add_Results(Imputed_birds_phylopars, TraitsBird.Predicts)


## Reptiles
TraitsReptileContinuous <- c("Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Body_mass_g")
Imputed_reptiles_phylopars <- Predictions_Molina_venegas(TraitsReptile.Predicts, Reptile_phylo, TraitsReptileContinuous)
TraitsReptile.Predicts_Imputed <- Add_Results(Imputed_reptiles_phylopars, TraitsReptile.Predicts)

# ## Prepare to export phylogenetic signal
# Phylo_signal_lambda_continuous <- data.frame(Imputed_mammals_phylopars[[1]][[1]], Imputed_mammals_phylopars[[2]][[1]],
#                                              Imputed_mammals_phylopars[[3]][[1]], Imputed_mammals_phylopars[[4]][[1]])
# colnames(Phylo_signal_lambda_continuous) <- c("Generation_length", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN")
# 
# Phylo_signal_lambda_continuous[2,1] <- Imputed_amphibians_phylopars[[2]][[1]]
# Phylo_signal_lambda_continuous[2,2] <- Imputed_amphibians_phylopars[[3]][[1]]
# Phylo_signal_lambda_continuous[2,3] <- Imputed_amphibians_phylopars[[4]][[1]]
# Phylo_signal_lambda_continuous[2,4] <- Imputed_amphibians_phylopars[[5]][[1]]
# Phylo_signal_lambda_continuous$Body_measure[2] <- Imputed_amphibians_phylopars[[1]][[1]]
# 
# row.names(Phylo_signal_lambda_continuous) <- c("Mammals", "Amphibians")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Save files
write.csv(TraitsMammal.Predicts_Imputed, "../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsMammal.Predicts_ContImputed.csv", row.names=F)
write.csv(TraitsMammal.Predicts_Imputed_Randomadd, "../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsMammal.Predicts_ContImputed_Randomadd.csv", row.names=F)

write.csv(TraitsAmphibian.Predicts_Imputed, "../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsAmphibians.Predicts_ContImputed.csv", row.names=F)
write.csv(TraitsAmphibian.Predicts_Imputed_Randomadd, "../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsAmphibians.Predicts_ContImputed_Randomadd.csv", row.names=F)

write.csv(TraitsBird.Predicts_Imputed, "../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsBird.Predicts_ContImputed.csv", row.names=F)

write.csv(TraitsReptile.Predicts_Imputed, "../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsReptile.Predicts_ContImputed.csv", row.names=F)


write.csv(Phylo_signal_lambda_continuous, "../Results/TraitsContinuous_after_imputations_Molina_Venegas/Phylogenetic_signals.csv", row.names=F)






# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# ##          Check for phylogenetic signal in continuous traits with package phytools (K, lambda)          ##
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# 
# 
# 
# 
# ## Use match.phylo.taxa to compare taxa present in phylogeny with trait data, 
# ## pruning and sorting the data to match one another for subsequent analysis
# 
# ## In Traits: species names as row names for subsequent analyses
# row.names(TraitsMammal.Predicts) <- TraitsMammal.Predicts$Best_guess_binomial
# 
# ## Match species by binomial names and prune when species do not match 
# ## Row names of trait data need to match phy$tip.label
# ## The trait datatframe returned (pruned) has rows in the same order as phy$tip.label
# Prune_Taxa <- match.phylo.data(Mammal_phylo, TraitsMammal.Predicts)
# 
# 
# Traits_match <- Prune_Taxa$data
# Cont_traits_match <- apply(Cont_traits_match, 2, as.character) %>% as.data.frame()
# Cont_traits_match[, c(2:4)] <- apply(Cont_traits_match[, c(2:4)], 2, as.numeric)
# row.names(Cont_traits_match) <- Cont_traits_match$Binomial_name
# 
# Mammal_phylo_match <- Prune_Taxa$phy
# 
# rm(Prune_Taxa)
# 
# 
# 
# 
# 
# # # Intersect by binomial names and subset
# # Y <- intersect(Mammal_phylo$tip.label, Traits$Binomial_name)
# # Body_mass_phylo <- Traits %>% subset(Binomial_name %in% Y) %>%
# #   select(c(Binomial_name, Body_mass_g))
# # 
# # # Add species that are only in the phylogeny and add NA for trait
# # To_add <- Mammal_phylo$tip.label[!(Mammal_phylo$tip.label %in% Y)] %>% as.data.frame() 
# # colnames(To_add) <- "tip.label"
# # To_add$Body_mass_g <- NA
# # 
# # Body_mass_phylo <- rbind(Body_mass_phylo, To_add)
# 
# 
# ## Measure phylogenetic signal using phylosig (package phytools): either "K" or "lambda"
# K.BM <- phylosig(Mammal_phylo_match, Cont_traits_match$Body_mass_g, method="K", test=T, nsim=1000, se=NULL, start=NULL, control=list())
# Lambda.BM <- phylosig(Mammal_phylo_match, Cont_traits_match$Body_mass_g, method="lambda", test=T, nsim=1000, se=NULL, start=NULL, control=list())
# 
# 
# 
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# ##               Check for phylogenetic signal in continuous traits (with package picante )               ##
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# 
# 
# 
# ## USING BLOMBERG'S K 
# 
# ## Test for phylogenetic signal in all continuous traits
# Phy_Signal <- multiPhylosignal(Prune_Taxa$data, Prune_Taxa$phy)
# 
# ## Test for phylogenetic signal in a given trait using the phylosignal function: 
# ## /!\ Need to order the trait vector as in phy$tip.label
# 
# Phy_Signal_BM <- phylosignal(Prune_Taxa$data["Body_mass_g"],  Prune_Taxa$phy)

