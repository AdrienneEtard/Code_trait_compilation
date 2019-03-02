## Script to assess the phylogenetic signal in Traits (continuous and categorical)

# TODO phylogenetic signal for categorical traits

## Code is parallelised
library(parallel)
library(dplyr)
library(phytools)

## START CLUSTER
Cluster <- makeCluster(detectCores())

## EXCECUTE ANY PRE PROCESSING CODE NECESSARY
clusterEvalQ(Cluster, {
  library(dplyr)
  library(phytools)
  library(picante)
  library(geiger)
  })

## Preamble

# Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}

# Function to apply for continuous traits: PhySignal
PhySignal <- function(Traitdata, Names, Phylo) {
    names(Traitdata) <- Names # names are binomial species names
    Signal <- phytools::phylosig(Phylo, Traitdata, method="lambda", test = TRUE) %>% 
      unlist
    return(Signal)
}

# # Function to apply for categorical traits: PhySignal_Cat 
# # TODO
# PhySignal_Cat <- function(Traitdata, Names, Phylo) {
#   
#   # Traitdata=Mammals[,"Trophic_level"] %>% as.vector()
#   # Names=Mammals$Best_guess_binomial
#   # Phylo=Phylo_Mammals
#   
#   ## For the current trait, match and prune phylogeny
#   names(Traitdata) <- Names
#   Match <- match.phylo.data(Phylo, Traitdata)
#   Phylo <- Match$phy
#   Trait <- Match$data
#   rm(Match)
#   
#   ## Fit different models of trait evolution and select model that best fits (model selection can be refined)
#   ER_model <- fitDiscrete(Phylo, Trait, type="discrete", model="ER")
#   ARD_model <- fitDiscrete(Phylo, Trait, type="discrete", model="ARD")
#   SYM_model <- fitDiscrete(Phylo, Trait, type="discrete", model="SYM")
#   
#   # Select model that best fits based on AICc and log likelihood
#   Model_selection <- as.data.frame(c(ER_model$opt$aicc, ARD_model$opt$aicc, SYM_model$opt$aicc))
#   colnames(Model_selection) <- "AICc"
#   row.names(Model_selection) <- c("ER", "ARD", "SYM")
# 
#   Model_selection$LogLk <- c(ER_model$opt$lnL, ARD_model$opt$lnL, SYM_model$opt$lnL)
# 
#   Model1 <- row.names(Model_selection[which.min(Model_selection$AICc),])
#   Model2 <- row.names(Model_selection[which.max(Model_selection$LogLk),])
#   
#   if(Model1!=Model2) {print("Selection on both AICc and logLk: clash!")}
#   
#   else{
#     
#     print(paste("Model selected was", Model1))
#     
#     ## Rescale tree: lambda for a "0" tree
#     Phylo_lambda_0 <- rescale(Phylo, model = "lambda", 0)
#        
#     # # To plot rescaled and non rescaled phylogenies:
#     # par(mfrow=c(1,2))
#     # plot(Phylo, edge.width = 1, show.tip.label = FALSE) #smaller tree
#     # add.scale.bar(cex = 0.7, font = 2, col = "red")
#     # plot(Phylo_lambda_0, edge.width = 1, show.tip.label = FALSE) #smaller tree
#     # add.scale.bar(cex = 0.7, font = 2, col = "red")
# 
#     ## Fit lambda to the tree
#     Lambda_trait <- fitDiscrete(Phylo, Trait, transform="lambda", type="discrete", model = Model1, niter=500)
# 
#     ## Fit lambda to the 0-rescaled tree
#     Lambda_trait_0 <- fitDiscrete(Phylo_lambda_0, Trait, transform="lambda", type="discrete", model = Model1, niter=500)
# 
#     ## Compare the outputs: compare the null model to the non null model - log likelihood and AICc difference
#     
#     # AICc difference and log likelihood chi-square test
#     AICc_diff <- abs(Lambda_trait_0$opt$aicc - Lambda_trait$opt$aicc)
#     d_logLk <- abs(2*(Lambda_trait_0$opt$lnL - Lambda_trait$opt$lnL))
#     p <- pchisq(d_logLk, 1, lower.tail=FALSE)
# 
#     if (p<0.05) {
#         cat("Model significantly different from null")
#     } else {cat("Model not significantly different from null")}
# 
#     Signal <- Lambda_trait$opt$lambda
# 
#     return(Signal)
#     }
# } 


# Read trait data
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

# Load phylogenies
Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk") %>% .Format_tiplabels() %>% compute.brlen()
Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk") %>% .Format_tiplabels() %>% compute.brlen() 
Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk") %>% .Format_tiplabels() %>% compute.brlen()
Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk") %>% .Format_tiplabels()  %>% compute.brlen()

# Traits
Continuous.Traits <- c("Body_mass_g",
                       "Longevity_d",
                       "Litter_size", 
                       "Diet_breadth",
                       "Range_size_m2", 
                       "Habitat_breadth_IUCN")

Categorical.Traits <- c("Specialisation",
                        "Trophic_level",
                        "Diel_activity",
                        "Primary_diet")

# Names
Names.Mammals <- Mammals$Best_guess_binomial
Names.Birds <- Birds$Best_guess_binomial
Names.Reptiles <- Reptiles$Best_guess_binomial
Names.Amphibians <- Amphibians$Best_guess_binomial


## Export variables in all clusters
clusterExport(cl=Cluster, list(".Format_tiplabels","PhySignal", #"PhySignal_Cat",
                               "Mammals", "Birds", "Reptiles", "Amphibians",
                               "Phylo_Mammals", "Phylo_Birds", "Phylo_Reptiles", "Phylo_Amphibians",
                               "Continuous.Traits", "Categorical.Traits",
                               "Names.Mammals", "Names.Birds", "Names.Reptiles", "Names.Amphibians"), envir=environment())


## PARALLEL CALCULATIONS WITH parApply

# Phylogenetic signal in continuous traits
Lambda_Mammals_continuous <- parApply(Cluster, Mammals[, Continuous.Traits], 2, PhySignal, Names=Names.Mammals, Phylo=Phylo_Mammals %>% as.data.frame())
Lambda_Birds_continuous <- parApply(Cluster, Birds[, Continuous.Traits], 2, PhySignal, Names=Names.Birds, Phylo=Phylo_Birds%>% as.data.frame())
Lambda_Reptiles_continuous <- parApply(Cluster, Reptiles[, Continuous.Traits], 2, PhySignal, Names=Names.Reptiles, Phylo=Phylo_Reptiles)
Lambda_Amphibians_continuous <- parApply(Cluster, Amphibians[, Continuous.Traits], 2, PhySignal, Names=Names.Amphibians, Phylo=Phylo_Amphibians)

# Phylogenetic signal in categorical traits
Lambda_Mammals_categorical <- parApply(Cluster, Mammals[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Mammals, Phylo=Phylo_Mammals)
Lambda_Birds_categorical <- parApply(Cluster, Birds[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Birds, Phylo=Phylo_Birds)
Lambda_Reptiles_categorical <- parApply(Cluster, Reptiles[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Reptiles, Phylo=Phylo_Reptiles)
Lambda_Amphibians_categorical <- parApply(Cluster, Amphibians[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Amphibians, Phylo=Phylo_Amphibians)


## Save files
write.csv(Lambda_Mammals_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousMammals.csv",row.names = FALSE)
write.csv(Lambda_Birds_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousBirds.csv",row.names = FALSE)
write.csv(Lambda_Reptiles_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousReptiles.csv",row.names = FALSE)
write.csv(Lambda_Amphibians_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousAmphibians.csv",row.names = FALSE)


## DESTROY CLUSTER
stopCluster(Cluster)



# Lambda_discrete <- function(TraitDB, Phylogeny, Discrete_trait) {
#   
#   row.names(TraitDB) <- TraitDB$Best_guess_binomial
#   Match <- match.phylo.data(Phylogeny, TraitDB)
#   Phylo <- Match$phy
#   Data <- Match$data
#   
#   Trait <- Data[, Discrete_trait]
#   
#   # Fit the different models (ARD, ER and SYM) -- model selection can be refined (Likelihood ratio test)
#   ER_model <- fitDiscrete(Phylo, Trait, type="discrete", model="ER")
#   ARD_model <- fitDiscrete(Phylo, Trait, type="discrete", model="ARD")
#   SYM_model <- fitDiscrete(Phylo, Trait, type="discrete", model="SYM")
#   
#   Model_select <- as.data.frame(c(ER_model$opt$aicc, ARD_model$opt$aicc, SYM_model$opt$aicc))
#   colnames(Model_select) <- "AICc"
#   row.names(Model_select) <- c("ER", "ARD", "SYM")
#   
#   Model_select$LogLk <- c(ER_model$opt$lnL, ARD_model$opt$lnL, SYM_model$opt$lnL)
#   
#   Model1 <- row.names(Model_select[which.min(Model_select$AICc),])
#   Model2 <- row.names(Model_select[which.max(Model_select$LogLk),])
#   
#   cat(Model1==Model2) 
#   
#   
#   # Rescale tree: lambda for a "0" tree
#   Phylo_lambda_0 <- rescale(Phylo, model = "lambda", 0)
#   
#   # par(mfrow=c(1,2))
#   # plot(Phylo, edge.width = 1, show.tip.label = FALSE) #smaller tree
#   # add.scale.bar(cex = 0.7, font = 2, col = "red")
#   # plot(Phylo_lambda_0, edge.width = 1, show.tip.label = FALSE) #smaller tree
#   # add.scale.bar(cex = 0.7, font = 2, col = "red")
#   
#   # Fit lambda to the tree
#   Lambda_trait <- fitDiscrete(Phylo, Trait, transform="lambda", type="discrete", model = Model1, niter=500)
#   
#   # Fit lambda to the 0-rescaled tree
#   Lambda_trait_0 <- fitDiscrete(Phylo_lambda_0, Trait, transform="lambda", type="discrete", model = Model1, niter=500)
#   
#   # Compare the outputs: compare the null model to the non null model - log likelihood and AICc difference
#   # AICc difference
#   AICc_diff <- abs(Lambda_trait_0$opt$aicc - Lambda_trait$opt$aicc)
#   d_logLk <- abs(2*(Lambda_trait_0$opt$lnL - Lambda_trait$opt$lnL))
#   p <- pchisq(d_logLk, 1, lower.tail=FALSE)
#   
#   if (p<0.05) {
#     cat("Model significantly different from null")
#     } else (cat("Model not significantly different from null"))
#   
#   Signal <- Lambda_trait$opt$lambda
#   
#   return(Signal)
#   
# }


# ### Phylogenetic signal on discrete traits
# 
# Trophic_level_mammal <- Lambda_discrete(TraitsMammal.Predicts, Mammal_phylo, "Trophic_level")
# Specialisation_mammal <- Lambda_discrete(TraitsMammal.Predicts, Mammal_phylo, "Specialisation")
# Diet_breadth_mammal <- Lambda_discrete(TraitsMammal.Predicts, Mammal_phylo, "Diet_breadth")
# 
# 
# Trophic_level_amphibian <- Lambda_discrete(TraitsAmphibian.Predicts_basic, Amphibian_phylo, "Trophic_level")
# Specialisation_amphibian <- Lambda_discrete(TraitsAmphibian.Predicts_basic, Amphibian_phylo, "Specialisation")
# Diet_breadth_amphibian <- Lambda_discrete(TraitsAmphibian.Predicts_basic, Amphibian_phylo, "Diet_breadth")
# 
# 
# Trophic_level_bird <- Lambda_discrete(TraitsBird.Predicts, Bird_phylo, "Trophic_level")
# Specialisation_bird <- Lambda_discrete(TraitsBird.Predicts, Bird_phylo, "Specialisation")
# Diet_breadth_bird <- Lambda_discrete(TraitsBird.Predicts, Bird_phylo, "Diet_breadth")
# 
# 
# Trophic_level_reptile <- Lambda_discrete(TraitsReptile.Predicts_basic, Reptile_phylo, "Trophic_level")
# Specialisation_reptile <- Lambda_discrete(TraitsReptile.Predicts_basic, Reptile_phylo, "Specialisation")

# Lambda_discrete.Reptile <- function(TraitDB, Phylogeny, Trophic_level, Specialisation) {
#   
#   Results <- as.data.frame(matrix(nrow=1, ncol=2))
#   colnames(Results) <- c(Trophic_level, Specialisation)
#   
#   row.names(TraitDB) <- TraitDB$Best_guess_binomial
#   Match <- match.phylo.data(Phylogeny, TraitDB)
#   Phylo <- Match$phy
#   Data <- Match$data
#   
#   Troph <- Data[, Trophic_level]
#   Lambda_trophic_level <- fitDiscrete(Phylo, Troph, transform="lambda")
#   Results[1, Trophic_level] <- Lambda_trophic_level$opt$lambda
#   
#   Spe <- Data[, Specialisation]
#   Lambda_Spe <- fitDiscrete(Phylo, Spe, transform="lambda")
#   Results[1, Specialisation] <- Lambda_Spe$opt$lambda
#   
#   return(Results)
# }
# 
# 
# Discrete_Mammals <- Lambda_discrete(TraitsMammal.Predicts, Mammal_phylo, "Trophic_level", "Specialisation", "Diet_breadth")
# Discrete_Amphibians <- Lambda_discrete(TraitsAmphibian.Predicts_basic, Amphibian_phylo, "Trophic_level", "Specialisation", "Diet_breadth")
# Discrete_Reptiles <- Lambda_discrete.Reptile(TraitsReptile.Predicts_basic, Reptile_phylo, "Trophic_level", "Specialisation")
# Discrete_Birds <- Lambda_discrete(TraitsBird.Predicts, Bird_phylo, "Trophic_level", "Specialisation", "Diet_breadth")


