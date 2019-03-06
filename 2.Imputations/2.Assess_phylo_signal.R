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
  library(ape)
  library(pbmcapply)
  library(pbapply)
  source("Functions_Borges_et_al_2018_delta_statistic.R")
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

# Function to apply for categorical traits: PhySignal_Cat. Based on Borges et al 2018, Bioinformatics
PhySignal_Cat <- function(Traitdata, Names, Phylo, n) {
  
  # n = number of simulations
  ## For the current trait, match and prune phylogeny
  names(Traitdata) <- Names
  Match <- match.phylo.data(Phylo, Traitdata)
  Phylo <- Match$phy
  Trait <- Match$data
  rm(Match)
  
  ## Priors and parameters
  lambda0 <- 0.1   #rate parameter of the proposal 
  se      <- 0.5   #standard deviation of the proposal
  sim     <- 10000 #number of iterations
  thin    <- 10    #we kept only each 10th iterate 
  burn    <- 100   #100 iterates are burned-in
  
  ## Run the function (Borges et al, 2018, Bioinformatics: Measuring phylogenetic signal between categorical traits and phylogenies)
  ## with trycatch to avoid the process crashing when encountering errors.
  Delta <- tryCatch(expr={delta(Trait, Phylo, lambda0, se, sim, thin, burn)}, error = function(e) {NA})
  
  ## If the signal is not NA, then calculate a null distribution of delta values for the trait
  if(!is.na(Delta)) {
    
    # Generate randomised trait values - n times to generate a null distribution of delta -- stored in a list
    Func <- function(Trait){
      N <- length(Trait)
      L <- levels(as.factor(Trait))
    return(sample(L, size=N, replace=TRUE))
    }
    ListRandom <- lapply(rep(list(Trait), n), Func)
    
    Func_delta_toapply <- function(trait, tree, lambda0, se, sim, thin, burn) {
      Result <- tryCatch(expr = {delta(trait, tree, lambda0, se, sim, thin, burn)}, 
                         error=function(e){NA})
      return(Result)
    }
   
    Random_Delta <- pbmapply(FUN=Func_delta_toapply,
                              trait=ListRandom,
                              tree=rep(list(Phylo), n),
                              lambda0=rep(list(lambda0), n),
                              se=rep(list(se), n),
                              sim=rep(list(sim), n),
                              thin=rep(list(thin), n),
                              burn=rep(list(burn), n)) %>%
      as.data.frame()
    return(list(Delta=Delta, Delta0=Random_Delta))
    
  }
  
  else{return(Delta)}
}


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
clusterExport(cl=Cluster, list(".Format_tiplabels","PhySignal", "PhySignal_Cat",
                               "Mammals", "Birds", "Reptiles", "Amphibians",
                               "Phylo_Mammals", "Phylo_Birds", "Phylo_Reptiles", "Phylo_Amphibians",
                               "Continuous.Traits", "Categorical.Traits",
                               "Names.Mammals", "Names.Birds", "Names.Reptiles", "Names.Amphibians"), envir=environment())


## PARALLEL CALCULATIONS WITH parApply

## 1. Phylogenetic signal in continuous traits (takes up to 6 hours for birds)

# start_time <- Sys.time()
# Lambda_Mammals_continuous <- parApply(Cluster, Mammals[, c(Continuous.Traits,"Generation_length_d", "Adult_svl_cm")], 2, PhySignal, Names=Names.Mammals, Phylo=Phylo_Mammals) %>% as.data.frame()
# end_time <- Sys.time()
# end_time-start_time
# 
# start_time <- Sys.time()
# Lambda_Amphibians_continuous <- parApply(Cluster, Amphibians[, c(Continuous.Traits, "Body_length_mm")], 2, PhySignal, Names=Names.Amphibians, Phylo=Phylo_Amphibians) %>% as.data.frame()
# end_time <- Sys.time()
# end_time-start_time
#
# RC <- Continuous.Traits[-which(Continuous.Traits=="Diet_breadth")]
# start_time <- Sys.time()
# Lambda_Reptiles_continuous <- parApply(Cluster, Reptiles[, c(RC, "Adult_svl_cm", "Maturity_d")], 2, PhySignal, Names=Names.Reptiles, Phylo=Phylo_Reptiles) %>% as.data.frame()
# end_time <- Sys.time()
# end_time-start_time
#
# start_time <- Sys.time()
# Lambda_Birds_continuous <- parApply(Cluster, Birds[, c(Continuous.Traits, "Generation_length_d")], 2, PhySignal, Names=Names.Birds, Phylo=Phylo_Birds) %>% as.data.frame()
# end_time <- Sys.time()
# end_time-start_time
#
## Save files
# write.csv(Lambda_Mammals_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousMammals.csv",row.names = FALSE)
# write.csv(Lambda_Birds_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousBirds.csv",row.names = FALSE)
# write.csv(Lambda_Reptiles_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousReptiles.csv",row.names = FALSE)
# write.csv(Lambda_Amphibians_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/ContinuousAmphibians.csv",row.names = FALSE)


## 2. Phylogenetic signal in categorical traits

start_time <- Sys.time()
delta_Mammals_categorical <- parApply(Cluster, Mammals[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Mammals, Phylo=Phylo_Mammals, n=50)
end_time <- Sys.time()
Mt <- end_time-start_time  # about 1.5hrs

start_time <- Sys.time()
delta_Birds_categorical <- parApply(Cluster, Birds[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Birds, Phylo=Phylo_Birds, n=50)
end_time <- Sys.time()
Bt <- end_time-start_time

start_time <- Sys.time()
delta_Reptiles_categorical <- parApply(Cluster, Reptiles[, Categorical.Traits[-which(Categorical.Traits=="Primary_diet")]], 2, PhySignal_Cat, Names=Names.Reptiles, Phylo=Phylo_Reptiles, n=50)
end_time <- Sys.time()
Rt <- end_time-start_time

start_time <- Sys.time()
delta_Amphibians_categorical <- parApply(Cluster, Amphibians[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Amphibians, Phylo=Phylo_Amphibians, n=50)
end_time <- Sys.time()
At <- end_time-start_time

## Save files
saveRDS(delta_Mammals_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalMammals.rds")
saveRDS(delta_Birds_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalBirds.rds")
saveRDS(delta_Reptiles_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalReptiles.rds")
saveRDS(delta_Amphibians_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/CategoricalAmphibians.rds")


## DESTROY CLUSTER
stopCluster(Cluster)





# # # Older script for categorical traits

# # Function to apply for categorical traits: PhySignal_Cat
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
#       cat("Model significantly different from null")
#     } else {cat("Model not significantly different from null")}
#     
#     Lambda <- Lambda_trait$opt$lambda
#     pvalue <- p
#     
#     List <- list(Lambda=Lambda, AICc_diff=AICc_diff, d_logLk=d_logLk, pvalue=pvalue)
#     
#     return(unlist(List))
#   }
# }

