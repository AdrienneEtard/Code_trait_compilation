## Script to assess the phylogenetic signal in Traits (continuous and categorical)

## Here with the unvorrected, original phylogenies

# TODO phylogenetic signal for categorical traits

## Code is parallelised
library(parallel)
library(dplyr)
library(phytools)

## START CLUSTER
n <- detectCores()-1
Cluster <- makeCluster(n)
rm(n)

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
  
  ## The Borges function does not work when branches have 0 length
  ## in that case, add a very small number to these branches
  Phylo$edge.length[Phylo$edge.length==0] <- 10e-10
  
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
Phylo_Mammals <- read.newick("../../Data/Phylogenies/TTOL_mammals_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels() 
Phylo_Amphibians <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()
Phylo_Reptiles <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()
Phylo_Birds <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels()  

# # Are there any polytomies in the trees?
# is.binary(Phylo_Birds)      # true
# is.binary(Phylo_Mammals)    # true
# is.binary(Phylo_Reptiles)   # true
# is.binary(Phylo_Amphibians) # false

Phylo_Amphibians <- multi2di(Phylo_Amphibians)

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

# log-transform continuous traits
Transf <- function(Data, Traits_log10) {
  Data[, Traits_log10] <- log10(Data[, Traits_log10])
  Data$Habitat_breadth_IUCN <- sqrt(Data$Habitat_breadth_IUCN)
  return(Data)
}

Mammals <- Transf(Mammals, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")],"Generation_length_d", "Adult_svl_cm"))
Birds <- Transf(Birds, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")],"Generation_length_d"))
Amphibians <- Transf(Amphibians, Traits_log10=c(Continuous.Traits[-which(Continuous.Traits=="Habitat_breadth_IUCN")],"Body_length_mm"))
Reptiles <- Transf(Reptiles, Traits_log10=c("Body_mass_g","Longevity_d","Litter_size","Adult_svl_cm", "Range_size_m2","Maturity_d"))


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

print("starting estimations on continuous traits")
Tstart <- Sys.time()
print(Tstart)

start_time <- Sys.time()
Lambda_Mammals_continuous <- parApply(Cluster, Mammals[, c(Continuous.Traits,"Generation_length_d", "Adult_svl_cm")], 2, PhySignal, Names=Names.Mammals, Phylo=Phylo_Mammals) %>% as.data.frame()
end_time <- Sys.time()
timeMcont <- end_time-start_time
print(timeMcont)
write.csv(Lambda_Mammals_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousMammals_log10.csv",row.names = FALSE)

print(Sys.time())
start_time <- Sys.time()
Lambda_Amphibians_continuous <- parApply(Cluster, Amphibians[, c(Continuous.Traits, "Body_length_mm")], 2, PhySignal, Names=Names.Amphibians, Phylo=Phylo_Amphibians) %>% as.data.frame()
end_time <- Sys.time()
timeAcont <- end_time-start_time
print(timeAcont)
write.csv(Lambda_Amphibians_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousAmphibians_log10.csv",row.names = FALSE)

print(Sys.time())
RC <- Continuous.Traits[-which(Continuous.Traits=="Diet_breadth")]
start_time <- Sys.time()
Lambda_Reptiles_continuous <- parApply(Cluster, Reptiles[, c(RC, "Adult_svl_cm", "Maturity_d")], 2, PhySignal, Names=Names.Reptiles, Phylo=Phylo_Reptiles) %>% as.data.frame()
end_time <- Sys.time()
timeRcont <- end_time-start_time
print(timeRcont)
write.csv(Lambda_Reptiles_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousReptiles_log10.csv",row.names = FALSE)

print(Sys.time())
start_time <- Sys.time()
Lambda_Birds_continuous <- parApply(Cluster, Birds[, c(Continuous.Traits, "Generation_length_d")], 2, PhySignal, Names=Names.Birds, Phylo=Phylo_Birds) %>% as.data.frame()
end_time <- Sys.time()
timeBcont <- end_time-start_time
print(timeBcont)
write.csv(Lambda_Birds_continuous, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousBirds_log10.csv",row.names = FALSE)


# ## 2. Phylogenetic signal in categorical traits
# 
# print("starting estimations on categorical traits")
# 
# start_time <- Sys.time()
# delta_Mammals_categorical <- parApply(Cluster, Mammals[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Mammals, Phylo=Phylo_Mammals, n=50)
# end_time <- Sys.time()
# timeMcat <- end_time-start_time 
# print(timeMcat)
# 
# start_time <- Sys.time()
# delta_Birds_categorical <- parApply(Cluster, Birds[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Birds, Phylo=Phylo_Birds, n=50)
# end_time <- Sys.time()
# timeBcat <- end_time-start_time
# print(timeBcat)
# 
# start_time <- Sys.time()
# delta_Reptiles_categorical <- parApply(Cluster, Reptiles[, Categorical.Traits[-which(Categorical.Traits=="Primary_diet")]], 2, PhySignal_Cat, Names=Names.Reptiles, Phylo=Phylo_Reptiles, n=50)
# end_time <- Sys.time()
# timeRcat <- end_time-start_time
# print(timeRcat)
# 
# start_time <- Sys.time()
# delta_Amphibians_categorical <- parApply(Cluster, Amphibians[, Categorical.Traits], 2, PhySignal_Cat, Names=Names.Amphibians, Phylo=Phylo_Amphibians, n=50)
# end_time <- Sys.time()
# timeAcat <- end_time-start_time
# print(timeAcat)
# 
# ## Save files
# saveRDS(delta_Mammals_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalMammals.rds")
# saveRDS(delta_Birds_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalBirds.rds")
# saveRDS(delta_Reptiles_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalReptiles.rds")
# saveRDS(delta_Amphibians_categorical, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/CategoricalAmphibians.rds")
# 
# TEnd <- Sys.time()
# print(Tstart-TEnd)

## DESTROY CLUSTER
stopCluster(Cluster)





