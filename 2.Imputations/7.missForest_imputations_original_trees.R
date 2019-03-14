# Imputations using missForest, with  phylogenetic eigenvectors as predictors
## Here for phylogenetic infomormation extracted from original phylogenies, untouched

# still to do:
# compare results from missForest imputations with results from phylogenetic imputations - at least for continuous traits?
# error of imputations / robustness of imputations.
# compare results with Rob Cooke's.

## Code is parallelised

# https://www.jottr.org/2018/06/23/future.apply_1.0.0/

library(parallel)

# # ## START CLUSTER
Cluster <- makeCluster(detectCores())

## EXCECUTE ANY PRE PROCESSING CODE NECESSARY
clusterEvalQ(Cluster, {
  library(dplyr)
  library(phytools)
  library(missForest)
  library(pbmcapply)
  library(pbapply)
})

# ## Preamble
`%nin%` <- Negate(`%in%`)

source("Functions_for_missForest_imputations.R")

# ## Load trait data: transformed and standardised, with phylogenetic imformation as eigenvectors
# Mammals_st <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Mammals.csv")
# Birds_st <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Birds.csv")
# Amphibians_st <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Amphibians.csv")
# Reptiles_st <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.transformed_traits/Reptiles.csv")

## Load trait data: NOT transformed and NOT standardised, with phylogenetic imformation as eigenvectors
## chnage here
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/corrected_datasets_but_uncorrected_phylogenies/Reptiles.csv")


## Define variables for imputations

Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
             "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")

Taxinfo <- "Order"

# Traits_cont_st <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", "sqrt_Diet_breadth",
#                   "Range_size_m2", "sqrt_Habitat_breadth_IUCN")

Traits_cont <-  c("Body_mass_g", "Longevity_d", "Litter_size", "Diet_breadth",
                  "Range_size_m2", "Habitat_breadth_IUCN")

Traits_cat <- c(Habitat, "Specialisation",      
                "Diel_activity","Trophic_level", Diet, "Primary_diet")


## Add some traits for each taxon
# MammalsCont_st <- c(Traits_cont_st, "log10_Generation_length_d", "log10_Adult_svl_cm")
# BirdsCont_st <- Traits_cont_st
# ReptilesCont_st <- c(Traits_cont_st, "log10_Adult_svl_cm", "log10_Maturity_d")
# AmphibiansCont_st <- c(Traits_cont_st, "log10_Body_length_mm")

MammalsCont <- c(Traits_cont, "Generation_length_d", "Adult_svl_cm")
BirdsCont <- c(Traits_cont, "Generation_length_d")
ReptilesCont <- c(Traits_cont, "Adult_svl_cm", "Maturity_d")
AmphibiansCont <- c(Traits_cont, "Body_length_mm")


## Function arguments as lists, nested into one bigger list - each of these list elements are agurments for the function Imputations_missForest
DF.TraitsList <- list(M=Mammals, B=Birds, R=Reptiles, A=Amphibians)
#DF.TraitsList_st <- list(M=Mammals_st, B=Birds_st, R=Reptiles_st, A=Amphibians_st)
Taxinfo.List <- list(M="Order", B="Order", R="Order", A="Order")
#Cont.TraitsList_st <- list(M=MammalsCont_st, B=BirdsCont_st, R=ReptilesCont_st, A=AmphibiansCont_st)
Cont.TraitsList <- list(M=MammalsCont, B=BirdsCont, R=ReptilesCont, A=AmphibiansCont)
Cat.TraitsList <- list(M=Traits_cat, B=Traits_cat, R=Traits_cat[Traits_cat %nin% Diet], A=Traits_cat)
EV.List <- list(M="EV_1", B="EV_1", R="EV_1", A="EV_1")
ErrorTrue.List <- list(M=TRUE, B=TRUE, R=TRUE, A=TRUE)
DietTRUE.List <- list(M=TRUE, B=TRUE, R=FALSE, A=TRUE)
#std.List_st <- list(M=TRUE, B=TRUE, R=TRUE, A=TRUE)
std.List <- list(M=FALSE, B=FALSE, R=FALSE, A=FALSE)

# List of function arguments. This list  will be replicated 8 times (number of cores) for parallel imputations.
# On each cluster, imputation of 4 datasets (one for each class). 
# ArgumentsList_st <- list(TraitDF=DF.TraitsList_st,
# 					            Taxinfo=Taxinfo.List,
# 					            Traits_cont=Cont.TraitsList_st,
# 					            Traits_cat=Cat.TraitsList,
# 					            EV=EV.List,
# 					            ErrorTrue=ErrorTrue.List,
# 					            DietTRUE=DietTRUE.List, 
# 					            std=std.List_st)

ArgumentsList <- list(TraitDF=DF.TraitsList,
                      Taxinfo=Taxinfo.List,
                      Traits_cont=Cont.TraitsList,
                      Traits_cat=Cat.TraitsList,
                      EV=EV.List,
                      ErrorTrue=ErrorTrue.List,
                      DietTRUE=DietTRUE.List, 
                      std=std.List)

# Replicate this list N times so that: list with N elements, each of these are ArgumentsLists
N <- 8
# To_impute_parallel_st <- rep(list(ArgumentsList_st), N)
To_impute_parallel <- rep(list(ArgumentsList), N)

rm(Mammals, Birds, Reptiles, Amphibians,
   Mammals_st, Birds_st, Reptiles_st, Amphibians_st, 
   Habitat, Diet, Taxinfo, Traits_cont, Traits_cat, Traits_cont_st, 
   MammalsCont, BirdsCont, ReptilesCont, AmphibiansCont,
   MammalsCont_st, BirdsCont_st, ReptilesCont_st, AmphibiansCont_st,
   DF.TraitsList, DF.TraitsList_st, Taxinfo.List, Cont.TraitsList, Cat.TraitsList, Cont.TraitsList_st,
   EV.List, ErrorTrue.List, DietTRUE.List, N, std.List, std.List_st,
   ArgumentsList, ArgumentsList_st)

## Export variables in all clusters
clusterExport(cl=Cluster, list("Imputations_missForest",
                               "To_apply_parallel_imputations",
                               #"To_impute_parallel_st",
                               "To_impute_parallel",
                               "%nin%"),
              envir=environment())

## Parallel imputations on 8 cores
# print("Imputations on tranformed and standardised data")
# system.time(Imputed_sets_st <- parLapply(cl=Cluster,
#                           X=To_impute_parallel_st,
#                           fun=To_apply_parallel_imputations))

print("Imputations on non-tranformed/standardised data")
system.time(Imputed_sets <- parLapply(cl=Cluster,
                                      X=To_impute_parallel,
                                      fun=To_apply_parallel_imputations))

## Save results
# saveRDS(Imputed_sets_st, "../../Results/2.imputed_trait_datasets/Imputed_on_standardised/List_of_8_sets.rds")
saveRDS(Imputed_sets, "../../Results/2.imputed_trait_datasets/Imputed_original_trees/List_of_8_sets.rds")

## DESTROY CLUSTER
stopCluster(Cluster)

# # # # # for quickly testing code
# Mammals <- Mammals[c(1:15),]
# Birds <- Birds[c(1:15),]
# Amphibians <- Amphibians[c(1:15),]
# Reptiles <- Reptiles[c(1:15),]
# DF.TraitsList <- list(M=Mammals, B=Birds, R=Reptiles, A=Amphibians)
# ArgumentsList <- list(TraitDF=DF.TraitsList,
#                       Taxinfo=Taxinfo.List,
#                       Traits_cont=Cont.TraitsList,
#                       Traits_cat=Cat.TraitsList,
#                       EV=EV.List,
#                       ErrorTrue=ErrorTrue.List,
#                       DietTRUE=DietTRUE.List,
#                       std=std.List)
# 
# Test <- To_apply_parallel_imputations(ArgumentsList)


