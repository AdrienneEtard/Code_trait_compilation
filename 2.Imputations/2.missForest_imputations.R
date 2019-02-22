# Imputations using missForest, with  phylogenetic eigenvectors as predictors -- imputation errors in another script.

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
# X <- c("dplyr", "phytools", "missForest", "pbapply")
# lapply(X, library, character.only=TRUE); rm(X)
`%nin%` <- Negate(`%in%`)

source("Functions_for_missForest_imputations.R")

## Load trait data, transformed, standardised, and with phylogenetic imformation
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv")

## Define variables for imputations

Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
             "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

Diet <- c("IN", "VE", "PL", "SE", "NE", "FR", "SCV")

Taxinfo <- "Order"

Traits_cont <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", "sqrt_Diet_breadth",
                  "Range_size_m2", "sqrt_Habitat_breadth_IUCN")

Traits_cat <- c(Habitat, "Specialisation",        # Include diet now or derive later from imputed primary diet?
               "Diel_activity","Trophic_level", Diet)
               
#EV <- c(); for(i in 1:10) {EV <- c(EV,paste("EV_",i, sep=""))}

## Add some traits for each taxon # TODO check trait names
MammalsCont <- c(Traits_cont, "log10_Generation_length_d", "log10_Adult_svl_cm")
BirdsCont <- Traits_cont
ReptilesCont <- c(Traits_cont, "log10_Adult_svl_cm", "log10_Maturity_d")
AmphibiansCont <- c(Traits_cont, "log10_Body_length_mm")


## Function arguments as lists, nested into one bigger list - each of these list elements are agurments for the function Imputations_missForest
DF.TraitsList <- list(M=Mammals, B=Birds, R=Reptiles, A=Amphibians)
Taxinfo.List <- list(M="Order", B="Order", R="Order", A="Order")
Cont.TraitsList <- list(M=MammalsCont, B=BirdsCont, R=ReptilesCont, A=AmphibiansCont)
Cat.TraitsList <- list(M=Traits_cat, B=Traits_cat, R=Traits_cat[Traits_cat %nin% Diet], A=Traits_cat)
EV.List <- list(M="EV_1", B="EV_1", R="EV_1", A="EV_1")
ErrorTrue.List <- list(M=TRUE, B=TRUE, R=TRUE, A=TRUE)
DietTRUE.List <- list(M=TRUE, B=TRUE, R=FALSE, A=TRUE)

# List of function arguments. This list  will be replicated 8 times (number of cores) for parallel imputations.
# On each cluster, imputation of 4 datasets (one for each class). 
ArgumentsList <- list(TraitDF=DF.TraitsList,
					            Taxinfo=Taxinfo.List,
					            Traits_cont=Cont.TraitsList,
					            Traits_cat=Cat.TraitsList,
					            EV=EV.List,
					            ErrorTrue=ErrorTrue.List,
					            DietTRUE=DietTRUE.List)

# Replicate this list N times so that: list with N elements, each of these are ArgumentsLists
N <- 8
To_impute_parallel <- rep(list(ArgumentsList), N)

rm(Mammals, Birds, Reptiles, Amphibians,
	Habitat, Diet, Taxinfo, Traits_cont, Traits_cat,
	MammalsCont, BirdsCont, ReptilesCont, AmphibiansCont,
	DF.TraitsList, Taxinfo.List, Cont.TraitsList, Cat.TraitsList,
	EV.List, ErrorTrue.List, DietTRUE.List, N)

## Export variables in all clusters
clusterExport(cl=Cluster, list("Imputations_missForest",
								"To_apply_parallel_imputations",
								"To_impute_parallel",
								"%nin%"), envir=environment())

## Parallel imputations on 8 cores
Imputed_sets <- parLapply(cl=Cluster,
							              X=To_impute_parallel,
							              fun=To_apply_parallel_imputations)

## Save results
saveRDS(Imputed_sets, "../../Results/2.imputed_trait_datasets/imputed_datasets/List_8_imputed_sets.rds")

## DESTROY CLUSTER
stopCluster(Cluster)


# ## to test: this function and also the apply version calling this one
# Test1 <- pbmapply (FUN=Imputations_missForest,
#                    TraitDF=DF.TraitsList,
#                    Taxinfo=Taxinfo,
#                    Traits_cont=Cont.TraitsList,
#                    Traits_cat=Cat.TraitsList,
#                    EV=EV.List,
#                    ErrorTrue=ErrorTrue.List,
#                    DietTRUE=DietTRUE.List)
# 

Imputed_set_1 <- To_apply_parallel_imputations(ArgumentsList)
saveRDS(Imputed_set_1, "../../Results/2.imputed_trait_datasets/imputed_datasets/Imputed_set_1.rds")

