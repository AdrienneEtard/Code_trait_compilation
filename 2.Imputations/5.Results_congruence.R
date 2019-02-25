## Imputation results congruence for each class, with and without phylogenetic information
library(reshape)
library(ggplot2)
source("Functions_for_results_congruence.R")

## Load imputed datasets (8)
Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_datasets/List_8_imputed_sets_v2.rds")

## Load data before imputations
Coll_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv")
Coll_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv")
Coll_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv")
Coll_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv")


## Amphibians
All_Amphibians <- Get_all_results(Imputed, "A")$Results
Amphibians_errors <- Get_all_results(Imputed, "A")$Errors
p <- Plot.Congruence.Continuous(All_Amphibians, Coll_Amphibians)
ggsave(p, file="../../Results/Plots/Congruence_imputations/Continuous/Amphibians.pdf", width=10, height=8)


## Reptiles
All_Reptiles <- Get_all_results(Imputed, "R")$Results
Reptiles_errors <- Get_all_results(Imputed, "R")$Errors
p <- Plot.Congruence.Continuous(All_Reptiles, Coll_Reptiles)
ggsave(p, file="../../Results/Plots/Congruence_imputations/Continuous/Reptiles.pdf", width=10, height=8)


## Mammals
All_Mammals <- Get_all_results(Imputed, "M")$Results
Mammals_errors <- Get_all_results(Imputed, "M")$Errors
p <- Plot.Congruence.Continuous(All_Mammals, Coll_Mammals)
ggsave(p, file="../../Results/Plots/Congruence_imputations/Continuous/Mammals.pdf", width=10, height=8)


## Birds
All_Birds <- Get_all_results(Imputed, "B")$Results
Birds_errors <- Get_all_results(Imputed, "B")$Errors
p <- Plot.Congruence.Continuous(All_Birds, Coll_Birds)
ggsave(p, file="../../Results/Plots/Congruence_imputations/Continuous/Birds.pdf", width=10, height=8)
