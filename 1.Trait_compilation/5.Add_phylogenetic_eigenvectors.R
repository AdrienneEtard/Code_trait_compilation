## Adding phylogenetic eigenvectors to trait datasets, using both corrected and uncorrected phylogenies -
# phylogenies that have been treated for pseudoreplication in (4.)

# # # # P r e a m b l e 

X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ape", "treeio", "ngram","phylobase")
lapply(X, library, character.only=TRUE); rm(X)
source("Functions_for_phylogenies.R")

## Load data

# No taxonomic correction
UN_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Amphibians.csv")
UN_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Birds.csv")
UN_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Mammals.csv")
UN_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Reptiles.csv")

# With taxonomic correction
C_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Amphibians.csv")
C_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Birds.csv")
C_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Mammals.csv")
C_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Reptiles.csv")

## Load PREDICTS with and without taxonomic correction
UN_Predicts <- readRDS("../../Data/PREDICTS_database.rds") %>%
  filter(Class %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"))
C_Predicts <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

## Load phylogenies

# Corrected
C_Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk")  %>% .Format_tiplabels() %>% compute.brlen()
C_Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk")  %>% .Format_tiplabels() %>% compute.brlen()
C_Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk")  %>% .Format_tiplabels() %>% compute.brlen()
C_Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk")  %>% .Format_tiplabels() %>% compute.brlen()

# Uncorrected
UN_Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Mammals.nwk")  %>% .Format_tiplabels() %>% compute.brlen()
UN_Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Birds.nwk")  %>% .Format_tiplabels() %>% compute.brlen()
UN_Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Amphibians.nwk")  %>% .Format_tiplabels() %>% compute.brlen()
UN_Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Reptiles.nwk")  %>% .Format_tiplabels() %>% compute.brlen()


## Phylogenetic eigenvectors added to each trait datasets, when possible (when the species is present in the phylogeny).
## I add the first 10 eigenvectors (enough to maximise the accuracy of further imputations.)

# For corrected datasets
C.EV_Amphibians <- Add_eigenvectors(C_Amphibians, C_Phylo_Amphibians, 10, TRUE)
C.EV_Reptiles <- Add_eigenvectors(C_Reptiles, C_Phylo_Reptiles, 10, TRUE)
C.EV_Birds <- Add_eigenvectors(C_Birds, C_Phylo_Birds, 10, TRUE)
C.EV_Mammals <- Add_eigenvectors(C_Mammals, C_Phylo_Mammals, 10, TRUE)


# For uncorrected datasets
UN.EV_Amphibians <- Add_eigenvectors(UN_Amphibians, UN_Phylo_Amphibians, 10, FALSE)
UN.EV_Mammals <- Add_eigenvectors(UN_Mammals, UN_Phylo_Mammals, 10, FALSE)
UN.EV_Reptiles <- Add_eigenvectors(UN_Reptiles, UN_Phylo_Reptiles, 10, FALSE)
UN.EV_Birds <- Add_eigenvectors(UN_Birds, UN_Phylo_Birds, 10, FALSE)

## Save datasets with added phylogenetic eigenvectors
write.csv(C.EV_Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv", row.names=FALSE)
write.csv(C.EV_Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv", row.names=FALSE)
write.csv(C.EV_Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv", row.names=FALSE)
write.csv(C.EV_Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv", row.names=FALSE)

write.csv(UN.EV_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv", row.names=FALSE)
write.csv(UN.EV_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv", row.names=FALSE)
write.csv(UN.EV_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv", row.names=FALSE)
write.csv(UN.EV_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv", row.names=FALSE)




