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


# ## Quick compilation (eigenvectors take a long time to be extracted) -- corrected
# EVC_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals_EV.csv")
# EVC_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds_EV.csv")
# EVC_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles_EV.csv")
# EVC_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians_EV.csv")
# 
# ## BindC.EV_Mammals <- cbind(C_Mammals, EVC_Mammals)
# C.EV_Birds <- cbind(C_Birds, EVC_Birds)
# C.EV_Reptiles <- cbind(C_Reptiles, EVC_Reptiles)
# C.EV_Amphibians <- cbind(C_Amphibians, EVC_Amphibians)
# 
# 
# 
# ## Quick compilation (eigenvectors take a long time to be extracted) -- uncorrected
# EVU_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals_EV.csv")
# EVU_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds_EV.csv")
# EVU_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles_EV.csv")
# EVU_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians_EV.csv")
# 
# ## Bind
# UN.EV_Mammals <- cbind(UN_Mammals, EVU_Mammals)
# UN.EV_Birds <- cbind(UN_Birds, EVU_Birds)
# UN.EV_Reptiles <- cbind(UN_Reptiles, EVU_Reptiles)
# UN.EV_Amphibians <- cbind(UN_Amphibians, EVU_Amphibians)
# 
# ## Save
# write.csv(C.EV_Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv", row.names=FALSE)
# write.csv(C.EV_Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv", row.names=FALSE)
# write.csv(C.EV_Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv", row.names=FALSE)
# write.csv(C.EV_Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv", row.names=FALSE)
# 
# write.csv(UN.EV_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv", row.names=FALSE)
# write.csv(UN.EV_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv", row.names=FALSE)
# write.csv(UN.EV_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv", row.names=FALSE)
# write.csv(UN.EV_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv", row.names=FALSE)



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


## For quicker compilations: saving just the eigenvectors

C_E_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)
C_E_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)
C_E_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)
C_E_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)

write.csv(C_E_Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals_EV.csv", row.names=F)
write.csv(C_E_Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds_EV.csv", row.names=F)
write.csv(C_E_Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles_EV.csv", row.names=F)
write.csv(C_E_Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians_EV.csv", row.names=F)

U_E_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)
U_E_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)
U_E_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)
U_E_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv") %>%
  select(EV_1, EV_2, EV_3, EV_4, EV_5, EV_6, EV_7, EV_8, EV_9, EV_10)

write.csv(U_E_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals_EV.csv", row.names=F)
write.csv(U_E_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds_EV.csv", row.names=F)
write.csv(U_E_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles_EV.csv", row.names=F)
write.csv(U_E_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians_EV.csv", row.names=F)



