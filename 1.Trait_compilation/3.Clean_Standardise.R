## Clean and standardise compiled datasets before imputations.
# rearrange columns
# merge longevity and max longevity columns
# standardise traits (log10 and sqrt transformations plus z-scoring)


# # PREAMBLE

X <-c("data.table", "plyr", "dplyr", "tidyr", "magrittr", "reshape", "reshape2", "stringr", "stringi", "lazyeval", "rlang", "PerformanceAnalytics") 
lapply(X, library, character.only=TRUE); rm(X)

Transform_zscore <- function(TraitDF, Trait, Transf) {
  
  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- log10(TraitDF[,Trait])
    TraitDF[, paste("log10", Trait, sep="_")] <- scale(TraitDF[, paste("log10", Trait, sep="_")], center=TRUE, scale=TRUE)
  }
  
  if(Transf=="sqrt") {
    TraitDF[,Trait]  <- sqrt(TraitDF[,Trait])
    TraitDF[, paste("sqrt", Trait, sep="_")] <- scale(TraitDF[,Trait] , center=TRUE, scale=TRUE)
  }
  
  return(TraitDF)
}

## Load trait data 
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered_taxinfo/All/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered_taxinfo/All/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered_taxinfo/All/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered_taxinfo/All/Reptiles.csv")
Reptiles$Maturity_d[Reptiles$Maturity_d==0] <- NA

## Predicts
PredictsVertebrates <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

# # # #


## Transforming and standardising traits

# 1. Merging together maturity and generation length, and longevity and max longevity

# Reptiles
Reptiles$Longevity_d <- apply(Reptiles[, c("Longevity_d","Max_longevity_d")], 1, mean, na.rm=TRUE)
Reptiles %<>% select(-Max_longevity_d)

# Amphibians
colnames(Amphibians)[7] <- c("Longevity_d")

# Birds
Birds$Longevity_d <- apply(Birds[, c("Max_longevity_d","Longevity_d")], 1, mean, na.rm=TRUE)
Birds %<>% select(-Max_longevity_d)

# Mammals
Mammals$Longevity_d <- apply(Mammals[, c("Max_longevity_d","Longevity_d")], 1, mean, na.rm=TRUE)
Mammals %<>% select(-Max_longevity_d)

# 2. Standardising traits

# Mammals
Mammals <- Transform_zscore(Mammals, "Body_mass_g", "log10")
Mammals <- Transform_zscore(Mammals, "Adult_svl_cm", "log10")
Mammals <- Transform_zscore(Mammals, "Forearm_length_mm", "log10")
Mammals <- Transform_zscore(Mammals, "Head_length_mm", "log10")
Mammals <- Transform_zscore(Mammals, "Generation_length_d", "log10")
Mammals <- Transform_zscore(Mammals, "Longevity_d", "log10")
Mammals <- Transform_zscore(Mammals, "Maturity_d", "log10")
Mammals <- Transform_zscore(Mammals, "AFR_d", "log10")
Mammals <- Transform_zscore(Mammals, "Litter_size", "log10")
Mammals <- Transform_zscore(Mammals, "Habitat_breadth_IUCN", "sqrt")

# Reptiles
Reptiles <- Transform_zscore(Reptiles, "Body_mass_g", "log10")
Reptiles <- Transform_zscore(Reptiles, "Adult_svl_cm", "log10")
Reptiles <- Transform_zscore(Reptiles, "Maturity_d", "log10")
Reptiles <- Transform_zscore(Reptiles, "Longevity_d", "log10")
Reptiles <- Transform_zscore(Reptiles, "Litter_size", "log10")
Reptiles <- Transform_zscore(Reptiles, "Habitat_breadth_IUCN", "sqrt")

# Birds
Birds <- Transform_zscore(Birds, "Body_mass_g", "log10")
Birds <- Transform_zscore(Birds, "Adult_svl_cm", "log10")
Birds <- Transform_zscore(Birds, "Maturity_d", "log10")
Birds <- Transform_zscore(Birds, "Longevity_d", "log10")
Birds <- Transform_zscore(Birds, "Litter_size", "log10")
Birds <- Transform_zscore(Birds, "Habitat_breadth_IUCN", "sqrt")

# Amphibians
Amphibians <- Transform_zscore(Amphibians, "Body_mass_g", "log10")
Amphibians <- Transform_zscore(Amphibians, "Body_length_mm", "log10")
Amphibians <- Transform_zscore(Amphibians, "Svl_length_mm", "log10")
Amphibians <- Transform_zscore(Amphibians, "Maturity_d", "log10")
Amphibians <- Transform_zscore(Amphibians, "Longevity_d", "log10")
Amphibians <- Transform_zscore(Amphibians, "Litter_size", "log10")
Amphibians <- Transform_zscore(Amphibians, "Habitat_breadth_IUCN", "sqrt")


## Reorganising columns
Diet <- c("IN", "VE", "PL", "SE", "NE", "FR", "SCV")
Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas",
             "Caves.and.subterranean","Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

Reptiles[, (ncol(Reptiles)+1):(ncol(Reptiles)+2)] <- NA
colnames(Reptiles)[c(34,35)] <- c("Primary_diet", "Diet_breadth")

Reptiles <- Reptiles[, c("Order", "Family", "Genus", "Best_guess_binomial",
                         "log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Maturity_d", "log10_Longevity_d",
                         "log10_Litter_size", "Range_size_m2", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", "Specialisation",
                         "sqrt_Habitat_breadth_IUCN", Habitat)]

Amphibians <- Amphibians[, c("Order", "Family", "Genus", "Best_guess_binomial",
                         "log10_Body_mass_g", "log10_Body_length_mm","log10_Svl_length_mm", "log10_Maturity_d", "log10_Longevity_d",
                         "log10_Litter_size", "Range_size_m2", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet",  Diet,"Specialisation",
                         "sqrt_Habitat_breadth_IUCN", Habitat)]

Birds <- Birds[, c("Order", "Family", "Genus", "Best_guess_binomial",
                             "log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Maturity_d", "log10_Longevity_d",
                             "log10_Litter_size", "Range_size_m2", "Diel_activity", "Trophic_level", "Diet_breadth","Primary_diet",Diet,"Specialisation",
                             "sqrt_Habitat_breadth_IUCN", Habitat)]


Mammals <- Mammals[, c("Order", "Family", "Genus", "Best_guess_binomial",
                   "log10_Body_mass_g", "log10_Adult_svl_cm", "log10_Forearm_length_mm","log10_Head_length_mm",
                   "log10_Generation_length_d","log10_Maturity_d", "log10_Longevity_d", "log10_AFR_d",
                   "log10_Litter_size", "Range_size_m2", "Diel_activity", "Trophic_level", "Diet_breadth", "Primary_diet", Diet, "Specialisation",
                   "sqrt_Habitat_breadth_IUCN", Habitat)]


# Saving results
write.csv(Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Reptiles.csv", row.names=FALSE)
write.csv(Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Mammals.csv", row.names=FALSE)
write.csv(Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Birds.csv", row.names=FALSE)
write.csv(Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Amphibians.csv", row.names=FALSE)

