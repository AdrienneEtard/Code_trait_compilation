## Add family / order information to compiled trait datasets.
## Filter out species that have no information for all traits (these species potentially are replicates); I have verified that none of these appear in PREDICTS.
## Filter out species that have no range information, unless they appear in the PREDICTS database?
## Filter out species that have completeness equal to zero


# # PREAMBLE
X <- c("dplyr", "stringr")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` = Negate(`%in%`)

## Traits
Traits <- c("Body_mass_g",
            "Longevity_d",
            "Litter_size", 
            "Range_size_m2", 
            "Habitat_breadth_IUCN",
            "Specialisation",
            "Trophic_level",
            "Diel_activity",
            "Primary_diet",
            "Diet_breadth")

TraitsReptiles <- c(Traits[1:8], "Adult_svl_cm", "Maturity_d")
TraitsAmphibians <- c(Traits[-which(Traits=="Longevity_d")], "Body_length_mm", "Max_longevity_d")
TraitsMammals <- c(Traits, "Adult_svl_cm", "Generation_length_d")
TraitsBirds <- c(Traits, "Generation_length_d")

Completeness <- function(TraitDF, Traits) {
  
  # completeness
  TraitDF$Percent <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y)))
  TraitDF$Percent <- TraitDF$Percent / length(Traits) * 100 
  
  Completeness_0 <- TraitDF$Best_guess_binomial[TraitDF$Percent==0]
  print(length(Completeness_0))

  
  return(Completeness_0)
  
}

## Load data 

# Corrected
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Reptiles.csv")

# Uncorrected
U_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Mammals.csv")
U_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Birds.csv")
U_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Amphibians.csv")
U_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Reptiles.csv")

Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals.csv")
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians.csv")
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv")
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles.csv")

## Predicts
Predicts <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

## Ranges
Range.Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Mammal_range.csv")
Range.Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Reptile_range.csv")
Range.Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bird_range.csv")
Range.Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibian_range.csv")

# Filter out of the trait compilation all species for which all information is NA (such as amphibian species "Eleutherodactylus fitzingereri")
FilterOut <- function(TraitDF) {

  TraitDF$Best_guess_binomial <- as.character(TraitDF$Best_guess_binomial)
  
  ToRemove <- c()
  for (i in 1:nrow(TraitDF)) { 
    
    X <- all(is.na(TraitDF[i, c(2:ncol(TraitDF))]))
    if(X==TRUE) {ToRemove <- c(ToRemove, TraitDF$Best_guess_binomial[i])}
  }
  return(ToRemove)
}

Species_Intersect <- function(TraitDF, Predicts, Range, Class) {
  
  Species_to_check <- FilterOut(TraitDF)
  print(length(Species_to_check))
  
  PSp <- unique(Predicts$Best_guess_binomial[Predicts$Class==Class])
  Y <- intersect(Species_to_check, PSp)
  print(paste(length(Y), "species intersect with PREDICTS species."))

  PRa <- unique(Range$Best_guess_binomial) %>% as.character()
  Z <- intersect(Species_to_check, PRa)
  print(paste(length(Z), "species intersect with Range dataset species."))
    
  Intersecting <- c(Y,Z) %>% unique()
  
  SpeciesToRemove <- setdiff(Species_to_check, Intersecting)
    
  return(SpeciesToRemove)
  
}

# Add family and order to the data
Add_family_order_genus <- function(TraitDF, SynDF) {
  
  # browser()
  
  TraitDF$Best_guess_binomial <- as.character(TraitDF$Best_guess_binomial)
  SynDF$Genus <- as.character(SynDF$Genus)
  SynDF$Family <- as.character(SynDF$Family)
  SynDF$Order <- as.character(SynDF$Order)
  
  TraitDF$Genus <- NA
  TraitDF$Family <- NA
  TraitDF$Order <- NA
  
  Sp <- intersect(TraitDF$Best_guess_binomial, SynDF$Accepted)
  
  for (i in 1:length(Sp)) {
    
    TraitDF$Genus[TraitDF$Best_guess_binomial==Sp[i]] <- SynDF$Genus[SynDF$Accepted==Sp[i]] %>% as.character %>% unique()
    TraitDF$Family[TraitDF$Best_guess_binomial==Sp[i]] <- SynDF$Family[SynDF$Accepted==Sp[i]] %>% as.character %>% unique()
    TraitDF$Order[TraitDF$Best_guess_binomial==Sp[i]] <- SynDF$Order[SynDF$Accepted==Sp[i]] %>% as.character %>% unique()
    
  }
  
  Sp2 <- setdiff(TraitDF$Best_guess_binomial, Sp)
  
  if (length(Sp2!=0)) {
    
    if (Sp2=="Desmodus draculae") {
      
      TraitDF$Genus[TraitDF$Best_guess_binomial=="Desmodus draculae"] <- "Desmodus"
      TraitDF$Family[TraitDF$Best_guess_binomial=="Desmodus draculae"] <- "PHYLLOSTOMIDAE"
      TraitDF$Order[TraitDF$Best_guess_binomial=="Desmodus draculae"] <- "CHIROPTERA"
      
    } 
    
    else {
      
      for (i in 1:length(Sp2)) {
        
        TraitDF$Genus[TraitDF$Best_guess_binomial==Sp2[i]] <- SynDF$Genus[SynDF$CorrectedTypos==Sp2[i]] %>% as.character %>% unique()
        TraitDF$Family[TraitDF$Best_guess_binomial==Sp2[i]] <- SynDF$Family[SynDF$CorrectedTypos==Sp2[i]] %>% as.character %>% unique()
        TraitDF$Order[TraitDF$Best_guess_binomial==Sp2[i]] <- SynDF$Order[SynDF$CorrectedTypos==Sp2[i]] %>% as.character %>% unique()
        
      }
    }
  }
  
  return(TraitDF)
}

# # Check if species that have no range information are also present in PREDICTS
# No_range <- function(TraitDF, Predicts, Class){
#   
#   # browser()
#   
#   SpNR <- TraitDF %>% filter(is.na(Range_size_m2))
#   SpNR <- unique(SpNR$Best_guess_binomial)
#   PSp <- unique(Predicts$Best_guess_binomial[Predicts$Class==Class]) 
#   
#   Y <- intersect(SpNR, PSp)
#   print(paste(length(Y), "species that have no range size information are also in PREDICTS."))
#   return(Y)
# }
# 
# Filter_No_range <- function(NTR, TraitDF) {
#   
#   ToRemove <- TraitDF$Best_guess_binomial[is.na(TraitDF$Range_size_m2)]
#   ToRemove <- setdiff(ToRemove, NTR)
#   
#   FTraitDF <- TraitDF %>% filter(Best_guess_binomial %nin% ToRemove)
#   Excluded <- TraitDF %>% filter(Best_guess_binomial %in% ToRemove)
#   
#   return(list(TraitDF=FTraitDF, Removed=Excluded))
# }


# # # #

# 1. Remove species that have 0 information from trait datasets (all traits are NA)

# Corrected
ToRemoveBirds <- Species_Intersect(Birds, Predicts, Range.Birds, "Aves")
ToRemoveMammals <- Species_Intersect(Mammals, Predicts, Range.Mammals, "Mammalia")
ToRemoveAmphibians <- Species_Intersect(Amphibians, Predicts, Range.Amphibians, "Amphibia")
ToRemoveReptiles <- Species_Intersect(Reptiles, Predicts, Range.Reptiles, "Reptilia")

Mammals <- Mammals %>% filter(Best_guess_binomial %nin% ToRemoveMammals)
Reptiles <- Reptiles %>% filter(Best_guess_binomial %nin% ToRemoveReptiles)
Amphibians <- Amphibians %>% filter(Best_guess_binomial %nin% ToRemoveAmphibians)


# Uncorrected
U_ToRemoveBirds <- Species_Intersect(U_Birds, Predicts, Range.Birds, "Aves")
U_ToRemoveMammals <- Species_Intersect(U_Mammals, Predicts, Range.Mammals, "Mammalia")
U_ToRemoveAmphibians <- Species_Intersect(U_Amphibians, Predicts, Range.Amphibians, "Amphibia")
U_ToRemoveReptiles <- Species_Intersect(U_Reptiles, Predicts, Range.Reptiles, "Reptilia")

U_Mammals <- U_Mammals %>% filter(Best_guess_binomial %nin% U_ToRemoveMammals)
U_Reptiles <- U_Reptiles %>% filter(Best_guess_binomial %nin% U_ToRemoveReptiles)
U_Amphibians <- U_Amphibians %>% filter(Best_guess_binomial %nin% U_ToRemoveAmphibians)


# # 2. Remove species that have no range size information, except if they appear in PREDICTS
# 
# # Check is.na(Range_size_m2) that intersect with PREDICTS
# NTRMammals <- No_range(Mammals, Predicts, "Mammalia")
# NTRAmphibians <- No_range(Amphibians, Predicts, "Amphibia")
# NTRReptiles <- No_range(Reptiles, Predicts, "Reptilia")
# NTRBirds <- No_range(Birds, Predicts, "Aves")
# 
# # Filter out species that have no RS information, except for those that appear in PREDICTS. Ideally, species that have no RS information should be checked manually
# # for taxonomic replication....
# TMammals <- Filter_No_range(NTR = NTRMammals, TraitDF = Mammals)
# TraitsMammals <- TMammals$TraitDF
# ExcludedMammals <- TMammals$Removed
# 
# TAmphibians <- Filter_No_range(NTR = NTRAmphibians, TraitDF = Amphibians)
# TraitsAmphibians <- TAmphibians$TraitDF
# ExcludedAmphibians <- TAmphibians$Removed
# 
# TBirds <- Filter_No_range(NTR = NTRBirds, TraitDF = Birds)
# TraitsBirds <- TBirds$TraitDF
# ExcludedBirds <- TBirds$Removed
# 
# TReptiles <- Filter_No_range(NTR = NTRReptiles, TraitDF = Reptiles)
# TraitsReptiles <- TReptiles$TraitDF
# ExcludedReptiles <- TReptiles$Removed


# Add family and order information to corrected datasets
Mammals <- Add_family_order_genus(Mammals, Syn_Mammals)
Amphibians <- Add_family_order_genus(Amphibians, Syn_Amphibians)
Reptiles <- Add_family_order_genus(Reptiles, Syn_Reptiles)
Birds <- Add_family_order_genus(Birds, Syn_Birds)


## Also remove species for which completeness = 0 for all predictor traits
## NB none of these species figure in PREDICTS
Comp.C_Mammals <- Completeness(Mammals, TraitsMammals)
Comp.C_Reptiles <- Completeness(Reptiles, TraitsReptiles)
Comp.C_Birds <- Completeness(Birds, TraitsBirds)
Comp.C_Amphibians <- Completeness(Amphibians, TraitsAmphibians)
Amphibians <- Amphibians %>% filter(Best_guess_binomial %nin% Comp.C_Amphibians) ## none of these species figure in PREDICTS


## write files
write.csv(Mammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Mammals.csv", row.names = FALSE)
write.csv(Amphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Amphibians.csv", row.names = FALSE)
write.csv(Reptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Reptiles.csv", row.names = FALSE)
write.csv(Birds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/2.filtered/Birds.csv", row.names=FALSE)

write.csv(U_Mammals, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Mammals.csv", row.names = FALSE)
write.csv(U_Amphibians, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Amphibians.csv", row.names = FALSE)
write.csv(U_Reptiles, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Reptiles.csv", row.names = FALSE)
write.csv(U_Birds, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/2.filtered/Birds.csv", row.names=FALSE)



# write.csv(ExcludedMammals, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/ExcludedMammals.csv", row.names = FALSE)
# write.csv(ExcludedAmphibians, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/ExcludedAmphibians.csv", row.names = FALSE)
# write.csv(ExcludedReptiles, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/ExcludedReptiles.csv", row.names = FALSE)
# write.csv(ExcludedBirds, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/ExcludedBirds.csv", row.names=FALSE)

