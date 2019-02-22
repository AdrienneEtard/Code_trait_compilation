## Performance of missForest VS phylopars to impute missing trait values -- continuous traits.

## Preamble
X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "data.table")
lapply(X, library, character.only=TRUE); rm(X)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}

## Load trait data 
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/MammalsComplete.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/BirdsComplete.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/AmphibiansComplete.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/ReptilesComplete.csv")

## Load phylogenies
Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/PhyloMammals.nwk") %>% .Format_tiplabels()
Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/PhyloAmphibians.nwk") %>% .Format_tiplabels() 
Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/PhyloReptiles.nwk") %>% .Format_tiplabels()
Phylo_Birds <- read.newick("../../Results/1.Phylogenies/PhyloBirds.nwk") %>% .Format_tiplabels() 



## Continuous traits: body mass, generation length, litter/clutch size, range size, habitat breadth


# ## TODO ## 
# Phylopars_error <- function(TraitDF, Phylo, Continuous_trait, N) {
#   
#   # 1. Assess the proportion of missing values for that trait.
#   SubTrait <- TraitDF[, c("Best_guess_binomial", Continuous_trait)]
#   SubTrait <-  subset(SubTrait, !is.na(SubTrait[, Continuous_trait]))
#   
#   Percent_NA <- (nrow(TraitDF)-(nrow(SubTrait)))/nrow(TraitDF)
#   print(paste("For the trait", Continuous_trait, "the percentage of missing values is", Percent_NA*100, "%."))
#   
#   # 2. Select all non missing values. Remove some at random.
#   ToRm <- sample(1:nrow(SubTrait), ceiling(Percent_NA*nrow(SubTrait)))
#   SubTrait$MAR_Trait <- SubTrait[, Continuous_trait]
#   SubTrait$MAR_Trait[ToRm] <- NA
#   
#   # 3. Impute missing values N times.
#   
#   # Are there species with missing values that are not in the phylogeny? If yes,
#   # attach these species to their genus in the phylogeny when possible.
#   NA_species <- SubTrait$Best_guess_binomial[is.na(SubTrait$MAR_Trait)]
#   ToAttach <- setdiff(NA_species, Phylo$tip.label)
#   print(paste("Attaching the following species randomly to their genus in the phylogeny:", ToAttach))
#   
#   ToAttach <- word(ToAttach, 1, 2, sep=" ")
#   
#   Phylo <- add.species.to.genus(Phylo, ToAttach, genus=NULL, "random")
#   
#   # Prune taxa
#   row.names(SubTrait) <- SubTrait$Best_guess_binomial
#   Match <- match.phylo.data(Phylo, SubTrait)
#   Phylo <- Match$phy
#   SubTrait <- Match$data
#   
#   browser()
#   
#   # Impute missing trait values, N times - with the function written by Molina-Venegas et al.
#   Names <- names(SubTrait[, "MAR_Trait"])
#   ToImpute <- SubTrait[, "MAR_Trait"] %>% 
#     as.character() %>%
#     as.numeric()
#   names(ToImpute) <- Names
# 
#   Imp1 <- tip_accuracy(Phylo, ToImpute, runs=10, method = "Rphylopars", threshold=0.75, range.lambda = 0.01)
# 
#   return(Imp1)
#   
# }



## TODO: species that do not have a match in the phylogeny should be imputed without phylogenetic eigenvectors 

missForest_error <- function (TraitDF, Phylo, Predictor_Traits, Trait_to_impute, N_eigen) {
  
  # 1. Assess the proportion of missing values for the trait to impute.
  SubTrait <-  subset(TraitDF, !is.na(TraitDF[, Trait_to_impute]))
  Percent_NA <- (nrow(TraitDF)-(nrow(SubTrait)))/nrow(TraitDF)
  print(paste("For the trait", Trait_to_impute, "the percentage of missing values is", Percent_NA*100, "%."))
  
  # 2. Select all non missing values. Remove some at random.
  
  # Create a new dataframe on which to impute
  DF_ToImpute <- SubTrait[, c("Best_guess_binomial", Predictor_Traits)]
  
  # Add to this DF the trait of interest with removed values
  ToRm <- sample(1:nrow(SubTrait), ceiling(Percent_NA*nrow(SubTrait)))
  MAR_Trait <- SubTrait[, Trait_to_impute]
  MAR_Trait[ToRm] <- NA
  
  DF_ToImpute$MAR_Trait <- MAR_Trait
  
  # 3. Impute the missing values with missForest, using all other traits as predictors + phylogenetic eigenvectors as predictors
  
  # Copy of original "DF to impute" for later use
  To_Impute_Original <- DF_ToImpute
  
  # Remove original trait values
  DF_ToImpute <- DF_ToImpute[, -which(colnames(DF_ToImpute)==Trait_to_impute)]
  
  # 3.1. Prune taxa: if phylogenetic information is included (if N_eigen != 0)
  
  if (N_eigen != 0) {
    
  row.names(DF_ToImpute) <- DF_ToImpute$Best_guess_binomial
  Match <- match.phylo.data(Phylo, DF_ToImpute)
  Phylo <- Match$phy
  DF_ToImpute <- Match$data
  
  # Get rid of replicated lines - should be resolved if no replicates in the phylogeny
  # DF_ToImpute <- DF_ToImpute[!duplicated(DF_ToImpute),]
  
  # Phylogenetic eigenvectors
  print("Eigenvector decomposition.")
  
  EigenV <- PVRdecomp(Phylo)
  Eigenvectors <- EigenV@Eigen$vectors
  Eigenvectors <- as.data.frame(Eigenvectors)
  Eigenvectors <- Eigenvectors[, 1:N_eigen]
  DF_ToImpute <- cbind(DF_ToImpute, Eigenvectors)
    
  }
  
  # Imputations (with or without phylogenetic eigenvectors as predictors)

  # Remove Best_guess_binomial
  rownames(DF_ToImpute) <- DF_ToImpute$Best_guess_binomial
  DF_ToImpute <- DF_ToImpute[, -which(colnames(DF_ToImpute)=="Best_guess_binomial")] 
  
  # Set classes 
  
  # Diet breadth as factor with 11 levels (1 to 11) 
  if (any(colnames(DF_ToImpute) %in% "Diet_breadth")) {
    DF_ToImpute[, "Diet_breadth"] <- as.factor(DF_ToImpute[, "Diet_breadth"])
    levels(DF_ToImpute[, "Diet_breadth"]) <- c(1:11)
  }
  
  # Continuous traits as numeric
  for (x in c("Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "MAR_Trait")) {
    if (any(colnames(DF_ToImpute)== x)) {
      DF_ToImpute[, x] <- DF_ToImpute[, x] %>%
        as.character() %>%
        as.numeric()
    }
  }

  # Impute
  browser()
  Imputed <- missForest(DF_ToImpute, variablewise = FALSE)
  Imputed <- Imputed$ximp
  Imputed$Best_guess_binomial <- gsub(".", " ", row.names(Imputed), fixed=TRUE)
  
  
  # 3.2. Imputations for the species that did not match with the phylogeny - for which no estimate of phylogenetic eigenvector
  
  if (N_eigen != 0) {
    
    Species <- setdiff(To_Impute_Original$Best_guess_binomial, Imputed$Best_guess_binomial)
  
    if (length(Species)!= 0){
      
      DF_ToImpute_2 <- subset(To_Impute_Original, Best_guess_binomial %in% Species)
      
      # Remove Best_guess_binomial and original trait
      DF_ToImpute_2 <- DF_ToImpute_2[, -c(which(colnames(DF_ToImpute_2)=="Best_guess_binomial"),
                                          which(colnames(DF_ToImpute_2)==Trait_to_impute))] 
      
      # Set classes 
      
      # Diet breadth as factor with 11 levels (1 to 11) 
      if (any(colnames(DF_ToImpute_2) %in% "Diet_breadth")) {
        DF_ToImpute_2[, "Diet_breadth"] <- as.factor(DF_ToImpute_2[, "Diet_breadth"])
        levels(DF_ToImpute_2[, "Diet_breadth"]) <- c(1:11)
      }
      
      # Continuous traits as numeric
      for (x in c("Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "MAR_Trait")) {
        if (any(colnames(DF_ToImpute_2)== x)) {
          DF_ToImpute_2[, x] <- DF_ToImpute_2[, x] %>%
            as.character() %>%
            as.numeric()
        }
      }
      
      # Impute
      Imputed_2 <- missForest(DF_ToImpute_2, variablewise = FALSE)
      Imputed_2 <- Imputed_2$ximp
      Imputed_2$Best_guess_binomial <- gsub(".", " ", row.names(Imputed_2), fixed=TRUE)
      
      
      # 3.3. rbind both datasets + new column for whether phylogenetic information was used or not
      
      Imputed <- Imputed %>% select(-Eigenvectors)
      
      # Imputed <- Imputed[, -which(colnames(DF_ToImpute_2)=="Eigenvectors")]
      Imputed$Phylo_info <- TRUE
      Imputed_2$Phylo_info <- FALSE
      Imputed_final <- rbind(Imputed, Imputed_2)
      
    }
  } 
  
  else {
    Imputed_final <- Imputed 
    Imputed_final$Phylo_info <- FALSE 
  }

  # 4. Add both original + predicted values in a DF to be returned
  Imputed_final <- Imputed_final[order(Imputed_final$Best_guess_binomial), ] 
  To_Impute_Original <- To_Impute_Original[order(To_Impute_Original$Best_guess_binomial), ] 
  
  To_return <- To_Impute_Original[, c("Best_guess_binomial", "MAR_Trait")]
  To_return$Original <- To_Impute_Original[, which(colnames(To_Impute_Original)==Trait_to_impute)]
  To_return$Imputed <- Imputed_final$MAR_Trait
  To_return$Phylo_info <- Imputed_final[, c("Phylo_info")]
  # Keep only the predictions for which NA were introduced randomly
  To_return <- subset(To_return, is.na(MAR_Trait))
  To_return <- To_return %>% select(-MAR_Trait)
  
  return(To_return)
  
}




## runs



Predictors <- c("Body_mass_g", 
                 "Generation_length_d", 
                 "Litter_size", 
                 "Range_size_m2", 
                 "Habitat_breadth_IUCN", 
                 "Trophic_level", 
                 "Specialisation",
                 "Diet_breadth")

Test_0 <- missForest_error(Mammals, Phylo_Mammals, Predictors, "Body_mass_g", 0)
Test_10 <- missForest_error(Mammals, Phylo_Mammals, Predictors, "Body_mass_g", 10)


# TODO:  

plot(Test$Original, Test$Imputed, pch=19)
abline(a=0, b=1)
Error <- Test$Original - Test$Imputed
Test$Error <- Error
Test <- Test[order(Test$Error),]
plot(Test$Error, c(1:84), pch=19)
abline(v=0, col="red")


