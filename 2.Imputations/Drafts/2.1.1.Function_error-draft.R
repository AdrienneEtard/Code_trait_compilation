missForest_error <- function (TraitDF, Phylo, Predictor_Traits, Trait_to_impute, Phylo_info, N_eigen, N_imputations) {
  
  # 1. Among known values, replace one at random by NA - so that the proportion of NAs for that trait is overall unchanged
  
  # Create a new dataframe on which to impute; replace the sampled value by NA in this dataframe
  DF_ToImpute <- TraitDF[, c("Best_guess_binomial", Predictor_Traits)]
  
  ToReplace <- sample(1:nrow(TraitDF[, !is.na(Trait_to_impute)]), 1)
  DF_ToImpute[ToReplace, Trait_to_impute] <- NA
  
  # 3. Impute the missing values with missForest, using all other traits as predictors + phylogenetic eigenvectors as predictors if Phylo_info true
  
  # 3.1. Without phylogenetic information
  if (!Phylo_info) {
    
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
    
    # Impute N times
    List_imputed <- list()
    for (i in 1:N) {
      Imputed <- missForest(DF_ToImpute, variablewise = FALSE)
      Imputed <- Imputed$ximp
      Imputed$Best_guess_binomial <- gsub(".", " ", row.names(Imputed), fixed=TRUE)
      Imputed <-  Imputed[order(Imputed$Best_guess_binomial), ] 
      Imputed$Original_trait_value <- TraitDF[, Trait_to_impute]
      List_imputed[[i]] <- Imputed
    }
    
    return(List_imputed)
    
  } # end without phylo info
  
  
  # 3.2. With phylogenetic information
  else {
    
    
  }
   
  
  
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
  
  # N imputations (with or without phylogenetic eigenvectors as predictors)
  
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
  
  # Impute N times
  List_imputed <- list()
  for (i in 1:N) {
  Imputed <- missForest(DF_ToImpute, variablewise = FALSE)
  Imputed <- Imputed$ximp
  Imputed$Best_guess_binomial <- gsub(".", " ", row.names(Imputed), fixed=TRUE)
  List_imputed[[i]] <- Imputed
  }
  
  # 3.2. Imputations for the species that did not match with the phylogeny - for which no estimate of phylogenetic eigenvector
  
  if (N_eigen != 0) {
    
    Species <- setdiff(TraitDF$Best_guess_binomial, Imputed$Best_guess_binomial)
    
    if (length(Species)!= 0){
      
      DF_ToImpute_2 <- subset(TraitDF, Best_guess_binomial %in% Species)
      
      # Remove Best_guess_binomial and original trait
      rownames(DF_ToImpute_2) <- DF_ToImpute_2$Best_guess_binomial
      DF_ToImpute_2 <- DF_ToImpute_2[, -which(colnames(DF_ToImpute_2)=="Best_guess_binomial")] 
      DF_ToImpute_2 <- DF_ToImpute_2[, Predictor_Traits)]
      
      
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
      List_imputed2 <- list()
      for (i in 1:N) {
      Imputed_2 <- missForest(DF_ToImpute_2, variablewise = FALSE)
      Imputed_2 <- Imputed_2$ximp
      Imputed_2$Best_guess_binomial <- gsub(".", " ", row.names(Imputed_2), fixed=TRUE)
      List_imputed2[[i]] <- Imputed
      }
      
      
      # 3.3. rbind datasets in the lists + new column for whether phylogenetic information was used or not
      List_imputation_results  <- 
      Imputed <- Imputed %>% select(-Eigenvectors)
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
  TraitDF <- TraitDF[order(TraitDF$Best_guess_binomial), ] 
  Imputed_final$Original_traitvalues <- TraitDF[, Trait_to_impute] 
  
  print(ToReplace)
  return(Imputed_final)
  
}


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

Diff <- Test_0$Body_mass_g[5300] - Test_0$Original_traitvalues[5300]

