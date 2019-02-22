## TODO : 1. adapt to include phylogenetic information 2. adapt to include several traits at the same time 3. return imputations not just errors

missForest_error <- function (TraitDF, Phylo, Predictor_Traits, Trait_to_impute, Phylo_info, N_eigen, N_imputations) {
  
  # 1. Among known values, replace one at random by NA - so that the proportion of NAs for that trait is overall unchanged
  
  # Create a new dataframe on which to impute; replace the sampled value by NA in this dataframe
  DF_ToImpute <- TraitDF[, c("Best_guess_binomial", Predictor_Traits)]
  
  # Number of samples: such as no more variation than 0.5%
  x <- nrow(TraitDF)/200
  
  ToSample <- subset(TraitDF, !is.na(TraitDF[, Trait_to_impute]))
  ToSample <- rownames(ToSample) %>% as.numeric
  ToReplace <- sample(ToSample, x)
  DF_ToImpute[ToReplace, Trait_to_impute] <- NA
  
  
  
  # 2. Impute the missing values with missForest, using all other traits as predictors + phylogenetic eigenvectors as predictors if Phylo_info true
  
  # 2.1. Without phylogenetic information
  if (!Phylo_info) {
    
    # Remove Best_guess_binomial
    # rownames(DF_ToImpute) <- DF_ToImpute$Best_guess_binomial
    DF_ToImpute <- DF_ToImpute[, -which(colnames(DF_ToImpute)=="Best_guess_binomial")] 
    
    # Set classes 
    
    # Diet breadth as factor with 11 levels (1 to 11) 
    if (any(colnames(DF_ToImpute) %in% "Diet_breadth")) {
      DF_ToImpute[, "Diet_breadth"] <- as.factor(DF_ToImpute[, "Diet_breadth"])
      levels(DF_ToImpute[, "Diet_breadth"]) <- c(1:11)
    }
    
    # Continuous traits as numeric # to refine to include all kinds of predictors
    for (x in c("Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Body_mass_g")) {
      if (any(colnames(DF_ToImpute)== x)) {
        DF_ToImpute[, x] <- DF_ToImpute[, x] %>%
          as.character() %>%
          as.numeric()
      }
    }
    
    # Impute N_imputations times, saving the imputed dataframes in a list
    List <- list()
    
    for (i in 1:N_imputations) {
      print(paste("Imputation", i, "on", N_imputations))
      Imputed <- missForest(DF_ToImpute, variablewise = FALSE)
      Imputed <- Imputed$ximp
      # Imputed$Best_guess_binomial <- gsub(".", " ", row.names(Imputed), fixed=TRUE)
      # Imputed <-  Imputed[order(Imputed$Best_guess_binomial), ] 
      # TraitDF <- TraitDF[order(TraitDF$Best_guess_binomial), ]
      Imputed$Best_guess_binomial <- TraitDF$Best_guess_binomial
      Imputed$Original_trait_value <- TraitDF[, Trait_to_impute]
      Imputed$OfInterest <- NA
      Imputed$OfInterest[ToReplace] <- "Replaced"
      List[[i]] <- Imputed
    }
    
      
    # List of imputed dataframes to process: extract lines "Replaced"
    Errors <- lapply(List, function(x) x[, c(Trait_to_impute, "Original_trait_value", "OfInterest", "Best_guess_binomial")] %>% 
                         filter(OfInterest=="Replaced"))
    Errors <- rbindlist(Errors) %>% as.data.frame()
    Errors$Error <- abs(Errors$Original_trait_value - Errors[, Trait_to_impute])/Errors$Original_trait_value*100

    return(Errors)
  } 
  
  # 2.2. With phylogenetic information
  
}



ToRun_error <- function(TraitDF, Phylo, Predictor_Traits, Trait_to_impute, Phylo_info, N_eigen, N_iterations, N_imputations) {
  
  List <- list()
  
  for (i in 1:N_iterations) {
    print(paste("Iteration", i, "on", N_iterations))
    List[[i]] <- missForest_error(TraitDF, Phylo, Predictor_Traits, Trait_to_impute, Phylo_info, N_eigen, N_imputations)
  }
  
  Results_error <- rbindlist(List)
  
  
  return(Results_error)
}





Predictors <- c("Body_mass_g", 
                "Generation_length_d", 
                "Litter_size", 
                "Range_size_m2", 
                "Habitat_breadth_IUCN", 
                "Trophic_level", 
                "Specialisation",
                "Diet_breadth")


Test <- ToRun_error(Mammals, Phylo_Mammals, Predictors, "Generation_length_d", FALSE, 0, 5, 10)

# list.save(Test, "../../Results/ListMissForest.rds")
Test2 <- readRDS("../../Results/ListMissForest.rds")

unique(Test$Best_guess_binomial) %>% length()
Errors <- Test2 %>% group_by(Best_guess_binomial) %>%
  summarise(MeanError=mean(Error)) %>%
  as.data.frame()
boxplot(Errors$MeanError)
plot(density(Errors$MeanError))
abline(v=mean(Errors$MeanError))
abline(v=median(Errors$MeanError), col="red")

