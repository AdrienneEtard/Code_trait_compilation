# AIMS:
# 1. Interpolation from correlations between variables (Function Corr_Longevity)
#  log(Body_mass) ~ log(Body_length) + Genus + log(Body_length):Genus 
#  log(Longevity) ~ log(Max_longevity) + Genus + log(Max_longevity):Genus 
#  Also trying with Family as factor with interaction (but it does not increase coverage) 
#
#  Actually, no matches between species known by binomial name for which a prediction is needed, and
#  species for which predictions are made. So no interpolation using correlation on body mass or longevity for mammals
#  Will use genus /family / order averages instead.
#
# 2. Interpolations using genus/family/order averages (Function Add_Averages)



# Function to drop unknown factor levels for lm predictions ---------------------------------------------------
# Found on https://stackoverflow.com/questions/4285214/predict-lm-with-an-unknown-factor-level-in-test-data
# Arg: output from lm; data: DF

.missingLevelsToNA<-function(object, data){

  # Obtain factor predictors in the model and their levels
  factors <- (gsub("[-^0-9]|as.factor|\\(|\\)", "", names(unlist(object$xlevels))))
  factorLevels <- unname(unlist(object$xlevels))
  modelFactors <- as.data.frame(cbind(factors,factorLevels))

  # Select column names in your data that are factor predictors in your model
  predictors <- names(data[names(data) %in% factors])

  #For each factor predictor in your data if the level is not in the model set the value to NA
  for (i in 1:length(predictors)){
    found <- data[,predictors[i]] %in% modelFactors[modelFactors$factors==predictors[i],]$factorLevels
    if (any(!found)) data[!found,predictors[i]]<-NA
  }

  return(data)

}



# Corr_Longevity function ------------------------------------------------------------------------------

# Predictions Longevity using Myhrvold dataset
# NB: not using Pacifici for mammals for correlations on Longevity 
# because it has already been done within Pacifici (won't increase coverage)

Corr_Longevity <- function(Myhrvold, Traits, VClass) {

  # Preparing the dataset to store predictions
  Myhrvold_VClass <- subset(Myhrvold, class == VClass)

  # Predicting Longevity from max longevity, Myhrvold
  Mod_GL <- lm(log(Generation_length_d)~log(maximum_longevity_y) + Genus + log(maximum_longevity_y):Genus, data=Myhrvold_VClass)
  # Mod_GL <- lm(log(Generation_length_d)~log(maximum_longevity_y) + family + log(maximum_longevity_y):family, data=Myhrvold_VClass)
  # summary(Mod_GL)$r.squared # 0.991 using genus as factor, 0.979 using family as factor

  MissingGenus <- .missingLevelsToNA(Mod_GL, Myhrvold_VClass)
  Myhrvold_VClass$Predicted_generation_length_d <- exp(predict(Mod_GL, newdata = MissingGenus))
  Myhrvold_VClass$Generation_length_d[is.na(Myhrvold_VClass$Generation_length_d)] <- Myhrvold_VClass$Predicted_generation_length_d[is.na(Myhrvold_VClass$Generation_length_d)]

  # Adding Longevity predictions to the main trait dataset:
  # See if any of these predictions match with missing values in the Traits dataset
  # then replace by the predictions

  # Extract predictions for species of interest
  Missing_GL <- subset(Traits, is.na(Generation_length_d))
  Myhrvold_sub <- select(Myhrvold_VClass, one_of("Binomial_name", "Generation_length_d")) %>% na.omit()
  Y <- intersect(Missing_GL$Best_guess_binomial, Myhrvold_sub$Binomial_name)

  Pred_GL <- subset(Myhrvold_sub, Binomial_name %in% Y) %>% na.omit()
  Pred_GL <- arrange(Pred_GL, Binomial_name)

  # Replace values in the Traits dataset
  Traits$Generation_length_d[match(Pred_GL$Binomial_name, Traits$Best_guess_binomial)] <- Pred_GL$Generation_length_d

  return(list(Traits=Traits, Mammals.Myhrvold=Myhrvold_VClass))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Body mass ~ Body length, mammals -----------------------------------------------------------------------------

# ## The following commented script tested if correlations on Body mass
# ## were interesting to increase coverage. It does not increase coverage
# ## as species for which BM is predicted are not the same species for which
# ## the data is lacking in Predicts.
# 
# ## WILL USE AVERAGES BY GENUS INSTEAD
# 
# ## Predicting body mass from adult svl, Myhrvold dataset --------------------------------------------------------
# VClass <- "Mammalia"
# Myhrvold_VClass <- subset(Myhrvold, class == VClass)
# 
# Mod <- lm(log(Body_mass_g)~log(Adult_svl_cm)+ family + family:log(Adult_svl_cm), data=Myhrvold_VClass)
# # Mod <- lm(log(Body_mass_g)~log(Adult_svl_cm)+ genus + genus:log(Adult_svl_cm), data=Myhrvold_VClass)
# # summary(Mod)$r.squared: 0.993 with genus; 0.98 with family
# MissingGenus <- .missingLevelsToNA(Mod, Myhrvold_VClass)
# Predictions_Myhrvold <- Myhrvold_VClass[, "Binomial_name"] %>% as.data.frame()
# colnames(Predictions_Myhrvold) <- "Binomial_name"
# Predictions_Myhrvold$Pred_BM_Myhrvold <- exp(predict(Mod, newdata = MissingGenus))
# rm(Mod, MissingGenus)
# 
# # # Predicting body mass using Pantheria -------------------------------------------------------------------------
# Predictions_Pantheria <- Pantheria[, "Binomial_name"] %>% as.data.frame()
# colnames(Predictions_Pantheria) <- "Binomial_name"
# 
# # # Predicting Body mass from Head length and Forearm legnth(Pantheria) ------------------------------------------------------------
# 
# # Mod_1a <- lm(log(Body_mass_g)~log(X8.1_AdultForearmLen_mm)+ Genus + Genus:log(X8.1_AdultForearmLen_mm), data=Pantheria)
# # Mod_1b <- lm(log(Body_mass_g)~log(Head_length_mm)+ Genus + Genus:log(Head_length_mm), data=Pantheria)
# 
# # MissingGenus_1a <- .missingLevelsToNA(Mod_1a, Pantheria)
# # MissingGenus_1b <- .missingLevelsToNA(Mod_1b, Pantheria)
# 
# # Predictions_Pantheria$Pred_BM_Pantheria_Arm <- exp(predict(Mod_1a, newdata = MissingGenus_1a))
# # Predictions_Pantheria$Pred_BM_Pantheria_Head <- exp(predict(Mod_1b, newdata = MissingGenus_1b))
# 
# 
# Mod_2a <- lm(log(Body_mass_g)~log(X8.1_AdultForearmLen_mm)+ Family + Family:log(X8.1_AdultForearmLen_mm), data=Pantheria)
# Mod_2b <- lm(log(Body_mass_g)~log(Head_length_mm)+ Family + Family:log(Head_length_mm), data=Pantheria)
# 
# MissingGenus_2a <- .missingLevelsToNA(Mod_2a, Pantheria)
# MissingGenus_2b <- .missingLevelsToNA(Mod_2b, Pantheria)
# 
# Predictions_Pantheria$Pred_BM_Pantheria_Arm <- exp(predict(Mod_2a, newdata = MissingGenus_2a))
# Predictions_Pantheria$Pred_BM_Pantheria_Head <- exp(predict(Mod_2b, newdata = MissingGenus_2b))
# 
# 
# # Averaging Body mass predictions from the two datasets -------------------
# Pred_BM <- merge(Predictions_Myhrvold, Predictions_Pantheria, all=TRUE)
# # Test <- Pred_BM[ , .(Mean_sources = rowMeans(.SD, na.rm = TRUE)), by = Binomial_name]
# Pred_BM$Mean <- apply(Pred_BM[,c(2,3,4)], 1, mean, na.rm=TRUE)
# 
# # # Adding the predictions to the Trait dataset -----------------------------
# Missing_BM <- subset(Traits, is.na(Body_mass_g)) %>% filter(!is.na(Best_guess_binomial))
# Y <- intersect(Pred_BM$Binomial_name, Missing_BM$Best_guess_binomial)
# # Y IS LENGTH 0: no matches!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Interpolations at Genus/Family/Order level for continuous traits --------

# This function averages on continuous traits at a given taxonomic level for a given trait dataset
# The output for a given group and trait can be used as a proxy for species that
# require trait information and for which no correlation with other traits is possible
# or for species that are unknown by binomial names

.Av_by_taxon <- function(Dataset, Traits, Taxon_levels) {

  # Taxons_levels: c("Order", "Family", "Genus")
  
  # Get averages per taxonomic group for each Taxon_levels
  Df1 <- select(Dataset, one_of(Traits, Taxon_levels[1]))
  AvOrder <- setDT(Df1)[ , lapply(.SD, mean, na.rm=TRUE), by=eval(Taxon_levels[1])]
  
  Df2 <- select(Dataset, one_of(Traits, Taxon_levels[2]))
  AvFam <- setDT(Df2)[ , lapply(.SD, mean, na.rm=T), by=eval(Taxon_levels[2])]
  
  Df3 <- select(Dataset, one_of(Traits, Taxon_levels[3]))
  AvGenus <- setDT(Df3)[ , lapply(.SD, mean, na.rm=T), by=eval(Taxon_levels[3])]
  
  return(list(Av_Order=AvOrder, Av_Family=AvFam, Av_Genus=AvGenus))

}

# NB: when datasets share information (ex: Myhrvold and Pantheria for both BM and litter size)
# I average estimates from the datasets at a given taxonomic level
.MergeAv <- function(Av1, Av2, Av3) {
  Av <- merge(Av1, Av2, all=T) %>% as.data.frame()
  Av <- merge(Av, as.data.frame(Av3), all=T)
  Av <- setDT(Av)[, lapply(.SD, mean, na.rm=T), by=eval(colnames(Av)[1])] %>%
    as.data.frame()
  return(Av)

}


# Add averages to Traits dataset ----------------------------------------------------------------

Add_Averages <- function(Dataset, AvGenus, AvFam, AvOrder, Cont_traits) {
  
  for (i in 1:length(Cont_traits)) { 
    
    ## Considered continuous trait
    Trait <- Cont_traits[i]
  
    ## If Genus is known, and Trait is NA, add genus averaged values
    # First subset by unknown Trait value and known genus 
    Data <- Dataset %>% subset(is.na(Dataset[,Trait])) %>%
      subset(Genus != "")
    # Genuses in the subset that match Genuses in the AvGenus dataset (get positions where to replace)
    Y <- match(Data$Genus, AvGenus$Genus)
    # Replace by values for genuses that match
    Data[Trait] <- AvGenus[Y, Trait]
    # Place trait values back in the trait dataset
    Index <- as.numeric(row.names(Data))
    Dataset[Index, Trait] <- Data[Trait]
    
    rm(Data, Y, Index)
    
    ## If Genus is unknown and Family is known, add family averaged values
    # First subset by unknown Trait value, unknown genuses and known families 
    Data <- Dataset %>% subset(is.na(Dataset[,Trait])) %>%
      subset(Genus == "") %>%
      subset(Family != "")
    # Genuses in the subset that match Genuses in the AvGenus dataset (get positions where to replace)
    Y <- match(Data$Family, AvFam$Family)
    # Replace by values for genuses that match
    Data[Trait] <- AvFam[Y, Trait]
    # Place trait values back in the trait dataset
    Index <- as.numeric(row.names(Data))
    Dataset[Index, Trait] <- Data[Trait]
     
    rm(Data, Y, Index)
    
    ## Finally, Genus and Family unknown: add order averaged values
    # First subset by unknown Trait value, unknown genuses, unknown families 
    Data <- Dataset %>% subset(is.na(Dataset[,Trait])) %>%
      subset(Genus == "") %>%
      subset(Family == "")
    # Genuses in the subset that match Genuses in the AvGenus dataset (get positions where to replace)
    Y <- match(Data$Order, AvOrder$Order)
    # Replace by values for genuses that match
    Data[Trait] <- AvOrder[Y, Trait]
    # Place trait values back in the trait dataset
    Index <- as.numeric(row.names(Data))
    Dataset[Index, Trait] <- Data[Trait]
    
  }
  
  return(Dataset)
}
