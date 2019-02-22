## Format phylogeny tip labels
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}


## Function to impute missing trait values using missForest

## ARG:
# Formatted phylogeny (class phylo)
# Formatted trait dataset (class dataframe)
# Taxinfo: genus, family, order
# Traits to impute (character string): continuous and categorical + eigenvectors to retain
# ErrorTrue: error format returned by missForest (variablewise: TRUE or FALSE)

## RETURNS:
# Imputed trait values (dataframe)
# OOB errors from the imputations

Imputations_missForest <- function (TraitDF, Taxinfo, Traits_cont, Traits_cat, EV, ErrorTrue, DietTRUE) {
  
# browser()
  
  ## Select traits of interest, to impute, and phylogenetic eigenvectors, and taxinfo 
  To_impute <- TraitDF[, colnames(TraitDF) %in% c(Taxinfo, Traits_cont, Traits_cat, EV)]
  rownames(To_impute) <- TraitDF$Best_guess_binomial
  
  ## Set habitat and diet variables as binary (TRUE for 1 and FALSE for 0)
  
  ## Set class as numeric for all continuous traits
  To_impute[, Traits_cont] <- apply(To_impute[, Traits_cont], 2, as.numeric)
  
  ## Set habitat and diet variables as factors. Otherwise, when they are numeric or logical, outputs are numeric
  
  if(DietTRUE){
    
    Diet <- c("IN", "VE", "PL", "SE", "NE", "FR", "SCV")
    for (i in Diet) {
      To_impute[,i] <- factor(To_impute[,i], levels=c(0,1))
      
    }
  }
  
  
  Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
               "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
               "Artificial","Introduced.vegetation","Other.Unknown")
  
  for (i in Habitat) {
    To_impute[,i] <- factor(To_impute[,i], levels=c(0,1))
    
  }
  
  
  ## Impute missing values
  print("Imputing missing values.")
  
  if (isTRUE(ErrorTrue)) { R.Imputed <- missForest(To_impute, variablewise = TRUE) } 
  else { R.Imputed <- missForest(To_impute, variablewise = FALSE) }
  
  Imputed <- R.Imputed$ximp
  Errors <- R.Imputed$OOBerror
  
  ## Select traits and variables of interest
  Continuous <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", 
                    "Range_size_m2", "sqrt_Habitat_breadth_IUCN")
  
  if(DietTRUE) {
  Categorical <- c(Habitat, "Specialisation", "Diel_activity","Trophic_level", Diet)}
  
  else{Categorical <- c(Habitat, "Specialisation", "Diel_activity","Trophic_level")}
  
  Imputed$Best_guess_binomial <- rownames(Imputed)
  Imputed <- Imputed[order(Imputed$Best_guess_binomial), c(Taxinfo, "Best_guess_binomial", Continuous, Categorical, "EV_1")]
  
  ## Add a column for phylogenetic information (yes or no)
  Imputed$Phylo_info <- TraitDF$EV_1
  Imputed$Phylo_info[!is.na(Imputed$Phylo_info)] <- "YES"
  Imputed$Phylo_info[is.na(Imputed$Phylo_info)] <- "NO"
  
  ## Return
  return(list(Imputed.Values=Imputed, Imputations.Errors=Errors))
  
}




## Function to impute missing values for species that are not represented in the phylogeny 

# Imputations_missForest_other <- function (OriginalTraitDF, ImputedDF, TraitsCont, TraitsCat, Habitat, ErrorTrue) {
#   
#   ## Get species that are not in the imputed dataset, select traits and add to the imputed dataset
#   ToAdd <- setdiff(OriginalTraitDF$Best_guess_binomial, rownames(ImputedDF))
#   ToAdd <- subset(OriginalTraitDF, Best_guess_binomial %in% ToAdd)
#   rownames(ToAdd) <- ToAdd$Best_guess_binomial
#   
#   ## select traits to impute
#   ToAdd <- ToAdd[, colnames(ToAdd) %in% c(TraitsCont, TraitsCat, Habitat)]
#   
#   # ## Set Diet breadth as factor 
#   # if (any(colnames(ToAdd) %in% "Diet_breadth")) {
#   #   Upper <- max(ImputedDF[, "Diet_breadth"], na.rm=T)
#   #   ToAdd[, "Diet_breadth"] <- as.factor(ToAdd[, "Diet_breadth"])
#   #   levels(ToAdd[, "Diet_breadth"]) <- c(1:Upper)
#   # }
#   
#   ToImpute <- rbind(ImputedDF, ToAdd)
#   
#   ## Impute with missForest
#   print("IMPUTATIONS.")
#   if (isTRUE(ErrorTrue)) {
#     Imputed <- missForest(ToImpute, variablewise = TRUE)
#   } else {Imputed <- missForest(ToImpute, variablewise = FALSE)}
#   
#   ImputedValues <- Imputed$ximp
#   ImputedErrors <- Imputed$OOBerror
#   
#   ## Species for which there was no phylogenetic information:
#   # rownames(ToAdd)
#   
#   ## RETURN
#   return(list(ImputedValues=ImputedValues, ImputedErrors=ImputedErrors, Sp_No_phylo_info=rownames(ToAdd)))
#   
# }
# 




# Imputations_missForest <- function (Phylo, TraitDF, Traits_cont, Trait_cat, Habitat, N, ErrorTrue) {
#   
#   browser()
#   
#   ## Prune species that do not intersect
#   row.names(TraitDF) <- TraitDF$Best_guess_binomial
#   Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
#   Phylo <- Prune_Taxa$phy
#   TraitDF <- Prune_Taxa$data
#   
#   ## Select traits of interest
#   TraitDF <- TraitDF[, colnames(TraitDF) %in% c(Traits_cont, Trait_cat, Habitat)]
#   
#   ## Set class as numeric for continuous
#   for (i in Traits_cont) {TraitDF[, i] <- as.numeric(as.character(TraitDF[, i]))}
#   
#   # ## Set Diet breadth as factor 
#   # if (any(colnames(TraitDF) %in% "Diet_breadth")) {
#   #   Upper <- max(TraitDF[, "Diet_breadth"], na.rm=T)
#   #   TraitDF[, "Diet_breadth"] <- as.factor(TraitDF[, "Diet_breadth"])
#   #   levels(TraitDF[, "Diet_breadth"]) <- c(1:Upper)
#   # }
#   
#   ## Get phylogenetic eigenvectors from the phylogeny and select N first eigenvectors
#   print("EIGENVECTOR DECOMPOSITION.")
#   EigenV <- PVRdecomp(Phylo)
#   Eigenvectors <- EigenV@Eigen$vectors
#   
#   Eigenvectors <- as.data.frame(Eigenvectors)
#   Eigenvectors <- Eigenvectors[, 1:N]
#   
#   ## Bind eigenvectors to trait dataframe, to use them as predictors
#   ToImpute <- cbind(TraitDF, Eigenvectors)
#   
#   ## Impute
#   print("IMPUTATIONS.")
#   
#   if (isTRUE(ErrorTrue)) {
#   Imputed <- missForest(ToImpute, variablewise = TRUE)
#   } else {
#   Imputed <- missForest(ToImpute, variablewise = FALSE)}
#   
#   ImputedValues <- Imputed$ximp
#   ImputedErrors <- Imputed$OOBerror
#   
#   ## Remove eigenvectors from imputed dataset
#   ImputedValues <- ImputedValues[, colnames(ImputedValues) %in% c(Traits_cont, Trait_cat, Habitat)]
# 
#   ## RETURN
#   return(list(ImputedValues=ImputedValues, ImputedErrors=ImputedErrors))
#   
# }
# 
# 
# ## Function to impute missing values for species that are not represented in the phylogeny 
# 
# Imputations_missForest_other <- function (OriginalTraitDF, ImputedDF, TraitsCont, TraitsCat, Habitat, ErrorTrue) {
# 
#   ## Get species that are not in the imputed dataset, select traits and add to the imputed dataset
#   ToAdd <- setdiff(OriginalTraitDF$Best_guess_binomial, rownames(ImputedDF))
#   ToAdd <- subset(OriginalTraitDF, Best_guess_binomial %in% ToAdd)
#   rownames(ToAdd) <- ToAdd$Best_guess_binomial
#   
#   ## select traits to impute
#   ToAdd <- ToAdd[, colnames(ToAdd) %in% c(TraitsCont, TraitsCat, Habitat)]
#   
#   # ## Set Diet breadth as factor 
#   # if (any(colnames(ToAdd) %in% "Diet_breadth")) {
#   #   Upper <- max(ImputedDF[, "Diet_breadth"], na.rm=T)
#   #   ToAdd[, "Diet_breadth"] <- as.factor(ToAdd[, "Diet_breadth"])
#   #   levels(ToAdd[, "Diet_breadth"]) <- c(1:Upper)
#   # }
#   
#   ToImpute <- rbind(ImputedDF, ToAdd)
#   
#   ## Impute with missForest
#   print("IMPUTATIONS.")
#   if (isTRUE(ErrorTrue)) {
#     Imputed <- missForest(ToImpute, variablewise = TRUE)
#   } else {Imputed <- missForest(ToImpute, variablewise = FALSE)}
#   
#   ImputedValues <- Imputed$ximp
#   ImputedErrors <- Imputed$OOBerror
#   
#   ## Species for which there was no phylogenetic information:
#   # rownames(ToAdd)
#   
#   ## RETURN
#   return(list(ImputedValues=ImputedValues, ImputedErrors=ImputedErrors, Sp_No_phylo_info=rownames(ToAdd)))
#   
# }


# NUMBER OF EIGENVECTORS --------------------------------------------------

## function to prune phylo and trait dataset and to get set of phylogenetic eigenvectors

# Phylo_eigenvectors <- function(Phylo, TraitDF) {
#   
#   ## Prune species that do not intersect
#   row.names(TraitDF) <- TraitDF$Best_guess_binomial
#   Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
#   Phylo <- Prune_Taxa$phy
#   TraitDF <- Prune_Taxa$data
#   
#   ## Get phylogenetic eigenvectors from the phylogeny
#   EigenV <- PVRdecomp(Phylo)
#   Eigenvectors <- EigenV@Eigen$vectors
#   
#   Eigenvectors <- as.data.frame(Eigenvectors)
#   
#   return(Eigenvectors)
# }


# Imputations_missForest_Error <- function (Phylo, Phylo_eigenv, TraitDF, Traits_cont, Trait_cat, ErrorTrue, Seq) {
#   
#   ## Prune species that do not intersect
#   row.names(TraitDF) <- TraitDF$Best_guess_binomial
#   Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
#   Phylo <- Prune_Taxa$phy
#   TraitDF <- Prune_Taxa$data
#   
#   ## Select traits of interest
#   TraitDF <- TraitDF[, colnames(TraitDF) %in% c(Traits_cont, Trait_cat)]
#   
#   ## Set class as numeric for continuous
#   for (i in Traits_cont) {TraitDF[, i] <- as.numeric(as.character(TraitDF[, i]))}
#   
#   ## Set Diet breadth as factor with 11 levels (1 to 11) 
#   if (any(colnames(TraitDF) %in% "Diet_breadth")) {
#     TraitDF[, "Diet_breadth"] <- as.factor(TraitDF[, "Diet_breadth"])
#     levels(TraitDF[, "Diet_breadth"]) <- c(1:11)
#   }
# 
#   
#   ## Dataset to store error results
#   ErrorDF <- as.data.frame(Seq)
#   
#   
#   for (i in 1:length(Seq)) {
#     
#       ## Bind eigenvectors to trait dataframe, to use them as predictors
#       ToImpute <- cbind(TraitDF, Phylo_eigenv[1:Seq[i]])
#       
#       
#       ## Impute
#       if (isTRUE(ErrorTrue)) {
#         Imputed <- missForest(ToImpute, variablewise = TRUE)
#       } else {Imputed <- missForest(ToImpute, variablewise = FALSE)}
#       
#      
#       ErrorDF$NRMSE[i] <- Imputed$OOBerror["NRMSE"]
#       ErrorDF$PFC[i] <- Imputed$OOBerror["PFC"]
#   
#       cat("\n \n", i, "\n \n of \n \n", length(Seq), "\n \n")
#       
#   }
#   
#   ## RETURN
#   return(ErrorDF)
#   
# }
# 
# ###################
# Imputations_missForest_ErrorSE <- function (Phylo, Phylo_eigenv, TraitDF, Traits_cont, Trait_cat, ErrorTrue, Seq, Rep) {
#   
#   ## Prune species that do not intersect
#   row.names(TraitDF) <- TraitDF$Best_guess_binomial
#   Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
#   Phylo <- Prune_Taxa$phy
#   TraitDF <- Prune_Taxa$data
#   
#   ## Select traits of interest
#   TraitDF <- TraitDF[, colnames(TraitDF) %in% c(Traits_cont, Trait_cat)]
#   
#   ## Set class as numeric for continuous
#   for (i in Traits_cont) {TraitDF[, i] <- as.numeric(as.character(TraitDF[, i]))}
#   
#   ## Set Diet breadth as factor with 11 levels (1 to 11) 
#   if (any(colnames(TraitDF) %in% "Diet_breadth")) {
#     TraitDF[, "Diet_breadth"] <- as.factor(TraitDF[, "Diet_breadth"])
#     levels(TraitDF[, "Diet_breadth"]) <- c(1:11)
#   }
#   
#   
#   ## Dataset to store error results
#   ErrorDF <- as.data.frame(Seq)
# 
#   
#   for (i in 1:length(Seq)) {
#     
#     ## Bind eigenvectors to trait dataframe, to use them as predictors
#     ToImpute <- cbind(TraitDF, Phylo_eigenv[1:Seq[i]])
#     
#     ErrorRep <- as.data.frame(c(1:Rep))
#     
#      for (j in 1:Rep) {
#            ## Impute
#       if (isTRUE(ErrorTrue)) {
#         Imputed <- missForest(ToImpute, variablewise = TRUE)
#       } else {Imputed <- missForest(ToImpute, variablewise = FALSE)}
#     
#        ErrorRep$NRMSE[j] <- Imputed$OOBerror["NRMSE"]
#        ErrorRep$PFC[j] <- Imputed$OOBerror["PFC"]
#      }
#     
#     ## Calculate mean and SE for each replicate and store in ErrorDF
#     ErrorDF$Mean_NRMSE[i] <- mean(ErrorRep$NRMSE)
#     ErrorDF$SE_NRMSE[i] <- std.error(ErrorRep$NRMSE)
#     
#     ErrorDF$Mean_PFC[i] <- mean(ErrorRep$PFC)
#     ErrorDF$SE_PFC[i] <- std.error(ErrorRep$PFC)
#     
#     cat("\n \n ", i, "\n \n of \n \n", length(Seq), "\n \n")
#     
#   }
#   
#   
#   ## RETURN
#   return(ErrorDF)
#   
# }
# 








