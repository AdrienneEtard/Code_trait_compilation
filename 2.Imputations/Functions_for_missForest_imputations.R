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
  if(DietTRUE) {
    Continuous <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", 
                    "Range_size_m2", "sqrt_Habitat_breadth_IUCN", "sqrt_Diet_breadth")
  
    Categorical <- c(Habitat, "Specialisation", "Diel_activity","Trophic_level", Diet)
  }
  
  else{
    Continuous <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", 
                     "Range_size_m2", "sqrt_Habitat_breadth_IUCN")
    Categorical <- c(Habitat, "Specialisation", "Diel_activity","Trophic_level")}
  
  Imputed$Best_guess_binomial <- rownames(Imputed)
  Imputed <- Imputed[order(Imputed$Best_guess_binomial), c(Taxinfo, "Best_guess_binomial", Continuous, Categorical, "EV_1")]
  
  ## Add a column for phylogenetic information (yes or no)
  Imputed$Phylo_info <- TraitDF$EV_1
  Imputed$Phylo_info[!is.na(Imputed$Phylo_info)] <- "YES"
  Imputed$Phylo_info[is.na(Imputed$Phylo_info)] <- "NO"
  
  ## Add taxonomic information
  Imputed$Order <- TraitDF$Order
  Imputed$Family <- TraitDF$Family
  Imputed$Genus <- TraitDF$Genus
  
  ## Reprocess Primary diet and diet breadth (sqrt + normalise) if diet is included
  ## Then reorder columns.
  
  if(DietTRUE){
    
    Func <- function(X) {
      names(X) <- Diet
      ToPaste <- names(X)[which(X==1)]
      return(paste(ToPaste, collapse = "|"))
    }
    
    Imputed$Primary_diet <- apply(Imputed[,Diet], 1, Func)
    
    Imputed[, Diet] <- apply(Imputed[, Diet], 2, as.numeric)
    Imputed$sqrt_Diet_breadth_reprocessed <- apply(Imputed[, Diet], 1, sum, na.rm=T) %>% sqrt()

    Imputed <- Imputed[, c("Order", "Family", "Genus", "Best_guess_binomial",
                         "log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size",
                         "Diel_activity", "Specialisation", 
                         "Trophic_level", "sqrt_Diet_breadth", "sqrt_Diet_breadth_reprocessed","Primary_diet", Diet,
                         "Range_size_m2", "sqrt_Habitat_breadth_IUCN", Habitat, "Phylo_info")] 
  
  }
  
  else{ 
    
    Imputed$Primary_diet <- NA
    Imputed$sqrt_Diet_breadth <- NA
    
    Imputed <- Imputed[, c("Order", "Family", "Genus", "Best_guess_binomial",
                           "log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size",
                           "Diel_activity", "Specialisation","Trophic_level", "sqrt_Diet_breadth", "Primary_diet",
                           "Range_size_m2", "sqrt_Habitat_breadth_IUCN", Habitat, "Phylo_info")]
  }
  
  rownames(Imputed) <- c(1:nrow(Imputed))
  
  ## Return
  
  ToReturn <- list(Imputed.Dataset=Imputed, Imputation.errors=Errors)
  ToReturn <- list(ToReturn)
  
  return(ToReturn)
  
}


## Function to apply in parallel (runs the aboce function Imputations_missForest)
To_apply_parallel_imputations <- function (List_of_arguments) {

  ## NB: which function to use here on windows to replace pbmapply (or mapply)?
  ## mapply only works with forking methods (unix or linux)

  Imputations_results <- pbmapply (FUN=Imputations_missForest,
                          TraitDF=List_of_arguments[["TraitDF"]],
                          Taxinfo=List_of_arguments[["Taxinfo"]],
                          Traits_cont=List_of_arguments[["Traits_cont"]],
                          Traits_cat=List_of_arguments[["Traits_cat"]],
                          EV=List_of_arguments[["EV"]],
                          ErrorTrue=List_of_arguments[["ErrorTrue"]],
                          DietTRUE=List_of_arguments[["DietTRUE"]])

  return (Imputations_results)
}








