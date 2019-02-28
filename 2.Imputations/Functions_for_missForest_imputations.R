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

Imputations_missForest <- function (TraitDF, Taxinfo, Traits_cont, Traits_cat, EV, ErrorTrue, DietTRUE, std) {

  browser()
  
  ## Select traits of interest, to impute, and phylogenetic eigenvectors, and taxinfo 
  To_impute <- TraitDF[, colnames(TraitDF) %in% c(Taxinfo, Traits_cont, Traits_cat, EV)]
  rownames(To_impute) <- TraitDF$Best_guess_binomial
  
  ## Set habitat and diet variables as binary (TRUE for 1 and FALSE for 0)
  
  ## Set class as numeric for all continuous traits
  To_impute[, Traits_cont] <- apply(To_impute[, Traits_cont], 2, as.numeric)
  
  ## Set habitat and diet variables as factors. Otherwise, when they are numeric or logical, outputs are numeric
  
  if(DietTRUE){
    
    Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")
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
  
  ## Select traits and variables of interest after imputations
  
  if(DietTRUE) {
    
    Continuous <- unique(c(colnames(Imputed[grepl("Body_mass_g", colnames(Imputed))]),
                   colnames(Imputed[grepl("Longevity_d", colnames(Imputed))]),
                   colnames(Imputed[grepl("Litter_size", colnames(Imputed))]),
                   colnames(Imputed[grepl("Range_size_m2", colnames(Imputed))]),
                   colnames(Imputed[grepl("Habitat_breadth_IUCN", colnames(Imputed))]),
                   colnames(Imputed[grepl("Diet_breadth", colnames(Imputed))])))
                            
    Categorical <- c(Habitat, "Specialisation", "Diel_activity","Trophic_level", "Primary_diet", Diet)
      
    }
  
  else{
    
    Continuous <- unique(c(colnames(Imputed[grepl("Body_mass_g", colnames(Imputed))]),
                   colnames(Imputed[grepl("Longevity_d", colnames(Imputed))]),
                   colnames(Imputed[grepl("Litter_size", colnames(Imputed))]),
                   colnames(Imputed[grepl("Range_size_m2", colnames(Imputed))]),
                   colnames(Imputed[grepl("Habitat_breadth_IUCN", colnames(Imputed))])))
    Categorical <- c(Habitat, "Specialisation", "Diel_activity","Trophic_level")
    
    }
  
  Imputed$Best_guess_binomial <- rownames(Imputed)
  Imputed <- Imputed[order(Imputed$Best_guess_binomial), c(Taxinfo, "Best_guess_binomial", Continuous, Categorical, "EV_1")]
  
  ## Add a column for phylogenetic information (yes or no)
  Imputed$Phylo_info <- TraitDF$EV_1
  Imputed$Phylo_info[!is.na(Imputed$Phylo_info)] <- "YES"
  Imputed$Phylo_info[is.na(Imputed$Phylo_info)] <- "NO"
  
  ## Add taxonomic information
  Imputed$Class <- TraitDF$Class
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
    
    # Reprocess primary diet, for comparison with imputed values
    Imputed$Primary_diet_reprocessed <- apply(Imputed[,Diet], 1, Func)
    # Imputed$Primary_diet_reprocessed[Imputed$Primary_diet_reprocessed==""] <- "OM"
    
    # Reprocess diet breadth, for comparison with imputed values
    Imputed[, Diet] <- apply(Imputed[, Diet], 2, as.numeric)
    Imputed$Diet_breadth_reprocessed <- apply(Imputed[, Diet], 1, sum, na.rm=T)
    
    if(std) {
      Imputed$Diet_breadth_reprocessed <- sqrt(Imputed$Diet_breadth_reprocessed)
      Imputed$Diet_breadth_reprocessed <- scale(Imputed$Diet_breadth_reprocessed, center = TRUE, scale = TRUE)
      colnames(Imputed)[colnames(Imputed)=="Diet_breadth_reprocessed"] <- "sqrt_Diet_breadth_reprocessed"
    }
      
    ## Rearrange columns
    Imputed <- Imputed[, unique(c("Class","Order", "Family", "Genus", "Best_guess_binomial",
                           colnames(Imputed[grepl("Body_mass_g", colnames(Imputed))]),
                           colnames(Imputed[grepl("Longevity_d", colnames(Imputed))]),
                           colnames(Imputed[grepl("Litter_size", colnames(Imputed))]),
                           colnames(Imputed[grepl("Range_size_m2", colnames(Imputed))]),
                           "Diel_activity", 
                           "Trophic_level", 
                           colnames(Imputed[grepl("Diet_breadth", colnames(Imputed))]),
                           colnames(Imputed[grepl("Diet_breadth_reprocessed", colnames(Imputed))]),
                           "Primary_diet", 
                           colnames(Imputed[grepl("Primary_diet_reprocessed", colnames(Imputed))]),
                           Diet,
                           "Specialisation", 
                           colnames(Imputed[grepl("Habitat_breadth_IUCN", colnames(Imputed))]),
                           Habitat,
                           "Phylo_info"))] 
  
  }
  
  else { 
    Diet <- c("IN", "VE", "PL", "SE", "NE", "FR")
    Imputed$Primary_diet <- NA
    Imputed$Primary_diet_reprocessed <- NA
    Imputed$IN <- NA
    Imputed$VE <- NA
    Imputed$SE <- NA
    Imputed$FR <- NA
    Imputed$NE <- NA
    Imputed$PL <- NA
    
    
    if (std) {
      Imputed$sqrt_Diet_breadth <- NA
      Imputed$sqrt_Diet_breadth_reprocessed <- NA
    } 
    
    else{
      Imputed$Diet_breadth <- NA
      Imputed$Diet_breadth_reprocessed <- NA
      }
 
    
    ## Rearrange columns
    Imputed <- Imputed[, unique(c("Class","Order", "Family", "Genus", "Best_guess_binomial",
                           colnames(Imputed[grepl("Body_mass_g", colnames(Imputed))]),
                           colnames(Imputed[grepl("Longevity_d", colnames(Imputed))]),
                           colnames(Imputed[grepl("Litter_size", colnames(Imputed))]),
                           colnames(Imputed[grepl("Range_size_m2", colnames(Imputed))]),
                           "Diel_activity", 
                           "Trophic_level", 
                           colnames(Imputed[grepl("Diet_breadth", colnames(Imputed))]),
                           colnames(Imputed[grepl("Diet_breadth_reprocessed", colnames(Imputed))]),
                           "Primary_diet", 
                           colnames(Imputed[grepl("Primary_diet_reprocessed", colnames(Imputed))]),
                           Diet,
                           "Specialisation", 
                           colnames(Imputed[grepl("Habitat_breadth_IUCN", colnames(Imputed))]),
                           Habitat,
                           "Phylo_info"))] 
  }
  
  rownames(Imputed) <- c(1:nrow(Imputed))
  
  ## Return
  
  ToReturn <- list(Imputed.Dataset=Imputed, Imputation.errors=Errors)
  ToReturn <- list(ToReturn)
  
  return(ToReturn)
  
}


## Function to apply in parallel (runs the above function Imputations_missForest)
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
                          DietTRUE=List_of_arguments[["DietTRUE"]],
                          std=List_of_arguments[["std"]])

  return (Imputations_results)
}








