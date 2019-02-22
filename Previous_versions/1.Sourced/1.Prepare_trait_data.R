## Functions to prepare the trait datasets and the data extraction:
## (necessary runs of these functions, sourced to the main script, are at the end atm)



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Normalise trait datasets -------------------------------------------------------------------------

## This function changes column names, and format trait datasets
## so that names are unique and transferable across datasets

.Normalise_TDB.Mammals <- function(Kissling, Pantheria, Pacifici, Myhrvold, Range) {
  
  # Kissling: function to format Family and Order names
  .Cap <- function(s) {
    paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2)),
          sep = "")}
  
  
  Kissling$Binomial_name <- paste(Kissling$Genus, Kissling$Species, sep= " ")
  Mammal_traits.Kissling <- c("Trophic_level")
  colnames(Kissling)[24] <- Mammal_traits.Kissling
  Kissling$Order <- lapply(Kissling$Order, .Cap)
  Kissling$Family <- lapply(Kissling$Family, .Cap)
  
  # Pantheria
  Pantheria[Pantheria==-999] <- NA
  Mammal_traits.Pantheria <- c("Body_mass_g", "Diet_breadth", "Habitat_breadth",
                               "Home_range_group_km2", "Home_range_ind_km2", "Litter_size", "Generation_length_d",
                               "Terrestriality", "Trophic_level")
  colnames(Pantheria)[c(1, 2, 3, 5, 7, 14, 17:19, 21, 23, 31, 32)] <- c("Order", "Family", "Genus", "Binomial_name", Mammal_traits.Pantheria)
  Pantheria$Trophic_level[Pantheria$Trophic_level==1] <- "Herbivore"
  Pantheria$Trophic_level[Pantheria$Trophic_level==2] <- "Omnivore"
  Pantheria$Trophic_level[Pantheria$Trophic_level==3] <- "Carnivore"
  Pantheria$Generation_length_d <- Pantheria$Generation_length_d*30.44
  
  # Pacifici
  Mammal_traits.Pacifici <- c("Body_mass_g", "Generation_length_d")
  colnames(Pacifici)[c(5, 6, 14)] <- c("Binomial_name", Mammal_traits.Pacifici) 
  Pacifici[Pacifici=="no information"] <- NA
  Pacifici$Max_longevity_d  %<>% as.character %>% as.numeric()  

  # Myhrvold
  Myhrvold[Myhrvold==-999] <- NA
  Mammal_traits.Myhrvold <- c("Litter_size", "Body_mass_g", "Generation_length_d")
  Myhrvold$Binomial_name <- paste(Myhrvold$genus, Myhrvold$species, sep= " ")
  colnames(Myhrvold)[c(9, 11, 20)] <- Mammal_traits.Myhrvold
  colnames(Myhrvold)[c(2:4)] <- c("Order", "Family", "Genus")
  Myhrvold$Generation_length_d <- Myhrvold$Generation_length_d * 365.25
  Myhrvold$maximum_longevity_y <- Myhrvold$maximum_longevity_y * 365.25
  
  # Format the Range size dataset and add Order/Family/Genus information to it
  colnames(Range) <- c("Binomial_name", "Range_size_m2")
  Range$Genus <- word(Range$Binomial_name, 1)
  # With information from Pantheria and from Kissling, all genuses covered except one to add manually
  Range$Family <- Pantheria$Family[match(Range$Genus, Pantheria$Genus)] %>% as.character()
  Range$Order <- Pantheria$Order[match(Range$Genus, Pantheria$Genus)] %>% as.character()
  # Get remaining genuses from Kissling
  Index <- as.numeric(row.names(Range[is.na(Range$Family),]))
  Y.match <- match(Range$Genus[Index], Kissling$Genus)
  Range$Family[Index] <- Kissling$Family[Y.match] %>% as.character()
  Range$Order[Index] <- Kissling$Order[Y.match] %>% as.character()
  # Add the last missing information manually
  Range$Family[Range$Genus=="Pennatomys"] <- "Cricetidae" 
  Range$Order[Range$Genus=="Pennatomys"] <- "Rodentia" 
  Range$Family <- as.factor(Range$Family)
 
  
  # OUTPUTS
  return(list(Kissling=Kissling, 
              Pantheria=Pantheria, 
              Pacifici=Pacifici, 
              Myhrvold=Myhrvold,
              Range=Range,
              Kissling=Kissling,
              Mammal_traits.Pantheria=Mammal_traits.Pantheria,
              Mammal_traits.Pacifici=Mammal_traits.Pacifici,
              Mammal_traits.Myhrvold=Mammal_traits.Myhrvold,
              Mammal_traits.Kissling=Mammal_traits.Kissling))

  }





#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Match by binomial species names ------------------------------------------------------------------

## This function matches species in Predicts against species 
## in the trait datasets by binomial names and extracts the relevant trait
## information
## ARG: Veterbrate class under consideration, Predicts DB, Trait DB, Trait names
## RETURNS: a dataset with all species known by binomial name in Predicts for the 
## considered class and the relevant trait information extracted from the trait DB

.Match_Species <- function (VClass, Predicts, TraitDB, Traits) {
  
  # Select species in Predicts, that are known by binomial name
  Sp <- Predicts %>% filter(Class==VClass) %>%  
    distinct(Order, Family, Genus, Best_guess_binomial) %>% 
    filter(Best_guess_binomial!="") %>%
    droplevels()
  
  # Match by species binomial name
  Y <- intersect(Sp$Best_guess_binomial, TraitDB$Binomial_name)
  TraitVal <- subset(TraitDB, Binomial_name %in% Y) %>%
    select(Binomial_name, Traits)
  
  # Merge the selected traits with Sp
  Sp <- merge(Sp, TraitVal, by.x = "Best_guess_binomial", by.y = "Binomial_name", all=TRUE)
  
  cat("Matching by best binomial guess on provides", length(Y), "matches for", VClass, 
      "on a total of", nrow(Sp)," binomial guesses for all traits\n")
  
  return (Sp)
} 



#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
# Get taxa not known by binomial names in the PREDICTS DB and add to the Traits dataset ------------

.GetOther <-function (Predicts, VClass, Traits) {
  
  Vert <- Predicts %>% filter(Class == VClass )
  
  # Get the other species, that are unknown by binomial name, extract Genus, Family and Order
  # Add them to the Traits dataset
  
  # Species known by genus only
  Other_genus <- Vert %>%
    filter(Best_guess_binomial=="") %>%
    filter(Genus != "") %>% distinct(Order, Family, Genus) %>% droplevels()
  
  # Species known by family only
  Other_family <- Vert %>% 
    filter(Best_guess_binomial=="") %>%
    filter(Genus == "") %>%
    filter(Family != "") %>%
    distinct(Order, Family, Genus) %>% droplevels()
  
  # Species known by Order only
  Other_order <- Vert %>% 
    filter(Best_guess_binomial=="") %>%
    filter(Genus == "") %>%
    filter(Family == "") %>%
    distinct(Order, Family, Genus) %>% droplevels()
  
  Other <- rbind(Other_genus, Other_family, Other_order) 
  rm(Other_family, Other_genus, Other_order)
  
  # Create same columns as in Traits and fil with NAs
  Df <- data.frame(matrix(nrow = nrow(Other), ncol = length(4:ncol(Traits))))
  colnames(Df) <- c("Best_guess_binomial", colnames(Traits)[5:ncol(Traits)])
  Df$Best_guess_binomial <- ""
  
  Other <- cbind(Other, Df)
  
  # Put the information together
  Other <- rbind(Traits, Other)
  
  return(Other)
  
}


