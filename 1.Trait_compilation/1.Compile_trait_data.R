## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                AIM: COMPILING TRAIT DATASETS BEFORE PERFORMING PHYLOGENETIC IMPUTATIONS                ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# TODO remove the "Meiri" dataset and use freely accessible data instead for reptiles

# Loading trait datasets
# Merging trait information from different datasets together (averaging on cont. traits when replicates)
# (Interpolating data to get a maximum of information through correlations among cont variables)
# Matching by species (/genus /family /order) name in Predicts using "Best_binomial_guess" as species name)
# Extracting trait information for matching species

# RETURNS: a complete trait dataset for each vertebrate class for which Binomial names are known
# with averages on lowest known taxonomic group when binomial name is unknown -> separate script


## Preamble ----------------------------------------------------------------
X <- c("data.table", "plyr", "dplyr", "tidyr", "magrittr", "reshape", "reshape2", "stringr", "stringi", "lazyeval", "rlang", "PerformanceAnalytics") 
invisible(lapply(X, library, character.only=TRUE)); rm(X)

source("Trait_data_compilation_functions.R")


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                     Load and prepare trait datasets                                    ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Load data ---------------------------------------------------------------------------------------------

Predicts <-  readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

# # # # # # # # # # # # # # # # # # #    A M N I O T E S     # # # # # # # # # # # # # # # # # # # 
Myhrvold <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Myhrvold.csv")

# # # # # # # # # # # # # # # # # # #    M A M M A L I A     # # # # # # # # # # # # # # # # # # # 
Pantheria <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Pantheria.csv") 
Pacifici <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Pacifici.csv") 
Kissling <- read.csv("../../Results/0.Processed_diet_datasets/MammalDIET_processed_diet.csv")
Elton_MD <- read.csv("../../Results/0.Processed_diet_datasets/Elton_mammals_processed.csv")  

# # # # # # # # # # # # # # # # # # #       B I R D S     # # # # # # # # # # # # # # # # # # # #  
Butchart_BM <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Butchart_BM.csv")
Butchart_GL <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Butchart_GL.csv")
# Sekercioglu_Diet <- read.csv("../../Results/0.Processed_diet_datasets/Sekercioglu_processed_diet.csv")
Elton_BD <- read.csv("../../Results/0.Processed_diet_datasets/Elton_birds_processed.csv")

# # # # # # # # # # # # # # # # # # #    R E P T I L E S     # # # # # # # # # # # # # # # # # # # 
Scharf <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Scharf.csv")

# This data has not been published. Using data from GARD (Global Assessments of Reptile Distributions)
# Meiri <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Meiri.csv") 
Vidan <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Vidan.csv")
Stark <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Stark.csv")
Schwarz <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Schwarz.csv")
Novosolov <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Novosolov.csv")
Novosolov_2 <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Novosolov_2.csv")
Slavenko <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Slavenko.csv") 
Meiri <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Meiri.csv")


# # # # # # # # # # # # # # # # # #    A M P H I B I A N S   # # # # # # # # # # # # # # # # # # # 
Amphibio <- read.csv("../../Results/0.Processed_diet_datasets/Amphibio_processed_diet.csv")
Cooper <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Cooper.csv")
Senior <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Senior.csv")
Bickford <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bickford.csv")

# # # # # # # # # # # # # # # # # # #      R A N G E S        # # # # # # # # # # # # # # # # # # # 
Mammal_range <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Mammal_range.csv")
Amphibian_range <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibian_range.csv")
Bird_range <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bird_range.csv")
Reptile_range <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Reptile_range.csv") 


# # # # # # # # # # # # # # # # # # #    I U C N    D A T A -- processed     # # # # # # # # # # # # 
Habitat_amphibian <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/Amphibians.csv")
Habitat_bird <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/Birds.csv")
Habitat_reptile <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/Reptiles.csv")
Habitat_mammal <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/Mammals.csv")


## Normalise all datasets ---------------------------------------------------------------------------------
X <- .Normalise_TDB(Kissling, Elton_MD, Pantheria, Pacifici, Myhrvold, 
                    Amphibio, Cooper, Senior, Bickford, 
                    Butchart_BM, Butchart_GL, Elton_BD,#Sekercioglu_Diet,
                    Scharf, Vidan, Stark, Schwarz, Novosolov, Novosolov_2, Slavenko, Meiri,
                    Mammal_range, Amphibian_range, Bird_range, TRUE)

Myhrvold <- X$Myhrvold
Pantheria <- X$Pantheria; Pacifici <- X$Pacifici; Kissling <- X$Kissling; Elton_MD <- X$Elton_MD
Amphibio <- X$Amphibio; Cooper <- X$Cooper; Senior <- X$Senior; Bickford <- X$Bickford
Butchart_BM <- X$Butchart_BM; Butchart_GL <- X$Butchart_GL; Elton_BD <- X$Elton_BD #Sekercioglu_Diet <- X$Sekercioglu_Diet
Scharf <- X$Scharf; Meiri <- X$Meiri; Vidan <- X$Vidan; Schwarz <- X$Schwarz; Stark <- X$Stark; Novosolov <- X$Novosolov; Novosolov_2 <- X$Novosolov_2; Slavenko <- X$Slavenko

Mammal_range <- X$Range_mammal; Bird_range <- X$Range_bird; Amphibian_range <- X$Range_amphibian



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                             Merge the trait datasets separatly for each class                          ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


## Merge Mammals -----------------------------------------------------------
TraitsMammal <- merge(Pacifici[, c("Best_guess_binomial", X$Mammal_traits.Pacifici)], 
                      Myhrvold[Myhrvold$class=="Mammalia", c("Best_guess_binomial", X$Traits.Myhrvold)], all=T)
TraitsMammal <- merge(TraitsMammal, Pantheria[, c("Best_guess_binomial", X$Mammal_traits.Pantheria)], all=T)
TraitsMammal <- merge(TraitsMammal, Mammal_range[, c("Best_guess_binomial", "Range_size_m2")], all=T)
TraitsMammal <- merge(TraitsMammal, Kissling[, c("Best_guess_binomial", "Trophic_level")], all=T)
TraitsMammal <- merge(TraitsMammal, Elton_MD[, c("Best_guess_binomial", X$Mammal_traits.Elton)], all=T)
TraitsMammal <- merge(TraitsMammal, Habitat_mammal, all=T)

## Merge Amphibians --------------------------------------------------------
TraitsAmphibian <- merge(Amphibio[, c("Best_guess_binomial", X$Amphibian_traits.Amphibio)],
                          Cooper[, c("Best_guess_binomial", X$Amphibian_traits.Cooper)], all=T)
TraitsAmphibian <- merge(TraitsAmphibian, Senior[, c("Best_guess_binomial", X$Amphibian_traits.Senior)], all=T)
TraitsAmphibian <- merge(TraitsAmphibian, Bickford[, c("Best_guess_binomial", X$Amphibian_traits.Bickford)], all=T)
TraitsAmphibian <- merge(TraitsAmphibian, Amphibian_range, all=T)
TraitsAmphibian <- merge(TraitsAmphibian, Habitat_amphibian, all=T)
# TraitsAmphibian <- merge(TraitsAmphibian, Amphibia_AnAge[, c("Best_guess_binomial", "Max_longevity_d", "Maturity_d")], all=T)

## Merge Reptiles ----------------------------------------------------------
TraitsReptile <- merge(Scharf[, c("Best_guess_binomial", X$Reptile_traits.Scharf)], 
                       Myhrvold[Myhrvold$class=="Reptilia", c("Best_guess_binomial", X$Traits.Myhrvold)], all=T) 
TraitsReptile <- merge(TraitsReptile, Reptile_range, all=T)
# TraitsReptile <- merge(TraitsReptile, Meiri[, c("Best_guess_binomial", "Body_mass_g")], all=T)
TraitsReptile <- merge(TraitsReptile, Meiri[, c("Best_guess_binomial", X$Traits.Meiri)], all=T)
TraitsReptile <- merge(TraitsReptile, Habitat_reptile, all=T)
TraitsReptile <- merge(TraitsReptile, Slavenko[,c("Best_guess_binomial", "Body_mass_g")], all=T)
TraitsReptile <- merge(TraitsReptile, Novosolov_2[,c("Best_guess_binomial", "Litter_size")], all=T)

TraitsReptile <- merge(TraitsReptile, Novosolov[,c("Best_guess_binomial", X$Traits.Novosolov)], all=T)
TraitsReptile <- merge(TraitsReptile, Stark[,c("Best_guess_binomial", X$Traits.Stark)], all=T)
TraitsReptile <- merge(TraitsReptile, Vidan[,c("Best_guess_binomial", "Diel_activity")], all=T)
TraitsReptile <- merge(TraitsReptile, Schwarz[,c("Best_guess_binomial", "Litter_size")], all=T)

# Merge Birds -------------------------------------------------------------
TraitsBird <- merge(Butchart_BM[, c("Best_guess_binomial", "Body_mass_g")], 
                    Myhrvold[Myhrvold$class=="Aves", c("Best_guess_binomial", X$Traits.Myhrvold)], all=T)
TraitsBird <- merge(TraitsBird, Butchart_GL[, c("Best_guess_binomial", "Longevity_d")], all=T)
TraitsBird <- merge(TraitsBird, Bird_range[, c("Best_guess_binomial", "Range_size_m2")], all=T)
TraitsBird <- merge(TraitsBird, Elton_BD[, c("Best_guess_binomial", X$Bird_traits.Elton)], all=T)
TraitsBird <- merge(TraitsBird, Habitat_bird, all=T)

rm(X)

cat(length(unique(TraitsMammal$Best_guess_binomial)), "unique species for mammals")
cat(length(unique(TraitsAmphibian$Best_guess_binomial)), "unique species for amphibians")
cat(length(unique(TraitsBird$Best_guess_binomial)), "unique species for birds")
cat(length(unique(TraitsReptile$Best_guess_binomial)), "unique species for reptiles")



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                              Reduce redundancy: average on continuous traits                           ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# To exclude non continuous traits or select categorical traits
Diet <- c("Trophic_level", "Primary_diet","IN", "VE", "SCV", "FR", "NE", "SE", "PL")
Habitat_IUCN_traits <- colnames(Habitat_amphibian)[3:16]


# # # # # # # # # # # # # # # # # # #    R E P T I L E S     # # # # # # # # # # # # # # # # # # # 
Continuous.Reptile <- TraitsReptile %>% select(-c(Trophic_level, Habitat_IUCN_traits, Binomial_name, Diel_activity))
Continuous.Reptile <- setDT(Continuous.Reptile)[, lapply(.SD, mean, na.rm=TRUE), by = Best_guess_binomial]  %>% as.data.frame()
Continuous.Reptile %<>% filter(Best_guess_binomial!="")

# # # # # # # # # # # # # # # # # # #      B I R D S      # # # # # # # # # # # # # # # # # # # # # 
Continuous.Bird <- TraitsBird %>% select(-Diet, -Habitat_IUCN_traits, -Diel_activity)
Continuous.Bird <- setDT(Continuous.Bird)[, lapply(.SD, mean, na.rm=TRUE), by = Best_guess_binomial]  %>% as.data.frame()


# # # # # # # # # # # # # # # # # # #      M A M M A L S      # # # # # # # # # # # # # # # # # # # 
Continuous.Mammal <- TraitsMammal %>% select(-c(Diet, Habitat_breadth, Terrestriality, Habitat_IUCN_traits, Diel_activity))
Continuous.Mammal <- setDT(Continuous.Mammal)[, lapply(.SD, mean, na.rm=TRUE), by = Best_guess_binomial]


# # # # # # # # # # # # # # # # # # #    A M P H I B I A N S    # # # # # # # # # # # # # # # # # # #
Continuous.Amphibians <- TraitsAmphibian %>% select(-c(Fos, Ter, Aqu, Arb, Habitat_breadth, Terrestriality,
                                                       Diet, Primary_diet, Habitat_IUCN_traits, Species,Diel_activity))
Continuous.Amphibians <- setDT(Continuous.Amphibians)[, lapply(.SD, mean, na.rm=TRUE), by = Best_guess_binomial]



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                              Reduce redundancy: for selected categorical traits                        ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##



# # # # # # # # # # # # # # # # # # #    R E P T I L E S     # # # # # # # # # # # # # # # # # # # 
TraitsReptile$Trophic_level %<>% as.character()
Categorical.Reptile <- TraitsReptile %>% select(Best_guess_binomial, Trophic_level, Habitat_IUCN_traits, Diel_activity)
Categorical.Reptile %<>% .MutHabIUCN()
Categorical.Reptile %<>% mutate(Trophic_level=
           ifelse(length(unique(Trophic_level))==1, unique(Trophic_level), unique(Trophic_level[!is.na(Trophic_level)])))
Categorical.Reptile %<>% mutate(Diel_activity=
                                  ifelse(length(unique(Diel_activity))==1, unique(Diel_activity), unique(Diel_activity[!is.na(Diel_activity)])))

Categorical.Reptile %<>% distinct() %>% as.data.frame()
Categorical.Reptile %<>% filter(Best_guess_binomial!="")

# Reptile lacking Diet


# # # # # # # # # # # # # # # # # # #      B I R D S      # # # # # # # # # # # # # # # # # # # # # 
Categorical.Bird <- TraitsBird %>% select(Best_guess_binomial, Diet, Habitat_IUCN_traits, Diel_activity)
Categorical.Bird %<>% .MutHabIUCN()
Categorical.Bird %<>% .MutDiet()
Categorical.Bird %<>% mutate(SCV= 
                                 ifelse(length(unique(SCV))==1, unique(SCV), unique(SCV[!is.na(SCV)])))
Categorical.Bird %<>% mutate(Diel_activity=
                                  ifelse(length(unique(Diel_activity))==1, unique(Diel_activity), unique(Diel_activity[!is.na(Diel_activity)])))

Categorical.Bird %<>% distinct() %>% as.data.frame()


# # # # # # # # # # # # # # # # # # #     M A M M A L S      # # # # # # # # # # # # # # # # # # # 
Categorical.Mammal <- TraitsMammal %>% select(Best_guess_binomial, Diet, Habitat_IUCN_traits, Habitat_breadth, Terrestriality, Diel_activity)
Categorical.Mammal %<>% .MutHabIUCN()
Categorical.Mammal %<>% .MutDiet()
Categorical.Mammal %<>% mutate(SCV= 
                                 ifelse(length(unique(SCV))==1, unique(SCV), unique(SCV[!is.na(SCV)])))
Categorical.Mammal %<>% mutate(Habitat_breadth=
                                 ifelse(length(unique(Habitat_breadth))==1, unique(Habitat_breadth), unique(Habitat_breadth[!is.na(Habitat_breadth)]))) %>%
                        mutate(Terrestriality=
                                 ifelse(length(unique(Terrestriality))==1, unique(Terrestriality), unique(Terrestriality[!is.na(Terrestriality)]))) %>%
                        mutate(Diel_activity=
                                ifelse(length(unique(Diel_activity))==1, unique(Diel_activity), unique(Diel_activity[!is.na(Diel_activity)])))
  

Categorical.Mammal %<>% distinct() %>% as.data.frame()



# # # # # # # # # # # # # # # # # # #    A M P H I B I A N S    # # # # # # # # # # # # # # # # # # #
Categorical.Amphibian <- TraitsAmphibian %>% select(Best_guess_binomial, Diet, Habitat_IUCN_traits,
                                                    Fos, Ter, Aqu, Arb, 
                                                    Terrestriality, Habitat_breadth, Diel_activity)

Categorical.Amphibian %<>% .MutHabIUCN()
Categorical.Amphibian %<>% .MutDiet()
Categorical.Amphibian %<>% 
  mutate(Fos= ifelse(length(unique(Fos))==1, unique(Fos), unique(Fos[!is.na(Fos)]))) %>%
  mutate(Ter= ifelse(length(unique(Ter))==1, unique(Ter), unique(Ter[!is.na(Ter)]))) %>%
  mutate(Aqu= ifelse(length(unique(Aqu))==1, unique(Aqu), unique(Aqu[!is.na(Aqu)]))) %>%
  mutate(Arb= ifelse(length(unique(Arb))==1, unique(Arb), unique(Arb[!is.na(Arb)]))) %>%
  mutate(Habitat_breadth=
           ifelse(length(unique(Habitat_breadth))==1, unique(Habitat_breadth), unique(Habitat_breadth[!is.na(Habitat_breadth)]))) %>%
  mutate(Terrestriality=
           ifelse(length(unique(Terrestriality))==1, unique(Terrestriality), unique(Terrestriality[!is.na(Terrestriality)]))) %>%
  mutate(Diel_activity=
           ifelse(length(unique(Diel_activity))==1, unique(Diel_activity), unique(Diel_activity[!is.na(Diel_activity)])))
  
Categorical.Amphibian %<>% distinct() %>% as.data.frame()


  


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                              Merge continuous and categorical                                          ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
TraitsReptile <- merge(Continuous.Reptile, Categorical.Reptile)
TraitsMammal <- merge(Continuous.Mammal, Categorical.Mammal)
TraitsAmphibian <- merge(Continuous.Amphibians, Categorical.Amphibian)
TraitsBird <- merge(Continuous.Bird, Categorical.Bird)


rm(Amphibian_range, Amphibio, Bickford, Bird_range, Butchart_BM, Butchart_GL,
   Categorical.Amphibian, Categorical.Mammal, Categorical.Reptile, Categorical.Bird,
   Continuous.Amphibians, Continuous.Mammal, Continuous.Reptile, Continuous.Bird,
   Kissling, Mammal_range, Meiri, Myhrvold, Pacifici, Pantheria, Reptile_range, Scharf,
   Sekercioglu_Diet, Senior, Cooper, Habitat_amphibian, Habitat_IUCN_traits, Habitat_mammal, Habitat_reptile, Habitat_bird,
   Diet)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                                     Manually adding some information                                   ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
TraitsMammal$Trophic_level[TraitsMammal$Best_guess_binomial=="Bos frontalis"] <- "Herbivore"
TraitsMammal$Trophic_level[TraitsMammal$Best_guess_binomial=="Bos taurus"] <- "Herbivore"
TraitsMammal$Trophic_level[TraitsMammal$Best_guess_binomial=="Capra hircus"] <- "Herbivore"
TraitsMammal$Trophic_level[TraitsMammal$Best_guess_binomial=="Equus caballus"] <- "Herbivore"
TraitsMammal$Trophic_level[TraitsMammal$Best_guess_binomial=="Felis catus"] <- "Carnivore"
TraitsMammal$Trophic_level[TraitsMammal$Best_guess_binomial=="Taurotragus oryx"] <- "Herbivore"


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                              Intersection with PREDICTS species                                        ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## What are the species that do not match?
# DiffMammals <- NoMatch(Predicts, TraitsMammal, "Mammalia")
# DiffAmphibians <- NoMatch(Predicts, TraitsAmphibian, "Amphibia")
# DiffReptiles <- NoMatch(Predicts, TraitsReptile, "Reptilia")
# DiffBirds <- NoMatch(Predicts, TraitsBird, "Aves")


## Adding the species that do not intersect to the Trait dataset for future phylogenetic imputation
AddM <- .Add_species(Predicts, TraitsMammal, "Mammalia")
AddR <- .Add_species(Predicts, TraitsReptile, "Reptilia")
AddA <- .Add_species(Predicts, TraitsAmphibian, "Amphibia")
AddB <- .Add_species(Predicts, TraitsBird, "Aves")

TraitsMammal <- AddM$TraitDB
TraitsAmphibian <- AddA$TraitDB
TraitsReptile <- AddR$TraitDB
TraitsBird <- AddB$TraitDB

## Match against PREDICTS species name and retrieve trait values before interpolation
TraitsMammal.Predicts <- .Match_Species("Mammalia", Predicts, TraitsMammal)
TraitsAmphibian.Predicts <- .Match_Species("Amphibia", Predicts, TraitsAmphibian)
TraitsReptile.Predicts <- .Match_Species("Reptilia", Predicts, TraitsReptile)
TraitsBird.Predicts <- .Match_Species("Aves", Predicts, TraitsBird)


cat(length(unique(TraitsMammal$Best_guess_binomial)), "unique species for mammals")
cat(length(unique(TraitsAmphibian$Best_guess_binomial)), "unique species for amphibians")
cat(length(unique(TraitsBird$Best_guess_binomial)), "unique species for birds")
cat(length(unique(TraitsReptile$Best_guess_binomial)), "unique species for reptiles")

cat(length(unique(TraitsMammal.Predicts$Best_guess_binomial)), "unique species for PREDICTS mammals")
cat(length(unique(TraitsAmphibian.Predicts$Best_guess_binomial)), "unique species for PREDICTS amphibians")
cat(length(unique(TraitsBird.Predicts$Best_guess_binomial)), "unique species for PREDICTS birds")
cat(length(unique(TraitsReptile.Predicts$Best_guess_binomial)), "unique species for PREDICTS reptiles")


## Saves files: basic trait compilation
write.csv(TraitsMammal, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Mammals.csv", row.names=F)
write.csv(TraitsAmphibian, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Amphibians.csv", row.names=F)
write.csv(TraitsBird, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Birds.csv", row.names=F)
write.csv(TraitsReptile, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/1.compiled/Reptiles.csv", row.names=F)

write.csv(TraitsMammal.Predicts, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/Predicts_subsets/PredictsMammals.csv", row.names=F)
write.csv(TraitsAmphibian.Predicts, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/Predicts_subsets/PredictsAmphibians.csv", row.names=F)
write.csv(TraitsBird.Predicts, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/Predicts_subsets/PredictsBirds.csv", row.names=F)
write.csv(TraitsReptile.Predicts, "../../Results/1.Traits_before_imputations/With_taxonomic_correction/Predicts_subsets/PredictsReptiles.csv", row.names=F)



