## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                AIM: COMPILING TRAIT DATASETS BEFORE PERFORMING PHYLOGENETIC IMPUTATIONS                ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

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

Predicts <-  readRDS("../../Data/PREDICTS_database.rds")
Predicts <- subset(Predicts, Class %in% c("Aves", "Amphibia", "Mammalia", "Reptilia"))

# # # # # # # # # # # # # # # # # # #    A M N I O T E S     # # # # # # # # # # # # # # # # # # # 
Myhrvold <- read.csv("../../Data/Amniotes_Myhrvold_2015/Amniote_Database_Aug_2015.csv")
Myhrvold$Best_guess_binomial <- paste(Myhrvold$genus, Myhrvold$species, sep=" ")

# # # # # # # # # # # # # # # # # # #    M A M M A L I A     # # # # # # # # # # # # # # # # # # # 
Pantheria <- read.csv("../../Data/Mammals/PanTHERIA/Pantheria_1_0_WR05_Aug2008.csv")
colnames(Pantheria)[5] <- "Best_guess_binomial"

Pacifici <- read.csv("../../Data/Mammals/PacificiMammals.csv")
colnames(Pacifici)[5] <- "Best_guess_binomial"

Kissling <- read.csv("../../Data/Mammals/MammalDIET_processed_diet.csv")
colnames(Kissling)[1] <- "Best_guess_binomial"

Elton_MD <- read.csv("../../Data/Mammals/Elton_mammals_processed.csv")
colnames(Elton_MD)[2] <- "Best_guess_binomial"

# # # # # # # # # # # # # # # # # # #       B I R D S     # # # # # # # # # # # # # # # # # # # #  
Butchart_BM <- read.csv("../../Data/Birds/Butchart_BM.csv")
Butchart_GL <- read.csv("../../Data/Birds/ButchartGenerationLength.csv")
colnames(Butchart_BM)[1] <- "Best_guess_binomial"
colnames(Butchart_GL)[1] <- "Best_guess_binomial"
Elton_BD <- read.csv("../../Data/Birds/Elton_birds_processed.csv")
colnames(Elton_BD)[8] <- "Best_guess_binomial"

# Sekercioglu_Diet <- read.csv("../../Data/Birds/Sekercioglu_processed_diet.csv")
# colnames(Sekercioglu_Diet)[2]<- "Best_guess_binomial"


# # # # # # # # # # # # # # # # # # #    R E P T I L E S     # # # # # # # # # # # # # # # # # # # 
Scharf <- read.csv("../../Data/Reptiles/Scharf.csv")
colnames(Scharf)[4] <- "Best_guess_binomial"
Vidan <- read.csv("../../Data/Reptiles/Vidan2017_Dielactivity.csv")
colnames(Vidan)[1] <- "Best_guess_binomial"
Stark <- read.csv("../../Data/Reptiles/Stark2018_GEB_longevity.csv")
colnames(Stark)[1] <- "Best_guess_binomial"
Schwarz <- read.csv("../../Data/Reptiles/Schwarz_Meiri_GEB_2017.csv")
colnames(Schwarz)[1] <- "Best_guess_binomial"
Novosolov <- read.csv("../../Data/Reptiles/Novosolov_2017_GEB.csv")%>%
  filter(Taxonomic.group!="Birds")%>%
  filter(Taxonomic.group!="Mammals")
colnames(Novosolov)[1] <- "Best_guess_binomial"
Novosolov_2 <- read.csv("../../Data/Reptiles/Novosolov_GEB_2013.csv")
colnames(Novosolov_2)[1] <- "Best_guess_binomial"
Slavenko <- read.csv("../../Data/Reptiles/Body_sizes_of_all_extant_reptiles_Slavenko_2016_GEB.csv") 
colnames(Slavenko)[2] <- "Best_guess_binomial"
Meiri <- read.csv("../../Data/Reptiles/Meiri_2015_Evolutionary_Biology.csv")
colnames(Meiri)[2] <- "Best_guess_binomial"

# # # # # # # # # # # # # # # # # #    A M P H I B I A N S   # # # # # # # # # # # # # # # # # # # 
Amphibio <- read.csv("../../Data/Amphibians/Amphibio_processed_diet.csv")
Amphibio$Species <- Amphibio$Best_guess_binomial

Cooper <- read.csv("../../Data/Amphibians/Cooper2008.csv")
Cooper <- subset(Cooper, Binomial!="")
colnames(Cooper)[3] <- "Best_guess_binomial"
Senior <- read.csv("../../Data/Amphibians/Senior_svl_data.csv")
Senior <- subset(Senior, Rank=="Species")
colnames(Senior)[1] <- "Best_guess_binomial"
Bickford <- read.csv("../../Data/Amphibians/Bickford.csv")
Bickford$Best_guess_binomial <-  paste(Bickford$Genus, Bickford$Species, sep=" ")

# # # # # # # # # # # # # # # # # # #      R A N G E S        # # # # # # # # # # # # # # # # # # # 
Mammal_range <- read.csv("../../Data/Range_sizes/mammal_range_areas.csv")
Mammal_range$Best_guess_binomial <- Mammal_range$Species

Amphibian_range <- read.csv("../../Data/Range_sizes/amphibian_range_areas.csv")
Amphibian_range$Best_guess_binomial <- Amphibian_range$Species

Bird_range <- read.csv("../../Data/Range_sizes/bird_range_areas.csv")
Bird_range$Best_guess_binomial <- Bird_range$Species

Reptile_range <- read.csv("../../Data/Range_sizes/reptile_range_areas.csv")
Reptile_range <- subset(Reptile_range, Binomial_name!="Chelonia mydas Hawaiian subpopulation")
colnames(Reptile_range)[1] <- "Best_guess_binomial"
Reptile_range$Species <- Reptile_range$Best_guess_binomial


# # # # # # # # # # # # # # # # # # #    I U C N    D A T A -- processed     # # # # # # # # # # # # 
Habitat_mammal <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Mammals.csv")
Habitat_amphibian <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Amphibians.csv")
Habitat_bird <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Birds.csv")
Habitat_reptile <- read.csv("../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Reptiles.csv")



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
Diet <- c("Trophic_level", "Primary_diet","IN", "VE", "FR", "NE", "SE", "PL")
Habitat_IUCN_traits <- colnames(Habitat_amphibian)[2:16]


# # # # # # # # # # # # # # # # # # #    R E P T I L E S     # # # # # # # # # # # # # # # # # # # 
Continuous.Reptile <- TraitsReptile %>% select(-c(Trophic_level, Habitat_IUCN_traits, Species, Diel_activity))
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
Categorical.Reptile %<>% mutate(Habitat_breadth_IUCN=
                                  ifelse(length(unique(Habitat_breadth_IUCN))==1, unique(Habitat_breadth_IUCN), unique(Habitat_breadth_IUCN[!is.na(Habitat_breadth_IUCN)])))
Categorical.Reptile %<>% distinct() %>% as.data.frame()
Categorical.Reptile %<>% filter(Best_guess_binomial!="")

# Reptile lacking Diet


# # # # # # # # # # # # # # # # # # #      B I R D S      # # # # # # # # # # # # # # # # # # # # # 
Categorical.Bird <- TraitsBird %>% select(Best_guess_binomial, Diet, Habitat_IUCN_traits, Diel_activity)
Categorical.Bird %<>% .MutHabIUCN()
Categorical.Bird %<>% .MutDiet()
Categorical.Bird %<>% mutate(Diel_activity=
                               ifelse(length(unique(Diel_activity))==1, unique(Diel_activity), unique(Diel_activity[!is.na(Diel_activity)])))

Categorical.Bird %<>% distinct() %>% as.data.frame()


# # # # # # # # # # # # # # # # # # #     M A M M A L S      # # # # # # # # # # # # # # # # # # # 
Categorical.Mammal <- TraitsMammal %>% select(Best_guess_binomial, Diet, Habitat_IUCN_traits, Habitat_breadth, Terrestriality, Diel_activity)
Categorical.Mammal %<>% .MutHabIUCN()
Categorical.Mammal %<>% .MutDiet()
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
   Senior, Cooper, Habitat_amphibian, Habitat_IUCN_traits, Habitat_mammal, Habitat_reptile, Habitat_bird,
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
TraitsMammal <- .Add_species(Predicts, TraitsMammal, "Mammalia")
TraitsReptile <- .Add_species(Predicts, TraitsReptile, "Reptilia")
TraitsAmphibian <- .Add_species(Predicts, TraitsAmphibian, "Amphibia")
TraitsBird <- .Add_species(Predicts, TraitsBird, "Aves")

TraitsMammal <- TraitsMammal$TraitDB
TraitsAmphibian <- TraitsAmphibian$TraitDB
TraitsBird <- TraitsBird$TraitDB
TraitsReptile <- TraitsReptile$TraitDB


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
write.csv(TraitsMammal, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Mammals.csv", row.names=F)
write.csv(TraitsAmphibian, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Amphibians.csv", row.names=F)
write.csv(TraitsBird, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Birds.csv", row.names=F)
write.csv(TraitsReptile, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/1.compiled/Reptiles.csv", row.names=F)

write.csv(TraitsMammal.Predicts, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/Predicts_subsets/PredictsMammals.csv", row.names=F)
write.csv(TraitsAmphibian.Predicts, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/Predicts_subsets/PredictsAmphibians.csv", row.names=F)
write.csv(TraitsBird.Predicts, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/Predicts_subsets/PredictsBirds.csv", row.names=F)
write.csv(TraitsReptile.Predicts, "../../Results/1.Traits_before_imputations/Without_taxonomic_correction/Predicts_subsets/PredictsReptiles.csv", row.names=F)





