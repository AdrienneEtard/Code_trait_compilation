## Standardising taxonomy: looking for synonyms, resolving typos and species entered by common name instead of scientific

# In this script:
# (1) Species lists are collated by merging species names from the trait datasets, the Predicts db and the phylogenies
# (2) A function is applied to each species "original name" to correct potential typos and mispellings
# (3) A function is applied to each corrected species name to find wether the name is accepted or synonymic, according to the Red List. If accepted, 
#      the function retrieves synonyms from the Red List; if synonymic, the function finds the accepted name and other synonyms

# (5) Species that are neither in RL or ITIS: use of the online GBIF name matching tool 
# (6) Manual checks for the remaining species identified as synonyms by GBIF ~~ not really..........
# (7) Finally replacing names by accepted names in the original datasets

# Because Elton trait dataset was considered in a second time (second round of data compilation), synonyms for Elton traits
# are gathered only for the species that synonyms had not been gathered for in the first place

## Because this script takes a long time to run, I have saved datasets throughout at different stages


# #  P R E A M B L E  # #  

X <- c("dplyr", "taxize", "phytools", "stringr", "rredlist", "stringdist", "plyr", "pbapply", "GlobalOptions", "data.table", "ngram")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` = Negate(`%in%`)

opt <- options(iucn_redlist_key="ba30954f38deda075bd9b52495d65092ccf1b220b0c7c67a41465646e50ef72c")

source("Resolve_taxonomy_functions.R")

## Load datasets
Predicts <-  readRDS("../../Data/Predicts_database.rds")
Predicts <- subset(Predicts, Class %in% c("Aves", "Amphibia", "Mammalia", "Reptilia"))
Predicts$Best_guess_binomial <- as.character(Predicts$Best_guess_binomial)

# subset Predicts species only known at genus (or other) level

# Non identifiable species
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Duellmanohyla eutisanota"] <- ""

# Just known at genus level
PredictsGenusLevel <- subset(Predicts, Best_guess_binomial == "")

# Phylogenies
PhyloMammal <- read.newick("../../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk") %>% .Format_tiplabels()
PhyloAmphibian <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()
PhyloBird <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels()
PhyloReptile <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk") %>% .Format_tiplabels()

# Traits
Myhrvold <- read.csv("../../Data/Amniotes_Myhrvold_2015/Amniote_Database_Aug_2015.csv", sep=",")
Myhrvold$Binomial_name <- paste(Myhrvold$genus, Myhrvold$species, sep=" ")
MyhrvoldMammal <- subset(Myhrvold, class=="Mammalia")
MyhrvoldBird <- subset(Myhrvold, class=="Aves")
MyhrvoldReptile <- subset(Myhrvold, class=="Reptilia")

## Mammal datsets
Pantheria <- read.csv("../../Data/Mammals/PanTHERIA/Pantheria_1_0_WR05_Aug2008.csv", sep=",")
Pacifici <- read.csv("../../Data/Mammals/PacificiMammals.csv", sep=",")
Elton_mammals <- read.csv("../../Data/Mammals/EltonTraits_Mammals.csv")

## Put mammalian diet information together
Kissling <- read.csv("../../Data/Mammals/Kissling_Mammal_diet_2014.csv", sep=",")
Kissling$Binomial_name <- paste(Kissling$Genus, Kissling$Species, sep=" ")
Kissling <- subset(Kissling, Binomial_name != "Mico sp. nov.")
Kissling <- Kissling %>% select(-TaxonID, -TaxonomicNote, -FillCode)
Kissling <- Kissling[, c(28, 1:27)]

MammalDIET2 <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET2.csv")[, c(1:28)]
MammalDIET2_supp <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET_2_supp.csv")
colnames(MammalDIET2)[c(1, 28)] <- c("Binomial","DataSource")
MammalDIET2 <- rbind(MammalDIET2, MammalDIET2_supp)
MammalDIET2$Binomial <- as.character(MammalDIET2$Binomial)
MammalDIET2$Binomial[MammalDIET2$Binomial=="Rhinolophus hildebrandtii"] <- "Rhinolophus hildebrandti"
MammalDIET2$Binomial[MammalDIET2$Binomial=="Tadarida bivittatus"] <- "Tadarida bivittata"
colnames(MammalDIET2)[1] <- "Binomial_name"

Y <- intersect(MammalDIET2$Binomial_name, Kissling$Binomial_name)
Kissling <- Kissling %>% filter(Binomial_name %nin% Y)
MammalDIET <- rbind(Kissling, MammalDIET2) # dataset to use

## Bird datasets
Butchart_BM <- read.csv("../../Data/Birds/Butchart_BM.csv", sep=",")
Butchart_GL <- read.csv("../../Data/Birds/ButchartGenerationLength.csv", sep=",")
Butchart_GL <- subset(Butchart_GL, Binomial !="")
Sekercioglu_Diet <- read.csv("../../Data/Birds/SekerciogluDiet.csv", sep=",")
Elton_birds <- read.csv("../../Data/Birds/EltonTraits_Birds.csv")

## Reptile datasets
Scharf <- read.csv("../../Data/Reptiles/Scharf.csv", sep=",")
Meiri_0 <- read.csv("../../Data/Reptiles/MeiriReptileMasses.csv", sep=",")
Meiri_0 <- subset(Meiri_0, Rank=="Species")
Vidan <- read.csv("../../Data/Reptiles/Vidan2017_Dielactivity.csv")

Stark <- read.csv("../../Data/Reptiles/Stark2018_GEB_longevity.csv")
colnames(Stark)[1] <- "species"
Schwarz <- read.csv("../../Data/Reptiles/Schwarz_Meiri_GEB_2017.csv")
Novosolov <- read.csv("../../Data/Reptiles/Novosolov_2017_GEB.csv") %>%
  filter(Taxonomic.group!="Birds")%>%
  filter(Taxonomic.group!="Mammals")
Novosolov_2 <- read.csv("../../Data/Reptiles/Novosolov_GEB_2013.csv")
Slavenko <- read.csv("../../Data/Reptiles/Body_sizes_of_all_extant_reptiles_Slavenko_2016_GEB.csv") 
Meiri <- read.csv("../../Data/Reptiles/Meiri_2015_Evolutionary_Biology.csv")

## Amphibian datasets
Amphibio <- read.csv("../../Data/Amphibians/AmphiBIO_v1.csv", sep=",")
Amphibio$Species <- as.character(Amphibio$Species)
Amphibio$Species[Amphibio$Species=="Duttaphrynus pariet alis"] <- "Duttaphrynus parietalis"
Cooper <- read.csv("../../Data/Amphibians/Cooper2008.csv", sep=",")
Cooper <- subset(Cooper, Binomial!="")
Senior <- read.csv("../../Data/Amphibians/Senior_svl_data.csv", sep=",")
Senior <- subset(Senior, Rank=="Species")
Bickford <- read.csv("../../Data/Amphibians/Bickford.csv", sep=",")
Bickford$Binomial_name <-  paste(Bickford$Genus, Bickford$Species, sep=" ")

## Range data
Mammal_range <- read.csv("../../Data/Range_sizes/mammal_range_areas.csv", sep=",")
Amphibian_range <- read.csv("../../Data/Range_sizes/amphibian_range_areas.csv", sep=",")
Bird_range <- read.csv("../../Data/Range_sizes/bird_range_areas.csv", sep=",")
Reptile_range <- read.csv("../../Data/Range_sizes/reptile_range_areas.csv", sep=",")
Reptile_range <- subset(Reptile_range, Binomial_name!="Chelonia mydas Hawaiian subpopulation")

## IUCN data
IUCN_mammal <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Mammals_20180726.csv", sep=",")
IUCN_amphibian <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Amphibians_20180726.csv", sep=",")
IUCN_bird <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Birds_20180726.csv", sep=",")
IUCN_reptile <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Reptiles_20180726.csv", sep=",")



# #  1.  Creating a list of species for which sysnonyms need to be extracted  # #  


## Create lists of species names with corrected typos

## 1.1 Replace vernacular names by scientific names

# Amphibians
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana cafe"] <- "Eleutherodactylus maurus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana hoja"] <- "Noblella lochites"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana hojarasca"] <- "Craugastor angelicus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Rana verde"] <- "Lithobates palmipes"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Arthroleptis sp"] <- ""

# Reptiles
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Serp cafe"] <- "Ninia sebae"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Serpiente cafe"] <- "Ninia sebae"

# Mammals
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Samango monkey"] <- "Cercopithecus albogularis"

# Birds
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Black backed starling"] <- "Acridotheres melanopterus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Black backed weaver"] <- "Ploceus bicolor"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Blue rock pigeon"] <- "Columba livia"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Blue whistling thrush"] <- "Myophonus caeruleus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Brown crested tit"] <- "Lophophanes cristatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Chestnut beilled nuthatch"] <- "Sitta cinnamoventris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Large yellow nape"] <- "Amazona auropalliata"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Common stone chat"] <- "Saxicola torquatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Rufous backed shrike"] <- "Lanius schach"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="White breasted forktail"] <- "Enicurus leschenaulti"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="White breasted kingfisher"] <- "Halcyon smyrnensis"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Jungle babbler"] <- "Turdoides striata"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Grey bulbul"] <- "Pycnonotus cyaniventris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Pied wagtail"] <- "Motacilla alba"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian rufous babbler"] <-"Argya subrufa"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Large wagtail"] <- "Motacilla madaraspatensis"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Common swallow"] <- "Hirundo rustica"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Eurasian tree creeper"] <- "Certhia familiaris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Red headed tit"] <- "Aegithalos iredalei"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Red jungle fowl"] <- "Gallus gallus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian myna"] <-"Acridotheres tristis"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian robin"] <- "Saxicoloides fulicatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Indian tree pie"] <-"Dendrocitta vagabunda"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Sparrow hawk"] <-"Accipiter nisus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Creasted bunting"] <-"Emberiza lathami"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Yellow eyed warbler"] <- "Chrysomma sinense"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Canary blue flycatcher"] <- "Eumyias thalassinus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Hawk booted eagle"] <- "Hieraaetus pennatus"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Yellow backed greenbul"] <- "Chlorocichla flaviventris"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="White piegon"] <- "Columba livia"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Simla crested tit"] <-"Lophophanes dichrous"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Cacomantis esculena"] <- "Collocalia esculenta"
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Bush lark"] <- ""
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Imperial pigeon"] <- ""
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Tailor bird"] <- ""
Predicts$Best_guess_binomial[Predicts$Parsed_name=="Paradise flycatcher"] <- ""


## 1.2. Correct a few noticeable taxonomic errors

# Amphibians
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Kaloula pleurostigma"] <- "Kalophrynus pleurostigma"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Pristimantis bransfordii"] <- "Craugastor bransfordii"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Pristimantis fitzingerii"] <-  "Craugastor fitzingeri"
# Reptiles
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Mabuya elegans"] <- "Lygosoma punctata"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Cryptoblepharus carnabyi"] <- "Cryptoblepharus australis"
# Birds
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Dicaeum maculatus"] <- "Prionochilus maculatus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Dicaeum percussus"] <- "Prionochilus percussus"
Predicts$Best_guess_binomial[Predicts$Best_guess_binomial=="Picoides leucotos"] <- "Dendrocopos leucotos"

# Select data where species known by binomial names
Predicts <- subset(Predicts, Best_guess_binomial != "")



## 1.3. Create lists of species names : check that all have two words (binomial), and return original species list

# for Predicts species
PredictsAmphibians <- SubClass(Predicts, "Amphibia")
PredictsReptiles <- SubClass(Predicts, "Reptilia")
PredictsMammals <- SubClass(Predicts, "Mammalia")
PredictsBirds <- SubClass(Predicts, "Aves")

# for phylogenies
PhyloMammal <- SpNames(PhyloMammal$tip.label)
PhyloAmphibian <- SpNames(PhyloAmphibian$tip.label)
PhyloReptile <- SpNames(PhyloReptile$tip.label)
PhyloBird <- SpNames(PhyloBird$tip.label)

# for trait datasets
MyhrvoldBird <- SpNames(MyhrvoldBird$Binomial_name)
MyhrvoldMammal <- SpNames(MyhrvoldMammal$Binomial_name)
MyhrvoldReptile <- SpNames(MyhrvoldReptile$Binomial_name)

Pantheria <- SpNames(Pantheria$MSW05_Binomial)
Pacifici <- SpNames(Pacifici$Scientific_name)
Kissling <- SpNames(Kissling$Binomial_name)
Elton_mammals <- SpNames(Elton_mammals$Scientific) %>%
  filter(Original!="")

Butchart_BM <- SpNames(Butchart_BM$Sci.name)
Butchart_GL <- SpNames(Butchart_GL$Binomial)
Sekercioglu_Diet <- SpNames(Sekercioglu_Diet$Latin)
Elton_birds <-  SpNames(Elton_birds$Scientific) %>%
  filter(Original!="")

Scharf <- SpNames(Scharf$species)
Meiri_0 <- SpNames(Meiri_0$Taxon)
Vidan <- SpNames(Vidan$Species)
Stark <- SpNames(Stark$species)  %>%
  filter(Original!="")
Schwarz <- SpNames(Schwarz$Species) %>%
  filter(Original!="")
Novosolov <- SpNames(Novosolov$Binomial) %>%
  filter(Original!="")
Novosolov_2 <- SpNames(Novosolov_2$species) %>%
  filter(Original!="")
Slavenko <- SpNames(Slavenko$Binomial) %>%
  filter(Original!="")
Meiri <- SpNames(Meiri$species) %>%
  filter(Original!="")

Amphibio <- SpNames(Amphibio$Species)
Cooper <- SpNames(Cooper$Binomial)
Senior <- SpNames(Senior$Taxon)
Bickford <- SpNames(Bickford$Binomial_name)

Mammal_range <- SpNames(Mammal_range$Species)
Bird_range <- SpNames(Bird_range$Species)
Reptile_range <- SpNames(Reptile_range$Binomial_name)
Amphibian_range <- SpNames(Amphibian_range$Species)

IUCN_amphibian <- SpNames(IUCN_amphibian$binomial)
IUCN_bird <- SpNames(IUCN_bird$binomial)
IUCN_mammal <- SpNames(IUCN_mammal$binomial)
IUCN_reptile <- SpNames(IUCN_reptile$binomial)


## 1.3. Merge species lists and correct typos

# Merge
Mammals <- Reduce(function(x,y) merge(x = x, y = y, by = "Original", all=TRUE),
       list(PredictsMammals, PhyloMammal, MyhrvoldMammal, Pacifici, Kissling, Pantheria, Elton_mammals, Mammal_range, IUCN_mammal))

Birds <- Reduce(function(x,y) merge(x = x, y = y, by = "Original", all=TRUE),
                  list(PredictsBirds, PhyloBird, MyhrvoldBird, Butchart_BM, Butchart_GL, Elton_birds,Sekercioglu_Diet, Bird_range, IUCN_bird))

Reptiles <- Reduce(function(x,y) merge(x = x, y = y, by = "Original", all=TRUE),
                   list(PredictsReptiles, PhyloReptile, MyhrvoldReptile, Meiri_0, Scharf, Vidan, Stark, Schwarz, Novosolov, Novosolov_2, Slavenko, Meiri, Reptile_range, IUCN_reptile))

Amphibians <- Reduce(function(x,y) merge(x = x, y = y, by = "Original", all=TRUE),
                   list(PredictsAmphibians, PhyloAmphibian, Amphibio, Cooper, Senior, Bickford, Amphibian_range, IUCN_amphibian))

rm(PredictsMammals, PhyloMammal, MyhrvoldMammal,
   PredictsBirds, PhyloBird, MyhrvoldBird,
   PredictsAmphibians, PhyloAmphibian,
   PredictsReptiles, PhyloReptile, MyhrvoldReptile)


## Correct for typos: when only one character in corrected names, take original name instead.

system.time({ListMammals <- CheckTypos(Mammals)})
system.time({ListAmphibians <- CheckTypos(Amphibians)})
system.time({ListReptiles <- CheckTypos(Reptiles)})
system.time({ListBirds <- ForBirds(Birds)})

# Manually correct a few entries for birds
ListBirds$CorrectedTypos[ListBirds$Original=="Oneillornis lunulatus"] <- "Oneillornis lunulatus"
ListBirds$CorrectedTypos[ListBirds$Original=="Oneillornis salvini"] <- "Oneillornis salvini"
ListBirds <- subset(ListBirds, Original!="'Scytalopus sp.") # this species is originally in the phylogeny


## Save these files
write.csv(ListBirds, "../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListBirds_V2.csv", row.names = FALSE)
write.csv(ListMammals, "../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListMammals_V2.csv", row.names = FALSE)
write.csv(ListAmphibians, "../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListAmphibians_V2.csv", row.names = FALSE)
write.csv(ListReptiles, "../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListReptiles_V2.csv", row.names = FALSE)



# # # # # # # # # #    S Y N O N Y M  E X T R A C T I O N    # # # # # # # # # # 
 
## Here, the extraction of synonyms from the red list starts.


# ## Load saved files
Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListMammals_V2.csv")
Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListAmphibians_V2.csv")
Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListBirds_V2.csv")
Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_before_resolving/ListReptiles_V2.csv")

## Correcting for typos:
cat("Correcting for typos accross Amphibians reduced number of species by", length(unique(Amphibians$Original)) - length(unique(Amphibians$CorrectedTypos)))
cat("Correcting for typos accross Reptiles reduced number of species by", length(unique(Reptiles$Original)) - length(unique(Reptiles$CorrectedTypos)))
cat("Correcting for typos accross Mammals reduced number of species by", length(unique(Mammals$Original)) - length(unique(Mammals$CorrectedTypos)))
cat("Correcting for typos accross Birds reduced number of species by", length(unique(Birds$Original)) - length(unique(Birds$CorrectedTypos)))


# # 2. Extracting synonyms - from the Red List  # #  

## Mammals
Mammals$CorrectedTypos <- as.character(Mammals$CorrectedTypos)
Synonym_mammals <- RunSyn(Mammals)
write.csv(Synonym_mammals, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_mammals_V2.csv", row.names=F)


## Reptiles
Reptiles$CorrectedTypos <- as.character(Reptiles$CorrectedTypos)
RunSynSave(Reptiles, Path="Reptiles_splits/Reptiles_")

files <- list.files(path = "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Reptiles_splits/", pattern = ".csv")
files <- paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Reptiles_splits/", files, sep="")
Synonym_reptiles <- do.call(rbind,lapply(files,read.csv))
Synonym_reptiles <- Synonym_reptiles[order(as.character(Synonym_reptiles$Original)),]
write.csv(Synonym_reptiles,"../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_reptiles_V2.csv", row.names = FALSE)


## Birds
Birds$CorrectedTypos <- as.character(Birds$CorrectedTypos)
Birds <- RunSynSave(Birds, Path="Birds_splits/Birds_")

files <- list.files(path = "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Birds_splits/", pattern = ".csv")
files <- paste("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Birds_splits/", files, sep="")
Synonym_birds <- do.call(rbind,lapply(files,read.csv))
Synonym_birds <- Synonym_birds[order(as.character(Synonym_birds$Original)),]
write.csv(Synonym_birds,"../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_birds_V2.csv", row.names = FALSE)


## Amphibians
Amphibians$CorrectedTypos <- as.character(Amphibians$CorrectedTypos)
Synonym_amphibians <- RunSyn(Amphibians)
write.csv(Synonym_amphibians, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/SYNONYMS_RL/Synonyms_amphibians.csv", row.names=F)





# #  4. Manual checks and other checks  # # 


# Load files
# Syn_Mammals <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Synonyms_mammals_ITIS.csv") %>% ToChar()
# Syn_Reptiles <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Synonyms_reptiles_ITIS.csv") %>% ToChar()
# Syn_Amphibians <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Synonyms_amphibians_ITIS.csv") %>% ToChar()
# Syn_Birds <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Synonyms_birds_ITIS.csv") %>% ToChar()

# Check that all species names have 2 words
# Check_binomial(Syn_Mammals)
# Check_binomial(Syn_Reptiles)
# Check_binomial(Syn_Birds)
# Check_binomial(Syn_Amphibians)


# decided not to do this bit in the end
# # # Subset and save species that have not been assessed by either Red List or ITIS
# # # NotAssessed_RL_ITIS(Syn_Amphibians, "Amphibians")
# # # NotAssessed_RL_ITIS(Syn_Birds, "Birds")
# # # NotAssessed_RL_ITIS(Syn_Mammals, "Mammals")
# # # NotAssessed_RL_ITIS(Syn_Reptiles, "Reptiles")
# # 
# # # # Load GBIF outputs
# # # GBIF_Mammals <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Not_assessed_RL_ITIS/GBIF_outputs/Mammals_GBIF.csv")
# # # GBIF_Birds <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Not_assessed_RL_ITIS/GBIF_outputs/Birds_GBIF.csv")
# # # GBIF_Reptiles <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Not_assessed_RL_ITIS/GBIF_outputs/Reptiles_GBIF.csv")
# # # GBIF_Amphibians <- read.csv("../../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Before_manual_checks/Not_assessed_RL_ITIS/GBIF_outputs/Amphibians_GBIF.csv")
# # 
# # # # Run function to subset and save names that require manual checks
# # # Require_manual_check(GBIF_Amphibians, "Amphibians")
# # # Require_manual_check(GBIF_Birds, "Birds")
# # # Require_manual_check(GBIF_Mammals, "Mammals")
# # # Require_manual_check(GBIF_Reptiles, "Reptiles")
# 
## NB: these species that have not been assessed by RL or ITIS are going to be matched 
## in GBIF using the name matching tool online. The output will tell whether names are 
## accepted or not so that manual checks can be done only on the names that are found
## to be not accepted.

# chose to go for the alternative instead
## Alternative: just consider that all the remaining names are accepted (function Remaining_species.)
# Species that are in the form Genus + sp. / cf. : known at the genus level only 
# Syn_Mammals <- Syn_Mammals %>% 
#   Remaining_species() %>%
#   ToGenus() %>%
#   GenusCol()
# 
# Syn_Birds <- Syn_Birds  %>% 
#   Remaining_species() %>%
#   ToGenus() %>%
#   GenusCol()
# 
# Syn_Reptiles <- Syn_Reptiles %>% 
#   Remaining_species() %>%
#   ToGenus() %>%
#   GenusCol()
# 
# Syn_Amphibians <- Syn_Amphibians %>% 
#   Remaining_species() %>%
#   ToGenus() %>%
#   GenusCol()


## Load files
# Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Mammals.csv")
# Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Birds.csv")
# Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Reptiles.csv")
# Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Amphibians.csv")
# 
# Elton_mammals_toadd$Genus_level <- "FALSE"
# Elton_mammals_toadd$Genus <- word(Elton_mammals_toadd$CorrectedTypos,1)
# 
# Elton_birds_toadd$Genus_level <- "FALSE"
# Elton_birds_toadd$Genus <- word(Elton_birds_toadd$CorrectedTypos,1)
# 
# Syn_Mammals <- rbind(Syn_Mammals, Elton_mammals_toadd)
# Syn_Mammals <- Syn_Mammals[order(Syn_Mammals$Original),]
# 
# Syn_Birds <- rbind(Syn_Birds, Elton_birds_toadd)
# Syn_Birds <- Syn_Birds[order(Syn_Birds$Original),]
# 
# write.csv(Syn_Birds,"../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Birds_with_elton.csv", row.names = FALSE)
# write.csv(Syn_Mammals,"../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Mammals_with_elton.csv", row.names = FALSE)


# # # # # # #        S Y N O N Y M     E X T R A C T I O N    F I N I S H E D       # # # # #  




# #  5. Matching original datasets and replacing by accepted names  # # 


## Need to run the relevant previous commented lines on Predicts, phylogenies and trait datasets before (correcting for names) (uncommented lines)

## Load files
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Mammals_with_elton.csv")
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Birds_with_elton.csv")
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Reptiles.csv")
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/After_manual_checks/Amphibians.csv")

length(unique(Syn_Mammals$Accepted)) - length(unique(Syn_Mammals$Original))
length(unique(Syn_Birds$Accepted)) - length(unique(Syn_Birds$Original))
length(unique(Syn_Reptiles$Accepted)) - length(unique(Syn_Reptiles$Original))
length(unique(Syn_Amphibians$Accepted)) - length(unique(Syn_Amphibians$Original))

# Trait datasets
MyhrvoldMammal <- Replace_by_accepted_name(Syn_Mammals, MyhrvoldMammal, "Binomial_name")
MyhrvoldBird <- Replace_by_accepted_name(Syn_Birds, MyhrvoldBird, "Binomial_name")
MyhrvoldReptile <- Replace_by_accepted_name(Syn_Reptiles, MyhrvoldReptile, "Binomial_name")

Pantheria <- Replace_by_accepted_name(Syn_Mammals, Pantheria, "MSW05_Binomial")
Pacifici <- Replace_by_accepted_name(Syn_Mammals, Pacifici, "Scientific_name")
#Kissling <- Replace_by_accepted_name(Syn_Mammals, Kissling, "Binomial_name")
MammalDIET <- Replace_by_accepted_name(Syn_Mammals, MammalDIET, "Binomial_name")

Butchart_BM <- Replace_by_accepted_name(Syn_Birds, Butchart_BM, "Sci.name")
Butchart_GL <- Replace_by_accepted_name(Syn_Birds, Butchart_GL, "Binomial")
Sekercioglu_Diet$Latin <- as.character(Sekercioglu_Diet$Latin)
Sekercioglu_Diet$Latin[Sekercioglu_Diet$Latin=="Zosterops xanthochrous"] <- "Zosterops xanthochroa"
Sekercioglu_Diet <- Replace_by_accepted_name(Syn_Birds, Sekercioglu_Diet, "Latin")

Scharf <- Replace_by_accepted_name(Syn_Reptiles, Scharf, "species")
Meiri_0 <- Replace_by_accepted_name(Syn_Reptiles, Meiri_0, "Taxon")

Amphibio <- Replace_by_accepted_name(Syn_Amphibians, Amphibio, "Species")
Cooper <- Replace_by_accepted_name(Syn_Amphibians, Cooper, "Binomial")
Senior <- Replace_by_accepted_name(Syn_Amphibians, Senior, "Taxon")
Bickford <- Replace_by_accepted_name(Syn_Amphibians, Bickford, "Binomial_name")

Mammal_range <- Replace_by_accepted_name(Syn_Mammals, Mammal_range, "Species")
Amphibian_range <- Replace_by_accepted_name(Syn_Amphibians, Amphibian_range, "Species")
Bird_range <- Replace_by_accepted_name(Syn_Birds, Bird_range, "Species")
Reptile_range <- Replace_by_accepted_name(Syn_Reptiles, Reptile_range, "Binomial_name")

IUCN_amphibian <-  Replace_by_accepted_name(Syn_Amphibians, IUCN_amphibian, "binomial")
IUCN_bird <-  Replace_by_accepted_name(Syn_Birds, IUCN_bird, "binomial")
IUCN_reptile <-  Replace_by_accepted_name(Syn_Reptiles, IUCN_reptile, "binomial")
IUCN_mammal <-  Replace_by_accepted_name(Syn_Mammals, IUCN_mammal, "binomial")

Elton_birds <- Replace_by_accepted_name(Syn_Birds, Elton_birds, "Scientific")
Elton_mammals <- Replace_by_accepted_name(Syn_Mammals, Elton_mammals, "Scientific")

# Predicts dataset
PredictsAmphibians <- subset(Predicts, Class=="Amphibia")
PredictsReptiles <- subset(Predicts, Class=="Reptilia")
PredictsMammals <- subset(Predicts, Class=="Mammalia")
PredictsBirds <- subset(Predicts, Class=="Aves")

PredictsAmphibians <- Replace_by_accepted_name(Syn_Amphibians, PredictsAmphibians, "Best_guess_binomial")
PredictsReptiles <- Replace_by_accepted_name(Syn_Reptiles, PredictsReptiles, "Best_guess_binomial")
PredictsMammals <- Replace_by_accepted_name(Syn_Mammals, PredictsMammals, "Best_guess_binomial")
PredictsBirds <- Replace_by_accepted_name(Syn_Birds, PredictsBirds, "Best_guess_binomial")

Predicts_Vertebrates <- rbind(PredictsAmphibians, PredictsBirds, PredictsMammals, PredictsReptiles, PredictsGenusLevel)
Predicts_Vertebrates <- Predicts_Vertebrates[order(Predicts_Vertebrates$SSBS),]


# Phylogenies
PhyloMammal <- Replace_by_accepted_name(Syn_Mammals, PhyloMammal)
PhyloBird <- Replace_by_accepted_name(Syn_Birds, PhyloBird)
PhyloReptile <- Replace_by_accepted_name(Syn_Reptiles, PhyloReptile)
PhyloAmphibian <- Replace_by_accepted_name(Syn_Amphibians, PhyloAmphibian)


## Saving processed datasets
Myhrvold <- rbind(MyhrvoldBird, MyhrvoldMammal, MyhrvoldReptile)
write.csv(Myhrvold, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Myhrvold.csv", row.names=FALSE)
write.csv(Pantheria, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Pantheria.csv", row.names=FALSE)
write.csv(Pacifici, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Pacifici.csv", row.names=FALSE)
write.csv(MammalDIET, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/MammalDIET.csv", row.names=FALSE)
write.csv(Butchart_BM, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Butchart_BM.csv", row.names=FALSE)
write.csv(Butchart_GL, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Butchart_GL.csv", row.names=FALSE)
write.csv(Sekercioglu_Diet, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Sekercioglu_Diet.csv", row.names=FALSE)
write.csv(Scharf, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Scharf.csv", row.names=FALSE)
write.csv(Meiri_0, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Meiri_0.csv", row.names=FALSE)
write.csv(Amphibio, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibio.csv", row.names=FALSE)
write.csv(Cooper, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Cooper.csv", row.names=FALSE)
write.csv(Senior, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Senior.csv", row.names=FALSE)
write.csv(Bickford, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bickford.csv", row.names=FALSE)
write.csv(Mammal_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Mammal_range.csv", row.names=FALSE)
write.csv(Amphibian_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Amphibian_range.csv", row.names=FALSE)
write.csv(Bird_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Bird_range.csv", row.names=FALSE)
write.csv(Reptile_range, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Reptile_range.csv", row.names=FALSE)
write.csv(IUCN_amphibian, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_amphibians.csv", row.names=FALSE)
write.csv(IUCN_bird, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_birds.csv", row.names=FALSE)
write.csv(IUCN_mammal, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_mammals.csv", row.names=FALSE)
write.csv(IUCN_reptile, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_reptiles.csv", row.names=FALSE)
write.csv(Elton_birds, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Elton_birds.csv", row.names=FALSE)
write.csv(Elton_mammals, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Elton_mammals.csv", row.names=FALSE)

saveRDS(Predicts_Vertebrates, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

write.tree(PhyloAmphibian, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloAmphibians.nwk")
write.tree(PhyloBird, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")
write.tree(PhyloMammal, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloMammals.nwk")
write.tree(PhyloReptile, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloReptiles.nwk")


