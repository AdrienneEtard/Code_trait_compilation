library(dplyr)

ModifiedAmphibians <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousAmphibians_log10.csv") 
Amphibians <- ModifiedAmphibians  %>% 
  t() %>% as.data.frame()%>% select(c(1,4)) %>% setNames(c("lambda", "p-value")) %>% mutate(Class="Amphibians")
Amphibians <- Amphibians %>% mutate(Trait=colnames(ModifiedAmphibians), Phylogeny="Original")

ModifiedBirds <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousBirds_log10.csv") 
Birds <- ModifiedBirds %>% 
  t() %>% as.data.frame()%>% select(c(1,4)) %>% setNames(c("lambda", "p-value"))%>% mutate(Class="Birds", Trait=colnames(ModifiedBirds), Phylogeny="Original")

ModifiedMammals <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousMammals_log10.csv") 
Mammals <- ModifiedMammals %>% 
  t() %>% as.data.frame()%>% select(c(1,4)) %>% setNames(c("lambda", "p-value"))%>% mutate(Class="Mammals", Trait=colnames(ModifiedMammals), Phylogeny="Original")


ModifiedReptiles <- read.csv("../../Results/1.Traits_before_imputations/Phylogenetic_signal/with_original_phylogenies/ContinuousReptiles_log10.csv") 
Reptiles <- ModifiedReptiles %>% 
  t() %>% as.data.frame()%>% select(c(1,4)) %>% setNames(c("lambda", "p-value"))%>% mutate(Class="Reptiles", Trait=colnames(ModifiedReptiles), Phylogeny="Original")


PhyContModified <- rbind(Birds, Reptiles, Mammals, Amphibians)
PhyContModified <- PhyContModified[, c("Class", "Trait", "lambda", "p-value", "Phylogeny")]
PhyContModified <- PhyContModified %>% mutate(Trait=ifelse(Trait=="Body_mass_g", "Body mass", 
                                                           ifelse(Trait=="Longevity_d", "Longevity",
                                                                  ifelse(Trait=="Litter_size", "Litter/clutch size", 
                                                                         ifelse(Trait=="Diet_breadth", "Diet breath", 
                                                                                ifelse(Trait=="Range_size_m2", "Range size", 
                                                                                       ifelse(Trait=="Habitat_breadth_IUCN", "Habitat breadth", 
                                                                                              ifelse(Trait=="Generation_length_d", "Generation length", 
                                                                                                     ifelse(Trait=="Maturity_d", "Maturity", 
                                                                                                            ifelse(Trait=="Adult_svl_cm", "Body length", 
                                                                                                                   ifelse(Trait=="Body_length_mm", "Body length", Trait)))))))))))


write.csv(PhyContModified, "../../Results/1.Traits_before_imputations/Phylogenetic_signal/Continuous_traits_significance_chisquared_original.csv", row.names = FALSE)

