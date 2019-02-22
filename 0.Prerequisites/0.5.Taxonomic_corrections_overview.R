## Taxonomic corrections: quick overview
library(dplyr)

## Load synonym datasets
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals.csv") 
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv")
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles.csv")
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians.csv")

## Differences in overall species number

print("Differences in species number across all compiled names after all taxonomic correction:")
print(paste(length(unique(Syn_Mammals$Original)) - length(unique(Syn_Mammals$Accepted)), "for mammals"))
print(paste(length(unique(Syn_Birds$Original)) - length(unique(Syn_Birds$Accepted)), "for birds"))
print(paste(length(unique(Syn_Amphibians$Original)) - length(unique(Syn_Amphibians$Accepted)), "for amphibians"))
print(paste(length(unique(Syn_Reptiles$Original)) - length(unique(Syn_Reptiles$Accepted)), "for reptiles"))


print("Differences in species number across all compiled names after correcting for typos:")
print(paste(length(unique(Syn_Mammals$Original)) - length(unique(Syn_Mammals$CorrectedTypos)), "for mammals"))
print(paste(length(unique(Syn_Birds$Original)) - length(unique(Syn_Birds$CorrectedTypos)), "for birds"))
print(paste(length(unique(Syn_Amphibians$Original)) - length(unique(Syn_Amphibians$CorrectedTypos)), "for amphibians"))
print(paste(length(unique(Syn_Reptiles$Original)) - length(unique(Syn_Reptiles$CorrectedTypos)), "for reptiles"))


print("Differences in species number across all compiled names after correcting for typos and extracting synonyms:")
print(paste(length(unique(Syn_Mammals$Accepted)) - length(unique(Syn_Mammals$CorrectedTypos)), "for mammals"))
print(paste(length(unique(Syn_Birds$Accepted)) - length(unique(Syn_Birds$CorrectedTypos)), "for birds"))
print(paste(length(unique(Syn_Amphibians$Accepted)) - length(unique(Syn_Amphibians$CorrectedTypos)), "for amphibians"))
print(paste(length(unique(Syn_Reptiles$Accepted)) - length(unique(Syn_Reptiles$CorrectedTypos)), "for reptiles"))


## Identified synonyms: which species was the most replicated across datasets? 

MaxSyn <- function(Syn) {
  Table <- Syn %>% 
    filter(Accepted!="") %>%
    group_by(Accepted) %>% 
    summarise(Count=n())
  Max <- Table$Accepted[Table$Count==max(Table$Count)]
  All <- Syn$CorrectedTypos[Syn$Accepted==Max]
  return(list(Accepted=Max, Replicates=All))
}

# For mammals: Tachyoryctes splendens appeared under 12 different names
Mammals <- MaxSyn(Syn_Mammals)
Birds <- MaxSyn(Syn_Birds)
Reptiles <- MaxSyn(Syn_Reptiles)
Amphibians <- MaxSyn(Syn_Amphibians)


