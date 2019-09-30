## quick fix (before imputing again)

Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

Mammals <- Mammals[order(Mammals$Best_guess_binomial), ]
Birds <- Birds[order(Birds$Best_guess_binomial), ]
Amphibians <- Amphibians[order(Amphibians$Best_guess_binomial), ]
Reptiles <- Reptiles[order(Reptiles$Best_guess_binomial), ]


Imputed <- readRDS("../../Results/2.imputed_trait_datasets/Imputed_corrected_trees/List_of_8_sets.rds")

for (i in 1:8) {
  
  Imputed[[i]][["M"]]$Imputed.Dataset$Order <- Mammals$Order
  Imputed[[i]][["M"]]$Imputed.Dataset$Family <- Mammals$Family
  Imputed[[i]][["M"]]$Imputed.Dataset$Genus <- Mammals$Genus
  
  Imputed[[i]][["B"]]$Imputed.Dataset$Order <- Birds$Order
  Imputed[[i]][["B"]]$Imputed.Dataset$Family <- Birds$Family
  Imputed[[i]][["B"]]$Imputed.Dataset$Genus <- Birds$Genus
  
  Imputed[[i]][["R"]]$Imputed.Dataset$Order <- Reptiles$Order
  Imputed[[i]][["R"]]$Imputed.Dataset$Family <- Reptiles$Family
  Imputed[[i]][["R"]]$Imputed.Dataset$Genus <- Reptiles$Genus
  
  Imputed[[i]][["A"]]$Imputed.Dataset$Order <- Amphibians$Order
  Imputed[[i]][["A"]]$Imputed.Dataset$Family <- Amphibians$Family
  Imputed[[i]][["A"]]$Imputed.Dataset$Genus <- Amphibians$Genus
  
}

saveRDS(Imputed, "../../Results/2.imputed_trait_datasets/Imputed_corrected_trees/List_of_8_sets_quick_fix.rds" )

Imputed_notfixed <- readRDS("../../Results/2.imputed_trait_datasets/Imputed_corrected_trees/List_of_8_sets.rds")

unique(Imputed_notfixed[[i]][["M"]]$Imputed.Dataset$Order==Imputed[[i]][["M"]]$Imputed.Dataset$Order)
unique(Imputed_notfixed[[i]][["B"]]$Imputed.Dataset$Order==Imputed[[i]][["B"]]$Imputed.Dataset$Order)
unique(Imputed_notfixed[[i]][["R"]]$Imputed.Dataset$Order==Imputed[[i]][["R"]]$Imputed.Dataset$Order)
unique(Imputed_notfixed[[i]][["A"]]$Imputed.Dataset$Order==Imputed[[i]][["A"]]$Imputed.Dataset$Order)

unique(Imputed_notfixed[[i]][["M"]]$Imputed.Dataset$Family==Imputed[[i]][["M"]]$Imputed.Dataset$Family)
unique(Imputed_notfixed[[i]][["B"]]$Imputed.Dataset$Family==Imputed[[i]][["B"]]$Imputed.Dataset$Family)
unique(Imputed_notfixed[[i]][["R"]]$Imputed.Dataset$Family==Imputed[[i]][["R"]]$Imputed.Dataset$Family)
unique(Imputed_notfixed[[i]][["A"]]$Imputed.Dataset$Family==Imputed[[i]][["A"]]$Imputed.Dataset$Family)

unique(Imputed_notfixed[[i]][["M"]]$Imputed.Dataset$Genus==Imputed[[i]][["M"]]$Imputed.Dataset$Genus)
unique(Imputed_notfixed[[i]][["B"]]$Imputed.Dataset$Genus==Imputed[[i]][["B"]]$Imputed.Dataset$Genus)
unique(Imputed_notfixed[[i]][["R"]]$Imputed.Dataset$Genus==Imputed[[i]][["R"]]$Imputed.Dataset$Genus)
unique(Imputed_notfixed[[i]][["A"]]$Imputed.Dataset$Genus==Imputed[[i]][["A"]]$Imputed.Dataset$Genus)


