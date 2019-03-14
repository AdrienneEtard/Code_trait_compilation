library(phytools)
library(dplyr)
library(ape)
library(picante)

## ftype="off"

## Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}


`%nin%` <- Negate(`%in%`)

PhyloAmphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk") %>% .Format_tiplabels()
PhyloMammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk")  %>% .Format_tiplabels()
PhyloBirds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk")  %>% .Format_tiplabels() 
PhyloReptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk")  %>% .Format_tiplabels() 

Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

PlotMissinfvalues_phylogeny <- function(Phylo, TraitDF, OorF, trait) {
  
  TraitDF[, OorF] <- tolower( TraitDF[, OorF])
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  TraitDF[, OorF] <- firstup(TraitDF[, OorF])
  
  # Prune the tree to keep species that will represent families or orders figuring in the trait dataset
  Taxa=TraitDF[, OorF] %>% as.vector()
  Names=TraitDF$Best_guess_binomial
  names(Taxa) <- Names
  Match = match.phylo.data(Phylo, Taxa)
  Phylo <- Match$phy
  Taxa <- Match$data
  
  # then for each order keep only one species
  Taxa <- as.data.frame(Taxa)
  Taxa$Species <- rownames(Taxa)
  Taxa <- Taxa %>% dplyr::group_by(Taxa) %>%
    dplyr::slice(1)
  
  Pattern <- paste(Taxa$Species, sep=" ")
  
  Phylo$tip.label[Phylo$tip.label %nin% Pattern] <- paste(Phylo$tip.label[Phylo$tip.label %nin% Pattern],
                                                          "to_drop", sep=" ")
  
  
  # drop tree tips that do not match the species names
  Phylo <- drop.tip(Phylo, Phylo$tip.label[grepl("to_drop", Phylo$tip.label)])
  
  
  # paste family or order next to the corresponding tree tip
  for (i in 1:nrow(Taxa)) {
    Phylo$tip.label[Phylo$tip.label==Taxa$Species[i]] <- Taxa$Taxa[i] %>% as.character()
  }
  
  # for these families or orders calculate percent missing values for the trait
  TraitDF[, trait] <- as.character(TraitDF[, trait])
  TraitDF[is.na(TraitDF[, trait]), trait] <- "Missing"
  TraitDF[ TraitDF[, trait]!="Missing", trait] <- "Non_missing"
  MissVal <- table(TraitDF[, OorF], TraitDF[, trait]) %>%
    as.data.frame() %>%
    setNames(., c("Taxa", "Trait","nspecies"))
  
  Missing <- MissVal %>% 
    filter(Trait=="Missing") %>%
    select(Taxa, nspecies) %>%
    setNames(., c("Taxa", "n_sp_missing"))
  NonMissing<- MissVal %>% filter(Trait=="Non_missing") %>%
    select("nspecies") %>%
    setNames(., "n_sp_not_missing")
  
  PercentNA <- cbind(Missing, NonMissing) %>%
    mutate(Percent=n_sp_missing/(n_sp_not_missing+n_sp_missing)*100)
  
  PercentNA <- PercentNA %>%
    filter(Taxa %in% Phylo$tip.label)

  
  trait_toplot <- PercentNA$Percent %>% as.vector()
 
  names(trait_toplot) <- PercentNA$Taxa
  
  # Use contMap from phytools
  tree.trait <- contMap(Phylo, trait_toplot, plot = FALSE, res=10)
  tree.trait <-setMap(tree.trait, colors=c("blue","cyan", "green","yellow","red"))
  # Carefull, we now plot using the plotting tool from phytools
  plot.contMap(tree.trait,type = "phylogram", cex=0.5, legend=FALSE)
  
  # plotSimmap(tree.trait$tree,tree.trait$cols, fsize=0, type = "phylogram")
  # add.color.bar(0.8,cols=tree.trait$cols,title="% missing values across species", lims=range(trait_toplot),digits=2)
  
  # Phylo$tip.label
  # plot(Phylo, type="fan", cex=0.5)
}

Traits_cont <-  c("Body_mass_g", "Longevity_d", "Litter_size", "Diet_breadth",
                  "Range_size_m2", "Habitat_breadth_IUCN")

Traits_cat <- c("Specialisation", "Diel_activity","Trophic_level", "Primary_diet")


## Amphibians
par(mfrow=c(5,3))

PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Range_size_m2")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Body_length_mm")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Habitat_breadth_IUCN")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Specialisation")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Diel_activity")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Litter_size")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Primary_diet")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Trophic_level")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Diet_breadth")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Body_mass_g")
PlotMissinfvalues_phylogeny(PhyloAmphibians, Amphibians, "Family", "Longevity_d")


# Reptiles
par(mfrow=c(5,2))

PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Body_mass_g")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Range_size_m2")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Litter_size")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Habitat_breadth_IUCN")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Specialisation")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Diel_activity")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Longevity_d")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Trophic_level")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Body_length_mm")
PlotMissinfvalues_phylogeny(PhyloReptiles, Reptiles, "Family", "Maturity_d")


# Birds
par(mfrow=c(5,3))
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Habitat_breadth_IUCN")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Specialisation")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Body_mass_g")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Range_size_m2")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Primary_diet")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Diel_activity")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Generation_length_d")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Trophic_level")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Diet_breadth")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Litter_size")
PlotMissinfvalues_phylogeny(PhyloBirds, Birds, "Family", "Longevity_d")

# Mammals
par(mfrow=c(5,3))
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Body_mass_g")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Generation_length_d")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Range_size_m2")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Trophic_level")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Diel_activity")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Primary_diet")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Habitat_breadth_IUCN")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Specialisation")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Diet_breadth")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Adult_svl_cm")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Litter_size")
PlotMissinfvalues_phylogeny(PhyloMammals, Mammals, "Family", "Longevity_d")

















