## Missing values at random?

## Preamble
X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree")
lapply(X, library, character.only=TRUE); rm(X)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}

## Load trait data 
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/MammalsComplete.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/BirdsComplete.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/AmphibiansComplete.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/ReptilesComplete.csv")

## Load phylogenies
Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/PhyloMammals.nwk") %>% .Format_tiplabels()
Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/PhyloAmphibians.nwk") %>% .Format_tiplabels() 
Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/PhyloReptiles.nwk") %>% .Format_tiplabels()
Phylo_Birds <- read.newick("../../Results/1.Phylogenies/PhyloBirds.nwk") %>% .Format_tiplabels() 



## 1. Are missing values missing at random?

MissingVal <- function (DF, Order_or_Family, Traits) {
  
  # browser()
  
  List <- list()
  
  for (i in Traits) {
    PercentNA <- table(DF[,Order_or_Family]) %>% 
      as.data.frame() %>%
      setNames(., c(Order_or_Family, "Total"))
    
    Missing <- subset(DF, is.na(DF[, i]))
    
    Missing <- table(Missing[, Order_or_Family])
    MissOrder <- names(Missing)
    PercentNA$Missing[PercentNA[, Order_or_Family] %in% MissOrder] <- Missing
    `%nin%` <- Negate(`%in%`)
    PercentNA$Missing[PercentNA[, Order_or_Family] %nin% MissOrder] <- 0
    
    
    ToBind <- c("TOTAL", sum(PercentNA$Total), sum(PercentNA$Missing))
    PercentNA[, Order_or_Family] <- as.character(PercentNA[, Order_or_Family])
    PercentNA <- rbind(PercentNA, ToBind)
    
    PercentNA$Total <- as.numeric(PercentNA$Total)
    PercentNA$Missing <- as.numeric(PercentNA$Missing)
    PercentNA$Percent <- PercentNA$Missing / PercentNA$Total * 100
    
    List[[i]] <- PercentNA
  }
  
  All <- Reduce(function(x, y) merge(x, y, by=Order_or_Family), List)
  ToPaste <- c()
  for (i in names(List)) { ToPaste <- c(ToPaste, rep(i,3)) }
  colnames(All) <- c(Order_or_Family, paste(colnames(All)[-1], ToPaste))
  
  All <- All %>% select(which(grepl(paste(Order_or_Family, "|Percent", sep=""), colnames(All))))
  colnames(All) <- c(Order_or_Family, names(List))
  
  del <- rownames(All[All[, Order_or_Family]=="TOTAL",]) %>% as.numeric
  All <- rbind(All[-del, ], All[All[, Order_or_Family]=="TOTAL", ])
  
  return(list(List=List, Percentages=All))
}

## Order level

NA_Amphibians <- MissingVal(Amphibians, "Order", c("Body_mass_g", "Maturity_d", "Litter_size", "Trophic_level", "Diet_breadth","Habitat_breadth_IUCN","Range_size_m2"))
NA_Mammals <- MissingVal(Mammals, "Order", c("Body_mass_g", "Generation_length_d", "Litter_size", "Trophic_level", "Diet_breadth", "Habitat_breadth_IUCN","Range_size_m2"))
NA_Reptiles <- MissingVal(Reptiles, "Order", c("Body_mass_g", "Maturity_d", "Litter_size", "Trophic_level","Habitat_breadth_IUCN","Range_size_m2"))
NA_Birds <- MissingVal(Birds, "Order", c("Body_mass_g", "Generation_length_d", "Litter_size", "Trophic_level", "Diet_breadth","Habitat_breadth_IUCN","Range_size_m2"))

## Plot the distribution of % of missing values - boxplot for each trait across orders/families
pdf(file="../../Results/Plots/Trait_missing_values/Per_order.pdf", width=10, height=7, family="Times", pointsize=11)
par(las=1, oma=c(4,5.5,0.5,0), mar=c(3,5,2,8))
par(mfrow=c(2,2), xpd=TRUE)
boxplot(NA_Amphibians$Percentages[1:nrow(NA_Amphibians$Percentages)-1,2:8], horizontal = TRUE, 
        main="Amphibians", ylim=c(0,100))
text(rep(120, 7), c(1:7), paste(round(NA_Amphibians$Percentages[4,2:8], 1), "%"), col="blue")

boxplot(NA_Reptiles$Percentages[1:nrow(NA_Reptiles$Percentages)-1,2:7], horizontal = TRUE, 
        main="Reptiles", ylim=c(0,100))
text(rep(120, 6), c(1:6), paste(round(NA_Reptiles$Percentages[5,2:7], 1), "%"), col="blue")

boxplot(NA_Birds$Percentages[1:nrow(NA_Birds$Percentages)-1,2:8], horizontal = TRUE, 
        main="Birds", ylim=c(0,100))
text(rep(120, 7), c(1:7), paste(round(NA_Birds$Percentages[38,2:8], 1), "%"), col="blue")

boxplot(NA_Mammals$Percentages[1:nrow(NA_Mammals$Percentages),2:8], horizontal = TRUE, 
        main="Mammals", ylim=c(0,100))
text(rep(120, 7), c(1:7), paste(round(NA_Mammals$Percentages[28,2:8], 1), "%"), col="blue")

mtext("Distribution of the proportion of missing values per orders", side=1, line=4, at=-50)
dev.off()


## Family level
NA_Amphibians_F <- MissingVal(Amphibians, "Family", c("Body_mass_g", "Maturity_d", "Litter_size", "Trophic_level", "Diet_breadth","Habitat_breadth_IUCN","Range_size_m2"))
NA_Mammals_F <- MissingVal(Mammals, "Family", c("Body_mass_g", "Generation_length_d", "Litter_size", "Trophic_level", "Diet_breadth", "Habitat_breadth_IUCN","Range_size_m2"))
NA_Reptiles_F <- MissingVal(Reptiles, "Family", c("Body_mass_g", "Maturity_d", "Litter_size", "Trophic_level","Habitat_breadth_IUCN","Range_size_m2"))
NA_Birds_F <- MissingVal(Birds, "Family", c("Body_mass_g", "Generation_length_d", "Litter_size", "Trophic_level", "Diet_breadth","Habitat_breadth_IUCN","Range_size_m2"))

pdf(file="../../Results/Plots/Trait_missing_values/Per_family.pdf", width=10, height=7, family="Times", pointsize=11)
par(las=1, oma=c(4,5.5,0.5,0), mar=c(3,5,2,8))
par(mfrow=c(2,2), xpd=TRUE)
boxplot(NA_Amphibians_F$Percentages[1:nrow(NA_Amphibians_F$Percentages)-1,2:8], horizontal = TRUE, 
        main="Amphibians", ylim=c(0,100))
text(rep(115, 7), c(1:7), paste(round(NA_Amphibians_F$Percentages[75,2:8], 1), "%"), col="blue")

boxplot(NA_Reptiles_F$Percentages[1:nrow(NA_Reptiles_F$Percentages)-1,2:7], horizontal = TRUE, 
        main="Reptiles", ylim=c(0,100))
text(rep(115, 6), c(1:6), paste(round(NA_Reptiles_F$Percentages[88,2:7], 1), "%"), col="blue")

boxplot(NA_Birds_F$Percentages[1:nrow(NA_Birds_F$Percentages)-1,2:8], horizontal = TRUE, 
        main="Birds", ylim=c(0,100))
text(rep(115, 7), c(1:7), paste(round(NA_Birds_F$Percentages[nrow(NA_Birds_F$Percentages),2:8], 1), "%"), col="blue")

boxplot(NA_Mammals_F$Percentages[1:nrow(NA_Mammals_F$Percentages),2:8], horizontal = TRUE, 
        main="Mammals", ylim=c(0,100))
text(rep(115, 7), c(1:7), paste(round(NA_Mammals_F$Percentages[nrow(NA_Mammals_F$Percentages),2:8], 1), "%"), col="blue")

mtext("Distribution of the proportion of missing values per families", side=1, line=4, at=-50)
dev.off()



## Plotting percent information across species

Percent_info <- function(TraitDF, Traits) {
  
  TraitDF <- TraitDF[, c("Best_guess_binomial", Traits)]
  Results <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y))) %>% as.data.frame() %>%
    setNames(., "Percent")
  Results$Percent <- Results$Percent/length(Traits)*100
  
  return(Results)
  
}

Test <- Percent_info(Mammals, TraitsMammals)
Test <- Percent_info(Birds, TraitsBirds)
Test <- Percent_info(Amphibians, TraitsAmphibians)
Test <- Percent_info(Reptiles, TraitsReptiles)


plot(density(Test$Percent))
hist(Test$Percent, breaks = 10)



