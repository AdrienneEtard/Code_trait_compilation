############################################################
# AIM: RETURNING TRAIT DATASETS
# Loading trait datasets
# Matching by species (/genus /family /order) name in Predicts using "Best_binomial_guess" as species name)
# Extracting trait information for matching species by trait dataset
# Interpolating data to get a maximum of information
# RETURNS: a complete trait dataset for each vertebrate class for which Binomial names are known
# + interpolated trait values for others

# Uses functions that are sourced

## /!\ UNITS TO BE CHECKED AND NORMALISED
## /!\ SCRIPT WRITTEN FOR MAMMALS (for now....)
## /!\ interpolations can be refined using phylogenetic imputations


## Preamble ----------------------------------------------------------------
X <- c("data.table", "plyr", "dplyr", "tidyr", "magrittr", "reshape", "reshape2", "stringr") 
invisible(lapply(X, library, character.only=TRUE)); rm(X)

Data_dir <- "../Data/"
Out_dir <- "../Results/Complete_Diversity_Trait_dataset"
source("./1.Sourced/1.Prepare_trait_data.R")
source("./1.Sourced/1.Trait_data_interpolation.R")

## Load data ---------------------------------------------------------------
Predicts <-  readRDS(paste(Data_dir,"PREDICTS/predicts_database.rds", sep=""))
Pantheria <- read.csv(paste(Data_dir, "Mammals/PanTHERIA/Pantheria_1_0_WR05_Aug2008.csv", sep=""), sep=",") 
Pacifici <- read.csv(paste(Data_dir, "Mammals/PacificiMammals.csv", sep=""), sep=",") 
Kissling <- read.csv(paste(Data_dir, "Mammals/Kissling_Mammal_diet_2014.csv", sep=""), sep=",")
Myhrvold <- read.csv(paste(Data_dir, "Amniotes_Myhrvold_2015/Amniote_Database_Aug_2015.csv", sep=""), sep=",")
Mammal_range <- read.csv(paste(Data_dir, "Range_sizes/mammal_range_areas.csv", sep=""), sep=",")

## Normalise datasets ------------------------------------------------------
X <- .Normalise_TDB.Mammals(Kissling, Pantheria, Pacifici, Myhrvold, Mammal_range)
Myhrvold <- X$Myhrvold; Pantheria <- X$Pantheria; Pacifici <- X$Pacifici; Mammal_range <- X$Range
Kissling <-  X$Kissling

## One bigger function? ##~~~~~~~~~
## Extract data from matches by binomial species names ---------------------
Mammals.Pantheria <- .Match_Species("Mammalia", Predicts, Pantheria, X$Mammal_traits.Pantheria)
Mammals.Pacifici <- .Match_Species("Mammalia", Predicts, Pacifici, X$Mammal_traits.Pacifici)
Mammals.Myhrvold <- .Match_Species("Mammalia", Predicts, Myhrvold, X$Mammal_traits.Myhrvold)
Mammals.Kissling <- .Match_Species("Mammalia", Predicts, Kissling, X$Mammal_traits.Kissling)
Mammals.range <- .Match_Species("Mammalia", Predicts, Mammal_range, "Range_size_m2")

# For Pantheria, remove Trophic level and Diet information; indeed, it clashes with Kissling for some species
# Diet information from Kissling.
Mammals.Pantheria %<>% select(-c(Trophic_level, Diet_breadth))

## Merge and reduce redundancy ---------------------------------------------
Traits <- merge(Mammals.Pacifici, Mammals.Myhrvold, all=TRUE)
Traits <- merge(Traits, Mammals.Pantheria, all=TRUE) 
Traits <- merge(Traits, Mammals.range, all=TRUE) 
Traits <- merge(Traits, Mammals.Kissling, all=TRUE) 
row.names(Traits) <- c(1:nrow(Traits))
rm(X, Mammals.Kissling, Mammals.Myhrvold, Mammals.Pacifici, Mammals.Pantheria, Mammals.range)

# Reduce redundancy (if) by averaging on each species

# Treat factors separately
T.1 <- Traits[, c("Order", "Family", "Genus", "Best_guess_binomial", "Trophic_level", "Habitat_breadth", "Terrestriality")] %>%
  group_by(Best_guess_binomial) %>%
  mutate(Habitat_breadth=
           ifelse(length(unique(Habitat_breadth))!=1, Habitat_breadth[!is.na(Habitat_breadth)], NA)) %>%
  mutate(Terrestriality=
           ifelse(length(unique(Terrestriality))!=1, Terrestriality[!is.na(Terrestriality)], NA)) %>%
  distinct(Order,Family, Genus, Best_guess_binomial, Trophic_level, Habitat_breadth, Terrestriality)
 
# Average continuous traits when replicates
T.2 <- Traits %>% select(-c(Habitat_breadth, Terrestriality, Trophic_level, Order, Family, Genus))
T.2 <- setDT(T.2)[, lapply(.SD, mean, na.rm=TRUE), by = Best_guess_binomial]

Traits <- merge(T.1, T.2)
rm(T.1, T.2)

## Get Mammal species that have unknown binomial names & add to Traits ------
Traits <- .GetOther(Predicts, "Mammalia", Traits)

## till here?? ##~~~~~~~~~~~~

## Interpolations for known binomial names  ---------------------------------
# For increase in coverage
Before_corr <- nrow(Traits[!is.na(Traits$Generation_length_d),])

T <- Corr_Longevity(Myhrvold, Traits, "Mammalia")
Traits <- T$Traits; Mammals.Myhrvold <- T$Mammals.Myhrvold; rm(T)

After_corr <- nrow(Traits[!is.na(Traits$Generation_length_d),])
cat("Correlation on Longevity increased taxonomic coverage for this trait by",
    (After_corr-Before_corr)/Before_corr*100, "%")
rm(Before_corr, After_corr)

# Interpolations by averaging at lowest known taxon level -------------------

# Getting Averages per Taxonomic groups for continuous traits (Body_mass_g, Generation_length_d, Litter_size)
Myhrvold.Averages <- .Av_by_taxon(Myhrvold, c("Body_mass_g", "Generation_length_d", "Litter_size"), c("Order", "Family", "Genus"))
Pacifici.Averages <-  .Av_by_taxon(Pacifici, c("Body_mass_g", "Generation_length_d"), c("Order", "Family", "Genus"))
Pantheria.Averages <- .Av_by_taxon(Pantheria, 
                                  c("Body_mass_g", "Generation_length_d", "Litter_size", "Home_range_group_km2", "Home_range_ind_km2"),
                                  c("Order", "Family", "Genus"))

Order.Av <- .MergeAv(Myhrvold.Averages$Av_Order,  Pacifici.Averages$Av_Order, Pantheria.Averages$Av_Order)
Family.Av <- .MergeAv(Myhrvold.Averages$Av_Family,  Pacifici.Averages$Av_Family, Pantheria.Averages$Av_Family)
Genus.Av <- .MergeAv(Myhrvold.Averages$Av_Genus,  Pacifici.Averages$Av_Genus, Pantheria.Averages$Av_Genus)

# Getting Averages per Taxonomic groups for range sizes
Range.Av <- .Av_by_taxon(Mammal_range, c("Range_size_m2"), c("Order", "Family", "Genus"))

rm(Myhrvold.Averages, Pacifici.Averages,Pantheria.Averages)


# Adding averaged estimates to the main trait dataset where necessary: Averages at Genus/Family/Order level
Traits <- Add_Averages(Traits, Genus.Av, Family.Av, Order.Av,
                       c("Body_mass_g", "Generation_length_d", "Litter_size", 
                         "Home_range_group_km2", "Home_range_ind_km2"))

Traits <- Add_Averages(Traits, as.data.frame(Range.Av$Av_Genus), 
                       as.data.frame(Range.Av$Av_Family),
                       as.data.frame(Range.Av$Av_Order), "Range_size_m2")



# Plots -------------------------------------------------------------------


# Coverage per trait, before interpolation
Completeness <- apply(Traits[, 5:ncol(Traits)], 2,  function(y) sum(!is.na(y))) 
Completeness <- as.data.frame(Completeness/nrow(Traits)*100)
colnames(Completeness) <- "Completeness"
Completeness <- Completeness[order(Completeness, decreasing=FALSE), , drop=FALSE]

pdf(file="../Results/Plots/Mammals_cov_no_interpolation.pdf", width=5, height=4)
par(family='serif', mar = c(3,12,3,2), tcl=0.4, cex.axis=1, mgp=c(1,0.2,0), oma=c(1,1,1,1))
barplot(Completeness$Completeness, horiz = TRUE, 
        xlim = c(0,100), names.arg=rownames(Completeness), las=1, main="Coverage before interpolation \n Mammals")
abline(v=100, lty="dotted")
mtext(1, text="% Mammal taxa in PREDICTS", line=2)
mtext(2, text="All traits used", line=12)
dev.off()

# Coverage after interpolation




