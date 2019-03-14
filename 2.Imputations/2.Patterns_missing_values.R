## Patterns in missing trait values with regards to the phylogenetic position of the species

# plot per traits (percent and number of species with missing values within families for each trait)
# plot completeness (within family median)

library(phytools)
library(dplyr)
library(ape)
library(picante)
`%nin%` <- Negate(`%in%`)
source("Functions_patterns_missing_values.R")

## Traits
Traits <- c("Body_mass_g",
            "Longevity_d",
            "Litter_size", 
            "Range_size_m2", 
            "Habitat_breadth_IUCN",
            "Specialisation",
            "Trophic_level",
            "Diel_activity",
            "Primary_diet",
            "Diet_breadth")

TraitsReptiles <- c(Traits[1:8], "Adult_svl_cm", "Maturity_d")
TraitsAmphibians <- c(Traits, "Body_length_mm")
TraitsMammals <- c(Traits, "Adult_svl_cm", "Generation_length_d")
TraitsBirds <- c(Traits, "Generation_length_d")

## Load phylogenies
PhyloAmphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk") %>% .Format_tiplabels()
PhyloMammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk")  %>% .Format_tiplabels()
PhyloBirds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk")  %>% .Format_tiplabels() 
PhyloReptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk")  %>% .Format_tiplabels() 

## Load traits
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")

Traits_cont <-  c("Body_mass_g", "Longevity_d", "Litter_size", "Diet_breadth",
                  "Range_size_m2", "Habitat_breadth_IUCN")

Traits_cat <- c("Specialisation", "Diel_activity","Trophic_level", "Primary_diet")

## 1. Plot within family completeness patterns (within family median of completeness across species)
Amphibians_completeness_results <- Completeness_families(Amphibians, TraitsAmphibians)
Reptiles_completeness_results <- Completeness_families(Reptiles, TraitsReptiles)
Birds_completeness_results <- Completeness_families(Birds, TraitsBirds)
Mammals_completeness_results <- Completeness_families(Mammals, TraitsMammals)

# Amphibians
pdf(file="../../Results/Plots/Phylogeny_missing_values/Amphibians_completeness.pdf", width=5, height=7, pointsize=9)
p <- Plot_NA_patterns(PhyloAmphibians, Amphibians_completeness_results, Amphibians, TRUE, TRUE)
add.color.bar(prompt=FALSE, 75, cols=p$cols, title="",  digits=2, lims = c(0,100), x=45, y=10, subtitle="\nwithin family completeness (%)")
box("outer")
dev.off()

# Reptiles
pdf(file="../../Results/Plots/Phylogeny_missing_values/Reptiles_completeness.pdf", width=5, height=7, pointsize=9)
Plot_NA_patterns(PhyloReptiles, Reptiles_completeness_results, Reptiles, TRUE, TRUE)
add.color.bar(prompt=FALSE, 75, cols=p$cols, title="",  digits=2, lims = c(0,100), x=30, y=10, subtitle="\nwithin family completeness (%)")
box("outer")
dev.off()

# # Birds
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Birds_completeness.pdf", width=4, height=4, pointsize=9)
# Plot_NA_patterns(PhyloBirds, Birds_completeness_results, Birds, FALSE, TRUE)
# add.color.bar(prompt=FALSE, 20, cols=p$cols, title="",  digits=2, lims = c(0,100), x=10, y=3, subtitle="within family completeness")
# box("outer")
# dev.off()
# 
# with tip labels
pdf(file="../../Results/Plots/Phylogeny_missing_values/Birds_completeness_labels.pdf", width=10, height=20, pointsize=9)
Birds_completeness_results <- Completeness_families(Birds, TraitsBirds)
Plot_NA_patterns(PhyloBirds, Birds_completeness_results, Birds, TRUE, TRUE)
add.color.bar(prompt=FALSE, 50, cols=p$cols, title="",  digits=2, lims = c(0,100), x=10, y=3, subtitle="\nwithin family completeness (%)", fsize=2)
box("outer")
dev.off()
# 
# # # Mammals
# # pdf(file="../../Results/Plots/Phylogeny_missing_values/Mammals_completeness.pdf", width=4, height=4, pointsize=9)
# # Plot_NA_patterns(PhyloMammals, Mammals_completeness_results, Mammals, FALSE, TRUE)
# # add.color.bar(prompt=FALSE, 20, cols=p$cols, title="",  digits=2, lims = c(0,100), x=20, y=5, subtitle="within family completeness")
# # box("outer")
# # dev.off()
# 
# with tip labels
pdf(file="../../Results/Plots/Phylogeny_missing_values/Mammals_completeness_labels.pdf", width=10, height=15, pointsize=9)
Mammals_completeness_results <- Completeness_families(Mammals, TraitsMammals)
Plot_NA_patterns(PhyloMammals, Mammals_completeness_results, Mammals, TRUE, TRUE)
add.color.bar(prompt=FALSE, 50, cols=p$cols, title="",  digits=2, lims = c(0,100), x=20, y=5, subtitle="\nwithin family completeness (%)", fsize=1.5)
box("outer")
dev.off()

## All together without the tip labels
pdf(file="../../Results/Plots/Phylogeny_missing_values/Completeness_all.pdf", width=7, height=5, pointsize=12)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(6,1,2,1))
par(xpd=NA)
par(mfrow=c(2,2))
Plot_NA_patterns(PhyloMammals, Mammals_completeness_results, Mammals, FALSE, TRUE); title("A", line=0.5, adj=0, font=2)
Plot_NA_patterns(PhyloBirds, Birds_completeness_results, Birds, FALSE, TRUE); title("B", line=0.5, adj=0, font=2)
Plot_NA_patterns(PhyloReptiles, Reptiles_completeness_results, Reptiles, FALSE, TRUE); title("C", line=0.5, adj=0, font=2)
p <- Plot_NA_patterns(PhyloAmphibians, Amphibians_completeness_results, Amphibians, FALSE, TRUE); title("D", line=0.5, adj=0, font=2)

add.color.bar(prompt=FALSE, 150, cols=p$cols, title="",  digits=2, lims = c(0,100), x=-40, y=-20, subtitle="\nwithin family completeness (%)", fsize=1.3)
box("outer")
dev.off()




## 2. Within family trait coverage - both as the percentages of missing trait value within the value and the number of species (log10?)

## Amphibians
X <- c("Range_size_m2", "Body_length_mm", "Habitat_breadth_IUCN","Specialisation","Diel_activity",'Litter_size',
       "Primary_diet", "Trophic_level", "Diet_breadth", "Body_mass_g", "Longevity_d")
Titles <- c("A. RS", "B. BL", "C. HB", "D. Sp", "E. DA", "F. LCS", "G. PD", "H. TL", "I. DB", "J. BM", "K. Longevity")

pdf(file="../../Results/Plots/Phylogeny_missing_values/Amphibians_coverage.pdf", width=6, height=4, pointsize=9)
par(mfrow=c(3,4)); par(xpd=NA)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(1,1,2,20))
for (i in 1:11) {
  Results <- PercentNA_families(Amphibians, X[i])
  p <- Plot_NA_patterns(PhyloAmphibians, Results, Amphibians, FALSE, FALSE)
  title(main=Titles[i], adj=0, line=0.5)
}
Results <- Family_rep(Amphibians, TRUE)
p2 <- Plot_NA_patterns(PhyloAmphibians, Results, Amphibians, FALSE, FALSE)
title(main="L. % rep", adj=0, line=0.5)

add.color.bar(prompt=FALSE, 300, cols=p$cols, title="",  digits=1, lims = c(0,100), x=600, y=100, subtitle="\nwithin family missingness (%): \nplots A - K", fsize=1.5)
add.color.bar(prompt=FALSE, 300, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=600, y=50, subtitle="\nrepresentation (log-10): \nplot L", fsize=1.5)
box("outer")
dev.off()


## Mammals
X <- c("Body_mass_g", "Generation_length_d", "Range_size_m2","Trophic_level",'Diel_activity',
       "Primary_diet", "Habitat_breadth_IUCN", "Specialisation", "Diet_breadth", "Adult_svl_cm", "Litter_size", "Longevity_d")
Titles <- c("A. BM", "B. GL", "C. RS", "D. TL", "E. DA", "F. PD", "G. HB", "H. Sp", "I. DB", "J. BL", "K. LCS", "L. Longevity")

pdf(file="../../Results/Plots/Phylogeny_missing_values/Mammals_coverage.pdf", width=6, height=4, pointsize=9)
par(mfrow=c(3,5)); par(xpd=NA)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(1,1,2,1))
for (i in 1:12) {
  Results <- PercentNA_families(Mammals, X[i])
  p <- Plot_NA_patterns(PhyloMammals, Results, Mammals, FALSE, FALSE)
  title(main=Titles[i], adj=0, line=0.5)
}
Results <- Family_rep(Mammals, TRUE)
p2 <- Plot_NA_patterns(PhyloMammals, Results, Mammals, FALSE, FALSE)
title(main="L. % rep", adj=0, line=0.5)

add.color.bar(prompt=FALSE, 200, cols=p$cols, title="",  digits=1, lims = c(0,100), x=300, y=130, subtitle="\nwithin family missingness (%):\nplots A - K", fsize=1.5)
add.color.bar(prompt=FALSE, 200, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=300, y=50, subtitle="\nrepresentation (log-10):\nplot L", fsize=1.5)
box("outer")
dev.off()


## Birds
X <- c("Habitat_breadth_IUCN", "Specialisation", "Body_mass_g","Range_size_m2","Primary_diet",'Diel_activity',
       "Generation_length_d", "Trophic_level", "Diet_breadth", "Litter_size", "Longevity_d")
Titles <- c("A. HB", "B. Sp", "C. BM", "D. RS", "E. PD", "F. DA", "G. GL", "H. TL", "I. DB", "J. LS", "K. Longevity")

pdf(file="../../Results/Plots/Phylogeny_missing_values/Birds_coverage.pdf", width=6, height=4, pointsize=9)
par(mfrow=c(3,4)); par(xpd=NA)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(1,1,2,20))
for (i in 1:11) {
  Results <- PercentNA_families(Birds, X[i])
  p <- Plot_NA_patterns(PhyloBirds, Results, Birds, FALSE, FALSE)
  title(main=Titles[i], adj=0, line=0.5)
}
Results <- Family_rep(Birds, TRUE)
p2 <- Plot_NA_patterns(PhyloBirds, Results, Birds, FALSE, FALSE)
title(main="L. % rep", adj=0, line=0.5)

add.color.bar(prompt=FALSE, 100, cols=p$cols, title="",  digits=1, lims = c(0,100), x=200, y=500, subtitle="\n within family missingness (%): \n plots A - K", fsize=1.5)
add.color.bar(prompt=FALSE, 100, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=200, y=200, subtitle="\n representation (log-10): \n plot L", fsize=1.5)
box("outer")
dev.off()


## Reptiles
X <- c("Body_mass_g", "Range_size_m2", "Litter_size","Habitat_breadth_IUCN","Specialisation",'Diel_activity',
       "Longevity_d", "Trophic_level", "Adult_svl_cm", "Maturity_d")
Titles <- c("A. BM", "B. RS", "C. LCS", "D. HB", "E. Sp", "F. DA", "G. Longevity", "H. TL", "J. BL", "K. SM")

pdf(file="../../Results/Plots/Phylogeny_missing_values/Reptiles_coverage.pdf", width=6, height=4, pointsize=9)
par(mfrow=c(3,4)); par(xpd=NA)
par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(5,1,5,1), oma=c(3,1,2,5))
for (i in 1:10) {
  Results <- PercentNA_families(Reptiles, X[i])
  p <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, FALSE, FALSE)
  title(main=Titles[i], adj=0, line=0.5)
}
Results <- Family_rep(Reptiles, TRUE)
p2 <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, FALSE, FALSE)
title(main="L. % rep", adj=0, line=0.5)

add.color.bar(prompt=FALSE, 150, cols=p$cols, title="",  digits=1, lims = c(0,100), x=330, y=60, subtitle="\nwithin family missingness (%): \nplots A - K", fsize=1.5)
add.color.bar(prompt=FALSE, 150, cols=p2$cols, title="",  digits=1, lims = range(Results$Percent), x=330, y=20, subtitle="\nrepresentation (log-10): \nplot L", fsize=1.5)
box("outer")
dev.off()







## Amphibians

# split.screen(c(1,2))
# screen(1)
# split.screen(c(3,4))
# screen(3)
# par(oma=c(4,0,2,4), las=1)
# for (i in 1:11) {
#   screen(i+2)
#   Results <- PercentNA_families(Amphibians, X[i])
#   p <- Plot_NA_patterns(PhyloAmphibians, Results, Amphibians, FALSE, FALSE)
#   mtext(text = Titles[i], line =-12, at=100, font = 2)
# }
# 
# screen(2)
# par(oma=c(4,2,1.5,0), las=1)
# Results <- Family_rep(Amphibians)
# p2 <- Plot_NA_patterns(PhyloAmphibians, Results, Amphibians, TRUE, FALSE)
# 
# par(xpd=NA)
# add.color.bar(prompt=TRUE, 400, cols=p$cols, title="",  digits=2, lims = c(0,100), x=20, y=5, subtitle="within family coverage (%)")
# add.color.bar(prompt=TRUE, 500, cols=p2$cols, title="",  digits=0, lims = range(Results$Percent), x=20, y=5, subtitle="% representation")
# 
# close.screen()
# box("outer")
# mtext("B", side=3, font=2, cex=1.5)
# mtext("A", side=3, font=2, cex=1.5, at=-650)
# 
# dev.off()


## Reptiles

# X <- c("Body_mass_g", "Range_size_m2", "Litter_size","Habitat_breadth_IUCN","Specialisation",'Diel_activity',
#        "Longevity_d", "Trophic_level", "Adult_svl_cm", "Maturity_d")
# Titles <- c("BM", "RS", "LCS", "HB", "Sp", "DA", "L", "TL", "BL", "M")
# 
# pdf(file="../../Results/Plots/Phylogeny_missing_values/Reptiles_coverage.pdf", width=7, height=8, pointsize=9)
# 
# split.screen(c(1,2))
# screen(1)
# split.screen(c(3,4))
# screen(3)
# par(oma=c(4,0,2,4), las=1)
# for (i in 1:11) {
#   screen(i+2)
#   Results <- PercentNA_families(Reptiles, X[i])
#   p <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, FALSE, FALSE)
#   mtext(text = Titles[i], line =-12, at=80, font = 2)
# }
# 
# screen(2)
# par(oma=c(4,2,1.5,0), las=1)
# Results <- Family_rep(Reptiles)
# p2 <- Plot_NA_patterns(PhyloReptiles, Results, Reptiles, TRUE, FALSE)
# 
# par(xpd=NA)
# add.color.bar(prompt=FALSE, 300, cols=p$cols, title="",  digits=2, lims = c(0,100), x=-500, y=-5, subtitle="within family missing values (%)")
# add.color.bar(prompt=FALSE, 300, cols=p2$cols, title="",  digits=0, lims = range(Results$Percent), x=150, y=-5, subtitle="% representation")
# 
# close.screen()
# box("outer")
# mtext("B", side=3, font=2, cex=1.5)
# mtext("A", side=3, font=2, cex=1.5, at=-400)
# 
# dev.off()






