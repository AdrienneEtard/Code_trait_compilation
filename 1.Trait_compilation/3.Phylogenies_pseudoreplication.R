## Treating pseudo-replication in phylogenetic tips

# TODO compute.brlen before plotting?


# # # # P r e a m b l e 

X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ape", "treeio", "ngram","phylobase")
lapply(X, library, character.only=TRUE); rm(X)
source("Functions_for_phylogenies.R")

## Load data
# No taxonomic correction
UN_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Amphibians.csv")
UN_Birds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Birds.csv")
UN_Mammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Mammals.csv")
UN_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.standardised/Reptiles.csv")

# With taxonomic correction
C_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Amphibians.csv")
C_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Birds.csv")
C_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Mammals.csv")
C_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.standardised/Reptiles.csv")

## Load PREDICTS with and without taxonomic correction
UN_Predicts <- readRDS("../../Data/PREDICTS_database.rds") %>% 
  filter(Class %in% c("Aves", "Mammalia", "Reptilia", "Amphibia"))
C_Predicts <- readRDS("../../Results/0.Data_resolved_taxonomy/Processed_datasets/PredictsVertebrates.rds")

## Load phylogenies

# # Either these files, to add species to genus
# # uncorrected
# PhyloMammal_UN <- read.newick("../../Data/Mammals/Phylogenies/TTOL_mammals_smoothed_interpolated.nwk")
# PhyloAmphibian_UN <- read.newick("../../Data/Phylogenies/TTOL_amphibians_unsmoothed_Hedges2015.nwk")
# PhyloBird_UN <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk")
# PhyloReptile_UN <- read.newick("../../Data/Phylogenies/TTOL_squamates_unsmoothed_Hedges2015.nwk")
# 
# # corrected
# Phylo_Mammals_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloMammals.nwk")
# Phylo_Amphibians_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloAmphibians.nwk")
# Phylo_Reptiles_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloReptiles.nwk")
# Phylo_Birds_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")
# 
# # ## If problem when opening bird phylogeny -- run these lines
# # source("../0.Prerequisites/Resolve_taxonomy_functions.R")
# # PhyloB <- read.newick("../../Data/Phylogenies/TTOL_birds_smoothed_interpolated_Hedges2015.nwk") %>% .Format_tiplabels()
# # Syn_birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv")
# # Phylo_Birds_C <- Replace_by_accepted_name(Syn_birds, PhyloB)
# # write.tree(Phylo_Birds_C, "../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")
# # Phylo_Birds_C <- read.newick("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Phylogenies/PhyloBirds.nwk")
# # rm(PhyloB, Syn_birds)


## Or these files, after species have been attached 

## Load data (from above commented script: phylogenies with ramdom addition of species at the genus level)
# corrected
Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Mammals.nwk") %>% .Format_tiplabels()
Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Birds.nwk") %>% .Format_tiplabels()
Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Amphibians.nwk") %>% .Format_tiplabels()
Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/1.Random_additions/Reptiles.nwk") %>% .Format_tiplabels()

# uncorrected
PhyloMammal_UN <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Mammals.nwk") %>% .Format_tiplabels()
PhyloAmphibian_UN <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Amphibians.nwk") %>% .Format_tiplabels()
PhyloBird_UN <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Birds.nwk") %>% .Format_tiplabels()
PhyloReptile_UN <- read.newick("../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Reptiles.nwk") %>% .Format_tiplabels()



# # # # #

# ## 1. Randomly attach species not present in the phylogenies to their genuses when applicable, for both corrected and uncorrected datasets
# 
# ## Corrected datasets
# C.ResMammals <- SpeciesPhylo(C_Predicts, C_Mammals, Phylo_Mammals_C, "Mammalia", "random")
# C.ResAmphibians <- SpeciesPhylo(C_Predicts, C_Amphibians, Phylo_Amphibians_C, "Amphibia","random")
# C.ResReptiles <- SpeciesPhylo(C_Predicts, C_Reptiles, Phylo_Reptiles_C, "Reptilia","random")
# C.ResBirds <- SpeciesPhylo(C_Predicts, C_Birds, Phylo_Birds_C, "Aves","random")
# 
# ## Uncorrected datasets
# UN.ResMammals <- SpeciesPhylo(UN_Predicts, UN_Mammals, PhyloMammal_UN, "Mammalia", "random")
# UN.ResAmphibians <- SpeciesPhylo(UN_Predicts,UN_Amphibians, PhyloAmphibian_UN, "Amphibia","random")
# UN.ResReptiles <- SpeciesPhylo(UN_Predicts, UN_Reptiles, PhyloReptile_UN, "Reptilia","random")
# UN.ResBirds <- SpeciesPhylo(UN_Predicts, UN_Birds, PhyloBird_UN, "Aves","random")
# 
# ## Save these phylogenies with attached genus
# C.Phylo_Mammals <- C.ResMammals$Phylo %>% .Format_tiplabels()
# C.Phylo_Amphibians <- C.ResAmphibians$Phylo %>% .Format_tiplabels()
# C.Phylo_Reptiles <- C.ResReptiles$Phylo %>% .Format_tiplabels()
# C.Phylo_Birds <- C.ResBirds$Phylo %>% .Format_tiplabels()
# 
# UN.Phylo_Mammals <- UN.ResMammals$Phylo %>% .Format_tiplabels()
# UN.Phylo_Amphibians <- UN.ResAmphibians$Phylo %>% .Format_tiplabels()
# UN.Phylo_Reptiles <- UN.ResReptiles$Phylo %>% .Format_tiplabels()
# UN.Phylo_Birds <- UN.ResBirds$Phylo %>% .Format_tiplabels()
# 
# write.tree(C.Phylo_Amphibians, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Amphibians.nwk")
# write.tree(C.Phylo_Reptiles, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Reptiles.nwk")
# write.tree(C.Phylo_Birds, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Birds.nwk")
# write.tree(C.Phylo_Mammals, "../../Results/1.Phylogenies/Corrected/1.Random_additions/Mammals.nwk")
# 
# write.tree(UN.Phylo_Amphibians, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Amphibians.nwk")
# write.tree(UN.Phylo_Reptiles, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Reptiles.nwk")
# write.tree(UN.Phylo_Birds, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Birds.nwk")
# write.tree(UN.Phylo_Mammals, "../../Results/1.Phylogenies/Uncorrected/1.Random_additions/Mammals.nwk")


## 2. Looking at replicated tips, and dropping replicated tips -- JUST FOR CORRECTED PHYLOGENIES. 


## 2.1. Drop tips 

# Mammals
Drop_Mammals <- DropTips(Phylo_Mammals, PhyloMammal_UN)
Phylo_Mammals.final <- Drop_Mammals$PhylogenyCor
# DroppedMammals <- Drop_Mammals$Replicated
# MTC <- Drop_Mammals$PbSpecies
# table(DroppedMammals$Count)


# Birds
Drop_Birds <- DropTips(Phylo_Birds, PhyloBird_UN)
Phylo_Birds.final <- Drop_Birds$PhylogenyCor
# DroppedBirds <- Drop_Birds$Replicated
# BTC <- Drop_Birds$PbSpecies
# table(DroppedBirds$Count)
# X <- Phylo_Birds.final$tip.label %>% as.data.frame() %>% setNames(., "Rep") %>% group_by(Rep) %>%
#   summarise(Count=n())

# Reptiles
Drop_Reptiles <- DropTips(Phylo_Reptiles, PhyloReptile_UN)
Phylo_Reptiles.final <- Drop_Reptiles$PhylogenyCor
DroppedReptiles <- Drop_Reptiles$Replicated
RTC <- Drop_Reptiles$PbSpecies
table(DroppedReptiles$Count)


# Amphibians
Drop_Amphibians <- DropTips(Phylo_Amphibians, PhyloAmphibian_UN)
Phylo_Amphibians.final <- Drop_Amphibians$PhylogenyCor
DroppedAmphibians <- Drop_Amphibians$Replicated
ATC <-Drop_Amphibians$PbSpecies
table(DroppedAmphibians$Count)


# ## Replicated species that are also in PREDICTS
# intersect(Rep_in_Predicts(C_Predicts, "Aves", DroppedBirds), BTC$Reps)
# intersect(Rep_in_Predicts(C_Predicts, "Mammals", DroppedMammals), MTC$Reps)
# intersect(Rep_in_Predicts(C_Predicts, "Amphibia", DroppedAmphibians), ATC$Reps)
# intersect(Rep_in_Predicts(C_Predicts, "Reptilia", DroppedReptiles), RTC$Reps)


## Save phylogenies with additions and dropped replicated tips
write.tree(Phylo_Amphibians.final, "../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk")
write.tree(Phylo_Reptiles.final, "../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk")
write.tree(Phylo_Birds.final, "../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk")
write.tree(Phylo_Mammals.final, "../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk")



# ## 2.3. Case studies (examples) 
# 
# ## 1: Rhinolophus pusillus - dropped tips will be the ones furthest from the original position
# # for which there are 8 identified synonyms, and 5 replicates
# # PhyloMammal_UN$tip.label <- paste(paste(substr(PhyloMammal_UN$tip.label,1,1), ".", sep=""), word(PhyloMammal_UN$))
# Reps <- c(PhyloMammal_UN$tip.label[grepl("Rhinolophus imaizumii|Rhinolophus monoceros|Rhinolophus blythi|Rhinolophus perditus|Rhinolophus minor|Rhinolophus gracilis|Rhinolophus pumilus|Rhinolophus cornutus",
#                                        PhyloMammal_UN$tip.label)], "Rhinolophus pusillus")
# 
# Reps <- paste(paste(substr(Reps, 1,1), ".", sep=""), word(Reps, 2), sep=" ")
# 
# par(mar=c(0,0,1,0), oma=c(0,1,1,1), family="sans")
# par(mfrow=c(1,3))
# 
# # For the uncorrected tree
# Rhinolophus <- drop.tip(PhyloMammal_UN, PhyloMammal_UN$tip.label[!grepl("Rhinolophus", PhyloMammal_UN$tip.label)])
# Rhinolophus$tip.label <- paste(paste(substr(Rhinolophus$tip.label, 1,1), ".", sep=""), word(Rhinolophus$tip.label, 2), sep=" ")
# plot(Rhinolophus,tip.color=ifelse(Rhinolophus$tip.label %in% Reps, ifelse(Rhinolophus$tip.label=="R. pusillus", "blue", "red"),"black"), 
#      cex=1, label.offset=1,
#      font=ifelse(Rhinolophus$tip.label %in% Reps, 4,3), main="A. Uncorrected",  node.depth=2) 
# 
# # For the corrected tree
# Rhinolophus <- drop.tip(Phylo_Mammals, Phylo_Mammals$tip.label[!grepl("Rhinolophus", Phylo_Mammals$tip.label)])
# Rhinolophus$tip.label <- paste(paste(substr(Rhinolophus$tip.label, 1,1), ".", sep=""), word(Rhinolophus$tip.label, 2), sep=" ")
# plot(Rhinolophus,tip.color=ifelse(Rhinolophus$tip.label %in% Reps, ifelse(Rhinolophus$tip.label=="R. pusillus", "blue", "red"),"black"),
#      cex=1, adj=0, label.offset=1,
#      font=ifelse(Rhinolophus$tip.label %in% "R. pusillus", 4,3), main="B. Corrected", node.depth=2) 
# 
# # For the corrected tree after reducing replication
# Rhinolophus <- drop.tip(Phylo_Mammals.final, Phylo_Mammals.final$tip.label[!grepl("Rhinolophus", Phylo_Mammals.final$tip.label)])
# Rhinolophus$tip.label <- paste(paste(substr(Rhinolophus$tip.label, 1,1), ".", sep=""), word(Rhinolophus$tip.label, 2), sep=" ")
# plot(Rhinolophus,tip.color=ifelse(Rhinolophus$tip.label %in% Reps, ifelse(Rhinolophus$tip.label=="R. pusillus", "blue", "red"),"black"),
#      cex=1, adj=0, label.offset=1,
#      font=ifelse(Rhinolophus$tip.label %in% "R. pusillus", 4,3), main="C. After removing replicated tips", node.depth=2) 
# 
# box("outer")
# 
# ## 2: Case study where the dropped tip is chosen randomly
# # Uperodon taprobanicus (Kaloula pulchra, Kaloula taprobanica) 
# Reps <- c(PhyloAmphibian_UN$tip.label[grepl("Kaloula pulchra|Kaloula taprobanica", PhyloAmphibian_UN$tip.label)],"Uperodon taprobanicus")
# Reps <- paste(paste(substr(Reps, 1,2), ".", sep=""), word(Reps, 2), sep=" ")
# 
# par(mar=c(0,0,1,0), oma=c(0,1,1,1), family="sans")
# par(mfrow=c(1,3))
# # For the uncorrected tree
# Uperodon <- drop.tip(PhyloAmphibian_UN, PhyloAmphibian_UN$tip.label[!grepl("Uperodon|Kaloula|Ramanella|Metaphrynella|Rhombophryne", PhyloAmphibian_UN$tip.label)])
# Uperodon$tip.label <- paste(paste(substr(Uperodon$tip.label, 1,2), ".", sep=""), word(Uperodon$tip.label, 2), sep=" ")
# plot(Uperodon,tip.color=ifelse(Uperodon$tip.label %in% Reps, "red","black"), cex=1, label.offset=1,
#      font=ifelse(Uperodon$tip.label %in% Reps, 4,3), main="A. Uncorrected",  node.depth=2) 
# 
# # For the corrected tree
# Uperodon <- drop.tip(Phylo_Amphibians, Phylo_Amphibians$tip.label[!grepl("Uperodon|Ramanella|Kaloula|Metaphrynella|Rhombophryne", Phylo_Amphibians$tip.label)])
# Uperodon$tip.label <- paste(paste(substr(Uperodon$tip.label, 1,2), ".", sep=""), word(Uperodon$tip.label, 2), sep=" ")
# plot(Uperodon,tip.color=ifelse(Uperodon$tip.label %in% "Up. taprobanicus", "blue","black"), cex=1, adj=0, label.offset=1,
#      font=ifelse(Uperodon$tip.label %in% "Up. taprobanicus", 4,3), main="B. Corrected", node.depth=2) 
# 
# # For the corrected tree after removing replicated tips
# Uperodon <- drop.tip(Phylo_Amphibians.final, Phylo_Amphibians.final$tip.label[!grepl("Uperodon|Ramanella|Kaloula|Metaphrynella|Rhombophryne", Phylo_Amphibians.final$tip.label)])
# Uperodon$tip.label <- paste(paste(substr(Uperodon$tip.label, 1,2), ".", sep=""), word(Uperodon$tip.label, 2), sep=" ")
# 
# plot(Uperodon,tip.color=ifelse(Uperodon$tip.label %in% "Up. taprobanicus", "blue","black"), cex=1, adj=0, label.offset=1,
#      font=ifelse(Uperodon$tip.label %in% "Up. taprobanicus", 4,3), main="C. After removing replicated tips", node.depth=2) 
# box("outer")





# # Plot distance range among replicated tips
# pdf(file="../../Results/Plots/Phylogenies/DropTipsDistances.pdf", width=10, height=7, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0))
# par(mfrow=c(2,2))
# plot(density(DroppedMammals$RangePosition), main="Mammals")
# abline(v=1, col="blue")
# plot(density(DroppedBirds$RangePosition), main="Birds")
# abline(v=1, col="blue")
# plot(density(DroppedReptiles$RangePosition), main="Reptiles")
# abline(v=1, col="blue")
# plot(density(DroppedAmphibians$RangePosition), main="Amphibians")
# abline(v=1, col="blue")
# dev.off()
# 
# ## Save phylogenies
# write.tree(Phylo_Amphibians, "../../Results/1.Phylogenies/PhyloAmphibians.nwk")
# write.tree(Phylo_Reptiles, "../../Results/1.Phylogenies/PhyloReptiles.nwk")
# write.tree(Phylo_Birds, "../../Results/1.Phylogenies/PhyloBirds.nwk")
# write.tree(Phylo_Mammals, "../../Results/1.Phylogenies/PhyloMammals.nwk")
# 
# 
# ## Plotting by family/order: Collapsing branches
# Phylo_order_family <- function (Phylo, TraitDF, Family_or_Order, Collapse) {
#   
#   # browser()
#   
#   TraitDF <- TraitDF %>% dplyr::select(Best_guess_binomial, Genus, Family, Order)
#   
#   # 1. Intersect species in the phylogeny with species in the trait dataset, and prune
#   row.names(TraitDF) <- TraitDF$Best_guess_binomial
#   Match <- match.phylo.data(Phylo, TraitDF)
#   Phylo_pruned <- Match$phy
#   
#   # 2. Group tree by family or order with ggtree::groupOTU
#   Groups <- subset(TraitDF, Best_guess_binomial %in% Phylo_pruned$tip.label)
#   Groups <- Groups[, c("Best_guess_binomial", Family_or_Order)]
#   Groups <- Groups %>% setNames(., c("Species", "TaxGroup"))
#   
#   List <- list()
#   
#   for (i in unique(Groups$TaxGroup)) {
#     List[[i]] <- as.character(Groups$Species[Groups$TaxGroup==i])
#   }
#   
#   if (Collapse) {
#     
#     browser()
#     
#     Nodes <- lapply(List, function(x) MRCA(Phylo_pruned, x))
#     # Nodes that are not null
#     Nodes <- Filter(Negate(is.null), Nodes)
#     
#     
#     Tree <- ggtree::groupOTU(Phylo_pruned, List)
#     Plot <- ggtree(Tree, aes(color=group)) + theme(legend.position="right")
#     for (i in Nodes){
#       Plot <- collapse(Plot, i)
#     }
#     
#     # Plot <- expand(Plot, Nodes[["ANSERIFORMES"]])
#     # Plot <- expand(Plot, Nodes[["PELECANIFORMES"]])
#     # Plot <- expand(Plot, Nodes[["STRUTHIONIFORMES"]])
#     
#     
#   } else {
#     
#     Tree <- ggtree::groupOTU(Phylo_pruned, List)
#     Plot <- ggtree(Tree, aes(color=group), layout = "circular") + theme(legend.position="right") 
#   }
#   
#   return(Plot)
# }
# 
# PlotAmphibians <- Phylo_order_family(Phylo_Amphibians, Amphibians, "Order", FALSE)
# Phylo_order_family(Phylo_Birds, Birds, "Order", FALSE)
# Phylo_order_family(Phylo_Mammals, Mammals, "Order", FALSE)
# # Reptiles: only Squamata... plot at family level
# Phylo_order_family(Phylo_Reptiles, Reptiles, "Family", FALSE)
# 
# 
# ## Collapsing branches and plotting tree by family or order
# Phylo_order_family(Phylo_Amphibians, Amphibians, "Order", TRUE)
# Phylo_order_family(Phylo_Birds, Birds, "Family", TRUE)
# Phylo_order_family(Phylo_Reptiles, Reptiles, "Family", TRUE)
# Phylo_order_family(Phylo_Mammals, Mammals, "Order", TRUE)
# 
