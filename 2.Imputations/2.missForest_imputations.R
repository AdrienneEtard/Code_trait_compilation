# Imputations using missForest, with  phylogenetic eigenvectors as predictors -- imputation errors in another script.

# still to do:
# compare results from missForest imputations with results from phylogenetic imputations - at least for continuous traits?
# error of imputations / robustness of imputations.
# compare results with Rob Cooke's.

## Code is parallelised

# https://www.jottr.org/2018/06/23/future.apply_1.0.0/

library(parallel)

# # ## START CLUSTER
# Cluster <- makeCluster(detectCores())
# 
# ## EXCECUTE ANY PRE PROCESSING CODE NECESSARY
# clusterEvalQ(Cluster, {
#   library(dplyr)
#   library(phytools)
#   library(missForest)
#   library(pbmcapply)
#   library(pbapply)
# })

## Preamble
X <- c("dplyr", "phytools", "missForest", "pbapply")
lapply(X, library, character.only=TRUE); rm(X)
`%nin%` <- Negate(`%in%`)

source("Functions_for_missForest_imputations.R")

## Load phylogenies; data obtained from www.biodiversitycenter.org/ttol; downloaded 06 July 2018 and then processed (resolved taxonomy, drop tips, etc)
# Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk") %>% .Format_tiplabels() %>% compute.brlen()
# Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk") %>% .Format_tiplabels() %>% compute.brlen() 
# Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk") %>% .Format_tiplabels() %>% compute.brlen()
# Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk") %>% .Format_tiplabels()  %>% compute.brlen()

## Load trait data, transformed, standardised, and with phylogenetic imformation
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv")
Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv")
Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Amphibians.csv")
Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Reptiles.csv")

## Define variables for imputations

Habitat <- c("Forest","Savanna","Shrubland","Grassland","Wetland","Rocky.areas","Caves.and.subterranean",
             "Desert","Marine","Marine.intertidal.or.coastal.supratidal",
             "Artificial","Introduced.vegetation","Other.Unknown")

Diet <- c("IN", "VE", "PL", "SE", "NE", "FR", "SCV")

Taxinfo <- "Order"

Traits_cont <-  c("log10_Body_mass_g", "log10_Longevity_d", "log10_Litter_size", 
                  "Range_size_m2", "sqrt_Habitat_breadth_IUCN")

Traits_cat <- c(Habitat, "Specialisation",        # Include diet now or derive later from imputed primary diet?
               "Diel_activity","Trophic_level", Diet)
               
#EV <- c(); for(i in 1:10) {EV <- c(EV,paste("EV_",i, sep=""))}

## Add some traits for each taxon # TODO check trait names
MammalsCont <- c(Traits_cont, "log10_Generation_length_d", "log10_Adult_svl_cm")
BirdsCont <- Traits_cont
ReptilesCont <- c(Traits_cont, "log10_Adult_svl_cm", "log10_Maturity_d")
AmphibiansCont <- c(Traits_cont, "log10_Body_length_mm")


## Function arguments as lists
#PhyloList <- list(M=Phylo_Mammals, B=Phylo_Birds, R=Phylo_Reptiles, A=Phylo_Amphibians)
Cont.TraitsList <- list(M=MammalsCont, B=BirdsCont, R=ReptilesCont, A=AmphibiansCont)
Cat.TraitsList <- list(M=Traits_cat, B=Traits_cat, R=Traits_cat[Traits_cat %nin% Diet], A=Traits_cat)
DF.TraitsList <- list(M=Mammals, B=Birds, R=Reptiles, A=Amphibians)
Taxinfo.List <- list(M="Order", B="Order", R="Order", A=c("Order", "Family"))
DietTRUE.List <- list(M=TRUE, B=TRUE, R=FALSE, A=TRUE)

# ## Export variables in all clusters
# clusterExport(cl=Cluster, list("Imputations_missForest", 
#                                "DF.TraitsList",
#                                "PhyloList",
#                                "Cont.TraitsList", 
#                                "Cat.TraitsList", 
#                                "Taxinfo", 
#                                "EV"), envir=environment())

## Parallelised calculations

## which function to use here? on windows? to parallelise mapply.....


# mcmapply only works with forking (Ubuntu) -- or pbmcmapply, only with forking

Imputed.Vertebrate_Traits <- pbmapply (FUN=Imputations_missForest,
                                       DF.TraitsList,
                                      # Phylo=PhyloList,
                                       Taxinfo=Taxinfo,
                                       Traits_cont=Cont.TraitsList,
                                       Traits_cat=Cat.TraitsList,
                                       EV="EV_1",
                                       ErrorTrue=TRUE,
                                       DietTRUE=DietTRUE.List)

saveRDS(Imputed.Vertebrate_Traits, "../../Results/2.imputed_trait_datasets/imputed_datasets/I1.rds")

## To add: prinary diet and taxinfo

## Save results

# ## DESTROY CLUSTER
# stopCluster(Cluster)


# Test <- Imputations_missForest(Mammals, "Order", MammalsCont, Traits_cat, EV = EV, ErrorTrue = TRUE, DietTRUE=TRUE)
# TestMammals <- Test$Imputed.Values
# Test$Imputations.Errors
# glimpse(TestMammals)

Test <- Imputations_missForest(Reptiles[c(1:50),], "Order", ReptilesCont, Traits_cat, EV = EV, ErrorTrue = TRUE, DietTRUE=FALSE)

Test <- Imputations_missForest(Birds[c(1:50),], "Order", BirdsCont, Traits_cat, EV = EV, ErrorTrue = TRUE, DietTRUE=TRUE)

Test <- Imputations_missForest(Amphibians[c(1:50),], "Order", AmphibiansCont, Traits_cat, EV = EV, ErrorTrue = TRUE, DietTRUE=TRUE)


















# # 4. Birds with phylogeny
# 
# # Birds with phylogeny
# Tr_cont <- c("Body_mass_g", "Adult_svl_cm","Generation_length_d","Max_longevity_d", "Longevity_d","Maturity_d","Litter_size", "Range_size_m2", "Habitat_breadth_IUCN")
# Tr_cat <- c("Trophic_level", "Specialisation", "Diet_breadth")
# 
# Birds <- Imputations_missForest(Phylo_Birds, Traits_Birds, Tr_cont, Tr_cat, 10, FALSE)
# Birds_imputed <- Birds$ImputedValues
# Birds_error <- Birds$ImputedErrors
# 
# # Other birds (without phylogeny)
# Birds2 <- Imputations_missForest_other(Traits_Birds, Birds_imputed, Tr_cont, Tr_cat, FALSE)
# Birds_imputed2 <- Birds2$ImputedValues
# Birds_error2 <- Birds2$ImputedErrors
# Birds_imputed2$Best_guess_binomial <- row.names(Birds_imputed2)
# Birds_imputed2 <- Birds_imputed2[, c(13,1:12)]


## Very rough imputations to use to prepare further analyses - temporary trait datasets.
# Traits <- c("Body_mass_g", "Generation_length_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation", "Diet_breadth")
# Traits2 <- c("Body_mass_g", "Maturity_d", "Litter_size", "Range_size_m2", "Habitat_breadth_IUCN", "Trophic_level", "Specialisation", "Diet_breadth")
# 
# rownames(TraitsMammal) <- TraitsMammal$Best_guess_binomial
# TraitsMammal <- missForest(TraitsMammal[,Traits])
# TraitsMammal <- TraitsMammal$ximp
# 
# rownames(TraitsAmphibian) <- TraitsAmphibian$Best_guess_binomial
# TraitsAmphibian <- missForest(TraitsAmphibian[,Traits2])
# TraitsAmphibian <- TraitsAmphibian$ximp
# 
# rownames(TraitsReptile) <- TraitsReptile$Best_guess_binomial
# TraitsReptile <- missForest(TraitsReptile[,Traits2[Traits2!="Diet_breadth"]])
# TraitsReptile <- TraitsReptile$ximp
# 
# rownames(TraitsBird) <- TraitsBird$Best_guess_binomial
# TraitsBird <- missForest(TraitsBird[,Traits])
# TraitsBird <- TraitsBird$ximp
# 
# write.csv(TraitsMammal, "../../Results/2.missForest_imputations_preliminary_results/Mammals.csv")
# write.csv(TraitsAmphibian, "../../Results/2.missForest_imputations_preliminary_results/Amphibians.csv")
# write.csv(TraitsReptile, "../../Results/2.missForest_imputations_preliminary_results/Reptiles.csv")
# write.csv(TraitsBird, "../../Results/2.missForest_imputations_preliminary_results/Birds.csv")


## Influence of number of eigenvectors on the NRMSE - selection of N eigenvectors
# MammalsPredictsEigenvectors <- Phylo_eigenvectors(Mammal_phylo, TraitsMammal.Predicts)
# Seq <- seq(0, ncol(MammalsPredictsEigenvectors), by=10)[2:length(Seq)]
# Seq <- c(1, Seq)
# 
# 
# # Without rep
# Test1 <- Imputations_missForest_Error(Mammal_phylo, MammalsPredictsEigenvectors, TraitsMammal.Predicts, Tr_cont, Tr_cat, FALSE, Seq)
# # With rep (calculation of mean and SE)
# Test2 <- Imputations_missForest_ErrorSE(Mammal_phylo, MammalsPredictsEigenvectors, TraitsMammal.Predicts, Tr_cont, Tr_cat, FALSE, Seq, 10)
# 
# p1 <- ggplot(data=Test1, aes(Seq, NRMSE))
# p1 <- p1 + geom_point(size =3) + geom_line()
# p1
# Test1$Seq[Test1$NRMSE==min(Test1$NRMSE)]
# 
# p2 <- ggplot(data=Test2, aes(Seq,Mean_NRMSE)) + 
#   geom_point(size =3) + geom_line() +
#   geom_errorbar(aes(ymax = Mean_NRMSE+SE_NRMSE, ymin = Mean_NRMSE-SE_NRMSE), size =1, width=0.3)
# p2
# Test2$Seq[Test2$Mean_NRMSE==min(Test2$Mean_NRMSE)]
# 
# 
# ## Conclusion: 10 eigenvectors are enough (?)


## Plottings
# 
# ToPlot <- function(OriginalDF, ImputedDF, Trait) {
#   ToPlot1 <- OriginalDF[, colnames(OriginalDF) %in% Trait] %>% as.data.frame()
#   ToPlot1$Which <- "Original"
#   ToPlot2 <- ImputedDF[,  colnames(ImputedDF) %in% Trait]  %>% as.data.frame()
#   ToPlot2$Which <- "Imputed"
#   ToPlot <- rbind(ToPlot1, ToPlot2)
#   ToPlot$. <- log(ToPlot$.)
#   return(ToPlot)
# }
# 
# Plot <- function(ToPlot, Trait) {
#   p <- ggplot(ToPlot, aes(x=., fill=Which), family="serif") +
#   geom_density(alpha=0.6) + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#   ylab("Density") + xlab(Trait) +
#   theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
#   guides(fill=guide_legend(title="Data"))
#   return(p)
# }
# 
# 
# ## Barplots for categorical variables
# ToPlotCat <- function(OriginalDF, ImputedDF, Trait) {
#   ToPlot1 <- OriginalDF[, colnames(OriginalDF) %in% Trait] %>% as.data.frame()
#   ToPlot1$Which <- "Original"
#   ToPlot2 <- ImputedDF[,  colnames(ImputedDF) %in% Trait]  %>% as.data.frame()
#   ToPlot2$Which <- "Imputed"
#   ToPlot <- rbind(ToPlot1, ToPlot2)
#   return(ToPlot)
# }
# 
# PlotCat <- function(ToPlot, Trait) {
#   
#   ToPlot <- ToPlot[!is.na(ToPlot$.), ]
#   
#   p <- ggplot(ToPlot, aes(x=., fill=Which), family="serif") +
#     geom_bar(alpha=0.6, color="black", stat="count", position=position_dodge()) +
#     scale_y_continuous(labels = scales::percent) +
#     theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#     ylab("Proportion") + xlab(Trait) +
#     theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
#     guides(fill=guide_legend(title="Data"))
#   return(p)
# }
# 
# ## Plot for all
# 
# PlotAll <- function(OriginalData, ImputedData, DBisTRUE, GenerationL_or_maturity) {
#   
#   browser()
#   
#   PlotLS <- ToPlot(OriginalData, ImputedData, "Litter_size")
#   PlotGL <- ToPlot(OriginalData, ImputedData, GenerationL_or_maturity)
#   PlotRS <- ToPlot(OriginalData, ImputedData, "Range_size_m2")
#   PlotHB <- ToPlot(OriginalData, ImputedData, "Habitat_breadth_IUCN")
#   PlotBM <- ToPlot(OriginalData, ImputedData, "Body_mass_g")
#   PlotTL <- ToPlotCat(OriginalData, ImputedData, "Trophic_level")
#   PlotSpe <- ToPlotCat(OriginalData, ImputedData, "Specialisation")
#   
#   # If Diet_breadth in taken into account
#   if (isTRUE(DBisTRUE)) {
#     PlotDB <- ToPlotCat(OriginalData, ImputedData, "Diet_breadth")
#       ggarrange(Plot(PlotLS, "log Litter size"),
#             Plot(PlotGL, "log Generation length (d)"),
#             Plot(PlotBM, "log Body mass (g)"),
#             Plot(PlotRS, "log Range size (m sq)"),
#             Plot(PlotHB, "log Habitat breadth"), 
#             PlotCat(PlotTL, "Trophic level"),
#             PlotCat(PlotDB, "Diet breadth"),
#             PlotCat(PlotSpe, "Specialisation"),
#             ncol=2, nrow=4, common.legend = TRUE)
#   } else {
#     ggarrange(Plot(PlotLS, "log Litter size"),
#               Plot(PlotGL, "log Generation length (d)"),
#               Plot(PlotBM, "log Body mass (g)"),
#               Plot(PlotRS, "log Range size (m sq)"),
#               Plot(PlotHB, "log Habitat breadth"), 
#               PlotCat(PlotTL, "Trophic level"),
#               PlotCat(PlotSpe, "Specialisation"),
#               ncol=2, nrow=4, common.legend = TRUE)
#   }
# }
# 
# 
# MammalsPlot <- PlotAll(TraitsMammal.Predicts, Mammals_imputed2, TRUE, "Generation_length_d")
# 
# AmphibiansPlot <- PlotAll(TraitsAmphibian.Predicts, Amphibians_imputed2, TRUE, "Maturity_d")
# ReptilesPlot <- PlotAll(TraitsReptile.Predicts, Reptiles_imputed2, FALSE, "Maturity_d")
# BirdsPlot <- PlotAll(TraitsBird.Predicts, Birds_imputed2, TRUE, "Generation_length_d")
# 
# ggarrange(MammalsPlot,
#           BirdsPlot,
#           ReptilesPlot,
#           AmphibiansPlot,
#           ncol=2, nrow=2,
#           common.legend=TRUE)
# 
# 
# ggsave("../Results/Plots/TraitDistMammalsPredicts.pdf", MammalsPlot, width=7, height =7)
# ggsave("../Results/Plots/TraitDistBirdsPredicts.pdf", BirdsPlot, width=7, height =7)
# ggsave("../Results/Plots/TraitDistReptilesPredicts.pdf", ReptilesPlot, width=7, height =7)
# ggsave("../Results/Plots/TraitDistAmphibiansPredicts.pdf", AmphibiansPlot, width=7, height =7)




