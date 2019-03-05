## Imputation results congruence for each class, with and without phylogenetic information
library(reshape)
library(ggplot2)
library(dplyr)
library(ggpubr)
source("Functions_for_results_congruence.R")

## Load imputed datasets (8)
Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_not_standardised/List_of_8_sets.rds")

## Load data before imputations
Coll_Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
Coll_Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")
Coll_Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
Coll_Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")


## Amphibians
All_Amphibians <- Get_all_results(Imputed, "A")$Results
Amphibians_errors <- Get_all_results(Imputed, "A")$Errors

## Reptiles
All_Reptiles <- Get_all_results(Imputed, "R")$Results
Reptiles_errors <- Get_all_results(Imputed, "R")$Errors

## Mammals
All_Mammals <- Get_all_results(Imputed, "M")$Results
Mammals_errors <- Get_all_results(Imputed, "M")$Errors

## Birds
All_Birds <- Get_all_results(Imputed, "B")$Results
Birds_errors <- Get_all_results(Imputed, "B")$Errors


## Congruence per trait - for continuous traits -- DB missing!
pA <- Plot.Congruence.Continuous(All_Amphibians, Coll_Amphibians, TRUE)
pR <- Plot.Congruence.Continuous(All_Reptiles, Coll_Reptiles, FALSE)
pM <- Plot.Congruence.Continuous(All_Mammals, Coll_Mammals, TRUE)
pB <- Plot.Congruence.Continuous(All_Birds, Coll_Birds, TRUE)

pBM <- ArrangePlots(pM$pBM, pB$pBM, pR$pBM, pA$pBM, posTagX = 0.3, posTagY = 0.9, TRUE)
ggsave(pBM, file="../../Results/Plots/Congruence_imputations/Continuous/BM.pdf", width=6, height = 6)

pLCS <- ArrangePlots(pM$pLCS, pB$pLCS, pR$pLCS, pA$pLCS, posTagX = 0.3, posTagY = 0.9, TRUE)
ggsave(pLCS, file="../../Results/Plots/Congruence_imputations/Continuous/LCS.pdf", width=6, height = 6)

pLG <- ArrangePlots(pM$pLG, pB$pLG, pR$pLG, pA$pLG, posTagX = 0.3, posTagY = 0.9, TRUE)
ggsave(pLG, file="../../Results/Plots/Congruence_imputations/Continuous/LG.pdf", width=6, height = 6)

pHB <- ArrangePlots(pM$pHB, pB$pHB, pR$pHB, pA$pHB, posTagX = 0.3, posTagY = 0.9, TRUE)
ggsave(pHB, file="../../Results/Plots/Congruence_imputations/Continuous/HB.pdf", width=6, height = 6)

pRS <- ArrangePlots(pM$pRS, pB$pRS, pR$pRS, pA$pRS, posTagX = 0.3, posTagY = 0.9, TRUE)
ggsave(pRS, file="../../Results/Plots/Congruence_imputations/Continuous/RS.pdf", width=6, height = 6)

pDB <- ArrangePlots(pM$pDB, pB$pDB, NULL, pA$pDB, posTagX = 0.3, posTagY = 0.9, FALSE)
ggsave(pDB, file="../../Results/Plots/Congruence_imputations/Continuous/DB.pdf", width=6, height = 6)

rm(pA,pR,pM,pB, pBM,pLCS,pLG,pHB,pRS,pDB)

## Congruence for categorical traits (DA, TL, PD, Sp)
pM <- Plot_cat(All_Mammals,Coll_Mammals, TRUE)
pB <- Plot_cat(All_Birds,Coll_Birds, TRUE)
pR <- Plot_cat(All_Reptiles,Coll_Reptiles, TRUE)
pA <- Plot_cat(All_Amphibians,Coll_Amphibians, TRUE)

pDA <- ArrangePlots(pM$pDA, pB$pDA, pR$pDA, pA$pDA, 0.3,0.9, TRUE)
ggsave(pDA, file="../../Results/Plots/Congruence_imputations/Categorical/DA.pdf", width=6, height = 6)

pSp <- ArrangePlots(pM$pSP, pB$pSP, pR$pSP, pA$pSP, 0.3,0.9, TRUE)
ggsave(pSp, file="../../Results/Plots/Congruence_imputations/Categorical/Sp.pdf", width=6, height = 6)

pTL <- ArrangePlots(pM$pTL, pB$pTL, pR$pTL, pA$pTL, 0.3,0.9, TRUE)
ggsave(pTL, file="../../Results/Plots/Congruence_imputations/Categorical/TL.pdf", width=6, height = 6)

pPD <- ArrangePlots(pM$pPD, pB$pPD, NULL, pA$pPD, 0.3,0.9, FALSE)
ggsave(pPD, file="../../Results/Plots/Congruence_imputations/Categorical/PD.pdf", width=6, height = 6)

rm(pA,pR,pM,pB, pDA,pSp,pTL,pPD)




