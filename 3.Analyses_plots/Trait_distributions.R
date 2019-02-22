## Comparing outputs in continuous trait distributions from missForest and phylogenetic imputations
## On PREDICTS subsets
## NB - Birds missing (need to run phylogenetic imputations)

X <- c("dplyr", "ggplot2", "ggpubr")
invisible(lapply(X, library, character.only=TRUE)); rm(X)


## Load data - results from missForest
MF_Amphibians <- read.csv("../../Results/missForest_imputations/AmphibiansPredicts.csv")
MF_Reptiles <- read.csv("../../Results/missForest_imputations/ReptilesPredicts.csv")
MF_Mammals <- read.csv("../../Results/missForest_imputations/MammalsPredicts.csv")


## Load data - results from phylogenetic imputations
PI_Amphibians <- read.csv("../../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsAmphibians.Predicts_ContImputed_Randomadd.csv")
PI_Reptiles <- read.csv("../../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsReptile.Predicts_ContImputed.csv")
PI_Mammals <- read.csv("../../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsMammal.Predicts_ContImputed_Randomadd.csv")



## Density
ToPlot <- function(mfdf, ipdf, Trait) {
  ToPlot1 <- mfdf[, colnames(mfdf) %in% Trait] %>% as.data.frame() %>% setNames(., "log")
  ToPlot1$Which <- "missForest"
  ToPlot2 <- ipdf[,  colnames(ipdf) %in% Trait]  %>% as.data.frame() %>% setNames(., "log")
  ToPlot2$Which <- "phylopars"
  ToPlot <- rbind(ToPlot1, ToPlot2)
  ToPlot$log <- log(ToPlot$log)
  return(ToPlot)
}


Plot <- function(ToPlot, Trait) {
  p <- ggplot(ToPlot, aes(x=log, fill=Which), family="Times") +
    geom_density(alpha=0.6) + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("Density") + xlab(Trait) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    guides(fill=guide_legend(title="Method"))
  return(p)
}


Plot.Boxplot <- function(ToPlot, Trait) {
  p <- ggplot(ToPlot, aes(y=log, x=Which, fill=Which), family="Times") +
    geom_boxplot(alpha=0.6, width=0.4) + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    ylab(Trait) + guides(fill=guide_legend(title="Method")) + xlab(NULL)
  return(p)
}


PlotAll <- function(mfdata, pidata, DBisTRUE, GenerationL_or_maturity, Boxplot) {
  
  PlotLS <- ToPlot(mfdata, pidata, "Litter_size")
  PlotGL <- ToPlot(mfdata, pidata, GenerationL_or_maturity)
  PlotRS <- ToPlot(mfdata, pidata, "Range_size_m2")
  PlotHB <- ToPlot(mfdata, pidata, "Habitat_breadth_IUCN")
  PlotBM <- ToPlot(mfdata, pidata, "Body_mass_g")
  
  
  if (isTRUE(Boxplot)) {
    # If Diet_breadth in taken into account
    if (isTRUE(DBisTRUE)) {
      PlotDB <- ToPlot(mfdata, pidata, "Diet_breadth")
      ggarrange(Plot.Boxplot(PlotLS, "log Litter size"),
                Plot.Boxplot(PlotGL, "log Generation length (d)"),
                Plot.Boxplot(PlotBM, "log Body mass (g)"),
                Plot.Boxplot(PlotRS, "log Range size (m sq)"),
                Plot.Boxplot(PlotHB, "log Habitat breadth"), 
                Plot.Boxplot(PlotDB, "log Diet breadth"),
                ncol=2, nrow=3, common.legend = TRUE)
    } else {
      ggarrange(Plot.Boxplot(PlotLS, "log Litter size"),
                Plot.Boxplot(PlotGL, "log Generation length (d)"),
                Plot.Boxplot(PlotBM, "log Body mass (g)"),
                Plot.Boxplot(PlotRS, "log Range size (m sq)"),
                Plot.Boxplot(PlotHB, "log Habitat breadth"), 
                ncol=2, nrow=3, common.legend = TRUE
      )
    }
  } else {
  
  
  # If Diet_breadth in taken into account
  if (isTRUE(DBisTRUE)) {
    PlotDB <- ToPlot(mfdata, pidata, "Diet_breadth")
    ggarrange(Plot(PlotLS, "log Litter size"),
              Plot(PlotGL, "log Generation length (d)"),
              Plot(PlotBM, "log Body mass (g)"),
              Plot(PlotRS, "log Range size (m sq)"),
              Plot(PlotHB, "log Habitat breadth"), 
              Plot(PlotDB, "log Diet breadth"),
              ncol=2, nrow=3, common.legend = TRUE)
  } else {
    ggarrange(Plot(PlotLS, "log Litter size"),
              Plot(PlotGL, "log Generation length (d)"),
              Plot(PlotBM, "log Body mass (g)"),
              Plot(PlotRS, "log Range size (m sq)"),
              Plot(PlotHB, "log Habitat breadth"), 
              ncol=2, nrow=3, common.legend = TRUE
              )
  }
  }
  }



## Densities
Amphibians <- annotate_figure(PlotAll(MF_Amphibians, PI_Amphibians, TRUE, "Maturity_d", FALSE),
                bottom = "PREDICTS amphibians: distribution of continuous trait values")

Mammals <- annotate_figure(PlotAll(MF_Mammals, PI_Mammals, TRUE, "Generation_length_d", FALSE), 
                           bottom = "PREDICTS mammals: distribution of continuous trait values")

Reptiles <- annotate_figure(PlotAll(MF_Reptiles, PI_Reptiles, FALSE, "Maturity_d", FALSE), 
                            bottom = "PREDICTS reptiles: distribution of continuous trait values")


## Boxplots
AmphibiansBoxplot <- annotate_figure(PlotAll(MF_Amphibians, PI_Amphibians, FALSE, "Maturity_d", TRUE),
                                     bottom = "PREDICTS amphibians: continuous trait values")

MammalsBoxplot <- annotate_figure(PlotAll(MF_Mammals, PI_Mammals, FALSE, "Generation_length_d", TRUE),
                                     bottom = "PREDICTS mammals: continuous trait values")

ReptilesBoxplot <- annotate_figure(PlotAll(MF_Reptiles, PI_Reptiles, FALSE, "Maturity_d", TRUE),
                                     bottom = "PREDICTS reptiles: continuous trait values")




## Save files

ggsave("../../Results/Plots/missForest_VS_phylopars_Amphibians_dist.pdf", Amphibians, width=6, height = 5)
ggsave("../../Results/Plots/missForest_VS_phylopars_Mammals_dist.pdf", Mammals, width=6, height = 5)
ggsave("../../Results/Plots/missForest_VS_phylopars_Reptiles_dist.pdf", Reptiles, width=6, height = 5)

ggsave("../../Results/Plots/missForest_VS_phylopars_Amphibians_bp.pdf", AmphibiansBoxplot, width=6, height = 5)
ggsave("../../Results/Plots/missForest_VS_phylopars_Mammals_bp.pdf", MammalsBoxplot, width=6, height = 5)
ggsave("../../Results/Plots/missForest_VS_phylopars_Reptiles_bp.pdf", ReptilesBoxplot, width=6, height = 5)


## Contrasting predictions

ToPlot_contrast <- function(mfdf, ipdf, Trait) {
  ToPlot1 <- mfdf[, colnames(mfdf) %in% Trait] %>% as.data.frame() %>% setNames(., "missForest")
  ToPlot2 <- ipdf[,  colnames(ipdf) %in% Trait]  %>% as.data.frame() %>% setNames(., "phylopars")
  ToPlot <- cbind(ToPlot1, ToPlot2)
  ToPlot <- log(ToPlot)
  return(ToPlot)
}

PlotContrast <- function(ToPlot, Trait) {
  p <- ggplot(ToPlot, aes(x=missForest, y=phylopars), family="Times") +
    geom_point() + theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
    ylab("phylopars predictions") + xlab("missForest predictions") + ggtitle(Trait) +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black")) +
    geom_abline(slope=1,intercept=0, alpha=0.7)
  return(p)
}

PlotAllContrast <- function(mfdata, pidata, DBisTRUE, GenerationL_or_maturity) {
  
  # browser()
  
  PlotLS <- ToPlot_contrast(mfdata, pidata, "Litter_size")
  PlotGL <- ToPlot_contrast(mfdata, pidata, GenerationL_or_maturity)
  PlotRS <- ToPlot_contrast(mfdata, pidata, "Range_size_m2")
  PlotHB <- ToPlot_contrast(mfdata, pidata, "Habitat_breadth_IUCN")
  PlotBM <- ToPlot_contrast(mfdata, pidata, "Body_mass_g")
  
    # If Diet_breadth in taken into account
    if (isTRUE(DBisTRUE)) {
      PlotDB <- ToPlot_contrast(mfdata, pidata, "Diet_breadth")
      ggarrange(PlotContrast(PlotLS, "log Litter size"),
                PlotContrast(PlotGL, "log Generation length (d)"),
                PlotContrast(PlotBM, "log Body mass (g)"),
                PlotContrast(PlotRS, "log Range size (m sq)"),
                PlotContrast(PlotHB, "log Habitat breadth"), 
                PlotContrast(PlotDB, "log Diet breadth"),
                ncol=2, nrow=3, common.legend = TRUE)
    } else {
      ggarrange(PlotContrast(PlotLS, "log Litter size"),
                PlotContrast(PlotGL, "log Generation length (d)"),
                PlotContrast(PlotBM, "log Body mass (g)"),
                PlotContrast(PlotRS, "log Range size (m sq)"),
                PlotContrast(PlotHB, "log Habitat breadth"), 
                ncol=2, nrow=3, common.legend = TRUE
      )
    }
}



MF_Amphibians <- MF_Amphibians[order(MF_Amphibians$Best_guess_binomial),]
MF_Mammals <- MF_Mammals[order(MF_Mammals$Best_guess_binomial),]
MF_Reptiles <- MF_Reptiles[order(MF_Reptiles$Best_guess_binomial),]

ContrastAmphibians <-  PlotAllContrast(MF_Amphibians, PI_Amphibians, FALSE, "Maturity_d")
ContrastMammals <-  PlotAllContrast(MF_Mammals, PI_Mammals, FALSE, "Generation_length_d")
ContrastReptiles <-  PlotAllContrast(MF_Reptiles, PI_Reptiles, FALSE, "Maturity_d")

ggsave("../../Results/Plots/missForest_VS_phylopars_contrast_amphibians.pdf", ContrastAmphibians, width = 6, height=5)
ggsave("../../Results/Plots/missForest_VS_phylopars_contrast_mammals.pdf", ContrastMammals, width = 6, height=5)
ggsave("../../Results/Plots/missForest_VS_phylopars_contrast_reptiles.pdf", ContrastReptiles, width = 6, height=5)



