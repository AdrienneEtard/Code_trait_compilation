# Function to get all imputed datasets for a class

Get_all_results <- function(List, Class=c("A", "B", "M", "R")) {
  
  Results <- list()
  Errors <- list()
  
  for (i in 1:length(List)) {
    Results[[i]] <- List[[i]][[Class]]$Imputed.Dataset  
    Errors[[i]] <- List[[i]][[Class]]$Imputation.errors
  }
   
  return(list(Results=Results, Errors=Errors))
} 


# Function to return all imputed values, with VS without phylogenetic information as column, for a trait, in a dataframe

Congruence_continuous <- function(List_results, Collected, TraitName) {
  
  # Get list of species for which values were imputed
  Sp <- Collected$Best_guess_binomial[is.na(Collected[, TraitName])] 
  
  DF <- Sp %>%
    as.data.frame() %>%
    setNames(., "Best_guess_binomial")
  
  DF[, c(2:(length(List_results)+1))] <- NA

  # Get all imputed values for these species
  for (i in 1:length(List_results)) {
    
    X <- List_results[[i]] %>%
      filter(Best_guess_binomial %in% Sp)
    
    DF[,i+1] <- X[, TraitName]
    colnames(DF)[i+1] <- paste("Imp", i, sep="_")
  }
  
  # Add whether phylogenetic info was present or not
  DF$Phylo_info <- List_results[[1]]$Phylo_info[List_results[[1]]$Best_guess_binomial %in% Sp]
  
  return(DF)
}


Plot_all_results <- function(All_results, Trait){
  
  X <-  melt(All_results[,c(2:(ncol(All_results)))], id.vars = c("Imp_1", "Phylo_info"))
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  p <- ggplot(X, aes(Imp_1, value, col=Phylo_info)) +
    geom_point() +
    xlab(paste(Trait, "1st imputed set)")) +
    ylab(paste(Trait, "all other imputed sets)")) +
    GGPoptions +
    scale_color_hue(labels = c("Without", "With"), name="Phylogenetic information") +
    geom_abline(slope=1, intercept=0)
  
  return(p)

}


Plot.Congruence.Continuous <- function(All_results, Collected) {
  
  LCS <- Congruence_continuous(All_results, Collected, "log10_Litter_size")
  LG <- Congruence_continuous(All_results, Collected, "log10_Longevity_d")
  HB <- Congruence_continuous(All_results, Collected, "sqrt_Habitat_breadth_IUCN")
  BM <- Congruence_continuous(All_results, Collected, "log10_Body_mass_g")
  RS <- Congruence_continuous(All_results, Collected, "Range_size_m2")
  RS[,c(2:(ncol(RS)-1))] <- log10(RS[,c(2:(ncol(RS)-1))])
  
  pBM <- Plot_all_results(BM, "BM (log 10,")
  pLCS <- Plot_all_results(LCS, "LCS (log 10,")
  pLG <- Plot_all_results(LG, "Longevity (log 10,")
  pHB <- Plot_all_results(HB, "HB (square-root,")
  pRS <- Plot_all_results(RS, "RS (log10,")
  
  p <- ggarrange(pBM, pLCS, pLG, pHB, pRS, common.legend = TRUE)
  
  return(p)
  
}


