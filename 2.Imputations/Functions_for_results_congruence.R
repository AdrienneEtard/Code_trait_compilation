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
Congruence <- function(List_results, Collected, TraitName) {
  
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

## For continuous traits:
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
    xlab(paste(Trait, "1st imputed set)", sep="")) +
    ylab(paste(Trait, "all other imputed sets)", sep="")) +
    GGPoptions +
    scale_color_hue(labels = c("Without", "With"), name="Phylogenetic information") +
    geom_abline(slope=1, intercept=0)
  
  return(p)

}

Plot.Congruence.Continuous <- function(All_results, Collected, isDB) {

  LCS <- Congruence(All_results, Collected, "Litter_size")
  LG <- Congruence(All_results, Collected, "Longevity_d")
  HB <- Congruence(All_results, Collected, "Habitat_breadth_IUCN")
  BM <- Congruence(All_results, Collected, "Body_mass_g")
  RS <- Congruence(All_results, Collected, "Range_size_m2")

  # Change units: BM -> kg; longevity -> years; RS -> Mm2
  BM[, c(2:9)] <-  BM[, c(2:9)]/1000
  RS[, c(2:9)] <-  RS[, c(2:9)]/1000000000000
  LG[, c(2:9)] <-  LG[, c(2:9)]/365.25
  
  pBM <- Plot_all_results(BM, "BM (kg, ")
  pLCS <- Plot_all_results(LCS, "LCS (")
  pLG <- Plot_all_results(LG, "L (yrs, ")
  pHB <- Plot_all_results(HB, "HB (")
  pRS <- Plot_all_results(RS, "RS (Mm-sq, ")
  
  if(isDB){
    DB <- Congruence(All_results, Collected, "Diet_breadth")
    pDB <- Plot_all_results(DB, "DB (")
    return(
      list(pBM=pBM, pLCS=pLCS, pLG=pLG, pHB=pHB, pRS=pRS, pDB=pDB)
    )
  }
  else{
  return(
    list(pBM=pBM, pLCS=pLCS, pLG=pLG, pHB=pHB, pRS=pRS)
  )}
  
}
ArrangePlots <- function(pM, pB, pR, pA, posTagX, posTagY, isDiet) {
  
  if(isDiet){
  p <- ggarrange(pM + xlab("") + labs(tag = "A") + theme(plot.tag.position = c(posTagX, posTagY)),
                   pB+ ylab("") + xlab("") + labs(tag = "B")+ theme(plot.tag.position = c(posTagX, posTagY)), 
                   pR+ labs(tag = "C")+theme(plot.tag.position = c(posTagX, posTagY)),
                   pA+ labs(tag = "D")+ylab("")+ theme(plot.tag.position = c(posTagX, posTagY)),
                   common.legend = TRUE)
  
  return(p)}
  
  else{
    p <- ggarrange(pM + xlab("") + labs(tag = "A") + theme(plot.tag.position = c(posTagX, posTagY)),
                   pB+ ylab("") + xlab("") + labs(tag = "B")+ theme(plot.tag.position = c(posTagX, posTagY)), 
                   pA+ labs(tag = "D")+ylab("")+ theme(plot.tag.position = c(posTagX, posTagY)),
                   common.legend = TRUE)
    
  }
}

## For categorical traits:

# Assess which agree and which disagree, and barplot
Congruence_cat <- function(List_results, Collected, TraitName, AxisX){
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  Func <- function(X) {
  L <- length(unique(as.character(X)))
  if(L==1) {return("similar")}
  else{return("contradicting")}
}
  DF <- Congruence(List_results, Collected, TraitName)
  DF$Outcome <- apply(DF[, c(2:9)], 1, Func)
  
  ## Barplot of outcomes
  ToPlot <-  with(DF, table(Outcome, Phylo_info)) %>%
    as.data.frame() %>%
    group_by(Phylo_info) %>%
    mutate(prop=Freq/sum(Freq)*100) %>%
    setNames(., c("outcome", "Phylo_info", "freq","prop")) %>%
    as.data.frame()
  
  ToPlot <- ToPlot[order(ToPlot$prop),]
  ToPlot$Phylo_info <- factor(ToPlot$Phylo_info, levels=c("NO", "YES"), labels=c("without", "with"))
  
  p <- ggplot(ToPlot, aes(outcome, prop, fill=Phylo_info)) +
    GGPoptions +
    geom_bar(stat="identity", position="dodge") +
    xlab(AxisX) + ylab("% species") +
    scale_fill_discrete(name = "Phylogenetic information") +
    geom_text(
      aes(label = ToPlot$freq, y = prop + 0.05),
      position = position_dodge(0.9),
      vjust = -0.5
    ) + ylim(0,105)
  
  return(p)
}

# Plot for all the categorical traits
Plot_cat <- function(List_results, Collected, isDiet){
  
  DA <- Congruence_cat(List_results, Collected, "Diel_activity", AxisX = "DA")
  SP <- Congruence_cat(List_results, Collected, "Specialisation", AxisX = "Sp")
  TL <- Congruence_cat(List_results, Collected, "Trophic_level", AxisX = "TL")
  if(isDiet){  
    PD <- Congruence_cat(List_results, Collected, "Primary_diet", AxisX = "DA")
    return(list(pDA=DA, pSP=SP, pTL=TL, pPD=PD))
  }
  else{return(list(pDA=DA, pSP=SP, pTL=TL))}
  
}







