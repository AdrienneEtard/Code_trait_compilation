## Compare Rob Cooke's data to mine
# diet?
# volancy: I have not collected the data

library(dplyr)
library(ggplot2)
library(ggpubr)

## Load data

# Robert Cooke's data, imputed and collected (need to be standardised within class and scaled to 0-mean and unit variance)
RC_Collected <- read.csv("../../Data/RCooke_Bird_and_Mammal_imputed_data/Cooke_et_al_2018_GEB_collected_trait_data.csv") %>%
  mutate(activity=ifelse(activity==1, "Other", ifelse(activity==2, "Nocturnal", NA)))

RC_Imputed <- read.csv("../../Data/RCooke_Bird_and_Mammal_imputed_data/Cooke_et_al_2018_GEB_single_imputed.csv") %>%
  mutate(activity=ifelse(activity==1, "Other", ifelse(activity==2, "Nocturnal", NA)))

# My data, collected and imputed (one dataframe selected randomly)

# Collected
X <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Mammals.csv") %>% 
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN) %>%
  mutate(Class="Mammalia") 

Y <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/4.with_phylo_eigenvectors/Birds.csv") %>%
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN) %>%
  mutate(Class="Aves") 

AE_Collected <- rbind(X,Y) %>%  mutate(Best_guess_binomial=as.character(Best_guess_binomial))
AE_Collected <- AE_Collected[order(AE_Collected$Best_guess_binomial),]
rm(X,Y)

# Imputed 
# select one randomly
i <- sample(1:10,1)
AE_Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_datasets/List_8_imputed_sets.rds")[[i]]

X <- AE_Imputed$M$Imputed.Dataset %>%
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN) %>%
  mutate(Class="Mammalia")

Y <- AE_Imputed$B$Imputed.Dataset %>%
  select(Best_guess_binomial, Diel_activity, log10_Body_mass_g, log10_Litter_size, sqrt_Habitat_breadth_IUCN) %>%
  mutate(Class="Aves")

AE_Imputed <- rbind(X,Y) %>%  mutate(Best_guess_binomial=as.character(Best_guess_binomial))
AE_Imputed <- AE_Imputed[order(AE_Imputed$Best_guess_binomial),] 
rm(X,Y)

## 1. Compare initial trait coverage


## 2. Comparison of collected data for mammals and birds
p1 <- Compare(RC_Collected, NULL, AE_Collected, NULL, "body_mass_median","log10_Body_mass_g", FALSE,"Body mass","collected", "collected", FALSE)
p2 <- Compare(RC_Collected, NULL, AE_Collected, NULL,  "litter_clutch_size","log10_Litter_size", FALSE,"Litter/clutch size", "collected", "collected",FALSE)
p3 <- Compare(RC_Collected, NULL, AE_Collected, NULL,"hab_breadth","sqrt_Habitat_breadth_IUCN", FALSE, NULL, "collected", "collected", FALSE)
p4 <- Compare(RC_Collected, NULL,  AE_Collected, NULL,  "activity","Diel_activity", TRUE, "Diel activity", "collected", "collected", FALSE)
p <- ggarrange(p1,p2,p3,p4,common.legend = TRUE)
p

## 3. Comparison of imputed data for mammals and birds
p1 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed, "body_mass_median","log10_Body_mass_g", FALSE,"Body mass","imputed","imputed", TRUE)
p2 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,  "litter_clutch_size","log10_Litter_size", FALSE,"Litter/clutch size", "imputed", "imputed",TRUE)
p3 <- Compare(RC_Collected, RC_Imputed, AE_Collected, AE_Imputed,"hab_breadth","sqrt_Habitat_breadth_IUCN", FALSE, NULL, "imputed", "imputed", TRUE)
p4 <- Compare(RC_Collected, RC_Imputed,  AE_Collected, AE_Imputed,  "activity","Diel_activity", TRUE, "Diel activity", "imputed","imputed", TRUE)
p <- ggarrange(p2,p1,p3,p4,common.legend = TRUE)
p


# #   F  U  N  C  T  I  O  N  S 
Transform_zscore <- function(TraitDF, Trait, Transf) {
  
  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- log10(TraitDF[,Trait])
    TraitDF[, paste("log10", Trait, sep="_")] <- scale(TraitDF[, paste("log10", Trait, sep="_")], center=TRUE, scale=TRUE)
  }
  
  if(Transf=="sqrt") {
    TraitDF[,Trait]  <- sqrt(TraitDF[,Trait])
    TraitDF[, paste("sqrt", Trait, sep="_")] <- scale(TraitDF[,Trait] , center=TRUE, scale=TRUE)
  }
  
  return(TraitDF)
}

Compare <- function(RC_data, RC_imputed, AE_data, AE_imputed, RC_TraitName, AE_TraitName, Categorical, Traitaxis, ImputedaxisX, ImputedaxisY, Imputed) {
  
  
  # ImputedaxisX <- ImputedaxisX
  # ImputedaxisY <- ImputedaxisY
  
  if(Imputed) {
    
    # filter species for which the trait value was imputed
    
    # AE data
    Sp_AE_I <- AE_data$Best_guess_binomial[is.na(AE_data[,AE_TraitName])]
    
    # if the coverage for that trait was 100% initially, compared my collected values to Rob's imputed values (if they exist)
    if (length(Sp_AE_I)!=0) {
      AE_Imputed <- AE_Imputed %>%
        filter(Best_guess_binomial %in% Sp_AE_I)
      AE_data <- AE_Imputed
    }
    else{ImputedaxisX <- "collected"}
    
    # RC data
    Sp_RC_I <- RC_data$binomial[is.na(RC_data[,RC_TraitName])]
    
    # if the coverage for that trait was 100% initially, compared my imputed values to Rob's collected values
    if (length(Sp_RC_I)!=0) {
      RC_Imputed <- RC_Imputed %>%
        filter(binomial %in% Sp_RC_I)
      RC_data <- RC_Imputed
    }
    else{ImputedaxisY <- "collected"}
    
  }
  
  # Standardise and scale RC data within each class
  if (RC_TraitName %in% c("body_mass_median", "litter_clutch_size")){
    RC_Birds <- subset(RC_data, class=="Aves")
    RC_Birds <- Transform_zscore(RC_Birds, RC_TraitName, "log10")
    RC_Mammals <- subset(RC_data, class=="Mammalia")
    RC_Mammals <- Transform_zscore(RC_Mammals, RC_TraitName, "log10")
    RC_data <- rbind(RC_Birds, RC_Mammals)
    RC_data <- RC_data[order(RC_data$binomial),]
    RC_TraitName <- paste("log10", RC_TraitName, sep="_")}
  
  if (RC_TraitName=="hab_breadth"){
    RC_Birds <- subset(RC_data, class=="Aves")
    RC_Birds <- Transform_zscore(RC_Birds, RC_TraitName, "sqrt")
    RC_Mammals <- subset(RC_data, class=="Mammalia")
    RC_Mammals <- Transform_zscore(RC_Mammals, RC_TraitName, "sqrt")
    RC_data <- rbind(RC_Birds, RC_Mammals)
    RC_data <- RC_data[order(RC_data$binomial),]
    RC_TraitName <- paste("sqrt", RC_TraitName, sep="_")}
  
  # intersect AE and RC species
  Y <- intersect(RC_data$binomial, AE_data$Best_guess_binomial)
  RC_data <- RC_data %>% filter(binomial %in% Y)
  AE_data <- AE_data %>% filter(Best_guess_binomial %in% Y)
  
  # plot
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=12, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,5,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=12))
  
  if(!Categorical) {
    
    ToPlot <- RC_data[, c("class","binomial", RC_TraitName)]
    colnames(ToPlot)[3] <- "RC"
    ToPlot$AE <- AE_data[, AE_TraitName]
    
    if(RC_TraitName %in% c("log10_litter_clutch_size", "log10_body_mass_median")){
      
      p <- ggplot(ToPlot, aes(AE, RC, color=ToPlot$class)) +
        GGPoptions +
        geom_point(alpha=0.7) +
        geom_abline(slope = 1, intercept = 0, alpha=0.8) +
        xlab(paste(Traitaxis, "(log10, ", ImputedaxisX," by AE)")) + 
        ylab(paste(Traitaxis,  "(log10, ", ImputedaxisY," by RC)")) + 
        scale_color_hue(name="Class", labels = c("Birds", "Mammals"))
      
      return(p)
    }
    
    else{
      
      p <- ggplot(ToPlot, aes(AE, RC, color=ToPlot$class)) +
        GGPoptions +
        geom_point(alpha=0.7) +
        geom_abline(slope = 1, intercept = 0, alpha=0.8) +
        xlab(paste("Habitat breadth (square-root, ", ImputedaxisX," by AE)")) + 
        ylab(paste("Habitat breadth (square-root, ",  ImputedaxisY, " by RC)")) + 
        scale_color_hue(name="Class", labels = c("Birds", "Mammals"))
      
      return(p)
    }
    
  }
  
  if(Categorical) {
    
    Outcome <- AE_data$Best_guess_binomial %>% as.data.frame()
    
    for (i in 1:nrow(AE_data)){
      
      x <- AE_data[i, AE_TraitName]
      y <- RC_data[i, RC_TraitName]
      
      if(is.na(x)|is.na(y)) {Outcome$Result[i] <- "unknown"}
      
      else{
        
        if(x==y) {Outcome$Result[i] <- "same"} 
        if(x!=y)  { Outcome$Result[i] <- "different"}
      }
    }
    
    ToPlot <- table(Outcome$Result) %>% 
      as.data.frame() %>%
      setNames(., c("outcome", "prop")) %>%
      mutate(prop=prop/nrow(Outcome)*100)
    
    ToPlot <- ToPlot[order(ToPlot$prop),]
    
    p <- ggplot(ToPlot, aes(outcome, prop)) +
      GGPoptions +
      geom_bar(stat="identity") +
      xlab(paste(Traitaxis, ImputedaxisX, sep = ", ")) + ylab("% species") +
      scale_x_discrete(limits=c("different", "unknown", "same"), labels=c("Contradicting", "Unkown", "Similar"))
    
    return(p)
  }
}
