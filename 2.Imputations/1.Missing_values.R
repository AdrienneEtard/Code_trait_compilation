## Missing values at random?

## Preamble
X <- c("Rphylopars", "dplyr", "phytools", "picante", "stringr", "PVR", "missForest", "colorspace", "ggtree", "ggplot2", "ggpubr", "plyr", "reshape", "reshape2", "mi", "cowplot")
lapply(X, library, character.only=TRUE); rm(X)

.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  return(Phylogeny)
}

## Load trait data, with and without taxonomic correction
Mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
UMammals <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")

Birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
UBirds <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")

Amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")
UAmphibians <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")

Reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")
UReptiles <- read.csv("../../Results/1.Traits_before_imputations/Without_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")


# ## Load phylogenies
# Phylo_Mammals <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Mammals.nwk") %>% .Format_tiplabels()
# Phylo_Amphibians <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Amphibians.nwk") %>% .Format_tiplabels() 
# Phylo_Reptiles <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Reptiles.nwk") %>% .Format_tiplabels()
# Phylo_Birds <- read.newick("../../Results/1.Phylogenies/Corrected/2.Dropped_tips/Birds.nwk") %>% .Format_tiplabels() 

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

TraitsReptiles <- Traits[1:8]

## 1. Plotting percent information across species ("trait filling" across species)

Percent_info_plot <- function(TraitDFC, TraitDFU, Traits, FontSize, BW) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=FontSize, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=FontSize), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=FontSize),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=FontSize))
  
  
  TraitDFC <- TraitDFC[, c("Best_guess_binomial", Traits)]
  TraitDFU <- TraitDFU[, c("Best_guess_binomial", Traits)]
  
  ResultsC <- apply(TraitDFC[,Traits], 1, function(y) sum(!is.na(y))) %>%
    as.data.frame() %>%
    setNames(., "Percent") %>%
    mutate(Percent=Percent/length(Traits)*100) %>%
    mutate(Taxonomy="Corrected")
  
  ResultsU <- apply(TraitDFU[,Traits], 1, function(y) sum(!is.na(y))) %>% as.data.frame() %>%
    setNames(., "Percent") %>%
    mutate(Percent=Percent/length(Traits)*100) %>%
    mutate(Taxonomy="Uncorrected")
  
  Results <- rbind(ResultsC, ResultsU)
  Results$Taxonomy <- factor(Results$Taxonomy, levels = c("Uncorrected", "Corrected"))
  
  # p <- ggplot(Results, aes(Percent, fill=Taxonomy)) +
  #   geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
  #                          ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
  #                  position='identity', binwidth=Bin) +
  #   ylab("Percentage") + xlab("Trait filling across species") + GGPoptions
  
  
  med <- ddply(Results, "Taxonomy", summarise, grp.median=median(Percent)); print(med)
  p <- ggplot(Results, aes(Percent, fill=Taxonomy)) +
    geom_density(alpha=0.5, adjust=BW) + GGPoptions +
    geom_vline(data=med, aes(xintercept=grp.median), linetype="dashed") +
    xlab("Trait filling (%)") + ylab("Density")
    
  return(p)
  
}

TraitFillMammals <- Percent_info_plot(Mammals, UMammals, Traits, 13, 2) + xlab(NULL)
TraitFillBirds <- Percent_info_plot(Birds, UBirds, Traits, 13, 5) + ylab(NULL) + xlab(NULL)
TraitFillReptiles <- Percent_info_plot(Reptiles, UReptiles, TraitsReptiles, 13, 2) 
TraitFillAmphibians <- Percent_info_plot(Amphibians, UAmphibians, Traits, 13, 2) + ylab(NULL)
p <- ggarrange(TraitFillMammals, TraitFillBirds, TraitFillReptiles, TraitFillAmphibians, common.legend = TRUE,
               labels=c("A", "B", "C", "D"), hjust=c(-6,-6,-22,-22), vjust=2, font.label = list(family="serif"))
ggsave(p, file="../../Results/Plots/Trait_missing_values/Traitfillingdist.pdf", width=6.5, height=5)

## 2. Are missing values missing at random? Association between NAs and taxonomy.

Xsquare_Tax_NA <- function(TraitDF, Traits, Family_or_Order) {
  
  TraitDF$Specialisation <- as.character(TraitDF$Specialisation)
  TraitDF$Trophic_level <- as.character(TraitDF$Trophic_level)
  TraitDF$Diel_activity <- as.character(TraitDF$Diel_activity)
  if(any(grepl("Primary_diet", colnames(TraitDF)))) { 
    TraitDF$Primary_diet <- as.character(TraitDF$Primary_diet)}
  
  List_Xsq <- list()
  List_cont_tbl <- list()
  
  for (i in 1:length(Traits)) {
    
    DF <- TraitDF[, c(Traits[i], Family_or_Order)]
    M <- subset(DF, is.na(DF[,1]))
    
    if(nrow(M)==0) {
      print("No missing values for that trait.")
      List_Xsq[[i]] <- "No missing values."
    }
    
    else{
      DF[is.na(DF[, Traits[i]]), Traits[i]] <- "Missing"
      DF[ DF[, Traits[i]]!="Missing", Traits[i]] <- "Non-missing"
      X <- table(DF[, Family_or_Order], DF[, Traits[i]])
      Y <- chisq.test(X) %>% as.list()
      
      if(Y[["p.value"]]>=0.05) {
        List_Xsq[[i]] <- "No significant association for this trait."
      }
      
      else{ 
        
        browser()
        
        List_Xsq[[i]] <- Y
        List_cont_tbl[[i]] <- as.data.frame(X)
        }
    }
  }
  
  names(List_Xsq) <- Traits
  names(List_cont_tbl) <- Traits
  
  return(list(Chi_squares=List_Xsq, Contingency_tables=List_cont_tbl))
}


Mammals_Order <- Xsquare_Tax_NA(Mammals, Traits, "Order")
XsqMammals_Order <- Mammals_Order$Chi_squares
Tbl <- Mammals_Order$Contingency_tables

Mammals_Family <- Xsquare_Tax_NA(Mammals, Traits, "Family")
XsqMammals_Family <- Mammals_Family$Chi_squares
TblMammals_Family <- Mammals_Family$Contingency_tables

Amphibians_Order <- Xsquare_Tax_NA(Amphibians, Traits, "Order")
XsqAmphibians_Order <- Amphibians_Order$Chi_squares
TblAmphibians_Order <- Amphibians_Order$Contingency_tables

Amphibians_Family <- Xsquare_Tax_NA(Amphibians, Traits, "Family")
XsqAmphibians_Family <- Amphibians_Family$Chi_squares
TblAmphibians_Family <- Amphibians_Family$Contingency_tables



X <- Tbl[[2]]
X <- 



## 3. Plotting median percent of missing values across family or order for each trait



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
  
  rownames(All) <- c(1:nrow(All))
  
  ## Add number of species in each family
  
  if(Order_or_Family=="Family") {
      Species <- DF %>% 
        dplyr::group_by(Family) %>%
        dplyr::summarise(Count=n()) %>% 
        as.data.frame()
      
      Species <- Species[order(Species$Family),]
  }
  
  if(Order_or_Family=="Order") {
    Species <- DF %>% 
      dplyr::group_by(Order) %>%
      dplyr::summarise(Count=n()) %>% 
      as.data.frame()
    Species <- Species[order(Species$Order),]
  }
  
  All$N <- c(as.vector(Species$Count), length(unique(DF$Best_guess_binomial)))
    
  return(list(List=List, Percentages=All))
}

# Calculate for each trait a weigthed median of percent missing information, with weigthed quantiles

WeightedMediansPlot <- function(NADF, FontSize, Ylab, TraitDF, Family, Xsub, Ysub) {

  NADF$N <- NADF$N / NADF$N[nrow(NADF)]
  
  x <- ncol(NADF)-1
  y <- nrow(NADF)-1
  
  R <- apply(NADF[c(1:y), c(2:x)], 2, spatstat::weighted.median, w=NADF$N) %>%
  as.data.frame() %>%
  setNames("median")
  R$qt.25 <-  apply(NADF[, c(2:x)], 2, spatstat::weighted.quantile, w=NADF$N, probs=0.25)
  R$qt.75 <-  apply(NADF[, c(2:x)], 2, spatstat::weighted.quantile, w=NADF$N, probs=0.75)
  R$Weighted <- "yes"
  R$Traits <- rownames(R)
  
  R2 <- apply(NADF[c(1:y), c(2:x)], 2, median, w=NADF$N) %>%
    as.data.frame() %>%
    setNames("median")
    R2$qt.25 <-  apply(NADF[, c(2:x)], 2, quantile, w=NADF$N, probs=0.25)
    R2$qt.75 <-  apply(NADF[, c(2:x)], 2, quantile, w=NADF$N, probs=0.75)
    R2$Weighted <- "no"
    R2$Traits <- rownames(R2)
  
  R <- rbind(R, R2)
  
  Names <- as.data.frame(c("Body_mass_g", "Adult_svl_cm", "Forearm_length_mm","Head_length_mm","Body_length_mm",
                           "Svl_length_mm","Generation_length_d", "Longevity_d", "Maturity_d", "AFR_d",
                           "Litter_size", "Range_size_m2", "Diel_activity", "Trophic_level","Primary_diet", "Diet_breadth", "Specialisation","Habitat_breadth_IUCN", "EV_1"))
  colnames(Names) <- "Original"
  Names$FP <- c("Body mass", "Svl length", "Forearm length", "Head length", "Body length", "Svl length", "Generation length", "Longevity", "Sexual maturity age",
                "Age 1st reproduction","Litter/clutch size", "Range size",
                "Diel activity", "Trophic level", "Primary diet","Diet breadth","Specialisation", "Habitat breadth", "Phylogenetic position")
  
  for (i in 1:nrow(R)) {R$Traits[i] <- Names$FP[Names$Original==R$Traits[i]]}
  
  
  GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=FontSize, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,5,0,"pt"), size=FontSize), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=FontSize),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=FontSize))

  
 p <-  ggplot(R, aes(reorder(Traits, median), median, group=Weighted, col=Weighted)) +
    geom_point(size=2, position=position_dodge(width=0.4)) +
    geom_errorbar(aes(ymax = qt.75, ymin = qt.25), size =0.5, width=0,position=position_dodge(width=0.4)) +
    scale_color_manual(values = c("black", "red")) +
    GGPoptions + xlab("") + ylab(Ylab) + coord_flip()

 
 ## Adding distribution of number of species per families as subplots
 if(Family) {
 NF <- table(TraitDF$Family) %>% as.data.frame()}
 else {
   NF <- table(TraitDF$Order) %>% as.data.frame()}
 
 q <- ggplot(NF, aes(Freq)) + geom_histogram(alpha=0.5, fill="black") + GGPoptions +
   xlab("freq") + ylab("") + 
   theme(text = element_text(size=FontSize-2, family="serif"),  
         axis.text.x = element_text(size=FontSize-2), 
         axis.text.y = element_text(size=FontSize-2))

 z <- ggdraw() +
   draw_plot(p) +
   draw_plot(q, x = Xsub, y = Ysub, width = .3, height = .3)
 
  return(z)
}


# Family level
Amphibians_F <- MissingVal(Amphibians, "Family", c(Traits, "Body_length_mm"))$Percentages
Reptiles_F <- MissingVal(Reptiles, "Family", c(TraitsReptiles, "Adult_svl_cm", "Maturity_d"))$Percentages
Birds_F <- MissingVal(Birds, "Family", Traits)$Percentages
Mammals_F <- MissingVal(Mammals, "Family", c(Traits, "Adult_svl_cm", "Generation_length_d"))$Percentages


pFA <- WeightedMediansPlot(Amphibians_F, 13, NULL, Family = TRUE, TraitDF = Amphibians, Xsub = 0.3, Ysub=0.5)
pFR <- WeightedMediansPlot(Reptiles_F, 13, NULL, Family=TRUE, TraitDF = Reptiles, Xsub = 0.3, Ysub=0.5 )
pFB <- WeightedMediansPlot(Birds_F, 13, NULL, Family=TRUE, TraitDF = Birds, Xsub = 0.3, Ysub=0.5)
pFM <- WeightedMediansPlot(Mammals_F, 13, NULL, Family=TRUE, TraitDF = Mammals, Xsub = 0.3, Ysub=0.5)

p <- ggarrange(pFM, pFB, pFR, pFA, common.legend = TRUE) #add corner labels, add distribtuion of number of species per family as inset plots
# ggsave

# How many species per families?
FM <- table(Mammals$Family) %>% as.data.frame()
FM <- ggplot(FM, aes(Freq)) +geom_histogram(alpha=0.5) + ylim(c(0,100)) + GGPoptions

FB <- table(Birds$Family) %>% as.data.frame()
FB <- ggplot(FB, aes(Freq))+geom_histogram(alpha=0.5, fill=I("black")) + ylim(c(0,100)) + GGPoptions

FR <- table(Reptiles$Family) %>% as.data.frame()
FR <- ggplot(FR, aes(Freq)) +geom_histogram(alpha=0.5, fill=I("black")) + ylim(c(0,100)) + GGPoptions

FA <- table(Amphibians$Family) %>% as.data.frame()
FA <- ggplot(FA, aes(Freq))+geom_histogram(alpha=0.5, fill=I("black")) + ylim(c(0,100)) + GGPoptions

ggarrange(FM, FB, FR, FA)

