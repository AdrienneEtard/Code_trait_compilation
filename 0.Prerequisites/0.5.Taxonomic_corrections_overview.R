## Taxonomic corrections: overview

#### Preamble

library(dplyr)
library(reshape)
library(ggplot2)
library(ggpubr)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))

#### Functions

## Overall differences in species number
Delta <- function(SynM, SynA, SynB, SynR) {
  DF <-as.data.frame(matrix(ncol=3, nrow=4)) %>%
    setNames(., c("original","corrected","accepted"))
  rownames(DF) <- c("mammals","birds","reptiles","amphibians")
  
  DF[1,] <- c(length(unique(SynM$Original)),
              length(unique(SynM$CorrectedTypos)),
              length(unique(SynM$Accepted)))
  
  DF[2,] <- c(length(unique(SynB$Original)),
              length(unique(SynB$CorrectedTypos)),
              length(unique(SynB$Accepted)))
  
  DF[3,] <- c(length(unique(SynR$Original)),
              length(unique(SynR$CorrectedTypos)),
              length(unique(SynR$Accepted)))
  
  DF[4,] <- c(length(unique(SynA$Original)),
              length(unique(SynA$CorrectedTypos)),
              length(unique(SynA$Accepted)))
  DF$Class <- rownames(DF)
  
  toplot <- DF %>% select(-corrected)
  toplot <- reshape::melt(toplot)
  toplot$Class <- factor(toplot$Class, levels=c("birds","reptiles","amphibians","mammals"), labels=c("Birds", "Reptiles",
                                                                                                     "Amphibians", "Mammals"))
  
  p <- ggplot(toplot, aes(Class,value, fill=variable)) + 
    geom_bar(stat="identity", position="dodge") + ylab("Names \nacross all datasets") +
    scale_fill_manual(values=c("cornflowerblue","coral"), name="", labels=c("original", "extracted")) + GGPoptions +
    theme(axis.text.x = element_text(angle = 17, hjust = 1))
  
  return(list(p=p, n=toplot))
  
}

## Species with the most number of replicates
MaxSyn <- function(Syn) {
  Table <- Syn %>% 
    filter(Accepted!="") %>%
    group_by(Accepted) %>% 
    summarise(Count=n())
  Max <- Table$Accepted[Table$Count==max(Table$Count)]
  if(length(Max)==1) {  All <- Syn$CorrectedTypos[Syn$Accepted==Max]}
  else{All1 <- Syn$CorrectedTypos[Syn$Accepted==Max[1]] %>% as.character()
  All2 <- Syn$CorrectedTypos[Syn$Accepted==Max[2]]%>% as.character()
  return(list(Accepted=Max, Replicates1=All1, Replicates2=All2))
  }
  return(list(Accepted=Max, Replicates=All))
}

## Function to plot distribution of number of names for each accepted name
Dist_NR <- function(Syn) {
  
  GGPoptions <- theme_bw() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  Table <- Syn %>% 
    filter(Accepted!="") %>%
    group_by(Class, Accepted) %>%
    summarise(Count=n())
  
  Table <- table(Table$Count, Table$Class) %>%
    as.data.frame() 
  Table_top <- Table %>% 
    #filter(Var1!=1) %>%
    filter(Freq!=0)
  
  ggplot(Table_top, aes(Var1, Freq, fill=Var2)) + geom_point(size=3,pch=21,alpha=0.5) + GGPoptions +
    scale_y_continuous(trans="log10", breaks=c(1,2,5, 10,25, 50, 100,250, 500, 1000,2000,5000)) +
    scale_fill_discrete(name="Class", labels=c("Amphibians","Birds","Mammals","Reptiles")) +
    ylab("Number of species") + xlab("Number of synonyms") %>%
    return()
  
}

#### Data

## Load synonym datasets
Syn_Mammals <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Mammals.csv") %>% mutate(Class="Mammalia")
Syn_Birds <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Birds.csv")%>% mutate(Class="Aves")
colnames(Syn_Birds)[15] <- "Manual_edits"
Syn_Reptiles <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Reptiles.csv") %>% mutate(Class="Reptilia")
colnames(Syn_Birds)[15]
Syn_Amphibians <- read.csv("../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/Amphibians.csv")%>% mutate(Manual_edits=NA, Notes=NA, Class="Amphibia")
Syn <- rbind(Syn_Mammals, Syn_Birds, Syn_Reptiles, Syn_Amphibians)
write.csv(Syn, "../../Results/0.Data_resolved_taxonomy/List_species_synonyms/Synonym_final_datasets/SingleFinalDataset.csv", row.names = FALSE)

## Load species in compiled trait datasets
M <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv")
B <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv")
R <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv")
A <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv")

SynM <- Syn_Mammals %>% filter(Accepted %in% M$Best_guess_binomial)
SynB <- Syn_Birds %>% filter(Accepted %in% B$Best_guess_binomial)
SynR <- Syn_Reptiles %>% filter(Accepted %in% R$Best_guess_binomial)
SynA <- Syn_Amphibians %>% filter(Accepted %in% A$Best_guess_binomial)


# # # # Script

# 1. Differences in species number before and after taxonomic corrections
ggarrange(
  Delta(Syn_Mammals, Syn_Amphibians, Syn_Birds, Syn_Reptiles)$p,
  Delta(SynM, SynA, SynB, SynR)$p, common.legend = TRUE
)

# 2. Identified synonyms: which species was the most replicated across datasets? 

# For mammals: Tachyoryctes splendens appeared under 12 different names
# across all datasets (and just for species in the trait datasets: same result)
Mammals <- MaxSyn(Syn_Mammals)
Birds <- MaxSyn(Syn_Birds)
Reptiles <- MaxSyn(Syn_Reptiles)
Amphibians <- MaxSyn(Syn_Amphibians)


# 3. Distribution of number of replicates across accepted names



## Plot
p1 <- Delta(SynM, SynA, SynB, SynR)$p + ylab("Binomial names") + xlab(NULL) + labs(tag = "A") + theme(plot.tag.position = c(0.7, 0.94))
p2 <- Dist_NR(Syn) + labs(tag = "B") + theme(plot.tag.position = c(0.7, 0.94))

p <- ggarrange(p1,p2)
p

