## Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}


## Function to calculate median trait completeness within families (trait DF is grouped by family, then median trait coverage is calculated)
Completeness_families <- function(TraitDF, Traits) {

  # completeness
  TraitDF$Percent <- apply(TraitDF[,Traits], 1, function(y) sum(!is.na(y)))
  TraitDF$Percent <- TraitDF$Percent / length(Traits) * 100 
  
  Completeness_0 <- TraitDF$Best_guess_binomial[TraitDF$Percent==0]
  print(length(Completeness_0))

  # group by families and summarise (median)
  Results <- TraitDF %>% 
    group_by(Family) %>%
    dplyr::summarise(Median=median(Percent)) %>% 
    as.data.frame()
  
  return(Results)

}

## Function to calculate coverage across species for a trait in each family (proportion of species within families with NA)

## % NA (Family) = n(sp NA in the family) / n(Sp total in the family)

PercentNA_families <- function(TraitDF, Trait) {
  
  TraitDF[, Trait] <- as.character(TraitDF[, Trait])
  TraitDF[is.na(TraitDF[, Trait]), Trait] <- "Missing"
  TraitDF[ TraitDF[, Trait]!="Missing", Trait] <- "Non_missing"
  MissVal <- table(TraitDF[, "Family"], TraitDF[, Trait]) %>%
    as.data.frame() %>%
    setNames(., c("Taxa", "Trait","nspecies"))
  
  Missing <- MissVal %>% 
    filter(Trait=="Missing") %>%
    select(Taxa, nspecies) %>%
    setNames(., c("Taxa", "n_sp_missing"))
  NonMissing<- MissVal %>% filter(Trait=="Non_missing") %>%
    select("nspecies") %>%
    setNames(., "n_sp_not_missing")
  
  PercentNA <- cbind(Missing, NonMissing) %>%
    mutate(n_tot_family=n_sp_missing+n_sp_not_missing) %>% 
    mutate(Percent_fam=n_sp_missing/n_tot_family*100)
  
  Results <- PercentNA %>%
    select(Taxa, Percent_fam) %>%
    setNames(., c("Family", "Percent"))
  
  return(Results)
  
}

## Proportion of species represented in each family (NspFAM/NspTOT)
Family_rep <- function(TraitDF, Log) {
  
  if(Log){
     TraitDF <- TraitDF %>% 
   group_by(Family) %>%
   summarise(Percent=n())%>%
   mutate(Percent=log10(Percent/nrow(TraitDF)*100))
  } 
  else{
    TraitDF <- TraitDF %>% 
      group_by(Family) %>%
      summarise(Percent=n())%>%
      mutate(Percent=Percent/nrow(TraitDF)*100)
  }

  return(TraitDF)
  
}


# ## Function to calculate coverage across species for a trait, then contribution of single species to the whole pattern
# 
# ## % NA  = n(sp NA in the family) / n(Sp total in the family) * (n(sp family)/n(Sp total)) * 100 
# 
# Coverage_families <- function(TraitDF, Trait) {
#   
#   # within families
#   TraitDF[, Trait] <- as.character(TraitDF[, Trait])
#   TraitDF[is.na(TraitDF[, Trait]), Trait] <- "Missing"
#   TraitDF[ TraitDF[, Trait]!="Missing", Trait] <- "Non_missing"
#   MissVal <- table(TraitDF[, "Family"], TraitDF[, Trait]) %>%
#     as.data.frame() %>%
#     setNames(., c("Taxa", "Trait","nspecies"))
#   
#   Missing <- MissVal %>% 
#     filter(Trait=="Missing") %>%
#     select(Taxa, nspecies) %>%
#     setNames(., c("Taxa", "n_sp_missing"))
#   NonMissing<- MissVal %>% filter(Trait=="Non_missing") %>%
#     select("nspecies") %>%
#     setNames(., "n_sp_not_missing")
#   
#   PercentNA <- cbind(Missing, NonMissing) %>%
#     mutate(n_tot_family=n_sp_missing+n_sp_not_missing) %>% 
#     mutate(n_tot=sum(n_sp_missing, n_sp_not_missing)) %>%
#     mutate(Percent_overall=n_sp_missing/n_tot*100)
#   
#   Results <- PercentNA %>%
#     select(Taxa, Percent_overall) %>%
#     setNames(., c("Family", "Percent"))
#   
#   return(Results)
#   
# }


## Function to plot the patterns

Plot_NA_patterns <- function(Phylo, ResultsDF, TraitDF, TipLabels, Comp) {
  
  colnames(ResultsDF)[2] <- "var"
  
  ## Not all capitalised
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  ResultsDF$Family <- tolower(ResultsDF$Family)
  ResultsDF$Family  <- firstup(ResultsDF$Family)
  TraitDF$Family <- tolower(TraitDF$Family)
  TraitDF$Family  <- firstup(TraitDF$Family)
  
  # Prune the tree to keep species that will represent families figuring in the result dataset
  
  # Families in the whole dataset
  Taxa <- TraitDF$Family %>% as.vector()
  Names <- TraitDF$Best_guess_binomial
  names(Taxa) <- Names
  Match <- match.phylo.data(Phylo, Taxa)
  Phylo <- Match$phy
  Taxa <- Match$data
  
  # For each family keep only one species in the tree
  Taxa <- as.data.frame(Taxa)
  Taxa$Species <- rownames(Taxa)
  Taxa <- Taxa %>%
    dplyr::group_by(Taxa) %>%
    dplyr::slice(1)
  
  # select families that are in the phylogenies in the result dataframe
  ResultsDF <- ResultsDF %>%
    filter(Family %in% Taxa$Taxa)
  
  # drop tips that are not the species that will be kept to represent each family
  Pattern <- paste(Taxa$Species, sep=" ")
  Phylo$tip.label[Phylo$tip.label %nin% Pattern] <- paste(Phylo$tip.label[Phylo$tip.label %nin% Pattern],"to_drop", sep=" ")
  Phylo <- drop.tip(Phylo, Phylo$tip.label[grepl("to_drop", Phylo$tip.label)])
  
  # paste family and order next to the corresponding tree tip
  for (i in 1:nrow(Taxa)) {
    Phylo$tip.label[Phylo$tip.label==Taxa$Species[i]] <- Taxa$Taxa[i] %>% as.character()
  }
  
  # Trait to plot is the variable in ResultsDF
  trait_toplot <- ResultsDF$var
  names(trait_toplot) <- ResultsDF$Family
  
  # add order information next to family information in both tip labels and names of trait_toplot
  
  for (i in 1:length(Phylo$tip.label)) {
    Y <- unique(as.character(TraitDF$Order[TraitDF$Family==Phylo$tip.label[i]]))
    Phylo$tip.label[i] <- paste(Y, Phylo$tip.label[i], sep=": ")}
  
  for (i in 1:length(trait_toplot)) {
    Y <- unique(as.character(TraitDF$Order[TraitDF$Family==names(trait_toplot)[i]]))
    names(trait_toplot)[i] <- paste(Y, names(trait_toplot)[i], sep=": ")
  }

  
  # Use contMap from phytools to plot
  tree.trait <- contMap(Phylo, trait_toplot, plot = FALSE, res=10)
  
  if(!Comp){
    tree.trait <-setMap(tree.trait, colors=c("blue","cyan", "green","yellow","red"))
  }
  
  else{
    tree.trait <-setMap(tree.trait, colors=c("red","yellow", "green","cyan","blue"))
  }
  
  if(!TipLabels) {
    plot.contMap(tree.trait,type = "phylogram", cex=0.5, legend=FALSE, ftype="off") 
    }
  
  else{
    plot.contMap(tree.trait,type = "phylogram", cex=0.5, legend=FALSE, fsize=0.7)
    }
  
  return(tree.trait)
}


