## Function to attach species that are not in the phylogeny (but that figure in the trait datasets) to their genus, at the root

Attach_Species_Phylo <- function(TraitDF, Phylo, Where) {

  
  N_original <- length(unique(Phylo$tip.label))
  
  TraitDF$Species <- paste(word(TraitDF$Best_guess_binomial,1), word(TraitDF$Best_guess_binomial,2), sep="_")
  Diff <- setdiff(TraitDF$Species, Phylo$tip.label)
  
  for (i in 1:length(Diff)) {
    print(paste("Iteration", i, "of", length(Diff)))
    Phylo <- tryCatch({expr=add.species.to.genus(Phylo, Diff[i], genus=NULL, where=Where)},
                      error=function(e) {Phylo}) 
  }
  
  NewDiff <- setdiff(TraitDF$Species, Phylo$tip.label)

  delta <-  length(unique(Phylo$tip.label))- N_original
  print(paste("Added", delta, "species on", length(Diff), "species that were not represented in the phylogeny"))
  
  return(Phylo)
}



## Function to resolve polytomies while adding branch length information to both edges that have 0 BL and their descendants
## Adapted from: Rangel TF, Colwell RK, Graves GR, Fucíková K, Rahbek C, & Diniz-Filho JAF (2015). Phylogenetic uncertainty revisited: 
#  Implications for ecological analyses. Evolution 69:1301-1312.

## Adapted to be ran with runs=1 only
bifurcatr <- function(phy,runs) {

  # trees <- vector("list", length=runs)
  
  # for(i in 1:runs) {
    
    tree <- phy
    resolves <- Ntip(tree)-Nnode(tree)-1
    for(j in 1:resolves) {
      descendent_counts <- rle(sort(tree$edge[,1]))
      polytomies <- descendent_counts$values[which(descendent_counts$lengths>2)] # these are parent nodes with more than 2 child nodes
      if(length(polytomies)>1) target_polytomy <- sample(polytomies,size=1) else target_polytomy <- polytomies
      polytomy_edges <- which(tree$edge[,1]==target_polytomy)
      target_edges <- sample(polytomy_edges,size=2)
      new_node <- max(tree$edge)+1
      tree$edge[target_edges,1] <- new_node
      new_edge <- c(target_polytomy,new_node)
      tree$edge <- rbind(tree$edge,new_edge)
      new_length <- runif(n=1,min=0,max= min(tree$edge.length[target_edges]))
      tree$edge.length <- c(tree$edge.length,new_length)
      tree$edge.length[target_edges] <- tree$edge.length[target_edges] - new_length
      tree$Nnode <- tree$Nnode+1
    }
    # trees[[i]] <- tree
  # }
  
  # if(runs==1) {
  #   trees <- trees[[1]]
  #   class(trees) <- "phylo"
  # }
  # else {
  #   class(trees) <- "multiPhylo"
  # }	
    
  return(tree)
}


## Function to look at how many species have been attached in the phylogeny
DeltaNspecies <- function(Phylo_Original, Phylo_Added, TraitDF) {
  
  Diff <- setdiff(TraitDF$Species, Phylo_Original$tip.label)
  
  n_original <- length(unique(Phylo_Original$tip.label))
  n_current <- length(unique(Phylo_Added$tip.label))
  
  delta_final <- setdiff(TraitDF$Species, Phylo_Added$tip.label)

  # percent species initially not represented in the tree
  print(paste(length(Diff)/nrow(TraitDF)*100, "species were initially not represented in tree (", length(Diff), "out of", nrow(TraitDF), ")."))
  print(paste(n_original-n_current, "species have been added on", length(Diff), "species that were not present in the phylogeny."))
  print(paste(length(delta_final)/nrow(TraitDF)*100, "species are finally not represented in tree (", length(delta_final), "out of", nrow(TraitDF), ")."))
  
}


## Function to drop replicated tips from the phylogeny
DropTips <- function(PhylogenyCor, PhylogenyUncor) {
  
  # convert phylogeny to a phylo4 object, to use with functions from phylobase package
  Phy <- as(PhylogenyCor,"phylo4") 
  
  # Duplicated tip labels
  Replicated <- PhylogenyCor$tip.label[duplicated(PhylogenyCor$tip.label)] %>% 
    as.data.frame() %>%
    setNames(., "Reps") %>%
    group_by(Reps) %>%
    summarise(Count=n() + 1) %>%
    as.data.frame() %>%
    mutate(Reps=as.character(Reps)) %>%
    mutate(Are_sister_species=NA) %>%
    mutate(RangePos=NA)
  
  ## For all replicates, check whether they are sister clades (sister species) in the phylogeny or not
  for(i in 1:nrow(Replicated)) {
    
    print(paste("Checking if sister species:", i, "on", nrow(Replicated)))
    
    if(wordcount(Replicated$Reps[i])>2) {
      Replicated$Reps[i] <- word(Replicated$Reps[i],1,2)
    }
    
    if(grepl("sp\\.|cf\\.|aff\\.", Replicated$Reps[i])) { next() } # those species will not match any species in the trait datasets; not need be considered
    
    else {
      
      # get the tip numbers of the replicates
      Nodes <- getNode(Phy, Replicated$Reps[i]) %>% 
        as.data.frame() %>%
        setNames(., "tip_number")
      
      Nodes$ancestral_node <- ancestor(Phy, Replicated$Reps[i])
    }
    
    Ancestral_nodes <- unique(Nodes$ancestral_node)
    
    if(length(Ancestral_nodes)==1) {Replicated$Are_sister_species[i] <- TRUE}
    else(Replicated$Are_sister_species[i] <- FALSE)
    
    # Range in replicate positions in the corrected tree -- just for information -- it is sensitive to the plotting
    # (that is, the difference between the most extreme positions)
    Replicated$RangePos[i] <- max(which(PhylogenyCor$tip.label==Replicated$Reps[i])) - min(which(PhylogenyCor$tip.label==Replicated$Reps[i]))
    
  }
  
  
  ## check if tip labels of replicates figure in the original, "uncorrected" tree for "other"
  Is_in_Uncor <- function(x) {
    
    if(any(grepl(x[1], PhylogenyUncor$tip.label))) {x[4] <- TRUE}
    else {x[4] <- FALSE}
    return(x[4])
    
  }
  
  Replicated$In_uncorrected_tree <- apply(Replicated, 1, Is_in_Uncor)
  
  Split <- split(Replicated, f=Replicated$Are_sister_species)
  Sistersp <- Split[["TRUE"]]
  Other <- Split[["FALSE"]]
  
  Other_ <- split(Other, f=Other$In_uncorrected_tree)
  Other_atrandom <- Other_[["FALSE"]]
  Drop_atrandom <- rbind(Sistersp, Other_atrandom)
  
  Drop_Not_atrandom <- Other_[["TRUE"]]
  
  
  ## for clades that are sister species, drop replicated tips randomly (choose randomly tips to drop) 
  ## + for species for which tip is not in the original tree
  
  for (j in Drop_atrandom$Reps) {
    
    print(paste("...at random:", j))
    
    # choose randomly the tip to drop among the 2 replicates, and drop
    ToDrop <- sample(which(PhylogenyCor$tip.label==j),1)
    # For these positions to drop, paste "to_drop" next to the species names.
    PhylogenyCor$tip.label[ToDrop] <- paste(PhylogenyCor$tip.label[ToDrop], "to_drop", sep="_")
    
  }
  
  
  ## for species that are not sister species, conserve the tip that is closest to the same tip in the original tree, drop other tips
  
  for (j in Drop_Not_atrandom$Reps) {
    
    print(paste("...according to position:", j))
    
    Cor_Pos <- which(PhylogenyCor$tip.label==j) %>%
      as.data.frame() %>%
      setNames(., "Pos.Cor") %>%
      mutate(Delta.Pos =
               abs(which(PhylogenyCor$tip.label==j) - which(PhylogenyUncor$tip.label==j)))
    
    # Get all positions that do not minimise the distance (those are the tips to drop)
    ToDrop <- Cor_Pos$Pos.Cor[Cor_Pos$Delta.Pos!=min(Cor_Pos$Delta.Pos)]
    
    if(length(ToDrop)==0) {
      ToDrop <- Cor_Pos$Pos.Cor[1] # random choice when similar distances
    }
    
    # For these positions to drop, paste "to_drop" next to the species names.
    PhylogenyCor$tip.label[ToDrop] <- paste(PhylogenyCor$tip.label[ToDrop], "to_drop", sep="_")
    
  }
  
  ## And finally drop
  print("Dropping tips")
  PhylogenyCor <- drop.tip(PhylogenyCor, tip=PhylogenyCor$tip.label[grepl("to_drop", PhylogenyCor$tip.label)])
  
  
  ## the "truly" problematic species are the ones for which sister species is false and in uncorrected is also false
  PbSpecies <- Replicated %>%
    filter(Are_sister_species==FALSE) %>%
    filter(In_uncorrected_tree==FALSE)
  
  return(list(PhylogenyCor=PhylogenyCor, Replicated=Replicated, PbSpecies=PbSpecies))
}


## Functions to extract and add phylogenetic eigenvectors from the phylogenies to the trait dataset as additional predictors
Extract_eigenvectors <- function(TraitDF, Phylo, N) {
  
  # N <- number of eigenvectors to extract
  
  ## Prune species that do not intersect
  row.names(TraitDF) <- TraitDF$Best_guess_binomial
  Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
  Phylo <- Prune_Taxa$phy
  
  ## Get phylogenetic eigenvectors from the phylogeny and select N first eigenvectors
  print("EIGENVECTOR DECOMPOSITION.")
  EigenV <- PVRdecomp(Phylo)
  
  Eigenvectors <- EigenV@Eigen$vectors
  Eigenvectors <- as.data.frame(Eigenvectors)
  Eigenvectors <- Eigenvectors[, 1:N]
  for (i in 1:10) {colnames(Eigenvectors)[i] <- paste("EV_",i, sep="")} 
  Eigenvectors$Best_guess_binomial <- Prune_Taxa$data$Best_guess_binomial
  Eigenvectors <- Eigenvectors[order(Eigenvectors$Best_guess_binomial), c(11, 1:10)]
  
  return(Eigenvectors)
  
}

Add_eigenvectors <- function (TraitDF, EV) {
  
  ColN <- vector()
  for (i in 1:10) {ColN <- c(ColN, paste("EV", i, sep="_")) }
  
  Species <- as.character(EV$Best_guess_binomial)
  x <- which(TraitDF$Best_guess_binomial %in% Species)
  TraitDF[, ColN] <- NA
  TraitDF[x, ColN] <- EV[, c(2:ncol(EV))]
  
  return(TraitDF)
  
}
