## Function to format phylogeny tip labels (from Genus_species to Genus species format)
.Format_tiplabels <- function (Phylogeny) {
  Phylogeny$tip.label <- gsub("_", " ", Phylogeny$tip.label)
  # Phylogeny$tip.label <- lapply(Phylogeny$tip.label, function(x) word(x, 1, 2)) %>% unlist()
  return(Phylogeny)
}


## Function to assess how many species in the trait dataset are represented in the phylogeny; 
## Randomly attach to their genus if possible. 
## How many species still left out of the phylogenies?
SpeciesPhylo <- function(Predicts, TraitDF, Phylo, Class, Where) {
  
  browser()
  
  PSp <- unique(Predicts$Best_guess_binomial[Predicts$Class==Class])
  PSp <- paste(word(PSp,1), word(PSp,2), sep="_")
  
  TraitDF$Species <- paste(word(TraitDF$Best_guess_binomial,1), word(TraitDF$Best_guess_binomial,2), sep="_")
  Diff <- setdiff(TraitDF$Species, Phylo$tip.label)
  
  PDiff <- setdiff(PSp, Phylo$tip.label)
  
  for (i in 1:length(Diff)) {
    print(paste("Iteration", i, "of", length(Diff)))
    Phylo <- tryCatch({expr=add.species.to.genus(Phylo, Diff[i], genus=NULL, where=Where)},
                      error=function(e) {Phylo}) 
  }
  
  PNewDiff <- setdiff(PSp, Phylo$tip.label)
  NewDiff <- setdiff(TraitDF$Species, Phylo$tip.label)
  
  print(paste("There were", length(PDiff), "PREDICTS species in", Class, "not represented in the phylogeny on ", length(PSp), "Predicts species."))
  print(paste("After attaching species, there are", length(PNewDiff), "PREDICTS species in", Class, "not represented in the phylogeny."))
  print(paste(length(PDiff)-length(PNewDiff),"species have been attached randomly."))
  print(paste("That is,",length(PNewDiff)/length(PSp)*100, "%  of PREDICTS species in the", Class, "trait dataset are finally not represented in the phylogeny."))
  
  
  print(paste("\n \nThere were", length(Diff), "species in the", Class, "trait dataset not represented in the phylogeny on", nrow(TraitDF), "species. Attaching these species to their genus if possible."))
  print(paste(length(Diff)-length(NewDiff), "have been attached randomly to their genus."))
  print(paste(length(NewDiff), "species in the", Class, "trait dataset are finally not represented in the phylogeny."))
  print(paste("That is,",length(NewDiff)/nrow(TraitDF)*100, "% species in the", Class, "trait dataset are finally not represented in the phylogeny."))
  
  return(list(Delta=NewDiff, Phylo=Phylo))
}

## Function to count the number of species replicated in corrected phylo that are also in PREDICTS
Rep_in_Predicts <- function(Predicts, Class, Dropped){
  Predicts.Sp <- unique(Predicts$Best_guess_binomial[Predicts$Class==Class])
  Reps <- gsub("_", " ", Dropped$Reps)
  print(paste(length(intersect(Predicts.Sp, Reps)), "species that have replicated tips are also in PREDICTS"))
  return(intersect(Predicts.Sp, Reps))
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
  
  
  ## for species that are sister species, drop all replicated tips randomly (choose randomly tips to drop) 
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

## Old version

# DropTips <- function(PhylogenyCor, PhylogenyUncor) {
#   
#   ToCheck <- vector()
#   
#   # Duplicated tip labels
#   Replicated <- PhylogenyCor$tip.label[duplicated(PhylogenyCor$tip.label)] %>% 
#     as.data.frame() %>%
#     setNames(., "Reps") %>%
#     group_by(Reps) %>%
#     summarise(Count=n() + 1) %>%
#     as.data.frame() %>%
#     mutate(Reps=as.character(Reps))
#   
#   Replicated$RangePos <- NA
#   
#   # Position of these tips on the trees (both corrected and uncorrected)
#   for(i in 1:nrow(Replicated)) {
#     
#     print(paste("dealing with replicate", i, "on", nrow(Replicated)))
#     
#     # For replicate i, get all positions in corrected and substract to position in uncorrected tree.
#     if(wordcount(Replicated$Reps[i])>2) {
#       Replicated$Reps[i] <- word(Replicated$Reps[i],1,2)
#     }
#     if(grepl("sp.|cf.|aff.", Replicated$Reps[i])) {next()} else{
#       
#       # First, assess the range in replicate positions in the corrected tree 
#       # (that is, the diffrence between the most extreme positions)
#       Replicated$RangePos[i] <- max(which(PhylogenyCor$tip.label==Replicated$Reps[i])) - min(which(PhylogenyCor$tip.label==Replicated$Reps[i]))
#       
#       # Then compare the positions, when possible, to position in the original tree
#       if(length(which(PhylogenyUncor$tip.label==Replicated$Reps[i]))==0){ # that is the case when no accepted name in the original phylogeny
#         ToCheck <- c(ToCheck, i)
#         
#         # # TODO all these have two replicates. what is the distance between the two replicates? hopefully all 1...
#         # Distances <- as.data.frame(ToCheck)
#         # colnames(Distances)
#         # which(PhylogenyCor$tip.label==Replicated$Reps[i])[1] - which(PhylogenyCor$tip.label==Replicated$Reps[i])[2]
#         
#         # choose randomly the tip to drop among the 2 replicates
#         ToDrop <- sample(which(PhylogenyCor$tip.label==Replicated$Reps[i]),1)
#         
#       } else {
#         
#         Cor_Pos <- which(PhylogenyCor$tip.label==Replicated$Reps[i]) %>%
#           as.data.frame() %>%
#           setNames(., "Pos.Cor") %>%
#           mutate(Delta.Pos = 
#                    abs(which(PhylogenyCor$tip.label==Replicated$Reps[i]) - which(PhylogenyUncor$tip.label==Replicated$Reps[i])))
#         
#         # Get all positions that do not minimise the distance (those are the tips to drop)
#         ToDrop <- Cor_Pos$Pos.Cor[Cor_Pos$Delta.Pos!=min(Cor_Pos$Delta.Pos)]
#       }
#     }
#     
#     # For these positions, paste "to_drop" next to the species names.
#     PhylogenyCor$tip.label[ToDrop] <- paste(PhylogenyCor$tip.label[ToDrop], "to_drop", sep="_")
#   }
#   
#   # Drop tips that have "to_drop"
#   PhylogenyCor <- drop.tip(PhylogenyCor, tip=PhylogenyCor$tip.label[grepl("to_drop", PhylogenyCor$tip.label)])
#   
#   return(list(PhylogenyCor=PhylogenyCor, ReplicatedTips=Replicated, ToCheck=ToCheck))
# }


# Add_eigenvectors <- function(TraitDF, Phylo, N, TaxInfo) {
#   
#   # N=number of eigenvectors
#   # Phylo=C_Phylo_Mammals
#   # TraitDF=C_Mammals
#   # N=10
#   
#   ## Prune species that do not intersect
#   row.names(TraitDF) <- TraitDF$Best_guess_binomial
#   Prune_Taxa <- match.phylo.data(Phylo, TraitDF)
#   Phylo <- Prune_Taxa$phy
#   ToBind1 <- Prune_Taxa$data
#   
#   ## Get remaining species, to be binded to species dataset with eigenvectors
#   Y <- setdiff(TraitDF$Best_guess_binomial, ToBind1$Best_guess_binomial)
#   ToBind2 <- TraitDF %>% filter(Best_guess_binomial %in% Y)
# 
#   if((nrow(ToBind2)+nrow(ToBind1))==nrow(TraitDF)){
#     
#     print("Good to start phylogenetic eigenvectors extraction.")
#   
#   ## Get phylogenetic eigenvectors from the phylogeny and select N first eigenvectors
#   print("EIGENVECTOR DECOMPOSITION.")
#   EigenV <- PVRdecomp(Phylo)
#   
#   Eigenvectors <- EigenV@Eigen$vectors
#   Eigenvectors <- as.data.frame(Eigenvectors)
#   Eigenvectors <- Eigenvectors[, 1:N]
#   for (i in 1:10) {colnames(Eigenvectors)[i] <- paste("EV_",i, sep="")} 
#   
#   ToBind1 <- cbind(ToBind1, Eigenvectors)
#   
#   # Set classes
#   if(TaxInfo) {
#     ToCharacter <- c("Order", "Family", "Genus", "Best_guess_binomial", "Diel_activity", "Trophic_level", "Specialisation", "Primary_diet") 
#   }
#   
#   else (ToCharacter <- c("Best_guess_binomial", "Diel_activity", "Trophic_level", "Specialisation", "Primary_diet"))
#   
#   ToBind1[, ToCharacter] <- apply(ToBind1[, ToCharacter], 2, as.character)
#   ToBind2[, ToCharacter] <- apply(ToBind2[, ToCharacter], 2, as.character)
#   
#   ToNumeric <- setdiff(colnames(ToBind2), ToCharacter)
#   ToBind1[, ToNumeric] <- apply(ToBind1[, ToNumeric], 2, as.numeric)
#   ToBind2[, ToNumeric] <- apply(ToBind2[, ToNumeric], 2, as.numeric)
#   
#   FinalDF <- dplyr::bind_rows(ToBind1, ToBind2)
#   FinalDF <- FinalDF[order(FinalDF$Best_guess_binomial),]
#   rownames(FinalDF) <- c(1:nrow(FinalDF))
#   
#   return(FinalDF)
#   
#   } else (print("Manual check for pseudoreplication in phylogenetic tips required."))
#   
# }


## Functions to extract and add eigenvectors to the phylogenies



