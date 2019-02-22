## Process Diet information from Kissling (Mammals) ("MammalDiet"), Amphibio (Amphibians) and Elton traits (Birds and mammals)
## Sekercioglu not used 

X <- c("data.table", "plyr", "dplyr", "tidyr", "magrittr", "reshape", "reshape2", "stringr", "stringi") 
invisible(lapply(X, library, character.only=TRUE)); rm(X)

`%nin%` = Negate(`%in%`)

# Sekercioglu <- read.csv("../../Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/Sekercioglu_Diet.csv")
Amphibio <- read.csv("../../Data/Amphibians/AmphiBIO_v1.csv")
Elton_birds <- read.csv("../../Data/Birds/EltonTraits_Birds.csv")
Elton_mammals <-  read.csv("../../Data/Mammals/EltonTraits_Mammals.csv")

## Put mammalian diet information together
Kissling <- read.csv("../../Data/Mammals/Kissling_Mammal_diet_2014.csv", sep=",")
Kissling$Binomial_name <- paste(Kissling$Genus, Kissling$Species, sep=" ")
Kissling <- subset(Kissling, Binomial_name != "Mico sp. nov.")
Kissling <- Kissling %>% select(-TaxonID, -TaxonomicNote, -FillCode)
Kissling <- Kissling[, c(28, 1:27)]

MammalDIET2 <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET2.csv")[, c(1:28)]
colnames(MammalDIET2)[28] <- "DataSource"

MammalDIET2_supp <- read.csv("../../Data/Mammals/MammalDIET_2/MammalDIET_2_supp.csv")
MammalDIET2_supp$TrophicLevel %<>% as.character()
MammalDIET2_supp <- MammalDIET2_supp %>%
  mutate(TrophicLevel=ifelse(TrophicLevel %in% c("Carnivore", "Omnivore", "Herbivore"), TrophicLevel, NA)) %>%
  filter(!is.na(TrophicLevel)) %>%
  mutate(Mammal=ifelse(Mammal %in% c(0:3), Mammal, NA)) %>%
  mutate(MammalEater=ifelse(MammalEater %in% c(0:3), MammalEater, NA)) %>%
  mutate(Insectivore=ifelse(Insectivore %in% c(0:3), Insectivore, NA))%>%
  mutate(Frugivore=ifelse(Frugivore %in% c(0:3), Frugivore, NA))%>%
  mutate(Granivore=ifelse(Granivore %in% c(0:3), Granivore, NA)) %>%
  mutate(Folivore=ifelse(Folivore %in% c(0:3), Folivore, NA))
MammalDIET2_supp[, c(6:27)] %<>% droplevels()

MammalDIET2 <- rbind(MammalDIET2, MammalDIET2_supp)
MammalDIET2$Binomial <- as.character(MammalDIET2$Binomial)
MammalDIET2$Binomial[MammalDIET2$Binomial=="Rhinolophus hildebrandtii"] <- "Rhinolophus hildebrandti"
MammalDIET2$Binomial[MammalDIET2$Binomial=="Tadarida bivittatus"] <- "Tadarida bivittata"
colnames(MammalDIET2)[1] <- "Binomial_name"

Y <- intersect(MammalDIET2$Binomial_name, Kissling$Binomial_name)
Kissling <- Kissling %>% filter(Binomial_name %nin% Y)
MammalDIET <- rbind(Kissling, MammalDIET2) # dataset to use


## Adopting a binary classification. The variables are: VE, IN, FR, SE, PL, NE


#---------------------------------------------------------------------------
# For Amphibio: just changing the column names to PL, NE, SE, FR, IN, VE 
colnames(Amphibio)[c(10:15)] <- c("PL", "NE", "SE", "FR", "IN", "VE") 
Amphibio$Diet_breadth <- apply(Amphibio[, c("PL", "NE", "SE", "FR", "IN", "VE")], 1, sum, na.rm=T)
Amphibio$Diet_breadth[Amphibio$Diet_breadth==0] <- NA

# Trophic levels
for (i in 1:nrow(Amphibio)) {
  
  if (any(Amphibio[i, c("PL", "NE", "SE", "FR")], na.rm=T)==T & any(Amphibio[i, c("IN", "VE")], na.rm=T)==F) {Amphibio$Trophic_level[i] <- "Herbivore"}
  if (any(Amphibio[i, c("PL", "NE", "SE", "FR")], na.rm=T)==F & any(Amphibio[i, c("IN", "VE")], na.rm=T)==T) {Amphibio$Trophic_level[i] <- "Carnivore"}
  if (any(Amphibio[i, c("PL", "NE", "SE", "FR")], na.rm=T)==T & any(Amphibio[i, c("IN", "VE")], na.rm=T)==T) {Amphibio$Trophic_level[i] <- "Omnivore"}
  if (any(Amphibio[i, c("PL", "NE", "SE", "FR")], na.rm=T)==F & any(Amphibio[i, c("IN", "VE")], na.rm=T)==F) {Amphibio$Trophic_level[i] <- NA}

  D <- c("PL", "NE", "SE", "FR", "IN", "VE")
  if (any(Amphibio[i, c("PL", "NE", "SE", "FR", "IN", "VE")], na.rm=T)==T) {
    ToPaste <- D[which(Amphibio[i, c("PL", "NE", "SE", "FR", "IN", "VE")]==1)]
    Amphibio$Primary_diet[i] <- paste(ToPaste, collapse = "|")} else {
      Amphibio$Primary_diet[i] <- NA  
    }

  }


colnames(Amphibio)[5] <- "Best_guess_binomial"
Amphibio$SCV <- NA

#-----------------------------------------------------------------------------
# # Sekercioglu
# # Add "or" in character strings where two elements are combined
# Sekercioglu$Primary_Diet %<>% as.character()
# ToChange <- Sekercioglu$Primary_Diet[nchar(Sekercioglu$Primary_Diet)==4]
# ToChange <- paste(substr(ToChange, 1, 2), "|", substr(ToChange, 3, 4), sep="") 
# Sekercioglu$Primary_Diet[nchar(Sekercioglu$Primary_Diet)==4] <- ToChange
# 
# # Diet as binary variables
# Diet <- c("FR", "NE", "SE", "PL", "IN", "VE", "SC", "OM")
# Sekercioglu[, Diet] <- NA
# 
# for (i in 1:nrow(Sekercioglu)) {
#   Item <- Sekercioglu$Primary_Diet[i]
#   Sekercioglu[i, c(which(grepl(Item, colnames(Sekercioglu))))] <- 1
# }
# 
# # For FI: in VE
# Sekercioglu$VE[grepl("FI", Sekercioglu$Primary_Diet)] <- 1
# Sekercioglu$VE[grepl("FI", Sekercioglu$Primary_Diet)] <- 1
# # For VG: in plant materials
# Sekercioglu$PL[grepl("VG", Sekercioglu$Primary_Diet)] <- 1
# 
# # Diet breadth as count
# Sekercioglu$Diet_breadth <- apply(Sekercioglu[, Diet], 1, sum, na.rm=T)
# Sekercioglu$Diet_breadth[Sekercioglu$Diet_breadth==0] <- NA
# 
# # Trophic level
# for (i in 1:nrow(Sekercioglu)) {
#   if (any(Sekercioglu[i, c("FR", "NE", "SE", "PL")], na.rm=T)==T & any(Sekercioglu[i, c("IN", "VE", "SC")], na.rm=T)==F) {Sekercioglu$Trophic_level[i] <- "Herbivore"}
#   if (any(Sekercioglu[i, c("FR", "NE", "SE", "PL")], na.rm=T)==F & any(Sekercioglu[i, c("IN", "VE","SC")], na.rm=T)==T) {Sekercioglu$Trophic_level[i] <- "Carnivore"}
#   if (any(Sekercioglu[i, c("FR", "NE", "SE", "PL")], na.rm=T)==T & any(Sekercioglu[i, c("IN", "VE", "SC")], na.rm=T)==T) {Sekercioglu$Trophic_level[i] <- "Omnivore"}
#   if (any(Sekercioglu[i, c("OM")], na.rm=T)==T) {Sekercioglu$Trophic_level[i] <- "Omnivore"}
#   if (any(Sekercioglu[i, c("FR", "NE", "SE", "PL","OM")], na.rm=T)==F & any(Sekercioglu[i, c("IN", "VE", "SC")], na.rm=T)==F) {Sekercioglu$Trophic_level[i] <- NA}
# }


#-----------------------------------------------------------------------------
# Kissling (MammalDIET)
Diet <- c("FR", "NE", "SE", "PL", "IN", "VE")

MammalDIET[, Diet] <- NA

for (i in 1:nrow(MammalDIET)) {

  if (any(MammalDIET[i, c("Vertebrate", "Mammal", "Bird", "Herptile","Fish")]==1, na.rm=T)) {MammalDIET$VE[i] <- 1}
  if (any(MammalDIET[i, "Invertebrate"]==1, na.rm=T)) {MammalDIET$IN[i] <- 1}
  if (any(MammalDIET[i, "Seed"]==1, na.rm=T)) {MammalDIET$SE[i] <- 1}
  if (any(MammalDIET[i, "Fruit"]==1, na.rm=T)) {MammalDIET$FR[i] <- 1}
  if (any(MammalDIET[i, "Nectar"]==1, na.rm=T)) {MammalDIET$NE[i] <- 1}
  if (any(MammalDIET[i, c("Root", "Leaf", "Woody", "Herbaceous")]==1, na.rm=T)) {MammalDIET$PL[i] <- 1}
  if (any(MammalDIET[i, "Plant"]==1 & all(MammalDIET[i,c("Seed", "Fruit", "Nectar")]==0), na.rm=T)) {MammalDIET$PL[i] <- 1} # this is an approximation
  
  if (any(MammalDIET[i, "MammalEater"]==1, na.rm=T)) {MammalDIET$VE[i] <- 1}
  if (any(MammalDIET[i, "Insectivore"]==1, na.rm=T)) {MammalDIET$IN[i] <- 1}
  if (any(MammalDIET[i, "Frugivore"]==1, na.rm=T)) {MammalDIET$FR[i] <- 1}
  if (any(MammalDIET[i, "Granivore"]==1, na.rm=T)) {MammalDIET$SE[i] <- 1}
  if (any(MammalDIET[i, "Folivore"]==1, na.rm=T)) {MammalDIET$PL[i] <- 1}
  

  if (any(MammalDIET[i, Diet], na.rm=TRUE)) {
    ToPaste <- D[which(MammalDIET[i, c("FR", "NE", "SE", "PL", "IN", "VE")]==1)]
    MammalDIET$Primary_diet[i] <- paste(ToPaste, collapse = "|")} else {
      MammalDIET$Primary_diet[i] <- NA  
    }
  }


# Diet breadth
MammalDIET$Diet_breadth <- apply(MammalDIET[, Diet], 1, sum, na.rm=T)
MammalDIET$Diet_breadth[MammalDIET$Diet_breadth==0] <- NA
unique(MammalDIET$Diet_breadth)

MammalDIET$SCV <- NA

#  -----------------------------------------------------------------------
# Elton trait diet
Process_Elton_diet <- function(Elton){
  
Diet <- c("FR", "NE", "SE", "PL", "IN", "VE", "SCV")

Elton[, Diet] <- NA

for (i in 1:nrow(Elton)) {
  
  if (any(Elton[i, c("Diet.Vunk", "Diet.Vend", "Diet.Vect", "Diet.Vfish")]!=0, na.rm=T)) {Elton$VE[i] <- 1}
  if (any(Elton[i, "Diet.Inv"]!=0, na.rm=T)) {Elton$IN[i] <- 1}
  if (any(Elton[i, "Diet.Scav"]!=0, na.rm=T)) {Elton$SCV[i] <- 1}
  if (any(Elton[i, "Diet.Fruit"]!=0, na.rm=T)) {Elton$FR[i] <- 1}
  if (any(Elton[i, "Diet.Nect"]!=0, na.rm=T)) {Elton$NE[i] <- 1}
  if (any(Elton[i, "Diet.Seed"]!=0, na.rm=T)) {Elton$SE[i] <- 1}
  if (any(Elton[i, "Diet.PlantO"]!=0, na.rm=T)) {Elton$PL[i] <- 1}
  
  if (any(Elton[i, Diet], na.rm=TRUE)) {
    ToPaste <- Diet[which(Elton[i, c("FR", "NE", "SE", "PL", "IN", "VE", "SCV")]==1)]
    Elton$Primary_diet[i] <- paste(ToPaste, collapse = "|")} else {
      Elton$Primary_diet[i] <- NA  
    }
  
  
Elton[, Diet] <- apply(Elton[, Diet], 2, as.numeric)
  
  # trophic levels
  if (any(Elton[i, c("PL", "NE", "SE", "FR")], na.rm=T)==T & any(Elton[i, c("IN", "VE", "SCV")], na.rm=T)==F) {Elton$Trophic_level[i] <- "Herbivore"}
  if (any(Elton[i, c("PL", "NE", "SE", "FR")], na.rm=T)==F & any(Elton[i, c("IN", "VE", "SCV")], na.rm=T)==T) {Elton$Trophic_level[i] <- "Carnivore"}
  if (any(Elton[i, c("PL", "NE", "SE", "FR")], na.rm=T)==T & any(Elton[i, c("IN", "VE", "SCV")], na.rm=T)==T) {Elton$Trophic_level[i] <- "Omnivore"}
  if (any(Elton[i, c("PL", "NE", "SE", "FR")], na.rm=T)==F & any(Elton[i, c("IN", "VE")], na.rm=T)==F) {Elton$Trophic_level[i] <- NA}
  
  print(i)
}

Elton$Diet_breadth <- apply(Elton[, Diet], 1, sum, na.rm=T)
Elton$Diet_breadth[Elton$Diet_breadth==0] <- NA

return(Elton)
}

Elton_birds <- Process_Elton_diet(Elton_birds)
Elton_mammals <- Process_Elton_diet(Elton_mammals)


# Function to transform diet variable into binary (used inside the Imputations_missForest function)
Diet_as_binary <- function(DF) {
  
  D <- DF %>% select(IN, VE, PL, SE, NE, FR, SCV)
  
  for (i in 1:nrow(D)) {
    
    if(any(D[i,]==1, na.rm = TRUE)) { 
      Col <- which(is.na(D[i,]))
      D[i, Col] <- 0}
    else{next()}
    print(i)
  }
  
  DF[, c("IN", "VE", "PL", "SE", "NE", "FR", "SCV")] <- D
  
  return(DF)
}

Amphibio <- Diet_as_binary(Amphibio)
MammalDIET <- Diet_as_binary(MammalDIET)
Elton_birds <- Diet_as_binary(Elton_birds)
Elton_mammals <- Diet_as_binary(Elton_mammals)


# Save files --------------------------------------------------------------
write.csv(Amphibio, "../../Data/Amphibians/Amphibio_processed_diet.csv", row.names = F)
#write.csv(Sekercioglu, "../../Results/0.Processed_diet_datasets/Sekercioglu_processed_diet.csv", row.names = F)
write.csv(MammalDIET, "../../Data/Mammals/MammalDIET_processed_diet.csv", row.names = F)
write.csv(Elton_birds, "../../Data/Birds/Elton_birds_processed.csv", row.names = F)
write.csv(Elton_mammals, "../../Data/Mammals/Elton_mammals_processed.csv", row.names = F)
