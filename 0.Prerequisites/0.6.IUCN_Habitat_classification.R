## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                AIM: PREPARING HABITAT DATA FROM IUCN FILES FOR FURTHER TRAIT COMPILATION               ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

X <- c("grid", "gridExtra", "ggpubr", "dplyr", "lattice", "stringr")
invisible(lapply(X, library, character.only=TRUE)); rm(X)

## Data corrected for taxonomy
IUCN_mammal_C <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_mammals.csv")
IUCN_amphibian_C <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_amphibians.csv")
IUCN_bird_C <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_birds.csv")
IUCN_reptile_C <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_reptiles.csv")

# all <- rbind(IUCN_amphibian_C, IUCN_bird_C, IUCN_mammal_C, IUCN_reptile_C)
# 
# all$Major.importance <- as.character(all$Major.importance)
# all$Major.importance[is.na(all$Major.importance)] <- "Unknown"
# 
# all$Suitability <- as.character(all$Suitability)
# all$Suitability[is.na(all$Suitability)] <- "Unknown"
# 
# table( all$Suitability, all$Major.importance)

IUCN_amphibian_C$Best_guess_binomial %<>% as.character()
IUCN_mammal_C$Best_guess_binomial %<>% as.character()
IUCN_bird_C$Best_guess_binomial %<>% as.character()
IUCN_reptile_C$Best_guess_binomial %<>% as.character()

## Data uncorrected for taxonomy
IUCN_mammal_UN <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Mammals_20180726.csv", sep=",")
IUCN_amphibian_UN <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Amphibians_20180726.csv", sep=",")
IUCN_bird_UN <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Birds_20180726.csv", sep=",")
IUCN_reptile_UN <- read.csv("../../Data/HabitatData_09_08_18_IUCN_unprocessed/API_HabitatLevel2_Reptiles_20180726.csv", sep=",")

colnames(IUCN_mammal_UN)[7] <- "Best_guess_binomial"
colnames(IUCN_amphibian_UN)[7] <- "Best_guess_binomial"
colnames(IUCN_reptile_UN)[7] <- "Best_guess_binomial"
colnames(IUCN_bird_UN)[7] <- "Best_guess_binomial"

IUCN_amphibian_UN$Best_guess_binomial %<>% as.character()
IUCN_mammal_UN$Best_guess_binomial %<>% as.character()
IUCN_bird_UN$Best_guess_binomial %<>% as.character()
IUCN_reptile_UN$Best_guess_binomial %<>% as.character()


# # # # # # # # # # # # # # # # # # # # # #    F U N C T I O N S   # # # # # # # # # # # # # # # # # # # # # # 

## HABITAT AFFINITY: pooling by "code" using IUCN data, and pooling into broader habitat categories

.Set_H_type <- function(IUCN_data) {
  
  Func <- function(X) { 
    x <- strsplit(as.character(X), "[.]") %>% unlist
    x <- x[1] %>% as.numeric(); return(x)}
  
  IUCN_data$Code <- sapply(IUCN_data$code, Func)
  
  IUCN_data <- IUCN_data %>% mutate(Affinity=ifelse(Code==1, "Forest",
                                                    ifelse(Code==2, "Savanna", 
                                                           ifelse(Code==3, "Shrubland", 
                                                                  ifelse(Code==4, "Grassland", 
                                                                         ifelse(Code==5, "Wetland", 
                                                                                ifelse(Code==6, "Rocky areas", 
                                                                                       ifelse(Code==7, "Caves and subterranean", 
                                                                                              ifelse(Code==8, "Desert", 
                                                                                                     ifelse(Code==9|Code==10|Code==11, "Marine", 
                                                                                                            ifelse(Code==12|Code==13, "Marine intertidal or coastal/supratidal", 
                                                                                                                   ifelse(Code==14|Code==15, "Artificial", 
                                                                                                                          ifelse(Code==16, "Introduced vegetation", 
                                                                                                                                 ifelse(Code==17|Code==18, "Other/unknown", NA))))))))))))))
  

  return(IUCN_data)
} 

## HABITAT BREADTH AND DEGREE OF SPECIALISATION
## using weights so that suitable habitats count for more

IUCN_Habitat_calc <- function(IUCN_Habitat, With_weights) {

  IUCN_Habitat <- IUCN_Habitat %>% 
    mutate(Suitability=as.character(Suitability)) %>%
    mutate(Suitability=ifelse(is.na(Suitability), "Unknown", Suitability))
  
  IUCN_Habitat <- IUCN_Habitat %>%
    mutate(Major.importance=as.character(Major.importance)) %>%
    mutate(Major.importance=ifelse(is.na(Major.importance), "Unknown", Major.importance))
  
  Species <- unique(IUCN_Habitat$Best_guess_binomial) %>%
    as.data.frame() %>%
    setNames(., "Best_guess_binomial") %>%
    mutate(Habitat_breadth_IUCN=NA) %>%
    mutate(Specialisation=NA)
    

  for (i in 1:nrow(Species)) {
    
  # Calculation of Habitat breadth
  s <- subset(IUCN_Habitat, Best_guess_binomial==Species$Best_guess_binomial[i])
  
  tbl <- table(s$Suitability, s$Major.importance) %>% 
    as.data.frame() %>%
    setNames(., c("Suitability", "Major.importance", "N"))
  
  # if(With_weights){
  #     tbl <- tbl %>% 
  #   mutate(Weight=ifelse((Suitability=="Suitable"| Suitability=="Unknown") & (Major.importance=="Yes"|Major.importance=="Unknown"), 1, 
  #                        ifelse(Suitability=="Suitable" & Major.importance=="No", 0.5, 0.3))) %>%
  #   mutate(Score=N*Weight)
  # }
  
  if(With_weights){

    tbl <- tbl %>% 
      mutate(Weight=ifelse((Suitability=="Suitable") & (Major.importance=="Yes"), 1, 
                           ifelse(Suitability=="Suitable" & (Major.importance=="No"|Major.importance=="Unknown"), 0.5,
                                  ifelse(Suitability=="Unknown" & Major.importance=="Unknown", 0.5, 0.3)))) %>%
      mutate(Score=N*Weight)
    

  }
  
  if(!With_weights) {
    tbl <- tbl %>% 
      mutate(Weight=1) %>%
      mutate(Score=N*Weight)
  }
  
  
  Species$Habitat_breadth_IUCN[i] <- sum(tbl$Score)
  

   
  # Assigning a specialisation on Natural habitats (if no artificial habitats are suitable)
  if(any(s$Affinity=="Artificial" & (s$Suitability=="Suitable"|s$Suitability=="Unknown"), na.rm=TRUE)) {
    Species$Specialisation[i] <- "Generalist"} 
  
  else{
      if(any(s$Affinity=="Other/unknown" & (s$Suitability=="Suitable"|s$Suitability=="Unknown"), na.rm=TRUE)) {
      Species$Specialisation[i] <- NA} 
    else {Species$Specialisation[i] <- "Natural habitat specialist"}
  }
}  
  return(Species)
}

## HABITAT AFFINITY AS A BINARY VARIABLE

Habitat_as_binary <-  function(Habitat_data, IUCN_data){
  
  Types <- c("Forest", "Savanna", "Shrubland", "Grassland", "Wetland", "Rocky areas", "Caves and subterranean", "Desert",
             "Marine", "Marine intertidal or coastal/supratidal", "Artificial", "Introduced vegetation", "Other/Unknown")
  
  Habitat_data[, Types] <- NA

  
  for (i in 1:nrow(Habitat_data)) {
    
    s <- subset(IUCN_data, Best_guess_binomial==Habitat_data$Best_guess_binomial[i])
    
    for (j in 1:length(Types)) {
     
      if (any(s$Affinity==Types[j])) {
        Habitat_data[i, Types[j]] <- 1 
      } 
      else {
        Habitat_data[i, Types[j]] <- 0
        }
    }
  
  }
  return(Habitat_data)
}


# # # # # # # # # # # # # # # # # # # # # # #      R U N S      # # # # # # # # # # # # # # # # # # # # 

## Amphibians
IUCN_amphibian_C %<>% .Set_H_type() 
Habitat_amphibian_C <- IUCN_Habitat_calc(IUCN_amphibian_C, TRUE)
Habitat_amphibian_C <- Habitat_as_binary(Habitat_amphibian_C, IUCN_amphibian_C)

IUCN_amphibian_UN %<>% .Set_H_type() 
Habitat_amphibian_UN <- IUCN_Habitat_calc(IUCN_amphibian_UN, TRUE)
Habitat_amphibian_UN <- Habitat_as_binary(Habitat_amphibian_UN, IUCN_amphibian_UN)


## Birds
IUCN_bird_C %<>% .Set_H_type()
Habitat_bird_C <- IUCN_Habitat_calc(IUCN_bird_C, TRUE)
Habitat_bird_C <- Habitat_as_binary(Habitat_bird_C, IUCN_bird_C)

IUCN_bird_UN %<>% .Set_H_type()
Habitat_bird_UN <- IUCN_Habitat_calc(IUCN_bird_UN, TRUE)
Habitat_bird_UN <- Habitat_as_binary(Habitat_bird_UN, IUCN_bird_UN)

## Reptiles
IUCN_reptile_C %<>% .Set_H_type()
Habitat_reptile_C <- IUCN_Habitat_calc(IUCN_reptile_C, TRUE)
Habitat_reptile_C <- Habitat_as_binary(Habitat_reptile_C, IUCN_reptile_C)

IUCN_reptile_UN %<>% .Set_H_type()
Habitat_reptile_UN <- IUCN_Habitat_calc(IUCN_reptile_UN, TRUE)
Habitat_reptile_UN <- Habitat_as_binary(Habitat_reptile_UN, IUCN_reptile_UN)


## Mammals
IUCN_mammal_C %<>% .Set_H_type()
Habitat_mammal_C <- IUCN_Habitat_calc(IUCN_mammal_C, TRUE)
Habitat_mammal_C <- Habitat_as_binary(Habitat_mammal_C, IUCN_mammal_C)

IUCN_mammal_UN %<>% .Set_H_type()
Habitat_mammal_UN <- IUCN_Habitat_calc(IUCN_mammal_UN, TRUE)
Habitat_mammal_UN <- Habitat_as_binary(Habitat_mammal_UN, IUCN_mammal_UN)


## Without weights
Habitat_amphibian_nw <- IUCN_Habitat_calc(IUCN_amphibian_C, FALSE)
Habitat_bird_nw <- IUCN_Habitat_calc(IUCN_bird_C, FALSE)
Habitat_mammal_nw <- IUCN_Habitat_calc(IUCN_mammal_C, FALSE)
Habitat_reptile_nw <- IUCN_Habitat_calc(IUCN_reptile_C, FALSE)


## Save files
write.csv(Habitat_amphibian_C, "../../Results/0.Processed_IUCN_Habitatdata/Amphibians_v2.csv", row.names = FALSE)
write.csv(Habitat_bird_C, "../../Results/0.Processed_IUCN_Habitatdata/Birds_v2.csv", row.names = FALSE)
write.csv(Habitat_reptile_C, "../../Results/0.Processed_IUCN_Habitatdata/Reptiles_v2.csv", row.names = FALSE)
write.csv(Habitat_mammal_C, "../../Results/0.Processed_IUCN_Habitatdata/Mammals_v2.csv", row.names = FALSE)

write.csv(Habitat_amphibian_nw, "../../Results/0.Processed_IUCN_Habitatdata/Amphibians_UNWEIGHTED.csv", row.names = FALSE)
write.csv(Habitat_bird_nw, "../../Results/0.Processed_IUCN_Habitatdata/Birds_UNWEIGHTED.csv", row.names = FALSE)
write.csv(Habitat_reptile_nw, "../../Results/0.Processed_IUCN_Habitatdata/Reptiles_UNWEIGHTED.csv", row.names = FALSE)
write.csv(Habitat_mammal_nw, "../../Results/0.Processed_IUCN_Habitatdata/Mammals_UNWEIGHTED.csv", row.names = FALSE)

write.csv(Habitat_amphibian_UN, "../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Amphibians_v2.csv", row.names = FALSE)
write.csv(Habitat_bird_UN, "../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Birds_v2.csv", row.names = FALSE)
write.csv(Habitat_reptile_UN, "../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Reptiles_v2.csv", row.names = FALSE)
write.csv(Habitat_mammal_UN, "../../Results/0.Processed_IUCN_Habitatdata/No_taxonomic_correction/Mammals_v2.csv", row.names = FALSE)


## Looking at the distributions with and without weights


CompareDist <- function(DataW, DataNW) {
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=13), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=13),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  
  DataW <- DataW[order(DataW$Best_guess_binomial), ]
  DataNW <- DataNW[order(DataNW$Best_guess_binomial), ]
  
  DataW <- DataW %>% 
    select(Habitat_breadth_IUCN) %>%
    setNames(., "weighted")
  
  DataNW <- DataNW %>% 
    select(Habitat_breadth_IUCN) %>%
    setNames(., "not weighted")
  
  Data <- cbind(DataW, DataNW) %>%
    reshape::melt()
  
  p <- ggplot(Data, aes(value, fill=variable)) +
  geom_density(alpha=0.5, adjust=1.5) + GGPoptions +
  xlab("Habitat breadth") + ylab("Density") +
    scale_fill_hue(name="Sum of habitats")
  
  return(p)

  
}

pAmphibians <- CompareDist(Habitat_amphibian_C, Habitat_amphibian_nw) + 
  labs(tag = "D") + theme(plot.tag.position = c(0.9, 0.92))
pBirds <- CompareDist(Habitat_bird_C, Habitat_bird_nw)+ 
  labs(tag = "B") + theme(plot.tag.position = c(0.9, 0.92))
pReptiles <- CompareDist(Habitat_reptile_C, Habitat_reptile_nw)+ 
  labs(tag = "C") + theme(plot.tag.position = c(0.9, 0.92))
pMammals <- CompareDist(Habitat_mammal_C, Habitat_mammal_nw)+ 
  labs(tag = "A") + theme(plot.tag.position = c(0.9, 0.92))
p <- ggarrange(pMammals, pBirds, pReptiles, pAmphibians, ncol=2, nrow=2, common.legend = TRUE)
ggsave(p, file="../../Results/Plots/Weighted_HB/distributions.pdf", width = 6, height = 6)



table(Habitat_amphibian_C$Specialisation)
table(Habitat_bird_C$Specialisation)
table(Habitat_mammal_C$Specialisation)
table(Habitat_reptile_C$Specialisation)


