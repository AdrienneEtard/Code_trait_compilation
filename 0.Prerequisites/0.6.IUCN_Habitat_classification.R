## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##                AIM: PREPARING HABITAT DATA FROM IUCN FILES FOR FURTHER TRAIT COMPILATION               ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

X <- c("grid", "gridExtra", "ggpubr", "dplyr", "lattice", "stringr")
invisible(lapply(X, library, character.only=TRUE)); rm(X)

IUCN_mammal <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_mammals.csv")
IUCN_amphibian <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_amphibians.csv")
IUCN_bird <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_birds.csv")
IUCN_reptile <- read.csv("../../../1.Trait_compilation/Results/0.Data_resolved_taxonomy/Processed_datasets/Traits/IUCN_habitat_reptiles.csv")


IUCN_amphibian$Best_guess_binomial %<>% as.character()
IUCN_mammal$Best_guess_binomial %<>% as.character()
IUCN_bird$Best_guess_binomial %<>% as.character()
IUCN_reptile$Best_guess_binomial %<>% as.character()


# # # # # # # # # # # # # # # # # # # # # #    F U N C T I O N S   # # # # # # # # # # # # # # # # # # # # # # 

## Functions to be ran successively in a script (further down)


# Habitat affinity: pooling by "code" using IUCN data, and pooling --------

.Set_H_type <- function(IUCN_data) {
  
  for (i in 1:nrow(IUCN_data)){
    IUCN_data$Code[i] <- unlist(strsplit(as.character(IUCN_data$code[i]), "[.]"))[1] %>% as.numeric() }
  
  Types <- c("Forest", "Savanna", "Shrubland", "Grassland", "Wetland", "Rocky areas", "Caves and subterranean", "Desert",
             rep("Marine", 3), rep("Marine intertidal or coastal/supratidal", 2), rep("Artificial", 2), "Introduced vegetation", rep("Other/Unknown", 2))
  
  for (i in 1:length(Types)) {
    if (!is.null(IUCN_data$Code[IUCN_data$Code==i])){
    IUCN_data$Affinity[IUCN_data$Code==i] <- Types[i]}
  }
  return(IUCN_data)
} 


# Calculating Habitat breadth and degree of specialisation  --------
# using weights so that suitable habitats count for more

IUCN_Habitat_calc <- function(IUCN_Habitat, w_suitable_important, w_suitable, w_marginal, w_unknown) {

Species <- unique(IUCN_Habitat$Best_guess_binomial) %>% as.data.frame()
colnames(Species) <- "Best_guess_binomial"

for (i in 1:nrow(Species)) {
  
  # Calculation of Habitat breadth
  s <- subset(IUCN_Habitat, Best_guess_binomial==Species$Best_guess_binomial[i])
  
  Breadth_suitable_important <- nrow(s[s$Suitability=="Suitable" & s$Major.importance=="Yes",]) * w_suitable_important
  Breadth_suitable <- nrow(s[s$Suitability=="Suitable" & (s$Major.importance=="No" | is.na(s$Major.importance)),]) * w_suitable
  
  Breadth_unknown <- nrow(s[s$Suitability=="Unknown",]) * w_unknown
  Breadth_NA <- nrow(s[is.na(s$Suitability),]) * w_unknown
  Breadth_empty <- nrow(s[s$Suitability=="",]) * w_unknown
  
  Breadth_marginal <- nrow(s[s$Suitability=="Marginal",]) * w_marginal
  
  Species$Habitat_breadth_IUCN[i] <- Breadth_suitable_important + Breadth_suitable + Breadth_unknown + Breadth_NA + Breadth_empty + Breadth_marginal 
  
  # Assigning a specialisation on Natural habitats (if no artificial habitats are suitable)
  if(any(s$Affinity=="Artificial" & s$Suitability=="Suitable", na.rm=T)==T) {
    Species$Specialisation[i] <- "Generalist"} else {
    Species$Specialisation[i] <- "Natural habitat specialist"}
}  
return(Species)
}


# Coding habitat affinity as a binary variable ----------------------------

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





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # #      R U N S      # # # # # # # # # # # # # # # # # # # # 


IUCN_amphibian %<>% .Set_H_type() 
Habitat_amphibian <- IUCN_Habitat_calc(IUCN_amphibian, 1, 0.8, 0.3, 0.8)
Habitat_amphibian <- Habitat_as_binary(Habitat_amphibian, IUCN_amphibian)

IUCN_bird %<>% .Set_H_type()
Habitat_bird <- IUCN_Habitat_calc(IUCN_bird, 1, 0.8, 0.3, 0.8)
Habitat_bird <- Habitat_as_binary(Habitat_bird, IUCN_bird)

IUCN_reptile %<>% .Set_H_type()
Habitat_reptile <- IUCN_Habitat_calc(IUCN_reptile, 1, 0.8, 0.3, 0.8)
Habitat_reptile <- Habitat_as_binary(Habitat_reptile, IUCN_reptile)

IUCN_mammal %<>% .Set_H_type()
Habitat_mammal <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0.3, 0.8)
Habitat_mammal <- Habitat_as_binary(Habitat_mammal, IUCN_mammal)

## Save files
write.csv(Habitat_amphibian, "../../Results/0.Processed_IUCN_Habitatdata/Amphibians.csv", row.names = FALSE)
write.csv(Habitat_bird, "../../Results/0.Processed_IUCN_Habitatdata/Birds.csv", row.names = FALSE)
write.csv(Habitat_reptile, "../../Results/0.Processed_IUCN_Habitatdata/Reptiles.csv", row.names = FALSE)
write.csv(Habitat_mammal, "../../Results/0.Processed_IUCN_Habitatdata/Mammals.csv", row.names = FALSE)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # #       S e n s i t i v i t y       # # # # # # # # # # # # # # # # # 

# ## Sensitivity of the calculation of habitat breadth to the weighting + plots
# 
# # With Amphibians; each set is a different set of weigths
# A_1 <- IUCN_Habitat_calc(IUCN_amphibian, 1, 0.8, 0.3, 0.8)
# A_2 <- IUCN_Habitat_calc(IUCN_amphibian, 1, 0.8, 0, 0)
# A_3 <- IUCN_Habitat_calc(IUCN_amphibian, 1, 1, 0, 0.5)
# A_4 <- IUCN_Habitat_calc(IUCN_amphibian, 1, 0.5, 0, 0)
# 
# 
# # Mammals
# M_1 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0.3, 0.8)
# M_2 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0, 0)
# M_3 <- IUCN_Habitat_calc(IUCN_mammal, 1, 1, 0, 0.5)
# M_4 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.5, 0, 0)
# 
# 
# # Birds
# B_1 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0.3, 0.8)
# B_2 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0, 0)
# B_3 <- IUCN_Habitat_calc(IUCN_mammal, 1, 1, 0, 0.5)
# B_4 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.5, 0, 0)
# 
# 
# # Reptiles
# R_1 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0.3, 0.8)
# R_2 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.8, 0, 0)
# R_3 <- IUCN_Habitat_calc(IUCN_mammal, 1, 1, 0, 0.5)
# R_4 <- IUCN_Habitat_calc(IUCN_mammal, 1, 0.5, 0, 0)
# 
# 
# 
# ## Plots
# pdf(file="../../Results/Plots/Sensitivity weighting/Sensitivity_breadth.pdf", width=6, height=4, family="Times", pointsize=11)
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,3), oma=c(1,1,1,1))
# 
# par(mfrow=c(2,3))
# 
# plot(density(A_1$Habitat_breadth_IUCN), ylim=c(0, 0.21), main="Amphibians", xlab="Habitat breadth calculated from IUCN", zero.line = FALSE)
# lines(density(A_2$Habitat_breadth_IUCN), col="red")
# lines(density(A_3$Habitat_breadth_IUCN), col="blue")
# lines(density(A_4$Habitat_breadth_IUCN), col="orange")
# 
# plot(density(B_1$Habitat_breadth_IUCN), ylim=c(0, 0.33), main="Birds", xlab="Habitat breadth calculated from IUCN", zero.line = FALSE)
# lines(density(B_2$Habitat_breadth_IUCN), col="red")
# lines(density(B_3$Habitat_breadth_IUCN), col="blue")
# lines(density(B_4$Habitat_breadth_IUCN), col="orange")
# 
# plot(density(M_1$Habitat_breadth_IUCN), ylim=c(0, 0.35), main="Mammals", xlab="Habitat breadth calculated from IUCN", zero.line = FALSE)
# lines(density(M_2$Habitat_breadth_IUCN), col="red")
# lines(density(M_3$Habitat_breadth_IUCN), col="blue")
# lines(density(M_4$Habitat_breadth_IUCN), col="orange")
# 
# plot(density(R_1$Habitat_breadth_IUCN), ylim=c(0, 0.35), main="Reptiles", xlab="Habitat breadth calculated from IUCN", zero.line = FALSE)
# lines(density(R_2$Habitat_breadth_IUCN), col="red")
# lines(density(R_3$Habitat_breadth_IUCN), col="blue")
# lines(density(R_4$Habitat_breadth_IUCN), col="orange")
# 
# plot(1, type="n", axes=F, xlab="", ylab="")
# par(xpd=T)
# legend(x=1, y=1.3, legend=c("w1", "w2", "w3", "w4"), col=c("black", "red", "blue", "orange"), lty = "solid")
# 
# dev.off()
# 
# 
# 
# ## Checking if species are shifting alongside the distribution
# 
# Get_sample <- function(Sample_size, D_1, D_2, D_3, D_4) {
#   
#   Sp <- sample(D_1$Best_guess_binomial, Sample_size)
#   Y <-  Sp %>% as.data.frame()
#   colnames(Y) <- "Binomial_name"
#   Y$Binomial_name %<>% as.character()
# 
#   Y$W1 <- D_1$Habitat_breadth_IUCN[D_1$Binomial_name %in% Sp]
#   Y$W2 <- D_2$Habitat_breadth_IUCN[D_2$Binomial_name %in% Sp]
#   Y$W3 <- D_3$Habitat_breadth_IUCN[D_3$Binomial_name %in% Sp]
#   Y$W4 <- D_4$Habitat_breadth_IUCN[D_4$Binomial_name %in% Sp]
# 
#   Y <- Y[order(Y$W1),]
#   
#   Y$Var <- apply(Y[, c("W1", "W2", "W3", "W4")], 1, var)
#   
#   return(Y)
#   
# }
# 
# Plot_sample <- function(Y, Main) {
# 
#   plot(Y$W1 ~ c(1:nrow(Y)), pch=19, las=1, type="n",
#        ylab = "Habitat breadth",
#        xlab="Random species sample", 
#        ylim=c(0, max(Y$W1, Y$W2, Y$W3, Y$W4)),
#        main=Main)
#   for (i in 1:nrow(Y)){abline(v=i, col="lightgrey")}
#   points(Y$W1 ~ c(1:nrow(Y)), pch=19)
#   points(Y$W2 ~ c(1:nrow(Y)), pch=19)
#   points(Y$W3 ~ c(1:nrow(Y)), pch=19)
#   points(Y$W4 ~ c(1:nrow(Y)), pch=19)
# 
# }
# 
# Plot_var <- function(Y, Main) {
#   plot(Y$Var ~ c(1:nrow(Y)), pch=19, las=1, type="n", xlab="Species sample", ylab="Variance of Habitat breadth")
#   for (i in 1:nrow(Y)){abline(v=i, col="lightgrey")}
#   points(Y$Var~ c(1:nrow(Y)), pch=19)
# }
# 
# 
# 
# Y_Amphibian <- Get_sample(100, A_1, A_2, A_3, A_4)
# Y_Mammal <- Get_sample(100, M_1, M_2, M_3, M_4)
# Y_Bird <- Get_sample(100, B_1, B_2, B_3, B_4)
# Y_Reptile <- Get_sample(100, R_1, R_2, R_3, R_4)
# 
# 
# pdf(file="../Results/Plots/Sensitivity weighting/Sensitivity_sample.pdf", width=6, height=4, family="Times", pointsize=11)
# 
# par(family='serif', tcl=0.2, cex.lab=1, mgp=c(1.5,0.2,0), mar = c(3,3,3,3), oma=c(1,1,1,1))
# par(mfrow=c(2, 2))
# 
# Plot_sample(Y_Amphibian, "Amphibian species")
# Plot_sample(Y_Mammal, "Mammal species")
# Plot_sample(Y_Bird, "Bird species")
# Plot_sample(Y_Reptile, "Reptile species")
# 
# dev.off()
# 
# 
# ## Plotting variances
# 
# Plot_sample_var <- function(Sample_size, D_1, D_2, D_3, D_4, Main) {
#   
#   Sp <- sample(D_1$Binomial_name, Sample_size)
#   Y <-  Sp %>% as.data.frame()
#   colnames(Y) <- "Binomial_name"
#   Y$Binomial_name %<>% as.character()
#   
#   Y$W1 <- D_1$Habitat_breadth_IUCN[D_1$Binomial_name %in% Sp]
#   Y$W2 <- D_2$Habitat_breadth_IUCN[D_2$Binomial_name %in% Sp]
#   Y$W3 <- D_3$Habitat_breadth_IUCN[D_3$Binomial_name %in% Sp]
#   Y$W4 <- D_4$Habitat_breadth_IUCN[D_4$Binomial_name %in% Sp]
#   
#   Y <- Y[order(Y$W1),]
#   
#   plot(Y$W1 ~ c(1:nrow(Y)), pch=19, las=1, type="n",
#        ylab = "Habitat breadth",
#        xlab="Random species sample", 
#        ylim=c(0, max(Y$W1, Y$W2, Y$W3, Y$W4)),
#        main=Main)
#   for (i in 1:nrow(Y)){abline(v=i, col="lightgrey")}
#   points(Y$W1 ~ c(1:nrow(Y)), pch=19)
#   points(Y$W2 ~ c(1:nrow(Y)), pch=19)
#   points(Y$W3 ~ c(1:nrow(Y)), pch=19)
#   points(Y$W4 ~ c(1:nrow(Y)), pch=19)
#   
# }

