## Correlation among variables (on Predicts subsets)

library(PerformanceAnalytics)
library(dplyr)

TraitsAmphibian <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsAmphibians.Predicts_ContImputed_Randomadd.csv")
TraitsMammal <- read.csv("../Results/2.TraitsContinuous_after_imputations_Molina_Venegas/TraitsMammal.Predicts_ContImputed_Randomadd.csv")
source("../Code/1.Sourced_v.2/1.Trait_data_interpolation_v.2.R")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Amphibians
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


## Looking for correlations among continuous variables
Cont_traits <- c("Body_mass_g", "Litter_size", "Maturity_d", "Range_size_m2", "Habitat_breadth_IUCN")
Continuous_Amphibians <- TraitsAmphibian[, c("Family", "Best_guess_binomial", Cont_traits)]

Continuous_Amphibians_log <- log(Continuous_Amphibians[, Cont_traits])
chart.Correlation(Continuous_Amphibians_log, method="pearson")
cor(Continuous_Amphibians_log, use="na.or.complete", method = "pearson")

## Body_length / Litter size, correlation on the log with family as interacting factor
plot(Continuous_Amphibians_log$Litter_size ~ Continuous_Amphibians_log$Body_length_mm, pch=19,
     ylab="Litter size (log)", xlab="Body length mm (log)", main="Amphibians")
Predictions <- Corr_predictions(TraitsAmphibian, "Litter_size", "Body_length_mm", "Family")
Predicted <- Predictions$Traits
Summary <- Predictions$Summary
points(log(Predicted$Body_length_mm), log(Predicted$Litter_size))

## Body_length / Generation length
plot(Continuous_Amphibians_log$Generation_length_d ~ Continuous_Amphibians_log$Body_length_mm, pch=19,
     ylab="Generation length d (log)", xlab="Body length mm (log)", main="Amphibians")
LM <- lm(Generation_length_d ~ Body_length_mm, data=Continuous_Amphibians_log)
summary(LM)
abline(LM)

ToPredict <- Continuous_Amphibians_log$Body_length_mm[is.na(Continuous_Amphibians_log$Generation_length_d)] %>% as.data.frame()
colnames(ToPredict) <- "Body_length_mm"
ToPredict$Pred_Generation_length_mm_log <- predict(LM, newdata = ToPredict)

points(ToPredict$Body_length_mm, ToPredict$Pred_Generation_length_mm_log, col="red", pch=19)


pred_interval <- predict(LM, newdata=ToPredict, interval="prediction",
                         level = 0.95)
lines(ToPredict$Body_length_mm, pred_interval[,2], col="orange", lty=2)
lines(ToPredict$Body_length_mm, pred_interval[,3], col="orange", lty=2)


# # With Family as an interacting factor
# plot(Continuous_Amphibians_log$Generation_length_d ~ Continuous_Amphibians_log$Body_length_mm, pch=19,
#      ylab="Generation length d (log)", xlab="Body length mm (log)", main="Amphibians", col=as.factor(Continuous_Amphibians$Family))
# Predictions <- Corr_predictions(TraitsAmphibian, "Generation_length_d", "Body_length_mm", "Family")
# Predicted <- Predictions$Traits
# Summary <- Predictions$Summary
# points(log(Predicted$Body_length_mm), log(Predicted$Generation_length_d), col=as.factor(Continuous_Amphibians$Family))


## Adding to the main Trait dataset
TraitsAmphibian$Litter_size <- Predicted$Litter_size
TraitsAmphibian$Generation_length_d[is.na(TraitsAmphibian$Generation_length_d)] <- exp(ToPredict$Pred_Generation_length_mm_log)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Mammals
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Looking for correlations among continuous variables
Cont_traits <- c("Body_mass_g", "Litter_size", "Generation_length_d", "Range_size_m2", "Habitat_breadth_IUCN")
Continuous_Mammals <- TraitsMammal[, c("Order", "Family", "Best_guess_binomial", Cont_traits)]

Continuous_Mammals_log <- log(Continuous_Mammals[, Cont_traits])
chart.Correlation(Continuous_Mammals_log, method="pearson")
cor(Continuous_Mammals_log, use="na.or.complete", method = "pearson")

## Litter size -- generation length (only one NA)
plot(Continuous_Mammals_log$Litter_size ~ Continuous_Mammals_log$Generation_length_d, pch=19, col=as.factor(Continuous_Mammals$Order))

Species_to_pred <- TraitsMammal$Best_guess_binomial[is.na(TraitsMammal$Litter_size)]

Predictions <- Corr_predictions(TraitsMammal, "Litter_size", "Generation_length_d", "Order")
Predicted <- Predictions$Traits
Summary <- Predictions$Summary
points(log(Predicted$Generation_length_d[Predicted$Best_guess_binomial==Species_to_pred]), 
       log(Predicted$Litter_size[Predicted$Best_guess_binomial==Species_to_pred]), pch=17, col="black")

## Adding to the trait dataset
TraitsMammal$Litter_size <- Predicted$Litter_size

## Habitat breadth IUCN 
plot(Continuous_Mammals_log$Habitat_breadth_IUCN ~ Continuous_Mammals_log$Range_size_m2, pch=19)
Predictions <- Corr_predictions(TraitsMammal, "Habitat_breadth_IUCN", "Range_size_m2", "Family")
Predicted <- Predictions$Traits
Summary <- Predictions$Summary

Species_to_pred <- TraitsMammal$Best_guess_binomial[is.na(TraitsMammal$Habitat_breadth_IUCN)]
points(log(Predicted$Habitat_breadth_IUCN[Predicted$Best_guess_binomial %in% Species_to_pred]) ~ log(Predicted$Range_size_m2[Predicted$Best_guess_binomial %in% Species_to_pred]), pch=17, col="blue")

## Adding to the trait dataset
TraitsMammal$Habitat_breadth_IUCN <- Predicted$Habitat_breadth_IUCN



## Saving files
write.csv(TraitsAmphibian, "../Results/3.TraitsContinuous_after_correlations/TraitsAmphibians.Predicts_imputed_correlated.csv", row.names=F)
write.csv(TraitsMammal, "../Results/3.TraitsContinuous_after_correlations/TraitsMammals.Predicts_imputed_correlated.csv", row.names=F)






