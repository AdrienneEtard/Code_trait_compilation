# Checking imputed diet breadth versus reprocessed diet breadth

# Rationale: Diet was imputed on using binary variables corresponding to food categories. At the same time, diet breadth was 
# imputed on. For species for which diet breadth was imputed on, I recalculated DB as the sum of imputed binary diet variables.
# I am checking whether imputed DB == sum of imputed diet variables
# I not, use the "reprocessed diet breadth" column

i <- sample(1:8,1)
AE_Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_datasets/List_8_imputed_sets_v2.rds")[[i]]

Diet <- c("IN", "VE", "PL", "SE", "NE", "FR", "SCV")

M <- AE_Imputed$M$Imputed.Dataset[, c("Best_guess_binomial", "sqrt_Diet_breadth", "sqrt_Diet_breadth_reprocessed", Diet)] %>%
  mutate(Class="Mammalia")

B <- AE_Imputed$B$Imputed.Dataset[, c("Best_guess_binomial", "sqrt_Diet_breadth", "sqrt_Diet_breadth_reprocessed", Diet)] %>%
  mutate(Class="Aves")

A <- AE_Imputed$B$Imputed.Dataset[, c("Best_guess_binomial", "sqrt_Diet_breadth", "sqrt_Diet_breadth_reprocessed", Diet)] %>%
  mutate(Class="Amphibia")

