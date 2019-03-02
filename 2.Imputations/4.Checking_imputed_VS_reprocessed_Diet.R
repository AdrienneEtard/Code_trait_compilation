# Checking imputed diet breadth versus reprocessed diet breadth

# Rationale: Diet was imputed on using binary variables corresponding to food categories. At the same time, diet breadth was 
# imputed on. For species for which diet breadth was imputed on, I recalculated DB as the sum of imputed binary diet variables.
# I am checking whether imputed DB == sum of imputed diet variables
# If not, use the "reprocessed diet breadth" column

# Preamble
source("Functions_for_comparison_with_RC_data.R")
library(dplyr)
library(ggplot2)

GGPoptions <- theme_classic() + theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=13, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=13))


# Load imputed data, selecting a dataset randomly
i <- sample(1:8,1)
AE_Imputed <- readRDS("../../Results/2.Imputed_trait_datasets/Imputed_not_standardised/List_of_8_sets.rds")[[i]]

M <- AE_Imputed$M$Imputed.Dataset[, c("Class","Best_guess_binomial", "Diet_breadth", "Diet_breadth_reprocessed")] %>%
  mutate(Class="Mammalia")
B <- AE_Imputed$B$Imputed.Dataset[, c("Class", "Best_guess_binomial", "Diet_breadth", "Diet_breadth_reprocessed")] %>%
  mutate(Class="Aves")
A <- AE_Imputed$B$Imputed.Dataset[, c("Class", "Best_guess_binomial", "Diet_breadth", "Diet_breadth_reprocessed")] %>%
  mutate(Class="Amphibia")
Imputed <- rbind(M, A, B)


# Load collected data
M <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv") %>% 
  select(Class, Best_guess_binomial, Diet_breadth)
B <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv") %>%
  select(Class, Best_guess_binomial, Diet_breadth)
A <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv") %>%
  select(Class, Best_guess_binomial, Diet_breadth)
Collected <- rbind(M, B,A)

## Get species for which diet breadth was inputed
X <- Collected$Best_guess_binomial[is.na(Collected$Diet_breadth)]
ImputedDB <- Imputed %>%
  filter(Best_guess_binomial %in% X)
ImputedDB$Class <- factor(ImputedDB$Class, levels = c("Aves", "Amphibia", "Mammalia"), labels = c("Birds", "Amphibians", "Mammals"))
p <- ggplot(ImputedDB, aes(Diet_breadth_reprocessed,Diet_breadth)) + 
         geom_point(size=2, alpha=0.3) +  
         GGPoptions + 
         xlab("Diet breadth, calculated from imputed food items \n(binary variables)") + 
         ylab("Diet breadth, imputed") +
  facet_wrap(~Class) + scale_x_continuous(breaks = c(0,1,2))
ggsave(p, file="../../Results/Plots/Congruence_imputations/Dietbreadth_imputed.pdf", width=4.5, height=2)




