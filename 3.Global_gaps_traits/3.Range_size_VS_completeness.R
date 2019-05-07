## Relationship between global range size of species and number of traits sampled
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("../../1.Trait_compilation/Code/3.Global_gaps_traits/")

GGPoptions <- theme_classic()+ theme(
  panel.border = element_rect(colour = "black", fill=NA),
  text = element_text(size=12, family="serif"), 
  axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,3,0,"pt"), size=12), 
  axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,5,"pt"), size=12),
  axis.ticks.length=unit(-0.1, "cm"),
  legend.text=element_text(size=12)) 

## Load range size data and completeness

RS_mammals <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Mammals.csv") %>%
  select(Best_guess_binomial, Range_size_m2) 
RS_birds <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Birds.csv") %>%
  select(Best_guess_binomial, Range_size_m2)
RS_reptiles <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Reptiles.csv") %>%
  select(Best_guess_binomial, Range_size_m2)
RS_amphibians <- read.csv("../../Results/1.Traits_before_imputations/With_taxonomic_correction/All_species/3.with_phylo_eigenvectors/Amphibians.csv") %>%
  select(Best_guess_binomial, Range_size_m2)

Completeness_mammals <- read.csv("../../Results/3.Global_gaps_traits/Completeness_mammals.csv")
Completeness_birds <- read.csv("../../Results/3.Global_gaps_traits/Completeness_birds.csv")
Completeness_reptiles <- read.csv("../../Results/3.Global_gaps_traits/Completeness_reptiles.csv")
Completeness_amphibians <- read.csv("../../Results/3.Global_gaps_traits/Completeness_amphibians.csv")

Mammals <- merge(RS_mammals, Completeness_mammals) %>% mutate(Class="Mammals")
Birds <- merge(RS_birds, Completeness_birds)%>% mutate(Class="Birds")
Reptiles <- merge(RS_reptiles, Completeness_reptiles)%>% mutate(Class="Reptiles")
Amphibians <- merge(RS_amphibians, Completeness_amphibians)%>% mutate(Class="Amphibians")

Data <- rbind(Mammals, Birds, Reptiles, Amphibians)
Data$log_RS <- log(Data$Range_size_m2)
Data$N_predictors <- as.factor(Data$N_predictors)
Data$N_target <- as.factor(Data$N_target)

p1 <- ggplot(Data, aes(x=N_predictors,y=log_RS, fill=Class)) + 
  geom_boxplot(outlier.size =0.5, outlier.colour = "lightblue") +
  facet_wrap(~Class) + coord_flip() + GGPoptions +
  ylab(expression("Range size (log m"^2*")")) + xlab("Number of sampled traits")
  

# ## Poisson models for each class
Data <- Data %>% filter(!is.na(log_RS))
# Summary <- summary(model <- glm(as.numeric(N_predictors) ~ log_RS + Class, family="poisson", data=Data))
Summary <- summary(model <- glm(as.numeric(N_predictors) ~ log_RS + Class + log_RS:Class, family="poisson", data=Data))
saveRDS(Summary, "../../Results/3.Global_gaps_traits/Poisson_model_results.rds")
write.csv(Summary$coefficients, "../../Results/3.Global_gaps_traits/Poisson_model_coefficients.csv")

Model_pred <- predict(model, type="response")
p2 <- ggplot(Data, aes(x = log_RS, y = Model_pred, group=Class, col=Class)) +
  geom_point(aes(y = N_predictors), alpha=.2, size=0.1, position=position_jitter(h=.1), col="grey") +
  geom_line(size = 1) + GGPoptions + xlab(expression("Range size (log m"^2*")")) + ylab("Predicted numbers of known traits")

# Test for goodness of fit: should not be sgnificant if the model fits the data well
# estimate of overdispersion residual deviance/residual df (oversdispersion if this is significantly greater than 1)
model$deviance / model$df.residual
with(model, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))

# # Test for effect of class
# ## update model for a  model dropping class
# m0 <- update(model, . ~ . - Class)
m0 <- glm(as.numeric(N_predictors) ~ log_RS, family="poisson", data=Data)
# ## test model differences with chi square test
anova(m0, model, test="Chisq")

p <- ggarrange(p1+ theme(legend.title = element_blank()),p2, widths = c(0.5,0.5), labels = c("A", "B"), common.legend = TRUE, legend = "right")
ggsave(p, filename="../../Results/Plots/Trait_missing_values/Poisson_model_predictions.pdf", width=10, height=5)

plot(model)


