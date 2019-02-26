Reprocess_diet <- function(Dataset) {
  
  Dataset <- Dataset %>% 
    mutate(Reprocessed_PD=ifelse(Primary_diet %in% c("PL", "SE", "PL|SE"), 1, 
                                 ifelse(Primary_diet %in% c("FR", "NE", "NE|FR"), 2, 
                                        ifelse(Primary_diet %in% c("VE", "SCV", "VE|SCV"), 3, 
                                               ifelse(Primary_diet %in% c("IN"), 4, 5)))))
  
  return(Dataset)
  
}

Transform_zscore <- function(TraitDF, Trait, Transf) {
  
  #browser()
  
  if (Transf=="log10"){
    TraitDF[,paste("log10", Trait, sep="_")] <- as.numeric(log10(TraitDF[,Trait]))
    TraitDF[, paste("log10", Trait, sep="_")] <- scale(TraitDF[, paste("log10", Trait, sep="_")], center=TRUE, scale=TRUE)
    TraitDF[, paste("log10", Trait, sep="_")] <- as.numeric(TraitDF[, paste("log10", Trait, sep="_")])
  }
  
  if(Transf=="sqrt") {
    TraitDF[,Trait]  <- sqrt(TraitDF[,Trait])
    TraitDF[, paste("sqrt", Trait, sep="_")] <- scale(TraitDF[,Trait] , center=TRUE, scale=TRUE)
    TraitDF[, paste("sqrt", Trait, sep="_")] <- as.numeric(TraitDF[, paste("sqrt", Trait, sep="_")])
  }
  
  return(TraitDF)
}

Compare <- function(RC_data, RC_imputed, AE_data, AE_imputed, 
                    RC_TraitName, AE_TraitName, Categorical, 
                    Traitaxis, ImputedaxisX, ImputedaxisY, Imputed,
                    Comparison, CompRC_collected, Diet) {

  # Standardise and scale RC data within each class
  if (RC_TraitName %in% c("body_mass_median", "litter_clutch_size")){
    RC_Birds <- subset(RC_data, class=="Aves")
    RC_Birds <- Transform_zscore(RC_Birds, RC_TraitName, "log10")
    RC_Mammals <- subset(RC_data, class=="Mammalia")
    RC_Mammals <- Transform_zscore(RC_Mammals, RC_TraitName, "log10")
    RC_data <- rbind(RC_Birds, RC_Mammals)
    RC_data <- RC_data[order(RC_data$binomial),]
    rm(RC_Birds, RC_Mammals)
    Original_RC_TraitName <- RC_TraitName
    RC_TraitName <- paste("log10", RC_TraitName, sep="_")
    
    if(!is.null(RC_imputed)){
    RC_Birds <- subset(RC_imputed, class=="Aves")
    RC_Birds <- Transform_zscore(RC_Birds, Original_RC_TraitName, "log10")
    RC_Mammals <- subset(RC_imputed, class=="Mammalia")
    RC_Mammals <- Transform_zscore(RC_Mammals, Original_RC_TraitName, "log10")
    RC_imputed <- rbind(RC_Birds, RC_Mammals)
    RC_imputed <- RC_imputed[order(RC_imputed$binomial),]
    rm(RC_Birds, RC_Mammals)}}
  
  if (RC_TraitName=="hab_breadth"){
    RC_Birds <- subset(RC_data, class=="Aves")
    RC_Birds <- Transform_zscore(RC_Birds, RC_TraitName, "sqrt")
    RC_Mammals <- subset(RC_data, class=="Mammalia")
    RC_Mammals <- Transform_zscore(RC_Mammals, RC_TraitName, "sqrt")
    RC_data <- rbind(RC_Birds, RC_Mammals)
    RC_data <- RC_data[order(RC_data$binomial),]
    
    if(!is.null(RC_imputed)){
    RC_Birds <- subset(RC_imputed, class=="Aves")
    RC_Birds <- Transform_zscore(RC_Birds, RC_TraitName, "sqrt")
    RC_Mammals <- subset(RC_imputed, class=="Mammalia")
    RC_Mammals <- Transform_zscore(RC_Mammals, RC_TraitName, "sqrt")
    RC_imputed <- rbind(RC_Birds, RC_Mammals)
    RC_imputed <- RC_imputed[order(RC_imputed$binomial),]
    
    rm(RC_Birds, RC_Mammals)
    RC_TraitName <- paste("sqrt", RC_TraitName, sep="_")}}
  
  if(Imputed) {
    
    # filter species for which the trait value was imputed
    
    # AE data
    Sp_AE_I <- AE_data$Best_guess_binomial[is.na(AE_data[,AE_TraitName])]
    
    # if the coverage for that trait was 100% initially, compared my collected values to Rob's imputed values (if they exist)
    if (length(Sp_AE_I)!=0) {
      AE_imputed <- AE_imputed %>%
        filter(Best_guess_binomial %in% Sp_AE_I)
      AE_data <- AE_imputed
    }
    else{ImputedaxisX <- "collected"}
    
    # RC data
    Sp_RC_I <- RC_data$binomial[is.na(RC_data[,RC_TraitName])]
    
    # if the coverage for that trait was 100% initially, compared my imputed values to Rob's collected values
    if (length(Sp_RC_I)!=0) {
      RC_imputed <- RC_imputed %>%
        filter(binomial %in% Sp_RC_I)
      RC_data <- RC_imputed
    }
    else{ImputedaxisY <- "collected"}
    
  }
  
  if(Comparison) {
    
    # RC collected VS AE imputed
    if(CompRC_collected){
      
      Sp_RC_C <- RC_data$binomial[!is.na(RC_data[,RC_TraitName])]
      Sp_AE_I <- AE_data$Best_guess_binomial[is.na(AE_data[,AE_TraitName])]
      Sp <- intersect(Sp_RC_C, Sp_AE_I)
      
      if(length(Sp)==0) {print("The data does not allow this comparison.")
        stop()} 
      
      else{
        RC_data <- RC_data %>% filter(binomial %in% Sp)
        RC_data <- RC_data[, c("class","binomial",RC_TraitName)]
        
        AE_data <- AE_imputed %>% filter(Best_guess_binomial %in% Sp)
        AE_data <- AE_data[, c("Class","Best_guess_binomial",AE_TraitName)]

        ImputedaxisX <- "imputed"
        ImputedaxisY <- "collected"
        }
      
    }
    
    # RC imputed VS AE collected
    if(!CompRC_collected){
      Sp_RC_I <- RC_data$binomial[is.na(RC_data[,RC_TraitName])]
      Sp_AE_C <- AE_data$Best_guess_binomial[!is.na(AE_data[,AE_TraitName])]
      Sp <- intersect(Sp_RC_I, Sp_AE_C)
      
      if(length(Sp)==0) {print("The data does not allow this comparison.")
        stop()} 
      
      else{
        
        RC_data <- RC_imputed %>% filter(binomial %in% Sp)
        RC_data <- RC_data[, c("class","binomial",RC_TraitName)]
        
        AE_data <- AE_data %>% filter(Best_guess_binomial %in% Sp)
        AE_data <- AE_data[, c("Class","Best_guess_binomial",AE_TraitName)]
        
        ImputedaxisX <- "collected"
        ImputedaxisY <- "imputed"
      }
      
      }
  }
  
  
  # intersect AE and RC species
  Y <- intersect(RC_data$binomial, AE_data$Best_guess_binomial)
  RC_data <- RC_data %>% filter(binomial %in% Y)
  AE_data <- AE_data %>% filter(Best_guess_binomial %in% Y)
  
  # plot
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  if(!Categorical) {
    
    ToPlot <- RC_data[, c("class","binomial", RC_TraitName)]
    colnames(ToPlot)[3] <- "RC"
    ToPlot$AE <- AE_data[, AE_TraitName]
    
    if(RC_TraitName %in% c("log10_litter_clutch_size", "log10_body_mass_median")){
      
      p <- ggplot(ToPlot, aes(AE, RC, color=ToPlot$class)) +
        GGPoptions +
        geom_point(alpha=0.7) +
        geom_abline(slope = 1, intercept = 0, alpha=0.8) +
        xlab(paste(Traitaxis, "(log10, ", ImputedaxisX," by AE)")) + 
        ylab(paste(Traitaxis,  "(log10, ", ImputedaxisY," by RC)")) + 
        scale_color_hue(name="Class", labels = c("Birds", "Mammals"))
      
      return(p)
    }
    
    else{
      
      p <- ggplot(ToPlot, aes(AE, RC, color=ToPlot$class)) +
        GGPoptions +
        geom_point(alpha=0.7) +
        geom_abline(slope = 1, intercept = 0, alpha=0.8) +
        xlab(paste("HB (square-root, ", ImputedaxisX," by AE)")) + 
        ylab(paste("HB (square-root, ",  ImputedaxisY, " by RC)")) + 
        scale_color_hue(name="Class", labels = c("Birds", "Mammals"))
      
      return(p)
    }
    
  }
  
  if(Categorical) {
    
    Outcome <- AE_data$Best_guess_binomial %>% as.data.frame()
    
    if(Diet) {
      x1 <- AE_data[,c("Best_guess_binomial",AE_TraitName)]
      x2 <- RC_data[, c("binomial",RC_TraitName)]
      xreturn <- cbind(x1,x2)
      xreturn <- xreturn %>%
        filter(Reprocessed_PD!=diet_5cat)
    }
    
    for (i in 1:nrow(AE_data)){
      
      x <- AE_data[i, AE_TraitName]
      y <- RC_data[i, RC_TraitName]
      
      if(is.na(x)|is.na(y)) {Outcome$Result[i] <- "unknown"}
      
      else{
        
        if(x==y) {Outcome$Result[i] <- "same"} 
        if(x!=y)  { Outcome$Result[i] <- "different"}
      }
    }
    
    ToPlot <- table(Outcome$Result) %>% 
      as.data.frame() %>%
      setNames(., c("outcome", "prop")) %>%
      mutate(prop=prop/nrow(Outcome)*100)
    
    ToPlot <- ToPlot[order(ToPlot$prop),]
    
    if(Comparison) {
      if(CompRC_collected) {ImputedaxisX <- "collected RC versus imputed AE"}
      if(!CompRC_collected) {ImputedaxisX <- "imputed RC versus collected AE"}
    }
    
    p <- ggplot(ToPlot, aes(outcome, prop)) +
      GGPoptions +
      geom_bar(stat="identity") +
      xlab(paste(Traitaxis, ImputedaxisX, sep = ", ")) + ylab("% species") +
      scale_x_discrete(limits=c("different", "unknown", "same"), labels=c("Contradicting", "Unkown", "Similar"))
    
    if(Diet) {return(list(p=p, outputs=xreturn))}
    else{return(p)}
  }
}

Plot.Cov <- function(TraitData, Traits, Main) {
  
  Names <- as.data.frame(c("log10_Body_mass_g", "body_mass_median" ,
                           "log10_Litter_size", "litter_clutch_size",
                           "sqrt_Habitat_breadth_IUCN","hab_breadth",
                           "Diel_activity", "activity"))
  colnames(Names) <- "Original"
  Names$FP <- c(rep("Body mass",2), rep("Litter/clutch size",2), rep("Habitat breadth",2),rep("Diel activity",2))
  
  Completeness <- apply(TraitData[, Traits], 2,  function(y) sum(!is.na(y))) 
  Completeness <- as.data.frame(Completeness/nrow(TraitData)*100)
  colnames(Completeness) <- "Completeness"
  Completeness <- Completeness[order(Completeness, decreasing=FALSE), , drop=FALSE]
  
  Names_plot <- as.data.frame(row.names(Completeness))
  colnames(Names_plot) <- "Or"
  Names$Original <- as.character(Names$Original)
  Names_plot$Or <- as.character(Names_plot$Or)
  for (i in 1:nrow(Names_plot)) {Names_plot$TP[i] <- Names$FP[Names$Original==Names_plot$Or[i]]}
  
  barplot(Completeness$Completeness, horiz = TRUE, 
          xlim = c(0,100), las=1,col="deepskyblue3", main = Main, names.arg = Names_plot$TP)
  
  abline(v=100, lty="dotted")
}





