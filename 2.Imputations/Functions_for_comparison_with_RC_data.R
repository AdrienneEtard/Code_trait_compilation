Reprocess_diet <- function(Dataset) {

  Dataset <- Dataset %>%
    mutate(Primary_diet=ifelse(Primary_diet %in% c("PL|SE"), 1,
                                 ifelse(Primary_diet %in% c("FR|NE"), 2,
                                        ifelse(Primary_diet %in% "VE", 3,
                                               ifelse(Primary_diet %in% "IN", 4, 
                                                      ifelse(Primary_diet %in% "OM", 5, NA)))))) %>%
    mutate(Primary_diet=as.character(Primary_diet))

  return(Dataset)

}



Compare <- function(RC_data, RC_imputed, AE_data, AE_imputed, Imputed,
                    TraitName, Categorical, 
                    Traitaxis, ImputedaxisX, ImputedaxisY,
                    Comparison, CompRC_collected, Diet) {
  
  
  if(Imputed) {

    
    # Filter species for which the trait value was imputed
    
    # AE data
    Sp_AE_I <- AE_data$Best_guess_binomial[is.na(AE_data[,TraitName])]
    
    # if the coverage for that trait was 100% initially, compared my collected values to Rob's imputed values (if they exist)
    if (length(Sp_AE_I)!=0) {
      AE_imputed <- AE_imputed %>%
        filter(Best_guess_binomial %in% Sp_AE_I)
      AE_data <- AE_imputed
    }
    else{ImputedaxisX <- "collected"}
    
    # RC data
    Sp_RC_I <- RC_data$binomial[is.na(RC_data[,TraitName])]
    
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
      
      Sp_RC_C <- RC_data$binomial[!is.na(RC_data[,TraitName])]
      Sp_AE_I <- AE_data$Best_guess_binomial[is.na(AE_data[,TraitName])]
      Sp <- intersect(Sp_RC_C, Sp_AE_I)
      
      if(length(Sp)==0) {print("The data does not allow this comparison.")
        stop()} 
      
      else{
        RC_data <- RC_data %>% filter(binomial %in% Sp)
        RC_data <- RC_data[, c("class","binomial",TraitName)]
        
        AE_data <- AE_imputed %>% filter(Best_guess_binomial %in% Sp)
        AE_data <- AE_data[, c("Class","Best_guess_binomial",TraitName)]

        ImputedaxisX <- "imputed"
        ImputedaxisY <- "collected"
        }
      
    }
    
    # RC imputed VS AE collected
    if(!CompRC_collected){
      Sp_RC_I <- RC_data$binomial[is.na(RC_data[,TraitName])]
      Sp_AE_C <- AE_data$Best_guess_binomial[!is.na(AE_data[,TraitName])]
      Sp <- intersect(Sp_RC_I, Sp_AE_C)
      
      if(length(Sp)==0) {print("The data does not allow this comparison.")
        stop()} 
      
      else{
        
        RC_data <- RC_imputed %>% filter(binomial %in% Sp)
        RC_data <- RC_data[, c("class","binomial",TraitName)]
        
        AE_data <- AE_data %>% filter(Best_guess_binomial %in% Sp)
        AE_data <- AE_data[, c("Class","Best_guess_binomial",TraitName)]
        
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
    
    ToPlot <- RC_data[, c("class","binomial", TraitName)]
    colnames(ToPlot)[3] <- "RC"
    ToPlot$AE <- AE_data[, TraitName]
    
    # Ssize <- paste("N=",length(unique(ToPlot$binomial)), sep="")
    # Xp <- max(ToPlot$AE, na.rm=TRUE)*(3/4)  
    # Yp <- max(ToPlot$RC, na.rm=TRUE)*(1/10)  
    
      p <- ggplot(ToPlot, aes(AE, RC, color=ToPlot$class)) +
        GGPoptions +
        geom_point(alpha=0.7) +
        geom_abline(slope = 1, intercept = 0, alpha=0.8) +
        xlab(paste(Traitaxis, "(", ImputedaxisX," by AE)", sep="")) + 
        ylab(paste(Traitaxis,  "(", ImputedaxisY," by RC)", sep="")) + 
        scale_color_manual(name="Class", labels = c("Birds", "Mammals"), values = c("#268bd2", "#dc322f")) 
        # annotate("text", x=Xp,y=Yp, label = Ssize)
      
      return(p)
    
  }
  
  if(Categorical) {

    
    Outcome <- AE_data[, c("Class", "Best_guess_binomial")]
    
    if(Diet) {
      
      x1 <- AE_data[,c("Best_guess_binomial",TraitName)]; colnames(x1)[2] <- "Primary_diet_AE"
      x2 <- RC_data[, c("binomial",TraitName)]; colnames(x2)[2] <- "Primary_diet_RC"
      xreturn <- cbind(x1,x2)
      xreturn <- xreturn %>%
        filter(Primary_diet_AE!=Primary_diet_RC)
    }
    
    for (i in 1:nrow(AE_data)){
      
      x <- AE_data[i, TraitName]
      y <- RC_data[i, TraitName]
      
      if(is.na(x)|is.na(y)) {Outcome$Result[i] <- "unknown"}
      
      else{
        
        if(x==y) {Outcome$Result[i] <- "same"} 
        if(x!=y)  { Outcome$Result[i] <- "different"}
      }
    }
    
    ToPlot <-  with(Outcome, table(Result, Class)) %>%
      as.data.frame() %>%
      group_by(Class) %>%
      mutate(Freq=Freq/sum(Freq)*100) %>%
      setNames(., c("outcome", "class", "prop")) %>%
      as.data.frame()
    
    # ToPlot <- table(ToPlot) %>% 
    #   as.data.frame() %>%
    #   setNames(., c("outcome", "prop")) %>%
    #   mutate(prop=prop/nrow(Outcome)*100)
    
    ToPlot <- ToPlot[order(ToPlot$prop),]
    ToPlot$class <- factor(ToPlot$class, levels=c("Aves", "Mammalia"))
    
    if(Comparison) {
      if(CompRC_collected) {ImputedaxisX <- "collected RC versus imputed AE"}
      if(!CompRC_collected) {ImputedaxisX <- "imputed RC versus collected AE"}
    }
    
    if(length(unique(ToPlot$outcome))==3){
      p <- ggplot(ToPlot, aes(outcome, prop, fill=class)) +
        GGPoptions +
        geom_bar(stat="identity", position="dodge") +
        xlab(paste(Traitaxis, ImputedaxisX, sep = ", ")) + ylab("% species") +
        scale_x_discrete(limits=c("different","unknown", "same"), labels=c("Contradicting","Unknown", "Similar")) +
        scale_fill_manual(name="Class", labels = c("Birds", "Mammals"), values = c("#268bd2", "#dc322f")) + 
        theme(legend.position='none') +
        ggtitle(paste("mammals:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Mammalia"])),
                      "\n birds:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Aves"])))) + 
        theme(plot.title = element_text(size =7.5,margin = margin(t = 0), hjust=1))
    }
    
    if(length(unique(ToPlot$outcome))==2){
      p <- ggplot(ToPlot, aes(outcome, prop, fill=class)) +
        GGPoptions +
        geom_bar(stat="identity", position="dodge") +
        xlab(paste(Traitaxis, ImputedaxisX, sep = ", ")) + ylab("% species") +
        scale_x_discrete(limits=c("different", "same"), labels=c("Contradicting", "Similar")) +
        scale_fill_manual(name="Class", labels = c("Birds", "Mammals"), values = c("#268bd2", "#dc322f")) + 
        theme(legend.position='none') +
        ggtitle(paste("mammals:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Mammalia"])),
                      "\n birds:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Aves"])))) + 
        theme(plot.title = element_text(size =7.5,margin = margin(t = 0), hjust=1))

    }
    

    
    # if(!Imputed){
    #   
    #   p <- ggplot(ToPlot, aes(outcome, prop)) +
    #     GGPoptions +
    #     geom_bar(stat="identity") +
    #     xlab(paste(Traitaxis, ImputedaxisX, sep = ", ")) + ylab("% species") +
    #     scale_x_discrete(limits=c("different", "unknown", "same"), labels=c("Contradicting", "Unkown", "Similar"))}
    
    if(Diet) {return(list(p=p, outputs=xreturn))}
    else{return(p)}
  }
}


Plot_all_values <- function(RC_data, AE_data, TraitName, Categorical, Diet, AxisX, AxisY) {
  
  X <- intersect(RC_data$binomial, AE_data$Best_guess_binomial)
  RC_data <- RC_data %>% filter(binomial %in% X)
  AE_data <- AE_data %>% filter(Best_guess_binomial %in% X)
  
  GGPoptions <- theme_classic() + theme(
    panel.border = element_rect(colour = "black", fill=NA),
    text = element_text(size=13, family="serif"), 
    axis.text.x = element_text(color="black", margin=ggplot2::margin(10,0,2,0,"pt"), size=12), 
    axis.text.y = element_text(color="black", margin=ggplot2::margin(0,10,0,0,"pt"), size=12),
    axis.ticks.length=unit(-0.1, "cm"),
    legend.text=element_text(size=13))
  
  ToPlot <- RC_data[, c("class","binomial", TraitName)]
  colnames(ToPlot)[3] <- "RC"
  ToPlot$AE <- AE_data[, TraitName]
  
  if(!Categorical) {  

    
    p <- ggplot(ToPlot, aes(AE, RC, color=ToPlot$class)) +
      GGPoptions +
      geom_point(alpha=0.7) +
      geom_abline(slope = 1, intercept = 0, alpha=0.5) +
      xlab(AxisX) + 
      ylab(AxisY) + 
      scale_color_manual(name="Class", labels = c("Birds", "Mammals"), values = c("#268bd2", "#dc322f"))
  
    
    return(p)}
    
    else {
      
    Outcome <- AE_data[, c("Class","Best_guess_binomial")]
    
    if(Diet) {
      
      x1 <- AE_data[,c("Best_guess_binomial",TraitName)]; colnames(x1)[2] <- "Primary_diet_AE"
      x2 <- RC_data[, c("binomial",TraitName)]; colnames(x2)[2] <- "Primary_diet_RC"
      xreturn <- cbind(x1,x2)
      xreturn <- xreturn %>%
        filter(Primary_diet_AE!=Primary_diet_RC)
    }
    
    for (i in 1:nrow(AE_data)){
      
      x <- AE_data[i, TraitName]
      y <- RC_data[i, TraitName]
      if(x==y) {Outcome$Result[i] <- "same"} 
      if(x!=y)  { Outcome$Result[i] <- "different"}
      }
  
    ToPlot <-  with(Outcome, table(Result, Class)) %>%
      as.data.frame() %>%
      group_by(Class) %>%
      mutate(Freq=Freq/sum(Freq)*100) %>%
      setNames(., c("outcome", "class", "prop")) %>%
      as.data.frame()
    
    
    ToPlot <- ToPlot[order(ToPlot$prop),]
    ToPlot$class <- factor(ToPlot$class, levels=c("Aves", "Mammalia"))
    
    
    if(length(unique(ToPlot$outcome))==3){
      p <- ggplot(ToPlot, aes(outcome, prop, fill=class)) +
        GGPoptions +
        geom_bar(stat="identity", position="dodge") +
        xlab(AxisX) + ylab("% species") +
        scale_x_discrete(limits=c("different","unknown", "same"), labels=c("Contradicting","Unknown", "Similar")) +
        scale_fill_manual(name="Class", labels = c("Birds", "Mammals"), values = c("#268bd2", "#dc322f")) + 
        theme(legend.position='none') +
        ggtitle(paste("mammals:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Mammalia"])),
                      "\n birds:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Aves"])))) + 
        theme(plot.title = element_text(size =7.5,margin = margin(t = 0), hjust=1)) +
        theme(legend.position="none")
        
    }
    
    if(length(unique(ToPlot$outcome))==2){
      p <- ggplot(ToPlot, aes(outcome, prop, fill=class)) +
        GGPoptions +
        geom_bar(stat="identity", position="dodge") +
        xlab(AxisX) + ylab("% species") +
        scale_x_discrete(limits=c("different", "same"), labels=c("Contradicting", "Similar")) +
        scale_fill_manual(name="Class", labels = c("Birds", "Mammals"), values = c("#268bd2", "#dc322f")) + 
        theme(legend.position='none') +
        ggtitle(paste("mammals:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Mammalia"])),
                      "\n birds:", length(unique(Outcome$Best_guess_binomial[Outcome$Class=="Aves"])))) + 
        theme(plot.title = element_text(size =7.5, margin = margin(t = 0), hjust=1)) +
        theme(legend.position="none")
    }
    
    if(Diet) {return(list(p=p, outputs=xreturn))}
    else{return(p)}
    
  }}


Plot.Cov <- function(TraitData, Traits, Main) {
  
  Names <- as.data.frame(c("Body_mass_g", "body_mass_median" ,
                           "Litter_size", "litter_clutch_size",
                           "Habitat_breadth_IUCN","hab_breadth",
                           "Diel_activity", "activity", "Primary_diet"))
  colnames(Names) <- "Original"
  Names$FP <- c(rep("Body mass",2), rep("Litter/clutch size",2), rep("Habitat breadth",2),rep("Diel activity",2), "Primary diet")
  
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





