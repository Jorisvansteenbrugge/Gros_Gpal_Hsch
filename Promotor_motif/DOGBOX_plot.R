library(magrittr)
library(ggplot2)
library(tidyverse)

data <- read.csv2("Downloads/data_dogbox/Gpal_Gros_Hs/promotor/", 
                 sep = ';',
                 stringsAsFactors = F)

data %>% head

data$DOGboxBinary <-   data$DOGbox %>% sapply( function(x) ifelse(x>0, 1, 0))  


data_percentages <- data.frame(matrix(ncol=3,nrow=0))

for(specie in unique(data$Species)){
  data_species <- data[data$Species == specie,]
  
  
  percentages_species <-  sapply(unique(data_species$ï..Family), function(family){
    data_species_family_dbbinary <- data_species[data_species$ï..Family == family, 'DOGboxBinary']
    total <- length(data_species_family_dbbinary)
    c <- data_species_family_dbbinary %>% sum
    
    perc <- (c*100/total) %>% as.numeric
    return(c(specie, family, perc))
    
  }) %>% t %>% rbind
  
  data_percentages <- rbind(data_percentages, percentages_species)
  
  
}

colnames(data_percentages) <- c("Species", "Family", "Percentage")
rownames(data_percentages) <- NULL

data_percentages$Percentage %<>% as.numeric
data_percentages$Family <-  factor(data_percentages$Family, levels = c("CLE", "VAP", "GLAND5", "GLAND6", "SPRYSEC", "GLAND13","GLAND4"))
data_percentages$Species <- factor(data_percentages$Species, levels = c("Gr-Line19", "Gp-D383", "Hs-IRS"))

data_percentages %>% ggplot(aes(y = Percentage, x = Family, group = Species)) + geom_bar(stat = 'identity') + facet_wrap(~Species) +
  theme_light()


data %>% ggplot(aes(y = DOGboxBinary, x = ï..Family, group = Species)) + geom_bar(stat = 'identity') + facet_wrap(~Species) +
  scale_y_continuous(breaks=seq(0,50,1))
