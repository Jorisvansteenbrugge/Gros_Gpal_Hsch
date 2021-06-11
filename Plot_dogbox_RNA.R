# This script illustrates how the DOG-box vs RNAseq relation was tested for Gpal D383 with SPRYSEC genes

library(magrittr)
library(dplyr)
library(ggplot2)

data <- read.csv2(
  "F:/Gpal_Gros_Hs/DOG_Box_effector_family_overview.csv", 
  sep = ';',
  stringsAsFactors = F
  )

colnames(data)[1] <- "Family"

tpm_counts <- read.csv(
  "F:/Gpal_Gros_Hs/TPMs/D_6dpi_1.txt",
  sep = '\t',
  stringsAsFactors = F,
  row.names = NULL
)

colnames(tpm_counts) <- c("Names", "TranscriptID", "TPM")
tpm_counts$TranscriptID %<>% trimws

Pallida <- data %>% 
  filter(
    Species=="Gp-D383",
    Family=="SPRYSEC"
    )

TPM <- sapply(Pallida$TranscriptID, function(x){
  tpm_counts %>% 
    filter(TranscriptID==x) %>% .$TPM
})

Pallida_counts <- cbind(Pallida, TPM)

model <- glm(TPM~DOGbox, data = Pallida_counts)
model
model %>% summary
plot(x = Pallida_counts$DOGbox, y = TPM)

abline(model)


Pallida_counts %>% 
  ggplot(aes(x = DOGbox, y = TPM)) + 
  geom_point(size = 3, colour = 'gray46', fill = 'black') +
  geom_abline(
    slope = model$coefficients["DOGbox"],
    intercept = model$coefficients['(Intercept)'] ) +
  xlab("No. of DOG-box motifs") +
  ylab("Transcripts per Million") + 
  ggtitle("G. pallida D383 - SPRYSEC") +
  theme_bw()

print(paste("Correlation :", cor(Pallida_counts$DOGbox, Pallida_counts$TPM)))
