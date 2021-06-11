# This illustrates how the plot was created for Gpal - SPRYSEC genes


library(magrittr)
library(dplyr)
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


cor(Pallida_counts$DOGbox, Pallida_counts$TPM)
