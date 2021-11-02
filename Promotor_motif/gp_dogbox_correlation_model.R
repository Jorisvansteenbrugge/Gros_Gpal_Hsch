library(dplyr)
library(magrittr)
dogbox_file <- read.csv2("Downloads/data_dogbox/Gpal_Gros_Hs/promotor/DOG_Box_effector_family_overview.csv", 
                         sep = ',',
                         stringsAsFactors = F)

tpm <- read.csv2("Downloads/data_dogbox/Gpal_Gros_Hs/TPMs/D_3dpi_1_clean.txt",
                      sep = '\t', header = T, stringsAsFactors = F, row.names = 1)
tpm %>% head
gp <- dogbox_file %>% filter(Species=="Gp-D383",
                               Family=="SPRYSEC")

gp <- cbind(gp, tpm[gp$TranscriptID,])

colnames(gp)[7] <- "TPM"
gp$TPM %<>% as.numeric
gp %>% head

cor(gp$TPM, gp$DOGbox)

lm(TPM~DOGbox, data = gp)
