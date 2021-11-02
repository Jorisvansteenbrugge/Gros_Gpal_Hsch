library(dplyr)
library(magrittr)
dogbox_file <- read.csv2("Downloads/data_dogbox/Gpal_Gros_Hs/promotor/DOG_Box_effector_family_overview.csv", 
                          sep = ',',
                          stringsAsFactors = F)

gr19_tpm <- read.csv2("/mnt/nemahomes/steen176/genome_comparisons/DOGbox/Gr19_RNAseq_v10_counts_clean.gtf",
                      sep = '\t', header = F, row.names = 1, stringsAsFactors = F)
colnames(gr19_tpm) <- c("TPM")
gr19_tpm %>% head

colnames(dogbox_file)

gr19 <- dogbox_file %>% filter(Species=="Gr-Line19",
                               Family=="SPRYSEC")

gr19 <- cbind(gr19, gr19_tpm[gr19$TranscriptID,])
colnames(gr19)[7] <- "TPM"
gr19$TPM %<>% as.numeric
gr19 %>% head

cor(gr19$TPM, gr19$DOGbox)



lm(TPM~DOGbox, data = gr19)
