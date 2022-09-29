library(data.table)
library(tidyr)


# input
file_gene1 <- "exons_gene1.txt"
file_gene2 <- "exons_gene2.txt"

exons_gene1 <- read.table(file_gene1, header=F, sep="\t")
exons_gene2 <- read.table(file_gene2, header=F, sep="\t")

exons_gene1 <- as.data.table(exons_gene1)
exons_gene2 <- as.data.table(exons_gene2)

colnames(exons_gene1) <- c("exon1","sequence1")
colnames(exons_gene2) <- c("exon2","sequence2")


# cross join between the two files
cross1 <- tidyr::crossing(exons_gene1, exons_gene2)
write.table(cross1,"exons_forward.txt", sep = "\t", row.names = F, col.names = F, quote = F, dec = ",")

cross2 <- tidyr::crossing(exons_gene2, exons_gene1)
write.table(cross2,"exons_reverse.txt", sep = "\t", row.names = F, col.names = F, quote = F, dec = ",")

