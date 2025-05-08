#load packages

library(data.table)
library(dplyr)
library("tidyr")

#Load miniprot results

miniprot_rslt <- read.table("Prot_vs_Genome_toparse.tsv", header=FALSE, sep="\t")
prot_length <- read.table("Target_protein_db.prot.fai", header=FALSE, sep="\t")


#rename exonerate result columns
colnames(miniprot_rslt) <- c("scaffold_region", "program", "mRNA", "start", "end", 
                          "misc", "strand", "score","desc","scaffold", "scaffold_start", "prot_start", "prot_stop", "query")

colnames(prot_length) <- c("query", "length", "misc1", "misc2", "misc3")


miniprot_rslt <-
	left_join(miniprot_rslt, prot_length, by="query")


miniprot_rslt <- 
	miniprot_rslt %>%
	mutate(true_start = scaffold_start + start) %>%
	mutate(true_stop = scaffold_start + end) %>%
	mutate(prot_len = prot_stop - prot_start + 1) 



miniprot_rslt <- 
	miniprot_rslt %>%
	dplyr::select(scaffold, program, mRNA, true_start, true_stop, misc, strand, score, desc, prot_len, length)



write.table(miniprot_rslt, "Prot_vs_Genome.miniprot.len", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
