#load packages

library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")
library("tidyr")
library(stringr)



#Load miniprot results

miniprot_rslt <- read.table("Prot_vs_Genome.miniprot.len", header=FALSE, sep="\t")


#rename exonerate result columns
colnames(miniprot_rslt) <- c("seqnames", "program", "mRNA", "start", "end", 
                          "misc", "strand", "score","desc","prot_len", "query_length")




miniprot_rslt <- as.data.frame(miniprot_rslt  %>% mutate(scaffold_length = end-start))



 
miniprot_rslt <- 
  miniprot_rslt %>%
  mutate(perc_aligned = prot_len/query_length)




Best_hits_filtered <- as.data.frame(NULL) 



for(putative_protein_length in c(0.7, 0.65, 0.6, 0.55)){ 

	miniprot_rslt_round <- miniprot_rslt %>% filter(perc_aligned >= putative_protein_length)

	
	if(nrow(Best_hits_filtered) > 0){
		for (row in 1:nrow(Best_hits_filtered)) {
			
			
			miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
			
			
			miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (	Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
			
			miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
			
			miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
			
			  
		}
	}
	
		
	while (nrow(miniprot_rslt_round) > 0){
	
	
		miniprot_rslt_round_irange <- miniprot_rslt_round %>% as_granges()
		miniprot_rslt_round_disjoin <- reduce(miniprot_rslt_round_irange,with.revmap=TRUE)
		list_revmap <- as.data.frame(mcols(miniprot_rslt_round_disjoin))
	
		filtered_data <- c()
		for(i in 1:nrow(list_revmap)){
	  	filtered_data <- c(filtered_data, (slice(list_revmap, i) %>% unlist(use.names=FALSE))[which.min(slice(miniprot_rslt_round, slice(list_revmap, i) %>% unlist(use.names=FALSE))$scaffold_length)])
		}
	
	
		Best_hits_filtered <- rbind(Best_hits_filtered, slice(miniprot_rslt_round, filtered_data))
	
	
		if(nrow(Best_hits_filtered) > 0){
			for (row in 1:nrow(Best_hits_filtered)) {
			
			
				miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "start"])))
			
			
				miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start >= (	Best_hits_filtered[row, "start"]) & end <=(Best_hits_filtered[row, "end"])))
			
				miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "end"]) & end >= (Best_hits_filtered[row, "end"])))
			
				miniprot_rslt_round <- miniprot_rslt_round %>% filter(!(seqnames == as.character(Best_hits_filtered[row, "seqnames"]) & start <= (	Best_hits_filtered[row, "start"]) & end >= (Best_hits_filtered[row, "end"])))
			
			  
			}
		
		}
	
	
	}
	

}


#write the result in a table
write.table(Best_hits_filtered, file="Parsed_miniprot_rslt.tsv", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)



