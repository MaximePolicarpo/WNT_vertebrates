#load packages
library(data.table)
library(dplyr)
library("plyranges")
library("GenomicRanges")

args = commandArgs(trailingOnly=TRUE)
extension_length <- as.numeric(args[1])


#load tblastn results
complete_coord <- read.table("Coordinates_completes.csv", header=FALSE, sep=",")
pseudo_coord <- read.table("Coordinates_pseudogenes.csv", header=FALSE, sep=",")

colnames(complete_coord) <- c("seqnames", "start", "end")
colnames(pseudo_coord) <- c("seqnames", "start", "end")


#Extend pseudogenes coordinates


new_coordinates_table <- as.data.frame(NULL)
for (i in seq(1:nrow(pseudo_coord))){

  my_seqname <- pseudo_coord[i,]$seqnames
  my_start <- pseudo_coord[i,]$start
  my_end <- pseudo_coord[i,]$end



  if (nrow(complete_coord %>% filter(seqnames == my_seqname) %>% filter(end >= (my_start-extension_length) & start < my_start & end < my_start)) > 0){
    
      extended_start = (tail(complete_coord %>% filter(seqnames == my_seqname) %>% filter(end >= (my_start-extension_length) & start < my_start & end < my_start) %>% dplyr::arrange(end), 1)$end) + 50
    
  } else {
  
      extended_start = my_start-extension_length
  
  }
    
    
  if (nrow(complete_coord %>% filter(seqnames == my_seqname) %>% filter(start <= (my_end+extension_length) & end > my_end & start > my_end)) > 0){
    
      extended_end = (head(complete_coord %>% filter(seqnames == my_seqname) %>% filter(start <= (my_end+extension_length) & end > my_end & start > my_end) %>% dplyr::arrange(start), 1)$start) - 50
    
  } else { 
  
      extended_end = my_end+extension_length
    
  }


  new_coordinates_table <- 
  as.data.frame(
    rbind(new_coordinates_table,
    cbind(my_seqname, extended_start, extended_end)))


}

colnames(new_coordinates_table) <- c("seqnames", "start", "end")
new_coordinates_table$start <- as.numeric(new_coordinates_table$start)
new_coordinates_table$end <- as.numeric(new_coordinates_table$end)


#Put the table as a grange object
new_coordinates_table_grange <- new_coordinates_table %>% as_granges()


#Reduce the table to merge overlapping results
new_coordinates_table_disjoin <- reduce(new_coordinates_table_grange,with.revmap=TRUE)

#transform the table in a data.frame

new_coordinates_table_filtered <- as.data.frame(new_coordinates_table_disjoin)
new_coordinates_table_filtered <- new_coordinates_table_filtered %>% dplyr::select(seqnames, start, end)

#Write regions in a text file

write.table(new_coordinates_table_filtered, file="new_coordinates_table_filtered.csv", quote=FALSE, sep=',', row.names = FALSE, col.names=FALSE)






