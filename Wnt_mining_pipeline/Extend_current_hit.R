library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
extension_length <- as.numeric(args[1])


#load coordinates 
current_coordinates <- read.table("current_line_coordinates.csv", header=FALSE, sep=",")
found_genes_coordinates <- read.table("Coordinates_alrdy_examined.csv", header=FALSE, sep=",")


colnames(current_coordinates) <- c("seqnames", "start", "end")
colnames(found_genes_coordinates) <- c("seqnames", "start", "end")




#Extend the hit -extension_length and +extension_length, except if there is a gene/pseudogene before those limits 


my_seqname <- current_coordinates$seqnames
my_start <- current_coordinates$start
my_end <- current_coordinates$end


if (nrow(found_genes_coordinates %>% filter(seqnames == my_seqname) %>% filter(end >= (my_start-extension_length) & start < my_start & end < my_start)) > 0){
  
    extended_start = (tail(found_genes_coordinates %>% filter(seqnames == my_seqname) %>% filter(end >= (my_start-extension_length) & start < my_start & end < my_start) %>% dplyr::arrange(end), 1)$end) + 50
  
} else {

    extended_start = my_start-extension_length

}
  
  
if (nrow(found_genes_coordinates %>% filter(seqnames == my_seqname) %>% filter(start <= (my_end+extension_length) & end > my_end & start > my_end)) > 0){
  
    extended_end = (head(found_genes_coordinates %>% filter(seqnames == my_seqname) %>% filter(start <= (my_end+extension_length) & end > my_end & start > my_end) %>% dplyr::arrange(start), 1)$start) - 50
  
} else { 

    extended_end = my_end+extension_length
  
}


mynewcoordinates <- as.data.frame(cbind(my_seqname, extended_start, extended_end))

#Write regions in a text file
write.table(mynewcoordinates, "mynewcoordinates.csv", sep=",", quote=FALSE, col.names=FALSE, row.names=FALSE)




