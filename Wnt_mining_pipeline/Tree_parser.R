#load packages

library("ape")
library("dplyr")
library("phytools")

#load the tree
mytree <- read.tree("Final_ALL_verification_alignment.aln.treefile")

#name internal nodes of the tree
mytree_label <- makeNodeLabel(mytree, method = "number", prefix="Node")

#import the name of known OR genes
known_Prot <- scan("Known_Prot_id.txt", what="character")
Outgroups_seqs <- scan("Known_outgroups_id.txt", what="character")


#Root the tree at the common ancestor of outgroup sequences
MRCA_outgroup <- findMRCA(mytree_label, tips=Outgroups_seqs, type="node")
mytree_rooted <- root(mytree_label, node=MRCA_outgroup, resolve.root= TRUE) #root

#Check the MRCA or V2R genes

MRCA_Prot <- findMRCA(mytree_rooted, tips=known_Prot, type="node")

#grep all tips from these MRCS
Prot_genes <- as.vector(extract.clade(mytree_rooted, MRCA_Prot)$tip)

#Remove alrdy known TAAR genes
Current_species_Prot <- setdiff(Prot_genes, known_Prot)


write(x=Current_species_Prot, file="Current_species_Prot.txt")








