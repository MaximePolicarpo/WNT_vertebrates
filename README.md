# WNT_vertebrates

## I - WNT gene mining procedure

The mining of WNT genes is done using the script [WNT_Pipeline.sh](Wnt_mining_pipeline/WNT_Pipeline.sh) 

Every file needed for the script are available, except for the uniprot database which can be downloaded here : https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue&facets=reviewed%3Atrue

The script has to be launched like this : 

./[WNT_Pipeline.sh](Wnt_mining_pipeline/WNT_Pipeline.sh)  Genome_fasta_file.fa [Candidate_WNT_genes_filtered_length.prot](Candidate_WNT_genes_filtered_length.prot) uniprot_sprot.fasta Folder_containing_scripts_of_this_repository max_intron_size number_of_thread


- Genome_fasta_file.fa = Your genome of interest, in fasta format  
-  max_intron_size = The maximum intron size you want. I used 50000 for all ray-finned fishes genomes.  
- number_of_thread = The number of threads to use. 
- Folder_containing_scripts_of_this_repository : A folder containing all the .sh, .py and .R scripts contained in the folder "Wnt_mining_pipeline".


## II - Manuscript Figures and analysis

Rscript to reproduce the various figures and analysis : [WNT_figures_analysis.R](Rscript_analysis/WNT_figures_analysis.R)  

Files necessary to run this script are located in the folder "Rscript_analysis"