# WNT_vertebrates
 Pipeline to mine WNT genes from vertebrate genomes

## I - WNT gene mining procedure

The mining of WNT genes is done using the script [Opsin_miner.sh](Opsin_miner.sh) 

Every file needed for the script are available in the "Database" folder, except for the uniprot database which can be downloaded here : https://www.uniprot.org/uniprotkb?query=reviewed%3Atrue&facets=reviewed%3Atrue

The script has to be launched like this : 

./[Opsin_miner.sh](Opsin_miner.sh)  Genome_fasta_file.fa [Database/all_opsins_cdhit80.prot](Database/all_opsins_cdhit80.prot) uniprot_sprot.fasta [Database/Scripts_opsins/](Database/Scripts_opsins/) max_intron_size number_of_thread

- Genome_fasta_file.fa = Your genome of interest, in fasta format  
-  max_intron_size = The maximum intron size you want. I used 50000 for all ray-finned fishes genomes.  
- number_of_thread = The number of threads to use. 

Note that the script is intented to be used in a SLURM environment. 

If you work in a non-slum environment, replace those two sbatch commands (line 127 and line 128) : 

`sbatch --job-name=ops_vs_scaff -W -c 2 --qos=6hours --mem=4G --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.exo.rslt ; sleep 10" &`

`sbatch --job-name=ops_vs_scaff -W -c 2 --qos=6hours --mem=4G --wrap="$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.noexhaustive.exo.rslt ; sleep 10" &`

by = 

`nohup $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.exo.rslt &`

`nohup $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --showtargetgff TRUE --model protein2genome --minintron 50 --maxintron $maximum_intron_length --ryo '%tcs' Gene_to_exo_per_scaff/$prot_to_exo $i > Exonerate_raw_results_folder/$file_name.noexhaustive.exo.rslt &`


If you work under a slurm environment, modify `--qos=6hours` by an existing partition in your environment.
