genome=$1 
target_DB=$2 
target_outgroup_DB=$3 
scripts_folder_location=$4 ; scripts_location=$( echo "$scripts_folder_location" | sed 's/\/$//' )
maximum_intron_length=$5
prot_name=$6
number_of_threads=$7
extension_lines=$8

previous_iteration_nb_sequences=`if test -f "Complete_genes.fa" ; then grep -c ">" Complete_genes.fa ; else echo "0" ;  fi`
current_nb_sequences=$((previous_iteration_nb_sequences + 1))
number_regions_blast=`grep "[0-9]" Potential_regions.tsv | wc -l`


while [ "$current_nb_sequences" -gt "$previous_iteration_nb_sequences" ] && [ "$number_regions_blast" -gt "0" ] ; do



	cat Complete_genes.fa Pseudogenes.fa > All_sequences.fa
	previous_iteration_nb_sequences=`if test -f "All_sequences.fa" ; then grep -c ">" All_sequences.fa ; else echo "0" ;  fi`


	xargs samtools faidx $genome < Potential_regions.tsv > Potential_regions.fa
	
	
	### MINIPROT SECTION
	
	#cp /scicore/home/salzburg/polica0000/Vertebrates_Taste_Receptors/Danio_rerio/$genome ./
	
	#miniprot -t $number_of_threads -d GCA_000002035.4_GRCz11_genomic.trimmed.mpi $genome
	rm Potential_regions.mpi
	miniprot -t $number_of_threads -d Potential_regions.mpi Potential_regions.fa
	
	
	
	#rm Splitted_db/*.miniprot
	#for protein_db in Splitted_db/* ; do 
	#	file_name=`echo "$protein_db" | sed "s/Splitted_db\///g"`
	#	miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 GCA_000002035.4_GRCz11_genomic.trimmed.mpi Splitted_db/$file_name > Splitted_db/$file_name.miniprot
	#done
	
	
	
	#miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 Potential_regions.mpi Database_2022_V2R_vertebrates_cdhit_80.prot > Raw_Prot.miniprot
	miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 Potential_regions.mpi Target_protein_db.prot > Raw_Prot.miniprot
	
	for curr_region in `cat Potential_regions.tsv` ; do 
		samtools faidx $genome $curr_region > current_region_miniprot.fa 
		miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 current_region_miniprot.fa Target_protein_db.prot >> Raw_Prot.miniprot
	done


	#Parse results
	
	#rm Prot_vs_Genome.miniprot
	#for protein_db in Splitted_db/*.miniprot ; do 
	#	file_name=`echo "$protein_db" | sed "s/Splitted_db\///g"`
	#	grep "miniprot.*mRNA" Splitted_db/$file_name >> Prot_vs_Genome.miniprot
	#
	#done
	
	
	grep "miniprot.*mRNA" Raw_Prot.miniprot > Prot_vs_Genome.miniprot
	
	cut -f1 Prot_vs_Genome.miniprot | sed 's/:/-/g' > scaffold_region.txt
	cut -f1 -d "-" scaffold_region.txt > scaffold.txt
	cut -f2 -d "-" scaffold_region.txt > start.txt
	cut -f9 Prot_vs_Genome.miniprot | cut -f2 -d " " | tr " " "\n" > prot_start.txt
	cut -f9 Prot_vs_Genome.miniprot | cut -f3 -d " " | tr " " "\n" > prot_stop.txt
	cut -f9 Prot_vs_Genome.miniprot | sed 's/.*Target=//g' | sed 's/ .*//g' > list_query.txt

	paste -d "\t" Prot_vs_Genome.miniprot scaffold.txt start.txt prot_start.txt prot_stop.txt list_query.txt > Prot_vs_Genome_toparse.tsv


	nb_prot_vs_genome=`wc -l < Prot_vs_Genome_toparse.tsv`

	if [ "$nb_prot_vs_genome" -gt 0 ] ; then
		Rscript $scripts_location/Correct_miniprot_table.R
		Rscript $scripts_location/Parse_miniprot_results_fifth.R

	else 

		echo "" > Parsed_miniprot_rslt.tsv
	fi



	
	IFS=$'\n'
	for line in `cat Parsed_miniprot_rslt.tsv` ; do 
	
		scaffold=`echo "$line" | cut -f1`
		start=`echo "$line" | cut -f4`
		stop=`echo "$line" | cut -f5`
	
		echo "$scaffold,$start,$stop" > current_line_coordinates.csv
		Rscript $scripts_location/Extend_current_hit.R $extension_lines
		start_extanded=`cut -f2 mynewcoordinates.csv -d ","`
		stop_extanded=`cut -f3 mynewcoordinates.csv -d ","`
		if [ "$start_extanded" -lt 1 ] ; then start_extanded=1; fi
	
		if [ "$stop_extanded" -le "$start_extanded" ] ; then start_extanded=$((start - 300)) ; stop_extanded=$((stop + 300)) ; fi
		#Extract region
	
		samtools faidx $genome $scaffold:$start-$stop > region_for_blast.fa
		samtools faidx $genome $scaffold:$start_extanded-$stop_extanded > region.fa

		makeblastdb -in region_for_blast.fa -dbtype nucl
		makeblastdb -in region.fa -dbtype nucl

		#Find the best query
	
		blastx -query region.fa -db Target_protein_db.prot -evalue 1e-5 -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads $number_of_threads
	 	best_query=`head -1 blastx_result | cut -f2`
	 	samtools faidx Target_protein_db.prot $best_query > best_query.prot
	
	 	#use miniprot with the best query against the region

	 	miniprot -t $number_of_threads --outn=1 --gff --aln --trans -L 5 -G $maximum_intron_length region.fa best_query.prot > curr_miniprot_result.txt
	
	 	#Extract the result protein sequence without frameshift
	 	fs_less_prot=`grep "#STA" curr_miniprot_result.txt | cut -f2`
	 	echo ">curr_seq" > curr_prot.fa
	 	echo "$fs_less_prot" >> curr_prot.fa
	 	blastp -query curr_prot.fa -db $target_outgroup_DB -evalue 1e-5 -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_threads
	
	
	 	#Continue if the best match is a V2R receptor
	 	if grep -q -i "$prot_name" blastp_result ; then
	
	 		grep "miniprot.*mRNA" curr_miniprot_result.txt > mRNA.txt
	 		strand=`cut -f7 mRNA.txt`
	 		gene_start=`cut -f4 mRNA.txt` 
			gene_stop=`cut -f5 mRNA.txt` 
			gene_stop_plus_1=$((gene_stop + 1))
			gene_start_minus_1=$((gene_start - 1))

			grep "miniprot.*CDS" curr_miniprot_result.txt > all_exons.txt
	 		nb_exons=`wc -l < all_exons.txt`
	 		grep "miniprot.*stop_codon" curr_miniprot_result.txt >> all_exons.txt

	 		sed 's/>.*/>myregion/g' region.fa > renamed_region.fa
	 		
	 		rm renamed_region.fa.fai
	 		rm curr.ORF


	 		if [ $strand == "+" ] ; then 

	 			strand_print="plus"
	 			echo ">curr_seq" > curr_nucl.fa
	 			echo ">curr_seq" > only_exons.fa
	
	 			extend_upstream_coord=$((gene_start_minus_1 - 300))
	 			extend_downstream_coord=$((gene_stop_plus_1 + 300))

	 			rm Extend_upstream.fa ; rm Extend_downstream.fa

				samtools faidx renamed_region.fa myregion:$extend_upstream_coord-$gene_start_minus_1 > Extend_upstream.fa
				samtools faidx renamed_region.fa myregion:$gene_stop_plus_1-$extend_downstream_coord > Extend_downstream.fa
				
				grep -v ">" Extend_upstream.fa | tr -d '\n' | sed 's/NNN//g' > Extend_upstream.txt
				grep -v ">" Extend_downstream.fa | tr -d '\n' | sed 's/NNN//g' > Extend_downstream.txt


				cat Extend_upstream.txt >> curr_nucl.fa


	 			for exon in `cat all_exons.txt` ; do
	 				exon_start=`echo "$exon" | cut -f4`
	 				exon_stop=`echo "$exon" | cut -f5`
	 				samtools faidx renamed_region.fa myregion:$exon_start-$exon_stop | grep -v "myregion" >> curr_nucl.fa
	 				samtools faidx renamed_region.fa myregion:$exon_start-$exon_stop | grep -v "myregion" >> only_exons.fa
	 			done

	 			cat Extend_downstream.txt >> curr_nucl.fa

	 			query_length=`grep "$best_query" Target_protein_db.prot.fai | cut -f2`
	 			perc80_query_length=$((query_length*75/100*3))

	 			getorf -sequence curr_nucl.fa -outseq curr.ORF -minsize $perc80_query_length -find 3 -reverse FALSE




	 		elif [ $strand == "-" ] ; then 

	 			strand_print="minus"
	 			echo ">curr_seq" > curr_nucl.fa
	 			echo ">curr_seq" > only_exons.fa

	 			extend_upstream_coord=$((gene_start_minus_1 - 300))
				extend_downstream_coord=$((gene_stop_plus_1 + 300))

				rm Extend_upstream.fa ; rm Extend_downstream.fa

				samtools faidx renamed_region.fa myregion:$extend_upstream_coord-$gene_start_minus_1 > Extend_upstream.fa
				samtools faidx renamed_region.fa myregion:$gene_stop_plus_1-$extend_downstream_coord > Extend_downstream.fa
				revseq Extend_upstream.fa Extend_upstream.rev
				revseq Extend_downstream.fa Extend_downstream.rev
				grep -v ">" Extend_upstream.rev | tr -d '\n' | sed 's/NNN//g' > Extend_upstream.txt
				grep -v ">" Extend_downstream.rev | tr -d '\n' | sed 's/NNN//g' > Extend_downstream.txt

				cat Extend_downstream.txt >> curr_nucl.fa
				for exon in `cat all_exons.txt` ; do
	 				
	 				exon_start=`echo "$exon" | cut -f4`
	 				exon_stop=`echo "$exon" | cut -f5`

	 				samtools faidx renamed_region.fa myregion:$exon_start-$exon_stop > temp.fa
	 				revseq temp.fa temp.rev 
	 				grep -v ">" temp.rev >> curr_nucl.fa
	 				grep -v ">" temp.rev >> only_exons.fa

	 			done
	 			cat Extend_upstream.txt >> curr_nucl.fa



	 			query_length=`grep "$best_query" Target_protein_db.prot.fai | cut -f2`
	 			perc80_query_length=$((query_length*75/100*3))

	 			getorf -sequence curr_nucl.fa -outseq curr.ORF -minsize $perc80_query_length -find 3 -reverse FALSE


	 		fi



	 		if [ `grep -c ">" curr.ORF` -ge 1 ] ; then

	 			sed -i "s/ - /-/g" curr.ORF
	 			sed -i "s/ /-/g" curr.ORF
	 			sed -i "s/\[//g" curr.ORF
	 			sed -i "s/]//g" curr.ORF
	 			rm curr.ORF.fai
	 			good_seq=`grep ">" curr.ORF | head -1 | sed "s/>//g"`
	 			samtools faidx curr.ORF $good_seq > temp ; mv temp curr.ORF ; rm curr.ORF.fai

	 			transeq curr.ORF curr.ORFP

	 			miniprot -t $number_of_threads --gff --aln --trans -L 5 -G $maximum_intron_length region.fa curr.ORFP > curr_miniprot_result.txt

	 			grep "miniprot.*mRNA" curr_miniprot_result.txt > mRNA.txt


				start_incr=`cut -f4 mRNA.txt`
	 			stop_incr=`cut -f5 mRNA.txt`
	 		

	 			true_start=$((start_extanded + start_incr))
	 			true_stop=$((start_extanded + stop_incr))
	 			gene_coordinate=`echo "$scaffold-$true_start-$true_stop"`

	 			gene_name=`echo "$gene_coordinate-$nb_exons-$strand_print"`

	 			samtools faidx $genome $scaffold:$true_start-$true_stop > gene_region.fa
	 			makeblastdb -in gene_region.fa -dbtype nucl
				number_blast_hit=`tblastn -query curr.ORFP -db gene_region.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`


				 


				echo "$scaffold,$true_start,$true_stop" >> Coordinates_alrdy_examined.csv
				sed "s/$good_seq/$gene_name/g" curr.ORF >> Complete_genes.fa


				


			else 



				stop_codon_state="FALSE"
				edge_state="FALSE"
				frameshift_state="FALSE"

				start_incr=`cut -f4 mRNA.txt`
				stop_incr=`cut -f5 mRNA.txt`
				true_start=$((start_extanded + start_incr))
				true_stop=$((start_extanded + stop_incr))
				samtools faidx $genome $scaffold:$true_start-$true_stop > gene_region.fa
				makeblastdb -in gene_region.fa -dbtype nucl
				number_blast_hit=`tblastn -query curr_prot.fa -db gene_region.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`


				gene_coordinate=`echo "$scaffold-$true_start-$true_stop"`


				 


				if grep -q "Frameshift" mRNA.txt ; then frameshift_state="TRUE" ; fi
				if grep -q "StopCodon" mRNA.txt ; then stop_codon_state="TRUE" ; fi
				echo "$fs_less_prot" > test.fa

				sed '$ s/.$//' test.fa > test_nolast.fa
				if grep -q "\*" test_nolast.fa ; then stop_codon_state="TRUE" ; fi

				start_check=`grep ">" region.fa | sed 's/.*://g' | cut -f1 -d "-"`
				stop_check=`grep ">" region.fa | sed 's/.*://g' | cut -f2 -d "-"`
				#First check if these coordinates are near the end of scaffolds (<5000 bp)
				if [ "$start_check" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
				scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
				diff_lengths=$((scaffold_length - stop_check))
				if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
				
				#Now check if there are consecutive N near the gene that could indicate conting end
				extanded_start_coord=$((start_check - 200))
				extanded_end_coord=$((stop_check + 200))
				
				#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
				consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
				if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
				if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				gene_name=`echo "$gene_coordinate-$nb_exons-$strand_print---$edge_state-$stop_codon_state-$frameshift_state"`

				fasta_formatter -i only_exons.fa -o temp.fa -w 60
				mv temp.fa only_exons.fa
				sed "s/>.*/>$gene_name/g" only_exons.fa >> Pseudogenes.fa

	 			echo ">curr_seq" > curr_prot.fa
	 			echo "$fs_less_prot" >> curr_prot.fa
				sed "s/>curr_seq/>$gene_name/g" curr_prot.fa >> Pseudogenes_wo_Fs.prot
				echo "$scaffold,$true_start,$true_stop" >> Coordinates_alrdy_examined.csv



				
	



	 		fi 

	 	else 



	 		echo "$scaffold,$start_extanded,$stop_extanded" >> Coordinates_alrdy_examined.csv

	 	fi

	done




	if [ `wc -l < Coordinates_alrdy_examined.csv` -lt 1 ] ; then echo "Simulated_scaffold,1,10" >> Coordinates_alrdy_examined.csv ; fi
	tr ',' '\t' < Coordinates_alrdy_examined.csv > Coordinates_already_examined.tsv


	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Complete_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa Complete_genes.fa ; rm temp.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa Pseudogenes.fa ; rm temp.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_wo_Fs.prot | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa Pseudogenes_wo_Fs.prot ; rm temp.fa


	cat Complete_genes.fa Pseudogenes.fa > All_sequences.fa
	current_nb_sequences=`if test -f "All_sequences.fa" ; then grep -c ">" All_sequences.fa ; else echo "0" ;  fi`
	
	#re-process blast result to find potential V2R regions exclusing already found genes
	
	Rscript $scripts_location/Rscript_merge_filter_extend_blast_hit_Second.R $maximum_intron_length

	if test -f "Potential_regions.tsv" ; then number_regions_blast=`grep "[0-9]" Potential_regions.tsv | wc -l` ; else number_regions_blast=0 ; fi
	rm Parsed_miniprot_rslt.tsv


	cp $target_DB ./Target_protein_db.prot 

	transeq Complete_genes.fa Complete_genes.prot
	

done



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Remove filter on the nb of exons ############################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

previous_iteration_nb_sequences=`if test -f "Complete_genes.fa" ; then grep -c ">" Complete_genes.fa ; else echo "0" ;  fi`
current_nb_sequences=$((previous_iteration_nb_sequences + 1))
number_regions_blast=`grep "[0-9]" Potential_regions.tsv | wc -l`

while [ "$current_nb_sequences" -gt "$previous_iteration_nb_sequences" ] && [ "$number_regions_blast" -gt "0" ] ; do



	cat Complete_genes.fa Pseudogenes.fa > All_sequences.fa
	previous_iteration_nb_sequences=`if test -f "All_sequences.fa" ; then grep -c ">" All_sequences.fa ; else echo "0" ;  fi`


	xargs samtools faidx $genome < Potential_regions.tsv > Potential_regions.fa
	
	
	### MINIPROT SECTION
	
	#cp /scicore/home/salzburg/polica0000/Vertebrates_Taste_Receptors/Danio_rerio/$genome ./
	
	#miniprot -t $number_of_threads -d GCA_000002035.4_GRCz11_genomic.trimmed.mpi $genome
	rm Potential_regions.mpi
	miniprot -t $number_of_threads -d Potential_regions.mpi Potential_regions.fa
	
	
	
	#rm Splitted_db/*.miniprot
	#for protein_db in Splitted_db/* ; do 
	#	file_name=`echo "$protein_db" | sed "s/Splitted_db\///g"`
	#	miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 GCA_000002035.4_GRCz11_genomic.trimmed.mpi Splitted_db/$file_name > Splitted_db/$file_name.miniprot
	#done
	
	
	
	#miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 Potential_regions.mpi Database_2022_V2R_vertebrates_cdhit_80.prot > Raw_Prot.miniprot
	miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 Potential_regions.mpi Target_protein_db.prot > Raw_Prot.miniprot
	
	for curr_region in `cat Potential_regions.tsv` ; do 
		samtools faidx $genome $curr_region > current_region_miniprot.fa 
		miniprot -t $number_of_threads --gff --aln -L 5 -G $maximum_intron_length --outs 0.2 current_region_miniprot.fa Target_protein_db.prot >> Raw_Prot.miniprot
	done


	#Parse results
	
	#rm Prot_vs_Genome.miniprot
	#for protein_db in Splitted_db/*.miniprot ; do 
	#	file_name=`echo "$protein_db" | sed "s/Splitted_db\///g"`
	#	grep "miniprot.*mRNA" Splitted_db/$file_name >> Prot_vs_Genome.miniprot
	#
	#done
	
	
	grep "miniprot.*mRNA" Raw_Prot.miniprot > Prot_vs_Genome.miniprot
	
	cut -f1 Prot_vs_Genome.miniprot | sed 's/:/-/g' > scaffold_region.txt
	cut -f1 -d "-" scaffold_region.txt > scaffold.txt
	cut -f2 -d "-" scaffold_region.txt > start.txt
	cut -f9 Prot_vs_Genome.miniprot | cut -f2 -d " " | tr " " "\n" > prot_start.txt
	cut -f9 Prot_vs_Genome.miniprot | cut -f3 -d " " | tr " " "\n" > prot_stop.txt
	cut -f9 Prot_vs_Genome.miniprot | sed 's/.*Target=//g' | sed 's/ .*//g' > list_query.txt

	paste -d "\t" Prot_vs_Genome.miniprot scaffold.txt start.txt prot_start.txt prot_stop.txt list_query.txt > Prot_vs_Genome_toparse.tsv



	nb_prot_vs_genome=`wc -l < Prot_vs_Genome_toparse.tsv`

	if [ "$nb_prot_vs_genome" -gt 0 ] ; then
		Rscript $scripts_location/Correct_miniprot_table.R
		Rscript $scripts_location/Parse_miniprot_results_fifth.R

	else 

		echo "" > Parsed_miniprot_rslt.tsv
	fi
	
	IFS=$'\n'
	for line in `cat Parsed_miniprot_rslt.tsv` ; do 
	
		scaffold=`echo "$line" | cut -f1`
		start=`echo "$line" | cut -f4`
		stop=`echo "$line" | cut -f5`
	
		echo "$scaffold,$start,$stop" > current_line_coordinates.csv
		Rscript $scripts_location/Extend_current_hit.R $extension_lines
		start_extanded=`cut -f2 mynewcoordinates.csv -d ","`
		stop_extanded=`cut -f3 mynewcoordinates.csv -d ","`
		if [ "$start_extanded" -lt 1 ] ; then start_extanded=1; fi

		if [ "$stop_extanded" -le "$start_extanded" ] ; then start_extanded=$((start - 300)) ; stop_extanded=$((stop + 300)) ; fi
	
		#Extract region
	
		samtools faidx $genome $scaffold:$start-$stop > region_for_blast.fa
		samtools faidx $genome $scaffold:$start_extanded-$stop_extanded > region.fa

		makeblastdb -in region_for_blast.fa -dbtype nucl
		makeblastdb -in region.fa -dbtype nucl

		#Find the best query
	
		blastx -query region.fa -db Target_protein_db.prot -evalue 1e-5 -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out blastx_result -max_target_seqs 1 -num_threads $number_of_threads
	 	best_query=`head -1 blastx_result | cut -f2`
	 	samtools faidx Target_protein_db.prot $best_query > best_query.prot
	
	 	#use miniprot with the best query against the region

	 	miniprot -t $number_of_threads --outn=1 --gff --aln --trans -L 5 -G $maximum_intron_length region.fa best_query.prot > curr_miniprot_result.txt
	
	 	#Extract the result protein sequence without frameshift
	 	fs_less_prot=`grep "#STA" curr_miniprot_result.txt | cut -f2`
	 	echo ">curr_seq" > curr_prot.fa
	 	echo "$fs_less_prot" >> curr_prot.fa
	 	blastp -query curr_prot.fa -db $target_outgroup_DB -evalue 1e-5 -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out blastp_result -max_target_seqs 1 -num_threads $number_of_threads
	
	
	 	#Continue if the best match is a V2R receptor
	 	if grep -q -i "$prot_name" blastp_result ; then
	
	 		grep "miniprot.*mRNA" curr_miniprot_result.txt > mRNA.txt
	 		strand=`cut -f7 mRNA.txt`
	 		gene_start=`cut -f4 mRNA.txt` 
			gene_stop=`cut -f5 mRNA.txt` 
			gene_stop_plus_1=$((gene_stop + 1))
			gene_start_minus_1=$((gene_start - 1))

			grep "miniprot.*CDS" curr_miniprot_result.txt > all_exons.txt
	 		nb_exons=`wc -l < all_exons.txt`
	 		grep "miniprot.*stop_codon" curr_miniprot_result.txt >> all_exons.txt
	 		
	 		sed 's/>.*/>myregion/g' region.fa > renamed_region.fa
	 		
	 		rm renamed_region.fa.fai
	 		rm curr.ORF


	 		if [ $strand == "+" ] ; then 

	 			strand_print="plus"
	 			echo ">curr_seq" > curr_nucl.fa
	 			echo ">curr_seq" > only_exons.fa
	
	 			extend_upstream_coord=$((gene_start_minus_1 - 300))
	 			extend_downstream_coord=$((gene_stop_plus_1 + 300))

	 			rm Extend_upstream.fa ; rm Extend_downstream.fa

				samtools faidx renamed_region.fa myregion:$extend_upstream_coord-$gene_start_minus_1 > Extend_upstream.fa
				samtools faidx renamed_region.fa myregion:$gene_stop_plus_1-$extend_downstream_coord > Extend_downstream.fa
				
				grep -v ">" Extend_upstream.fa | tr -d '\n' | sed 's/NNN//g' > Extend_upstream.txt
				grep -v ">" Extend_downstream.fa | tr -d '\n' | sed 's/NNN//g' > Extend_downstream.txt


				cat Extend_upstream.txt >> curr_nucl.fa


	 			for exon in `cat all_exons.txt` ; do
	 				exon_start=`echo "$exon" | cut -f4`
	 				exon_stop=`echo "$exon" | cut -f5`
	 				samtools faidx renamed_region.fa myregion:$exon_start-$exon_stop | grep -v "myregion" >> curr_nucl.fa
	 				samtools faidx renamed_region.fa myregion:$exon_start-$exon_stop | grep -v "myregion" >> only_exons.fa
	 			done

	 			cat Extend_downstream.txt >> curr_nucl.fa

	 			query_length=`grep "$best_query" Target_protein_db.prot.fai | cut -f2`
	 			perc80_query_length=$((query_length*75/100*3))

	 			getorf -sequence curr_nucl.fa -outseq curr.ORF -minsize $perc80_query_length -find 3 -reverse FALSE




	 		elif [ $strand == "-" ] ; then 

	 			strand_print="minus"
	 			echo ">curr_seq" > curr_nucl.fa
	 			echo ">curr_seq" > only_exons.fa

	 			extend_upstream_coord=$((gene_start_minus_1 - 300))
				extend_downstream_coord=$((gene_stop_plus_1 + 300))

				rm Extend_upstream.fa ; rm Extend_downstream.fa

				samtools faidx renamed_region.fa myregion:$extend_upstream_coord-$gene_start_minus_1 > Extend_upstream.fa
				samtools faidx renamed_region.fa myregion:$gene_stop_plus_1-$extend_downstream_coord > Extend_downstream.fa
				revseq Extend_upstream.fa Extend_upstream.rev
				revseq Extend_downstream.fa Extend_downstream.rev
				grep -v ">" Extend_upstream.rev | tr -d '\n' | sed 's/NNN//g' > Extend_upstream.txt
				grep -v ">" Extend_downstream.rev | tr -d '\n' | sed 's/NNN//g' > Extend_downstream.txt

				cat Extend_downstream.txt >> curr_nucl.fa
				for exon in `cat all_exons.txt` ; do
	 				
	 				exon_start=`echo "$exon" | cut -f4`
	 				exon_stop=`echo "$exon" | cut -f5`

	 				samtools faidx renamed_region.fa myregion:$exon_start-$exon_stop > temp.fa
	 				revseq temp.fa temp.rev 
	 				grep -v ">" temp.rev >> curr_nucl.fa
	 				grep -v ">" temp.rev >> only_exons.fa

	 			done
	 			cat Extend_upstream.txt >> curr_nucl.fa



	 			query_length=`grep "$best_query" Target_protein_db.prot.fai | cut -f2`
	 			perc80_query_length=$((query_length*75/100*3))

	 			getorf -sequence curr_nucl.fa -outseq curr.ORF -minsize $perc80_query_length -find 3 -reverse FALSE


	 		fi



	 		if [ `grep -c ">" curr.ORF` -ge 1 ] ; then

	 			sed -i "s/ - /-/g" curr.ORF
	 			sed -i "s/ /-/g" curr.ORF
	 			sed -i "s/\[//g" curr.ORF
	 			sed -i "s/]//g" curr.ORF
	 			rm curr.ORF.fai
	 			good_seq=`grep ">" curr.ORF | head -1 | sed "s/>//g"`
	 			samtools faidx curr.ORF $good_seq > temp ; mv temp curr.ORF ; rm curr.ORF.fai

	 			transeq curr.ORF curr.ORFP

	 			miniprot -t $number_of_threads --gff --aln --trans -L 5 -G $maximum_intron_length region.fa curr.ORFP > curr_miniprot_result.txt

	 			grep "miniprot.*mRNA" curr_miniprot_result.txt > mRNA.txt


				start_incr=`cut -f4 mRNA.txt`
	 			stop_incr=`cut -f5 mRNA.txt`
	 		

	 			true_start=$((start_extanded + start_incr))
	 			true_stop=$((start_extanded + stop_incr))
	 			gene_coordinate=`echo "$scaffold-$true_start-$true_stop"`

	 			gene_name=`echo "$gene_coordinate-$nb_exons-$strand_print"`

	 			samtools faidx $genome $scaffold:$true_start-$true_stop > gene_region.fa
	 			makeblastdb -in gene_region.fa -dbtype nucl
				number_blast_hit=`tblastn -query curr.ORFP -db gene_region.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`


				


				echo "$scaffold,$true_start,$true_stop" >> Coordinates_alrdy_examined.csv
				sed "s/$good_seq/$gene_name/g" curr.ORF >> Complete_genes.fa


				


			else 



				stop_codon_state="FALSE"
				edge_state="FALSE"
				frameshift_state="FALSE"

				start_incr=`cut -f4 mRNA.txt`
				stop_incr=`cut -f5 mRNA.txt`
				true_start=$((start_extanded + start_incr))
				true_stop=$((start_extanded + stop_incr))
				samtools faidx $genome $scaffold:$true_start-$true_stop > gene_region.fa
				makeblastdb -in gene_region.fa -dbtype nucl
				number_blast_hit=`tblastn -query curr_prot.fa -db gene_region.fa -evalue 1e-10 -outfmt 6 | awk '{ if ($4 >= 30) { print } }' | wc -l`


				gene_coordinate=`echo "$scaffold-$true_start-$true_stop"`


				

				if grep -q "Frameshift" mRNA.txt ; then frameshift_state="TRUE" ; fi
				if grep -q "StopCodon" mRNA.txt ; then stop_codon_state="TRUE" ; fi
				echo "$fs_less_prot" > test.fa

				sed '$ s/.$//' test.fa > test_nolast.fa
				if grep -q "\*" test_nolast.fa ; then stop_codon_state="TRUE" ; fi

				start_check=`grep ">" region.fa | sed 's/.*://g' | cut -f1 -d "-"`
				stop_check=`grep ">" region.fa | sed 's/.*://g' | cut -f2 -d "-"`
				#First check if these coordinates are near the end of scaffolds (<5000 bp)
				if [ "$start_check" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the start of scaffold
				scaffold_length=`grep -m1 "^$scaffold	" $genome.fai | cut -f2` #extract scaffold length from .fai file
				diff_lengths=$((scaffold_length - stop_check))
				if [ "$diff_lengths" -le '5000' ] ; then edge_state="TRUE" ; fi #check if its near the end of scaffold
				
				#Now check if there are consecutive N near the gene that could indicate conting end
				extanded_start_coord=$((start_check - 200))
				extanded_end_coord=$((stop_check + 200))
				
				#Command below extract the region, put in a single line and count the consecutive number of N (only take the greatest number)
				consecutive_N_nb=`samtools faidx $genome $scaffold:$extanded_start_coord-$extanded_end_coord | sed 's/n/N/g' | grep -v ">" | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | grep N | awk -F '[^N]+' '{for (i=1; i<=NF; i++) if ($i != "") print length($i)}' | sort -n | tail -1`
				if [[ $consecutive_N_nb == "" ]] ; then consecutive_N_nb=0 ; fi
				if [ "$consecutive_N_nb" -ge '50' ] ; then edge_state="TRUE" ; fi 
				
				gene_name=`echo "$gene_coordinate-$nb_exons-$strand_print---$edge_state-$stop_codon_state-$frameshift_state"`

				fasta_formatter -i only_exons.fa -o temp.fa -w 60
				mv temp.fa only_exons.fa
				sed "s/>.*/>$gene_name/g" only_exons.fa >> Pseudogenes.fa

	 			echo ">curr_seq" > curr_prot.fa
	 			echo "$fs_less_prot" >> curr_prot.fa
				sed "s/>curr_seq/>$gene_name/g" curr_prot.fa >> Pseudogenes_wo_Fs.prot
				echo "$scaffold,$true_start,$true_stop" >> Coordinates_alrdy_examined.csv



				
	



	 		fi 

	 	else 



	 		echo "$scaffold,$start_extanded,$stop_extanded" >> Coordinates_alrdy_examined.csv

	 	fi

	done




	if [ `wc -l < Coordinates_alrdy_examined.csv` -lt 1 ] ; then echo "Simulated_scaffold,1,10" >> Coordinates_alrdy_examined.csv ; fi
	tr ',' '\t' < Coordinates_alrdy_examined.csv > Coordinates_already_examined.tsv


	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Complete_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa Complete_genes.fa ; rm temp.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa Pseudogenes.fa ; rm temp.fa
	awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_wo_Fs.prot | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa Pseudogenes_wo_Fs.prot ; rm temp.fa


	cat Complete_genes.fa Pseudogenes.fa > All_sequences.fa
	current_nb_sequences=`if test -f "All_sequences.fa" ; then grep -c ">" All_sequences.fa ; else echo "0" ;  fi`
	
	#re-process blast result to find potential V2R regions exclusing already found genes
	
	Rscript $scripts_location/Rscript_merge_filter_extend_blast_hit_Second.R $maximum_intron_length

	if test -f "Potential_regions.tsv" ; then number_regions_blast=`grep "[0-9]" Potential_regions.tsv | wc -l` ; else number_regions_blast=0 ; fi
	rm Parsed_miniprot_rslt.tsv


	cp $target_DB ./Target_protein_db.prot 

	transeq Complete_genes.fa Complete_genes.prot
	

done


