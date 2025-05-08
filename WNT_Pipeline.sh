#!/bin/bash


#SBATCH --job-name=Wnt_mining   # Job name



#conda = miniprot ; cd-hit ; seqkit

eval "$(conda shell.bash hook)"
conda activate miniprot


LMOD_DISABLE_SAME_NAME_AUTOSWAP=no

module purge
module load R/4.2.0-foss-2021a
module load BLAST/2.12.0-Linux_x86_64
module load EMBOSS/6.2.0-goolf-1.7.20
module load SAMtools/1.15-GCC-10.3.0
module load MAFFT/7.467-GCCcore-7.3.0-with-extensions
module load IQ-TREE/2.0-rc1-foss-2018b
module load Python/3.9.5-GCCcore-10.3.0
module load FASTX-Toolkit/0.0.14-goolf-1.7.20



dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"



genome=$1 
target_DB=$2 
target_outgroup_DB=$3 
scripts_folder_location=$4 ; scripts_location=$( echo "$scripts_folder_location" | sed 's/\/$//' )
maximum_intron_length=$5
prot_name=$6
number_of_threads=$7



#genome=/scicore/home/salzburg/polica0000/Vertebrates_Taste_Receptors/Danio_rerio/GCA_000002035.4_GRCz11_genomic.trimmed.fna
##target_DB=/scicore/home/salzburg/polica0000/WNT_Database/Wnt_DB.fa
#target_DB=/scicore/home/salzburg/polica0000/WNT_Database/Database_creation/Candidate_WNT_genes_filtered_length.prot
#target_outgroup_DB=/scicore/home/salzburg/polica0000/Uniprot_db/uniprot_sprot.fasta
#scripts_folder_location=/scicore/home/salzburg/polica0000/WNT_Database/ ; scripts_location=$( echo "$scripts_folder_location" | sed 's/\/$//' )
#maximum_intron_length=60000
#number_of_threads=8
#prot_name=Wnt



#Makeblastdb so we can blast genes against the genome

if test -f "$genome.ndb" ; then echo "Genome blast database already exist" ; else makeblastdb -in $genome -dbtype nucl ; fi 
if test -f "$genome.fai" ; then echo "Genome fai file already exist" ; else samtools faidx $genome ; fi 


#Perform tblastn using known V2R genes against the genome with an evalue of 1e-5

cp $target_DB ./Target_protein_db.prot
makeblastdb -in Target_protein_db.prot -dbtype prot
samtools faidx Target_protein_db.prot 



cd-hit -i Target_protein_db.prot -o Target_protein_db_c80.prot -c 0.9
tblastn -query Target_protein_db_c80.prot -db $genome -evalue 1e-3 -outfmt 6 -out Prot_vs_Genome.blastn -num_threads $number_of_threads

#tblastn -query Target_protein_db.prot -db $genome -evalue 1e-3 -outfmt 6 -out Prot_vs_Genome.blastn -num_threads $number_of_threads

#Lets launch a Rscript that will merge all blast hits 

Rscript $scripts_location/Rscript_merge_blast_hits.R

xargs samtools faidx $genome < Blast_nonoverlapping.tsv > Blast_nonoverlapping.fasta
cp Blast_nonoverlapping.tsv Prot_Regions.tsv


#Extend all best hits by upstream and downstream . Result file : Potential_regions.tsv
Rscript $scripts_location/Rscript_merge_filter_extend_blast_hit.R $maximum_intron_length


rm Coordinates_alrdy_examined.csv
rm Coordinates_already_examined.tsv
rm Pseudogenes.fa
rm Pseudogenes_wo_Fs.prot
rm Complete_genes.fa

echo "Simulated_scaffold,1,10" > Coordinates_alrdy_examined.csv


echo "Launch very strict miniprot loop + First relaxed loop"
$scripts_location/Launch_strict_miniprot_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 70000 &> error.txt
$scripts_location/Launch_relaxed_miniprot_first_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 70000 &>> error.txt


$scripts_location/Launch_strict_miniprot_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 40000  &>> error.txt
$scripts_location/Launch_relaxed_miniprot_first_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 40000 &>> error.txt

$scripts_location/Launch_strict_miniprot_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 20000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_first_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 20000 &>> error.txt


$scripts_location/Launch_strict_miniprot_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 10000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_first_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 10000 &>> error.txt

$scripts_location/Launch_strict_miniprot_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_first_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt


$scripts_location/Launch_strict_miniprot_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_first_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt

echo "Launch relaxed loop - Second protein filter"
$scripts_location/Launch_relaxed_miniprot_second_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 70000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_second_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 40000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_second_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 20000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_second_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 10000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_second_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_second_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt

echo "Launch relaxed loop - Third protein filter"

$scripts_location/Launch_relaxed_miniprot_third_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 40000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_third_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 20000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_third_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 10000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_third_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_third_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt

echo "Launch relaxed loop - Fourth protein filter"

$scripts_location/Launch_relaxed_miniprot_fourth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 20000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_fourth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 10000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_fourth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_fourth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt

echo "Launch relaxed loop - Fifth protein filter"

$scripts_location/Launch_relaxed_miniprot_fifth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 10000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_fifth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_fifth_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt

echo "Launch relaxed loop - Six protein filter"

$scripts_location/Launch_relaxed_miniprot_six_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 5000 &>> error.txt
$scripts_location/Launch_relaxed_miniprot_six_loop.sh $genome $target_DB $target_outgroup_DB $scripts_folder_location $maximum_intron_length $prot_name $number_of_threads 1000 &>> error.txt



echo "Lets apply a last filter sequences"


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Complete_genes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Complete_genes_uniq.fa 
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_wo_Fs.prot | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > Pseudogenes_wo_Fs_uniq.prot


#perform a BLASTP for functional genes and pseudogenes

blastx -query Complete_genes_uniq.fa -db $target_outgroup_DB -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out functional_verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads
grep -i "$prot_name" functional_verif_blastx_result | cut -f1 | sort | uniq > good_blast_complete.txt
xargs samtools faidx Complete_genes_uniq.fa < good_blast_complete.txt > Complete_genes_uniq_filtered.fa
rm *.fai


blastx -query Pseudogenes_uniq.fa -db $target_outgroup_DB -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads
grep -i "$prot_name" verif_blastx_result | cut -f1 | sort | uniq > good_blast_pseudos.txt
xargs samtools faidx Pseudogenes_uniq.fa < good_blast_pseudos.txt > Pseudogenes_uniq_filtered.fa
xargs samtools faidx Pseudogenes_wo_Fs.fa < good_blast_pseudos.txt > Pseudogenes_wo_Fs_filtered.fa




## Now lets filter these files : remove ambigous sequences 


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Complete_genes_uniq_filtered.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Complete_genes_uniq_filtered.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_uniq_filtered.fa | sed 's/\*$//g' | awk -F '\t'  '!($2 ~ /\N/)' | tr "\t" "\n" > clear_Pseudogenes_uniq_filtered.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Complete_genes_uniq_filtered.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Complete_genes_uniq_filtered.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Pseudogenes_uniq_filtered.fa | sed 's/\*$//g' | awk -F '\t'  '($2 ~ /\N/)' | tr "\t" "\n" > unclear_Pseudogenes_uniq_filtered.fa

grep ">" clear_Pseudogenes_uniq_filtered.fa | sed 's/>//g' > clear_id 
grep ">" unclear_Pseudogenes_uniq_filtered.fa | sed 's/>//g' > unclear_id 

xargs samtools faidx Pseudogenes_wo_Fs_filtered.fa < clear_id > clear_Pseudogenes_wo_Fs_filtered.fa
xargs samtools faidx Pseudogenes_wo_Fs_filtered.fa < unclear_id > unclear_Pseudogenes_wo_Fs_filtered.fa

cat unclear_Complete_genes_uniq_filtered.fa unclear_Pseudogenes_uniq_filtered.fa > FINAL_Ambigous.fa


# If two genes/pseudogenes are overlapping due to the multiple extensions, then keep only the longest found gene/pseudogene


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' FINAL_Ambigous.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > FINAL_Ambigous_uniq.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Pseudogenes_uniq_filtered.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa clear_Pseudogenes_uniq_filtered.fa ; rm temp.fa
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' clear_Complete_genes_uniq_filtered.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > temp.fa ; mv temp.fa clear_Complete_genes_uniq_filtered.fa



nb_seq=`grep -c ">" FINAL_Ambigous_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" FINAL_Ambigous_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambigous_final.tsv
fi

nb_seq=`grep -c ">" clear_Pseudogenes_uniq_filtered.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" clear_Pseudogenes_uniq_filtered.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
fi

nb_seq=`grep -c ">" clear_Complete_genes_uniq_filtered.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" clear_Complete_genes_uniq_filtered.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
fi


Rscript $scripts_location/Remove_redundancy.R



IFS=$'\n'

for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" clear_Complete_genes_uniq_filtered.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" FINAL_Ambigous_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" clear_Pseudogenes_uniq_filtered.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done


if test -f "functionnal_to_keep.txt" ; then xargs samtools faidx clear_Complete_genes_uniq_filtered.fa < functionnal_to_keep.txt > FINAL_Complete.fa ; else echo "" > FINAL_Complete.fa ; fi 
if test -f "ambigous_to_keep.txt" ; then xargs samtools faidx FINAL_Ambigous_uniq.fa < ambigous_to_keep.txt > FINAL_Ambigous.fa ; else echo "" > FINAL_Ambigous.fa ; fi 
if test -f "pseudogenes_to_keep.txt" ; then xargs samtools faidx clear_Pseudogenes_uniq_filtered.fa < pseudogenes_to_keep.txt > FINAL_Pseudogene.fa ; else echo "" > FINAL_Pseudogene.fa ; fi 
if test -f "pseudogenes_to_keep.txt" ; then xargs samtools faidx clear_Pseudogenes_wo_Fs_filtered.prot < pseudogenes_to_keep.txt > FINAL_Pseudogene.prot ; else echo "" > FINAL_Pseudogene.prot ; fi 

nb_complete_genes=`grep -c ">" FINAL_Complete.fa`
nb_ambigous_genes=`grep -c ">" FINAL_Ambigous.fa`
nb_pseudogenes=`grep -c ">" FINAL_Pseudogene.fa`

echo "Results of the pipeline = "

echo "$nb_complete_genes complete genes"
echo "$nb_pseudogenes incomplete genes"
echo "$nb_ambigous_genes ambiguous genes"



dt=$(date '+%d/%m/%Y %H:%M:%S');
echo "$dt"




################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################## Correction of miniprot pseuogene predictions with exonerate  ################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################


cp $target_DB ./Target_protein_db.prot
makeblastdb -in Target_protein_db.prot -dbtype prot
samtools faidx Target_protein_db.prot 


samtools faidx FINAL_Pseudogene.fa
grep ">" FINAL_Pseudogene.fa | sed "s/>//g" >  FINAL_Pseudogene.id
cut -f1,2,3 -d "-" FINAL_Pseudogene.id | sed 's/-/,/g' > Coordinates_pseudogenes.csv

grep ">" FINAL_Complete.fa | sed 's/>//g' | cut -f1,2,3 -d "-" | sed 's/-/,/g' > Coordinates_completes.csv
grep ">" FINAL_Ambigous.fa | sed 's/>//g' | cut -f1,2,3 -d "-" | sed 's/-/,/g' >> Coordinates_completes.csv

Rscript $scripts_location/Coordinates_pseudo_refine.R $maximum_intron_length



rm Exonerate_complete.fa
IFS=$'\n'

for new_coordinates in `cat new_coordinates_table_filtered.csv` ;  do 

	scaffold=`echo "$new_coordinates" | cut -f1 -d ","`
	start=`echo "$new_coordinates" | cut -f2 -d ","`
	stop=`echo "$new_coordinates" | cut -f3 -d ","`


	rm Region.fa.fai
	rm Region.renamed.fa.fai
	samtools faidx $genome $scaffold:$start-$stop > Region.fa

	blastx -query Region.fa -db Target_protein_db.prot -evalue 1e-5 -outfmt "6 qseqid sseqid" -out blastx_result -max_target_seqs 1 -num_threads $number_of_threads
	best_query=`head -1 blastx_result | cut -f2`
	samtools faidx Target_protein_db.prot $best_query > best_query.prot

	$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE best_query.prot Region.fa > Verif_prediction.exo


	if [ `grep -c "Query: " Verif_prediction.exo` -ge 2 ] ; then
		sed '/^C4 Alignment:/,/^# --- END OF GFF DUMP ---/!d;/^# --- END OF GFF DUMP ---/q'  Verif_prediction.exo > first_result_infos
		sed '/^# --- END OF GFF DUMP ---/,/^C4 Alignment:/!d;/^C4 Alignment:/q'  Verif_prediction.exo > first_result_sequence
		cat first_result_infos first_result_sequence > Verif_prediction.exo
	fi


	first_hit_range=`grep -m1 "Target range:" Verif_prediction.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
	second_hit_range=`grep -m1 "Target range:" Verif_prediction.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			

	strand=`grep "	similarity	" Verif_prediction.exo  | cut -f7`

	if [ $strand == "-" ] ; then 
			

		strand_print="minus"
	
		target_end=$((first_hit_range + 1))
		target_extanded_end=$((first_hit_range + 200))
		sed 's/>.*/>myregion/g' Region.fa > Region.renamed.fa


		rm Region.renamed.fa.fai
		samtools faidx Region.renamed.fa myregion:$target_end-$target_extanded_end > Extend_three_prime.fa
		samtools faidx Region.renamed.fa myregion:1-$second_hit_range > Extend_five_prime.fa
			
		grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
		grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
		grep "	exon	"  Verif_prediction.exo | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
		for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Region.renamed.fa myregion:$begin_exon-$end_exon >> Correct_cds.fa ; done
		grep -v ">" Correct_cds.fa > predicted_cds_rev.txt ; rm Correct_cds.fa 
	
		cat Extend_five_prime.txt predicted_cds_rev.txt Extend_three_prime.txt > Complete_extanded_sequence.fa
		sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
		rm Exonerate_pred.ORF*

		query_length=`grep "$best_query" Target_protein_db.prot.fai | cut -f2`
	 	perc80_query_length=$((query_length*75/100*3))

		getorf -sequence Complete_extanded_sequence.fa -outseq Exonerate_pred.ORF -minsize $perc80_query_length -find 3
		rm Exonerate_pred.ORF.fai ; if grep -q -i "reverse" Exonerate_pred.ORF ; then sequence_to_grep=`grep -i "reverse" Exonerate_pred.ORF | sed 's/>//g' | sed 's/ .*//g'`  ; samtools faidx Exonerate_pred.ORF $sequence_to_grep > temporary ; mv temporary Exonerate_pred.ORF ; rm Exonerate_pred.ORF.fai ; else echo "no_gene" > Exonerate_pred.ORF ; fi
			
	
		if [ `grep -c ">" Exonerate_pred.ORF` -ge 1 ] ; then
			
			rm Exonerate_pred.ORF.fai
			samtools faidx Exonerate_pred.ORF

			blastx -query Exonerate_pred.ORF -db $target_outgroup_DB -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out functional_verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads
			grep -i "$prot_name" functional_verif_blastx_result | cut -f1 | sort | uniq > good_blast_complete.txt
			xargs samtools faidx Exonerate_pred.ORF < good_blast_complete.txt > temp.fa ; mv temp.fa Exonerate_pred.ORF
			rm *.fai ; samtools faidx Exonerate_pred.ORF


			exonerate_gene_length=`head -1 Exonerate_pred.ORF.fai | cut -f2`
			exonerate_gene_name=`head -1 Exonerate_pred.ORF.fai | cut -f1`
			samtools faidx Exonerate_pred.ORF $exonerate_gene_name > temp.fa ; mv temp.fa Exonerate_pred.ORF ; rm Exonerate_pred.ORF.fai
			transeq Exonerate_pred.ORF Exonerate_pred.ORFP

			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Exonerate_pred.ORFP Region.fa > verif_coord.exo
			if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Exonerate_pred.ORFP Region.fa > verif_coord.exo ; fi

			extracted_scaffold_start=$start

			cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			cds_coord_start=$((extracted_scaffold_start + cds_start_extract - 50))
			cds_coord_end=$((extracted_scaffold_start + cds_end_extract - 1 + 50))
			exon_number=`grep "	exon	" verif_coord.exo | wc -l`
			sed "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end-$exon_number-$strand_print/g" Exonerate_pred.ORF >> Exonerate_complete.fa

	
		fi


	
	else
			

		strand_print="plus"

		target_end=$((second_hit_range + 1))
		target_extanded_end=$((second_hit_range + 200))
		sed 's/>.*/>myregion/g' Region.fa > Region.renamed.fa
			

		rm Region.renamed.fa.fai
		samtools faidx Region.renamed.fa myregion:1-$first_hit_range > Extend_three_prime.fa
		samtools faidx Region.renamed.fa myregion:$target_end-$target_extanded_end > Extend_five_prime.fa
			
		grep -v ">" Extend_three_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_three_prime.txt
		grep -v ">" Extend_five_prime.fa | tr -d '\n' | sed 's/NNN//g' > Extend_five_prime.txt
				
		grep "	exon	"  Verif_prediction.exo | cut -f3,4,5 | sort -n -k2 > target_seq.tsv
		for exons in `cat target_seq.tsv` ; do begin_exon=`echo "$exons" | cut -f2` ; end_exon=`echo "$exons" | cut -f3` ; samtools faidx Region.renamed.fa myregion:$begin_exon-$end_exon >> Correct_cds.fa ; done
		grep -v ">" Correct_cds.fa > predicted_cds.txt ; rm Correct_cds.fa
			
		cat Extend_three_prime.txt predicted_cds.txt Extend_five_prime.txt > Complete_extanded_sequence.fa
		sed -i '1 i\>Complete_seq' Complete_extanded_sequence.fa
		rm Exonerate_pred.ORF*

		query_length=`grep "$best_query" Target_protein_db.prot.fai | cut -f2`
	 	perc80_query_length=$((query_length*75/100*3))


		getorf -sequence Complete_extanded_sequence.fa -outseq Exonerate_pred.ORF -minsize $perc80_query_length -find 3 -reverse FALSE
			
		

		if [ `grep -c ">" Exonerate_pred.ORF` -ge 1 ] ; then
			
			rm Exonerate_pred.ORF.fai
			samtools faidx Exonerate_pred.ORF


			blastx -query Exonerate_pred.ORF -db $target_outgroup_DB -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out functional_verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads
			grep -i "$prot_name" functional_verif_blastx_result | cut -f1 | sort | uniq > good_blast_complete.txt
			xargs samtools faidx Exonerate_pred.ORF < good_blast_complete.txt > temp.fa ; mv temp.fa Exonerate_pred.ORF
			rm *.fai ; samtools faidx Exonerate_pred.ORF

			exonerate_gene_length=`head -1 Exonerate_pred.ORF.fai | cut -f2`
			exonerate_gene_name=`head -1 Exonerate_pred.ORF.fai | cut -f1`

			samtools faidx Exonerate_pred.ORF $exonerate_gene_name > temp.fa ; mv temp.fa Exonerate_pred.ORF ; rm Exonerate_pred.ORF.fai
			transeq Exonerate_pred.ORF Exonerate_pred.ORFP

			$scripts_location/exonerate-2.2.0-x86_64/bin/exonerate -E True --model protein2genome:bestfit --bestn 1 --showtargetgff TRUE Exonerate_pred.ORFP Region.fa > verif_coord.exo
			if grep -q "Query range:" verif_coord.exo ; then echo "No segmentation default" ; else $scripts_location/exonerate-2.2.0-x86_64/bin/exonerate --model protein2genome --bestn 1 --showtargetgff TRUE Exonerate_pred.ORFP Region.fa > verif_coord.exo ; fi

			extracted_scaffold_start=$start
			cds_start_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f1 -d " "`
			cds_end_extract=`grep -m1 "Target range:" verif_coord.exo | sed 's/^ *//g' | sed 's/Target range://g' | sed 's/ //g' | sed 's/->/ /g' | cut -f2 -d " "`
			cds_coord_start=$((extracted_scaffold_start + cds_start_extract - 50))
			cds_coord_end=$((extracted_scaffold_start + cds_end_extract + 50))
			exon_number=`grep "	exon	" verif_coord.exo | wc -l`
			sed "s/>.*/>$scaffold-$cds_coord_start-$cds_coord_end-$exon_number-$strand_print/g" Exonerate_pred.ORF >> Exonerate_complete.fa



		fi


	fi

done

#Remove empty genes

awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' Exonerate_complete.fa > temp.fa ; mv temp.fa Exonerate_complete.fa



## Verify exonerate predicted genes

blastx -query Exonerate_complete.fa -db $target_outgroup_DB -outfmt "6 qseqid sseqid stitle sblastnames sgi sacc" -out functional_verif_blastx_result -max_target_seqs 1 -num_threads $number_of_threads
grep -i "$prot_name" functional_verif_blastx_result | cut -f1 | sort | uniq > good_blast_complete.txt
xargs samtools faidx Exonerate_complete.fa < good_blast_complete.txt > Exonerate_complete_filtered.fa
rm *.fai



## Now lets correct final coordinates


awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Exonerate_complete_filtered.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > FINAL_Exonerate_complete.fa

cat FINAL_Exonerate_complete.fa FINAL_Complete.fa > Exonerate_plus_miniprot_complete.fa

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' Exonerate_plus_miniprot_complete.fa | sort -t $'\t' -k1,1 -u | tr "\t" "\n" > FINAL_Exonerate_plus_miniprot_complete.fa




nb_seq=`grep -c ">" FINAL_Ambigous_uniq.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" FINAL_Ambigous_uniq.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_ambigous_final.tsv
fi

nb_seq=`grep -c ">" clear_Pseudogenes_uniq_filtered.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" clear_Pseudogenes_uniq_filtered.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_pseudogenes_final.tsv
fi

nb_seq=`grep -c ">" FINAL_Exonerate_plus_miniprot_complete.fa`
if [ "$nb_seq" -gt "0" ] ; then
	grep ">" FINAL_Exonerate_plus_miniprot_complete.fa | sed 's/>//g' | sed 's/-/	/g' | cut -f1,2,3 > Coordinates_genes_final.tsv
fi


Rscript $scripts_location/Remove_redundancy.R




IFS=$'\n'

rm functionnal_to_keep.txt ; rm ambigous_to_keep.txt ; rm pseudogenes_to_keep.txt 
for line in `cat best_genes_functionnal.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" FINAL_Exonerate_plus_miniprot_complete.fa | sed 's/>//g' >> functionnal_to_keep.txt ; done
for line in `cat best_genes_ambigous.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" FINAL_Ambigous_uniq.fa | sed 's/>//g' >> ambigous_to_keep.txt ; done
for line in `cat best_genes_pseudogenes.tsv` ; do scaff=`echo "$line" | cut -f1` ; start=`echo "$line" | cut -f2` ; end=`echo "$line" | cut -f3` ; grep -m1 "$scaff.*$start.*$end" clear_Pseudogenes_uniq_filtered.fa | sed 's/>//g' >> pseudogenes_to_keep.txt ; done


if test -f "functionnal_to_keep.txt" ; then xargs samtools faidx FINAL_Exonerate_plus_miniprot_complete.fa < functionnal_to_keep.txt > FINAL_Complete.fa ; else echo "" > FINAL_Complete.fa ; fi 
if test -f "ambigous_to_keep.txt" ; then xargs samtools faidx FINAL_Ambigous_uniq.fa < ambigous_to_keep.txt > FINAL_Ambigous.fa ; else echo "" > FINAL_Ambigous.fa ; fi 
if test -f "pseudogenes_to_keep.txt" ; then xargs samtools faidx clear_Pseudogenes_uniq_filtered.fa < pseudogenes_to_keep.txt > FINAL_Pseudogene.fa ; else echo "" > FINAL_Pseudogene.fa ; fi 
if test -f "pseudogenes_to_keep.txt" ; then xargs samtools faidx clear_Pseudogenes_wo_Fs_filtered.prot < pseudogenes_to_keep.txt > FINAL_Pseudogene.prot ; else echo "" > FINAL_Pseudogene.prot ; fi 

nb_complete_genes=`grep -c ">" FINAL_Complete.fa`
nb_ambigous_genes=`grep -c ">" FINAL_Ambigous.fa`
nb_pseudogenes=`grep -c ">" FINAL_Pseudogene.fa`

echo "Results of the pipeline = "

echo "$nb_complete_genes complete genes"
echo "$nb_pseudogenes incomplete genes"
echo "$nb_ambigous_genes ambiguous genes"



