#!/bin/bash

# Enter the fastq file pathways
 
Input_DIR=


	# Subsampling sequencing data to 5 million reads per sample

cd $Input_DIR

mkdir $Input_DIR/Subset_files

for files in *.fastq.gz

do seqtk sample -s100 $files 5000000 > $Input_DIR/Subset_files/${files}_sub.fastq

done

cd $Input_DIR/Subset_files

for Sub in *.fastq

do seqtk seq -a $Sub > ${Sub}.fa

done


# DIAMOND DATABASE ANALYSIS OF SUBSAMPLED READS

# Enter Hydrogenase database pathway

HydroDB_DIR=

	#make diamond database
	
cd $HydroDB_DIR

diamond makedb --in $HydroDB_DIR/Hydrogenase_Updated_Fixed_GadMinus.fasta -d Hydrogenase_Latest --threads 8



mkdir $Input_DIR/Hydrogenase_blast
mkdir $Input_DIR/Hydrogenase_blast/Blast_output

cd $Input_DIR/Subset_files


for Subfa in *.fa

do nice diamond blastx -d $HydroDB_DIR/Hydrogenase_Latest.dmnd -q $Subfa -a $Input_DIR/Hydrogenase_blast/Blast_output/${Subfa}_reads -t $Input_DIR/Hydrogenase_blast/Blast_output --max-target-seqs 1 --threads 8

done

cd $Input_DIR/Hydrogenase_blast/Blast_output

for Subread in *.daa

do nice diamond view -a $Subread -o ${Subread}.m8 --threads 8
done

##FILTERING AND SUMMARISING DIAMOND OUTPUT
#       >= 65% sequence identity
#       >= 40 amino acid length

for daafile in *.daa
do diamond view -a $daafile | awk '{if ($3 >= 65 && $4 >= 40) print $0}' | cut -f 2,3 | sort | uniq -c | sort -r -k 1 > $Output_DIR/Final/${daafile}_summary.txt
done

cd $Output_DIR/Final

#Total reads summary
# cat *.m8 | awk '{if ($3 >= 65 && $4 >= 40) print $0}' | cut -f 2,3 | sort | uniq -c | sort -n -k 1 > $Output_DIR/Final/total_reads_summary.txt
awk '{print $0,"\t", FILENAME}' Rank*.txt >Total_reads_sum.txt
#the further low vs high emission summary will need to type the file names manually and change *.m8








