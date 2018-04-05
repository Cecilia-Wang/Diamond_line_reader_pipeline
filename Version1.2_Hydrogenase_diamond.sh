#!/bin/bash

# Enter the fastq file pathways
echo Hello, please enter the full pathway to the input files

read Input_DIR

 if test -d $Input_DIR
then echo Input files will be load from $Input_DIR 

else
        echo $Input_DIR does not exist, please rerun the bash file and enter the complete pathway to the input files 
        exit
fi


# Enter Hydrogenase database pathway
echo Please enter the full pathway to the database reference

read HydroDB_DIR

 if test -d $HydroDB_DIR
then echo Database reference files will be load from $HydroDB_DIR

else
        echo $HydroDB_DIR does not exist, please rerun the bash file and enter the complete pathway to the database reference file 
        exit
fi




        # Subsampling sequencing data to 5 million reads per sample

cd $Input_DIR

if test -d $Input_DIR/Subset_files
then
        read -p "$Input_DIR/Subset_files exists, do you want to remove the directory?" yn
    case $yn in
        [Yy]* ) rm -rf $Input_DIR/Subset_files; mkdir $Input_DIR/Subset_files;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac

else
        mkdir $Input_DIR/Subset_files
fi



for files in *.fastq #.gz

do seqtk sample -s100 $files 5000000 > $Input_DIR/Subset_files/${files}_sub.fastq

done

cd $Input_DIR/Subset_files

for Sub in *.fastq

do seqtk seq -a $Sub > ${Sub}.fa

done


# DIAMOND DATABASE ANALYSIS OF SUBSAMPLED READS

        #make diamond database

cd $HydroDB_DIR

diamond makedb --in $HydroDB_DIR/Hydrogenase_Updated_Fixed_GadMinus.fasta -d Hydrogenase_Latest --threads 8

if test -d $Input_DIR/Hydrogenase_blast
then
        read -p "$Input_DIR/Hydrogenase_blast exists, do you want to remove the directory?" yn
    case $yn in
        [Yy]* ) rm -rf $Input_DIR/Hydrogenase_blast; mkdir $Input_DIR/Hydrogenase_blast;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac

else
        mkdir $Input_DIR/Hydrogenase_blast
fi


cd $Input_DIR/Subset_files


for Subfa in *.fa

do nice diamond blastx -d $HydroDB_DIR/Hydrogenase_Latest.dmnd -q $Subfa -a $Input_DIR/Hydrogenase_blast/${Subfa}_reads -t $Input_DIR/Hydrogenase_blast --max-target-seqs 1 --threads 8

done

cd $Input_DIR/Hydrogenase_blast

for Subread in *.daa

do nice diamond view -a $Subread -o ${Subread}.m8 --threads 8
done

# FILTERING AND SUMMARISING DIAMOND OUTPUT 
        #       ≥60% sequence identity
        #       ≥40 amino acid length

for daafile in *.daa

do diamond view -a $daafile | awk '{if ($3 >= 60 && $4 >=40) print $0}' | cut -f 2 | sort | uniq -c | sort -r -k 1 > ${daafile}_summary.txt

done

        # total reads summary
cat *.m8 | awk '{if ($3 >= 60 && $4 >=40) print $0}' | cut -f 2 | sort | uniq -c | sort -n -k 1 > total_reads_summary.txt

        # the further low vs high emission summary will need to type the file names manually and change *.m8

