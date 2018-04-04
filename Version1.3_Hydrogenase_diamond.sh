#!/bin/bash

# using a yesno function for giving an default answer and time out option to some of the folder remove questions 

# ============================ Here is the function ===========================#
function yesno()
{
    local ans
    local ok=0
    local timeout=0
    local default
    local t

    while [[ "$1" ]]
    do
        case "$1" in
        --default)
            shift
            default=$1
            if [[ ! "$default" ]]; then error "Missing default value"; fi
            t=$(tr '[:upper:]' '[:lower:]' <<<$default)

            if [[ "$t" != 'y'  &&  "$t" != 'yes'  &&  "$t" != 'n'  &&  "$t" != 'no' ]]; then
                error "Illegal default answer: $default"
            fi
            default=$t
            shift
            ;;

        --timeout)
            shift
            timeout=$1
            if [[ ! "$timeout" ]]; then error "Missing timeout value"; fi
            if [[ ! "$timeout" =~ ^[0-9][0-9]*$ ]]; then error "Illegal timeout value: $timeout"; fi
            shift
            ;;

        -*)
            error "Unrecognized option: $1"
            ;;

        *)
            break
            ;;
        esac
    done

    if [[ $timeout -ne 0  &&  ! "$default" ]]; then
        error "Non-zero timeout requires a default answer"
    fi

    if [[ ! "$*" ]]; then error "Missing question"; fi

    while [[ $ok -eq 0 ]]
    do
        if [[ $timeout -ne 0 ]]; then
            if ! read -t $timeout -p "$*" ans; then
                ans=$default
            else
                # Turn off timeout if answer entered.
                timeout=0
                if [[ ! "$ans" ]]; then ans=$default; fi
            fi
        else
            read -p "$*" ans
            if [[ ! "$ans" ]]; then
                ans=$default
            else
                ans=$(tr '[:upper:]' '[:lower:]' <<<$ans)
            fi 
        fi

        if [[ "$ans" == 'y'  ||  "$ans" == 'yes'  ||  "$ans" == 'n'  ||  "$ans" == 'no' ]]; then
            ok=1
        fi

        if [[ $ok -eq 0 ]]; then warning "Valid answers are: yes y no n"; fi
    done
    [[ "$ans" = "y" || "$ans" == "yes" ]]
}

# ============================ Function finishes here ===========================#


# Enter the fastq file pathways
echo Hello, please enter the full pathway to the input files 

read Input_DIR

 if test -d $Input_DIR
then echo Input files will be load from $Input_DIR 
	echo ""

else
	echo $Input_DIR does not exist, please rerun the bash file and enter the complete pathway to the input files
	exit
fi


# Enter Hydrogenase database pathway
echo Please enter the full pathway to the database reference 

read HydroDB_DIR

 if test -d $HydroDB_DIR
then echo Database reference files will be load from $HydroDB_DIR 
	echo ""

else
	echo $HydroDB_DIR does not exist, please rerun the bash file and enter the complete pathway to the database reference file 
	exit
fi




	# Subsampling sequencing data to 5 million reads per sample

cd $Input_DIR

if test -d $Input_DIR/Subset_files
then 
	if yesno --timeout 5 --default yes "$Input_DIR/Subset_files exists, do you want to remove the directory?Yes or no (timeout 5, default yes) ? "; then
        
        rm -rf $Input_DIR/Subset_files
        mkdir $Input_DIR/Subset_files
        echo ""
        echo " You answered yes, existing directory removed, new directory made"
    else
        echo " You answered no, you will exit the program"
        exit
    fi
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
	if yesno --timeout 5 --default yes "$Input_DIR/Hydrogenase_blast, do you want to remove the directory?Yes or no (timeout 5, default yes) ? "; then
        
        rm -rf $Input_DIR/Hydrogenase_blast
        mkdir $Input_DIR/Hydrogenase_blast
        echo ""
        echo " You answered yes, existing directory removed, new directory made"
    else
        echo " You answered no, you will exit the program"
        exit
    fi  

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
	#	≥60% sequence identity
	#	≥40 amino acid length

for daafile in *.daa

do diamond view -a $daafile | awk '{if ($3 >= 60 && $4 >=40) print $0}' | cut -f 2 | sort | uniq -c | sort -r -k 1 > ${daafile}_summary.txt

done

	# total reads summary
cat *.m8 | awk '{if ($3 >= 60 && $4 >=40) print $0}' | cut -f 2 | sort | uniq -c | sort -n -k 1 > total_reads_summary.txt

	# the further low vs high emission summary will need to type the file names manually and change *.m8









