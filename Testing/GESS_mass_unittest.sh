#!/bin/bash

#H5 file location is defined
h5_file_location="FCA_H5s"

all_h5_files=( $( ls $h5_file_location ) )

#Outdir is where the resultant testing logs will be stored
outdir="GESS_BATCH_OUTPUTS"

#Each supplied H5 file is processed seperately
for h5file in ${all_h5_files[@]}; do

    parsed_name=(${h5file//_/ })

    h5tissue=${parsed_name[3]}

    full_file_name=$h5file
    generic_h5_file="ALL_H5_GENES.txt"
    
    #All genes in the supplied H5 file are extracted
    python GESS/Utils/GetH5genes.py -h5 $full_file_name

    readarray -t allgenes < $generic_h5_file
    array_length=${#allgenes[@]}

    exptypes=("expression" "prevalence")

    for runtype in ${exptypes[@]}; do

        logfile="$outdir/"$h5tissue"_"$runtype"_2_.log"

        for ((i = 1; i <= 100; i++)); do
            #We gather random genes from the supplied H5 files, then carry out GESS  using those genes within that file
            shuffled_genes=($(shuf -e "${allgenes[@]}"))

            index1=$((RANDOM % array_length))
            index2=$((RANDOM % array_length))

            echo "--- Iteration #$i ---"
            RANDOM=$$$(date +%s)
            gene1=${shuffled_genes[index1]}
            
            RANDOM=$$$(date +%s*2)
            gene2=${shuffled_genes[index2]}

            echo $gene1
            echo $gene2
	    
	    #GESS is run and run time stored in log file
            time python Scripts/GESS_Runner_Mar24UPDATE.py -qg $gene1 -tg $gene2 -qd $full_file_name -td $full_file_name -anno annotation annotation_broad -h5mode $runtype
            echo

        done 2>&1 | tee $logfile
    done

    rm $generic_h5_file 

done
