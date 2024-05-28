#!bin/bash

generic_h5_file="ALL_H5_GENES.txt"

#All genes in the FlyCellAtlas "Whole fly" dataset are gathered into an array
fca_file="r_fca_biohub_all_wo_blood_10x.h5ad"
python GESS/Utils/GetH5genes.py -h5 $fca_file
readarray -t fcagenes < $generic_h5_file

#The same is done for Aging Fly Cell Atlas, using the 5 days dataset (for simplicity, genes don't differ between datasets unlike FlyCellAtlas)
python GESS/Utils/GetH5genes.py -h5 AgingAtlas_H5s/_5days_ForGESS.h5ad
readarray -t aginggenes < $generic_h5_file

#Only genes which can be found in both datasets are analysed. Thus, we find gene names present in both arrays
findable_all=$( echo ${fcagenes[@]} ${aginggenes[@]} | tr " " "\n" | sort | uniq -d )
findable=(${findable_all})
echo ${#findable[@]}

#Aging Fly Cell Atlas files are found within the AgingAtlas_H5s folder
h5_file_location='AgingAtlas_H5s'
all_aging_files=( $( ls $h5_file_location ) )

i=1

#For each individual gene, we find its GESS between FlyCellAtlas and each individual Aging Atlas age group
for gene in ${findable[@]}; do
    echo "Gene " $i " of " ${#findable[@]}
    let "i++"

    generesults=($gene)
    
    #Each individual Aging Fly Cell Atlas age group is processed seperately
    for agingfile in ${all_aging_files[@]}; do

        parsed_name=(${agingfile//_/ })
        actual_age=${parsed_name[-2]}
	
	#Calculates the each gene's GESS between the FlyCellAtlas/Aging Atlas datasets  across annotation and annotation_broad 
	#	- Tissue is not useful within Aging Fly Cell Atlas, in which is is only defined as "Head" and "Whole"
        gess_results=$( python GESS/GESSfinder.py -qg $gene -tg $gene -qd $fca_file -td $agingfile -anno annotation annotation_broad -h5mode expression )
	
	#Results are stored to a results array
        gess_parse_1=(${gess_results})
        gess_percent=${gess_parse_1[-1]}
        gess_figure=${gess_percent::-1}

        generesults+=($gess_figure)

    done
    
    #results are printed in tab-separated format, one row per gene 
    printf '%s\t' "${generesults[@]}">> 'agingdata.tsv'
    echo '' >>agingdata.csv

done
