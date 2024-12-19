#!/bin/bash

# For use downstream of RNASeq processing pipeline

# Get info
type=${1}
strand=${2}

if [[ "${type}" == "Gene" ]]
then

    header_row=4
    file_pattern=_ReadsPerGene.out.tab
    element_type=GeneID
    output_file=MergedGeneCounts.tsv
    element_id_column=1

    if [[ "${strand}" == "unstranded" ]]
    then

        counts_column=2

    elif [[ "${strand}" == "stranded" ]]
    then

        counts_column=3

    else

        counts_column=4
            
    fi

elif [[ "${type}" == "Transcript" ]]
then

    header_row=2
    file_pattern=_abundance.tsv
    element_type=TranscriptID
    output_file=MergedTranscriptCounts.tsv
    element_id_column=1
    counts_column=4

elif [[ "${type}" == "Exon" ]]
then

    header_row=1
    file_pattern=_exon_counts.tsv
    element_type=ExonID
    output_file=MergedExonCounts.tsv
    element_id_column=1
    counts_column=7

else

    header_row=1
    file_pattern=_miRNA_counts.tsv
    element_type=miRNA
    output_file=MergedmiRNACounts.tsv
    element_id_column=1
    counts_column=7

fi

# N.B. This script assumes that all files processed have the same elements (genes or wathever)
# and that they are in the same order.

# Extract element ids
echo "Extracting element ids"
ls -1 *${file_pattern} | head -n 1 | xargs awk -v header_row=${header_row} '(NR > header_row) { print $1 }' > merged_temp.txt

# Extracting counts and merging files
echo "Merging count files"
for file in *${file_pattern}
do

    awk -v header_row=${header_row} -v counts_column=${counts_column} '(NR > header_row) { print $counts_column }' $file > counts_temp.txt

    paste merged_temp.txt counts_temp.txt > new_merged_temp.txt

    rm counts_temp.txt
    rm merged_temp.txt
    mv new_merged_temp.txt merged_temp.txt

done

# Removing duplicate lines. Duplicates happen if the GTF used for featureCounts had duplicated entries (common for exons, apparently)
echo "Removing duplicated lines"
sort merged_temp.txt | uniq >> merged_temp_unique.txt
rm merged_temp.txt

# Add header
echo "Adding header"
header=${element_type}

for file in *${file_pattern}
do

    sample_name=${file%$file_pattern*}
    header="${header}\t${sample_name}"

done

touch ${output_file}
echo -e "${header}" >> ${output_file}
cat merged_temp_unique.txt >> ${output_file}
rm merged_temp_unique.txt