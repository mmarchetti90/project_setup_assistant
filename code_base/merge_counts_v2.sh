#!/bin/bash

# N.B. This script assumes that all files processed have the same elements (genes or wathever)
# and that they are in the same order. If these assumptions are wrong, use version 1

header_row=1
file_pattern=_exon_counts.tsv
element_type=ExonID
output_file=MergedExonCounts.tsv
element_id_column=1
counts_column=7

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
