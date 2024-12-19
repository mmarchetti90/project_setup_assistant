#!/bin/bash

# N.B. very slow

header_row=1
file_pattern=_exon_counts.tsv
element_type=ExonID
output_file=MergedExonCounts.tsv
element_id_column=1
counts_column=7

# Extract unique element ids
echo "Extracting unique element ids"
touch element_ids.txt
for file in *${file_pattern}
do

	awk -v header_row=${header_row} '(NR > header_row) { print $1 }' $file | sort | uniq -d >> element_ids.txt

done

sort element_ids.txt | uniq -d > unique_element_ids.txt

rm element_ids.txt

# Touch output file
touch ${output_file}

# Add header
echo "Creating header"
header=${element_type}

for file in *${file_pattern}
do

	sample_name=${file%$file_pattern*}
	header="${header}\t${sample_name}"

done

echo -e "${header}" >> ${output_file}

# Extract counts for each unique element
echo "Extracting counts"
cat unique_element_ids.txt | while read element || [[ -n $element ]]
do
	
	new_line="${element}"

	# Extract count for each sample
	for file in *${file_pattern}
	do

		if (( $(grep -c $element $file) >= 1 ))
		then

			# Get element count from file
			new_count=$(grep $element $file | head -n 1 | awk -v counts_column=${counts_column} '{ print $counts_column }')

		else

			# Element missing, adding 0
			new_count=0

		fi

		new_line="${new_line}\t${new_count}"

	done

	# Add line to output file
	echo -e "${new_line}" >> ${output_file}

done
