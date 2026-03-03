#!/bin/bash

: '
This script parses the read identifiers of mate fastq files and removes non-paired reads.
If desired, files are also sorted by identifier after filtering.
'

### Args

r1=${1}

r2=${2}

sort_toggle=${3:-false}

### File name to be used for output

out_prefix=$(basename ${r1} | sed "s/_R1.fastq//g" | sed "s/_R1.fq//g" | sed "s/.gz//g")

### Edit read identifiers to avoid grep mistakes (e.g. grep "@12345" would find a match in @12345, @123456, @1234567, etc...)

if [[ ${r1: -2} == "gz" ]]
then

	zcat ${r1} | awk '{ if($0 ~ /@/) { print $0":A" } else { print $0 } }' > ${out_prefix}_R1_mod.fq

	zcat ${r2} | awk '{ if($0 ~ /@/) { print $0":A" } else { print $0 } }' > ${out_prefix}_R2_mod.fq

else

	awk '{ if($0 ~ /@/) { print $0":A" } else { print $0 } }' ${r1} > ${out_prefix}_R1_mod.fq

	awk '{ if($0 ~ /@/) { print $0":A" } else { print $0 } }' ${r2} > ${out_prefix}_R2_mod.fq

fi

### Filter fastq pair by only keeping common identifiers

# Extract identifiers

cat ${out_prefix}_R1_mod.fq | grep "@" > r1_identifiers.txt
cat ${out_prefix}_R2_mod.fq | grep "@" > r2_identifiers.txt

# Filter

cat ${out_prefix}_R1_mod.fq | grep -f r2_identifiers.txt -A 3 | sed '/^--$/d' | gzip > ${out_prefix}_R1_fixed.fq.gz
cat ${out_prefix}_R2_mod.fq | grep -f r1_identifiers.txt -A 3 | sed '/^--$/d' | gzip > ${out_prefix}_R2_fixed.fq.gz

# Cleanup

rm ${out_prefix}_R1_mod.fq ${out_prefix}_R2_mod.fq

rm r1_identifiers.txt r2_identifiers.txt

### Sort reads (if desired)

if [ ${sort_toggle} = true ]
then

	zcat ${out_prefix}_R1_fixed.fq.gz | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip > ${out_prefix}_R1_fixed_sorted.fq.gz

	zcat ${out_prefix}_R2_fixed.fq.gz | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" | gzip > ${out_prefix}_R2_fixed_sorted.fq.gz

fi
