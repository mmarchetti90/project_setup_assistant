#!/bin/bash

fasta_dir=$(cd ${1} && pwd)
work_dir=$(pwd)

# Run HOMER
for fasta in ${fasta_dir}/*.fa
do

	# Making output directory
	out_dir=$(basename -s .fa ${fasta})
	mkdir ${out_dir}

	# HOMER run
	echo -e "\n\n### Running HOMER for $(basename ${fasta})"
	docker run \
	--rm \
	-v ${fasta_dir}:${fasta_dir} \
	-v ${work_dir}:${work_dir} \
	findMotifs.pl \
	${fasta} \
	fasta \
	${work_dir}/${out_dir}
	echo -e "Done!\n\n"

done