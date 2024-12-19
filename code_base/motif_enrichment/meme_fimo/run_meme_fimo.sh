#!/bin/bash

fasta_dir=$(cd ${1} && pwd)
motif_file_dir="$(cd ${2} && pwd)"
motif_file="${motif_file_dir}/${3}"
work_dir=$(pwd)

for fasta in ${fasta_dir}/*.fa
do

	# Making output directory
	out_dir=$(basename -s .fa ${fasta})
	mkdir ${out_dir}

	# MEME FIMO to find EcR motifs
	echo -e "\n\n### Running MEME FIMO for $(basename ${fasta})"
	docker run \
	--rm \
	-v ${fasta_dir}:${fasta_dir} \
	-v ${motif_file}:${motif_file} \
	-v ${work_dir}:${work_dir} \
	memesuite/memesuite \
	fimo --oc ${work_dir}/${out_dir} ${motif_file} ${fasta}
	echo -e "Done!"

done
