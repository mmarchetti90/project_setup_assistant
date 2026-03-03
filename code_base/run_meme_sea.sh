#!/bin/bash


fasta_dir=$(cd ${1} && pwd)
motif_db_dir=$(cd ${2} && pwd)
work_dir=$(pwd)

# Build list of motif databases to use
motif_db_list=""
for meme_db in ${motif_db_dir}/*.meme
do 

	motif_db_list+=" --m ${meme_db}"

done

# Run SEA
for fasta in ${fasta_dir}/*.fa
do

	# Making output directory
	out_dir=$(basename -s .fa ${fasta})
	mkdir ${out_dir}

	# SEA run
	echo -e "\n\n### Running MEME SEA for $(basename ${fasta})"
	docker run \
	--rm \
	-v ${fasta_dir}:${fasta_dir} \
	-v ${motif_db_dir}:${motif_db_dir} \
	-v ${work_dir}:${work_dir} \
	memesuite/memesuite \
	sea	--p ${fasta} ${motif_db_list} --oc ${work_dir}/${out_dir} --seed 42
	echo -e "Done!\n\n"

done