#!/bin/bash

genes_of_interest="$(pwd)/${1}"
promoter_set=${2}
work_dir=$(pwd)

# Making output directory
out_dir=$(basename -s .txt ${genes_of_interest})
mkdir ${out_dir}

# HOMER run
echo -e "\n\n### Running HOMER for $(basename ${genes_of_interest})"
docker run \
--rm \
-v ${genes_of_interest}:${genes_of_interest} \
-v ${work_dir}:${work_dir} \
findMotifs.pl \
${genes_of_interest} \
${promoter_set} \
${work_dir}/${out_dir}
echo -e "Done!\n\n"
