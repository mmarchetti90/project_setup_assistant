#!/bin/bash

# Inputs

input_file=${1}
output_format=${2:-"html"}

run_pandoc () {

	in_file=${1}
	out_fmt=${2}
	out_file=$(echo "${in_file}" | sed "s/.md/.${out_fmt}/g")

	docker run \
	--rm \
	--volume "$(pwd):/data" \
	--user $(id -u):$(id -g) \
	pandoc/latex:3-ubuntu \
	${in_file} -o ${out_file} \
	--syntax-highlighting tango

}

run_pandoc ${input_file} ${output_format}
