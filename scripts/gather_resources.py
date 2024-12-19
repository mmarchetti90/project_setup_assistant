#!/usr/bin/env python3

### ---------------------------------------- ###

def parse_args():

    # Config file path
    config_path = argv[argv.index('--config') + 1]

    # Notes file path
    notes_path = argv[argv.index('--notes') + 1]

    return config_path, notes_path

### ---------------------------------------- ###

def init_model(model_checkpoint):

	# Tokenizer

	tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)

	# Init model

	#model_checkpoint = config["zeroshot_model_path"]
	model = BartForSequenceClassification.from_pretrained(model_checkpoint)

	# Init summarizer

	zeroshot = pipeline(task='zero-shot-classification', model=model, tokenizer=tokenizer, framework="pt")

	return zeroshot

### ---------------------------------------- ###

def refine_labels(tags, possible_tags):

	possible_tags = [pt.replace(" ", "_") for pt in possible_tags]

	print("\n" + "#" * 20)
	print("Project has been tagged as follows:")
	print(tags)

	print("Modify? (y/n)")
	answer = input().lower()

	if answer in ("y", "yes", "yay", "yeah", "hell yeah!"):

		print("\n" + "#" * 20)
		print("Choose among the following tags:")
		print(", ".join(possible_tags))

		toggle = True
		while toggle:
		
			updated_tags = input().replace(" ", ",").split(",")

			unrecognized_tags = [ut for ut in updated_tags if ut not in possible_tags]

			if len(unrecognized_tags):

				print("ERROR: unrecognized tags...")
				print(unrecognized_tags)

			else:

				toggle = False

	else:

		updated_tags = tags

	return updated_tags

### ---------------------------------------- ###

def gathering_code_generation(code_dir, labs, cdman):

	# Get code within scope of work
	useful_code = [c for c in cdman['code'].values() if c['main_tag'] in labs]

	# Define directories and subdirectories to be created and transfer instructions
	create_dir_commands, copy_code_commands = [], []
	for uc in useful_code:

		# Dir command
		target_dir = f'{uc["main_tag"]}/{uc["secondary_tag"]}'
		new_dir_command = f'mkdir -p useful_code/{target_dir}'

		if new_dir_command not in create_dir_commands:

			create_dir_commands.append(new_dir_command)

		# Copy command
		if uc["path"].startswith("git@github.com"):

			new_copy_command = f'cd useful_code/{target_dir}\n' + f'git clone {uc["path"]}' + '\ncd ${work_dir}'

		else:

			new_copy_command = f'cp {code_dir}/{uc["path"]} useful_code/{target_dir}/'

		copy_code_commands.append(new_copy_command)

	# Assemble shell script
	shell_script = (["#!/bin/bash",
					 "work_dir=$(pwd)",
					 "rm -r useful_code"] +
					 [""] +
					 create_dir_commands +
					 [""] +
					 copy_code_commands +
					 [""])

	shell_script = "\n".join(shell_script)

	return shell_script

### ------------------MAIN------------------ ###

import json

from transformers import AutoTokenizer, BartForSequenceClassification, pipeline
from sys import argv

### Read args

config_path, notes_path = parse_args()

### Read config file

with open(config_path) as config_open:

    config = json.load(config_open)

### Load notes

notes = open(notes_path, "r").read().lower()

truncation_idx = notes.index("### Notes") if "### Notes" in notes else 0

notes = notes[truncation_idx:]

### Load code manifest

with open(config["code_manifest_path"]) as code_manifest_open:

    code_manifest = json.load(code_manifest_open)

### Define keywords

general_labels = ['utils', 'visualizations']

candidate_labels = [l.replace('_', ' ') for l in code_manifest['scope_of_work'] if l not in general_labels]

### Init zeroshot model

zeroshot = init_model(config["zeroshot_model_path"])

### Classify

classification = zeroshot(notes, candidate_labels=candidate_labels)

threshold = 0.25

hit_labels = list(set([label for label,score in zip(classification['labels'], classification['scores']) if score >= threshold]))

### Refinement via user input

hit_labels = refine_labels(hit_labels, candidate_labels)

hit_labels += general_labels

### Gather resources

gathering_code = gathering_code_generation(config["code_dir_path"], hit_labels, code_manifest)

with open("gather_useful_code.sh", "w") as gathering_code_out:

	gathering_code_out.write(gathering_code)
