#!/usr/bin/env python3

### ---------------------------------------- ###

def parse_args():

    # Config file path
    config_path = argv[argv.index("--config") + 1]

    return config_path

### ------------------MAIN------------------ ###

import json

from huggingface_hub import snapshot_download
from os import getcwd
from sys import argv

### Read config file

config_path = parse_args()

with open(config_path) as config_open:

    config = json.load(config_open)

### Download models

for model_type in ["summarization", "zeroshot"]:

    # Gather info

    model_name = config[f'{model_type}_model']

    if config[f'{model_type}_model_path'] == "":

        model_local_dir = model_name.split("/")[-1]

    else:

        model_local_dir = config[f'{model_type}_model_path']

    # Download

    print(f'Downloading {model_name}')

    snapshot_download(repo_id=model_name,
                      local_dir=model_local_dir,
                      force_download=True,
                      ignore_patterns=config["ignored_patters"])

    print("All done\n")

    # Update config

    config[f'{model_type}_model_path'] = f'{getcwd()}/{model_local_dir}'

### Save updated config

with open(config_path, "w") as config_out:
        
    json.dump(config, config_out, indent=4)
