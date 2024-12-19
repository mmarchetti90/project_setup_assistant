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

    model = BartForConditionalGeneration.from_pretrained(model_checkpoint)

    # Init summarizer

    summarizer = SummarizationPipeline(model, tokenizer, framework="pt", task="summarization")

    return summarizer

### ---------------------------------------- ###

def summarize_text(model, text):

    summary = ['### Notes']

    for t in text:

        if not len(t):

            continue

        words = len(t.split(" "))

        min_length = min(200, int(words * 0.60))
        max_length = min(300, int(words * 0.75))

        if words <= 125:

            summary.append(t)

        else:

            summary.append(model(t, min_length=min_length, max_length=max_length)[0]["summary_text"])

    summary = '\n\n'.join(summary)

    return summary

### ------------------MAIN------------------ ###

import json

from transformers import AutoTokenizer, BartForConditionalGeneration, SummarizationPipeline
from sys import argv

### Read args

config_path, notes_path = parse_args()

### Read config file

with open(config_path) as config_open:

    config = json.load(config_open)

### Load notes

notes = open(notes_path, "r").read().lower().split('\n\n')

### Init summarizer model

summarizer = init_model(config["summarization_model_path"])

### Generate summary

summary = summarize_text(summarizer, notes)

with open("ReadME.txt", "w") as summary_out:

	summary_out.write(summary)
