#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import json

from collections.abc import Callable
from datetime import datetime
from gc import collect as gc_collect
from .manager.manager import Manager
from .models import *
from .tools import *
from sys import argv

### CLASSES AND FUNCTIONS ------------------ ###

def parse_args():

    # Config file path

    config_path = argv[argv.index('--config') + 1]

    # Notes file path

    notes_path = argv[argv.index('--notes') + 1]

    # Model

    mode = argv[argv.index('--mode') + 1]

    # Code filtering

    code_filtering = ('--code_filtering' in argv)

    return config_path, notes_path, mode, code_filtering

### ---------------------------------------- ###

def parse_config(config_path: str) -> dict:

    """
    Parsing config file containing global vars
    """

    with open(config_path) as opened_config_file:
    
        config_data = json.load(opened_config_file)

    return config_data

### ---------------------------------------- ###

if __name__ == '__main__':

    # Parsing CLI args

    config_path, notes_path, mode, code_filtering = parse_args()

    # Load data

    config_vars = parse_config(config_path)

    notes = open(notes_path, 'r').read()

    if "### Notes" in notes:

        notes = notes.replace('### Notes:', '### Notes').split('### Notes')[-1]

    # Init agent

    manager_instance = Manager(
        code_base_path = config_vars['CODE_BASE_PATH'],
        code_manifest_path = config_vars['CODE_MANIFEST_PATH'],
        zero_shot_model = zero_shot_model.init_zero_shot_model,
        zero_shot_model_checkpoint = config_vars['ZERO_SHOT_MODEL'],
        generative_model = generation_model.init_text_generation_model,
        generative_model_checkpoint = config_vars['TEXT_GENERATION_MODEL'],
        embedding_model = embedding_model.init_text_embedding_model,
        embedding_model_checkpoint = config_vars['TEXT_EMBEDDING_MODEL'],
        device_map = config_vars['DEVICE_MAP'],
        max_new_tokens = config_vars['MAX_NEW_TOKENS'],
        zero_shot_score_threshold = config_vars['ZERO_SHOT_SCORE_THRESHOLD'],
        rag_similarity_function = config_vars['RAG_SIMILARITY_FUNCTION'],
        rag_score_threshold = config_vars['RAG_SCORE_THRESHOLD'],
        code_filtering = code_filtering
    )

    # Add tools

    tool_function = notes_summarizer.notes_summarizer
    manager_instance.add_tool(tool_function(), is_notes_summarizer=True)

    tool_function = zero_shot_filter.zero_shot_filter
    manager_instance.add_tool(tool_function(), is_zero_shot_filter=True)

    tool_function = rag_filter.rag_filter
    manager_instance.add_tool(tool_function(), is_rag_filter=True)

    # Run agent

    manager_instance.forward(mode=mode, notes=notes)

    # Cleanup

    del manager_instance

    gc_collect()