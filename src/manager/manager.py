#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import json

from collections.abc import Callable
from datetime import datetime

### CLASSES AND FUNCTIONS ------------------ ###

class Manager:

    """
    Class for AI agent responsible for managing tasks execution
    """

    def __init__(
        self,
        code_base_path: str,
        code_manifest_path: str,
        generative_model: Callable,
        generative_model_checkpoint: str,
        zero_shot_model: Callable,
        zero_shot_model_checkpoint: str,
        embedding_model: Callable,
        embedding_model_checkpoint: str,
        device_map: str='cpu',
        max_new_tokens: int=1024,
        zero_shot_score_threshold: float=0.25,
        rag_similarity_function: str='cosine',
        rag_score_threshold: float=0.75,
        code_filtering=False

    ) -> None:

        """
        Initializes the AI Agent

        Parameters
        ----------
        code_base_path
            Path to local code base
        code_manifest_path
            Path to code manifest
        generative_model
            The LLM used to run the text generation tools
        generative_model_checkpoint
            Name or path of the checkpoint for the text generation tools
        zero_shot_model
            The model used for zero-shot classification
        zero_shot_model_checkpoint
            Name or path of the checkpoint for the zero-shot model
        embedding_model
            The model to be used for RAG
        embedding_model_checkpoint
            Name or path of the checkpoint for the RAG model
        device_map
            Device to use to run the models
            Default = 'cpu'
        max_new_tokens
            Maximum number of new tokens for the text generation tools
            Default = 1024
        zero_shot_score_threshold
            Threshold for matching tags to project descriptors using zero-shot
            Default = 0.25
        rag_similarity_function
            Similarity function for the RAG model
            Default = 'cosine'
        rag_score_threshold
            Threshold for matching tags to project descriptors using RAG
            Default = 0.75
        code_filtering
            Set to True if useful code should be selected using both scope of work and description
            Set to False if useful code should be selected using scope of work only
            Default = False
        """

        self.code_base_path = code_base_path

        self.code_manifest = [line.split('\t') for line in open(code_manifest_path, 'r').read().split('\n')[1:] if len(line)]

        self.init_timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S").replace(' ', '_')

        self.text_generation_model = generative_model(generative_model_checkpoint, device_map)

        #self.zero_shot_model = zero_shot_model(zero_shot_model_checkpoint, device_map) # cpu not working on Mac M1
        self.zero_shot_model = zero_shot_model(zero_shot_model_checkpoint, 'auto')

        self.text_embedding_model = embedding_model(embedding_model_checkpoint, device_map, rag_similarity_function)

        self.zero_shot_score_threshold = zero_shot_score_threshold

        self.max_new_tokens = max_new_tokens

        self.rag_score_threshold = rag_score_threshold

        self.code_filtering = code_filtering

    ### ------------------------------------ ###
    ### AGENT HISTORY                        ###
    ### ------------------------------------ ###

    def log_trace(self, trace_message: str) -> None:

        trace_timestamp = datetime.now().strftime("%H:%M:%S")

        timestamped_trace_message = f"[{trace_timestamp}] {trace_message}"

        print(timestamped_trace_message)

    ### ------------------------------------ ###
    ### TOOLS MANAGEMENT                     ###
    ### ------------------------------------ ###

    class empty_tool:

        """
        Mock tool that returns the query when run
        """

        def __init__(self):

            pass

        ### -------------------------------- ###

        def forward(self, query: str, *args, **kwargs) -> str:

            return query

    ### ------------------------------------ ###

    def add_tool(
        self,
        tool_call: Callable,
        is_notes_summarizer: bool=False,
        is_zero_shot_filter: bool=False,
        is_rag_filter: bool=False,
        log: bool=True
    ) -> None:

        """
        Registers a tool to the agent

        Parameters
        ----------
        tool_call : Callable
            A class with the following attributes: name, description, inputs, output_type
            Must also have the following method: forward (runs the tool) 
        is_notes_summarizer : bool
            Set to True to identify the tool for summarizing project notes
        is_zero_shot_filter : bool
            Set to True to identify the tool for filtering retrieved tags using zero shot classification
        is_rag_filter : bool
            Set to True to identify the tool for filtering retrieved tags using RAG
        log : bool
            Set to True to log the entry
        """

        if is_notes_summarizer:

            self.notes_summarizer = tool_call

            if log:

                self.log_trace(trace_message='Added tool: notes_summarizer')

        elif is_zero_shot_filter:

            self.zero_shot_filter = tool_call

            if log:

                self.log_trace(trace_message='Added tool: zero_shot_filter')

        elif is_rag_filter:

            self.rag_filter = tool_call

            if log:

                self.log_trace(trace_message='Added tool: rag_filter')

        else:

            self.tools[tool_call.name] = tool_call

            if log:

                self.log_trace(trace_message=f'Added tool: {tool_call.name}')

    ### ------------------------------------ ###
    ### RUN TASK                             ###
    ### ------------------------------------ ###

    def forward(self, mode: str, notes: str) -> None:

        accepted_modes = {
            'summarize_notes' : self.summarize_notes,
            'gather_code' : self.gather_code,
            'summarize_notes_and_gather_code' : self.summarize_and_gather
        }

        if mode not in accepted_modes.keys():

            print(f'ERROR: invalid run mode "{mode}"')

            return None

        else:

            _ = accepted_modes[mode](notes)

    ### ------------------------------------ ###

    def summarize_notes(self, notes: str) -> list[str]:

        # Extract most important info

        self.log_trace(trace_message='Reading notes and extracting most relevant info')

        project_descriptors = []

        while not len(project_descriptors):

            project_descriptors = self.notes_summarizer.forward(
                model=self.text_generation_model,
                query=notes,
                context='',
                max_new_tokens=self.max_new_tokens
            )

        self.log_trace(trace_message='  \u2714 Task completed, here are the identified project descriptors:')

        for pd in project_descriptors:

            self.log_trace(trace_message=f'    * {pd}')

        # Save to file

        with open('project_descriptors.txt', 'w') as output:

            output.write('\n'.join(project_descriptors))

        return project_descriptors

    ### ------------------------------------ ###

    def gather_code(self, notes: str) -> str:

        # Split notes into list of sentences

        notes = [n for n in notes.split('\n') if len(n)]

        # Extract unique scopes of work

        scopes = list(set('; '.join([cm[1] for cm in self.code_manifest]).split('; ')))

        # Find useful code based on scope of work

        self.log_trace(trace_message='Extracting most relevant scopes of work')

        selected_scopes = []

        while not len(selected_scopes):
            
            selected_scopes = self.zero_shot_filter.forward(
                model=self.zero_shot_model,
                query=notes,
                tags=scopes,
                score_threshold=self.zero_shot_score_threshold
            )

            selected_scopes = list(set([v[0] for val in selected_scopes.values() for v in val]))

        self.log_trace(trace_message='  \u2714 Task completed. Project has been tagged as follows:')

        for s in selected_scopes:

            self.log_trace(trace_message=f'    * {s}')

        # Manual edit
            
        self.log_trace(trace_message='Modify? (y/n)')

        answer = input().lower()

        if answer in ("y", "yes", "yay", "yeah", "hell yeah!"):

            toggle = True

            while toggle:
                
                self.log_trace(trace_message="Choose among the following scopes of work:")

                self.log_trace(trace_message=', '.join(scopes))

                selected_scopes = [s.strip(' ') for s in input().split(",")]

                unrecognized_scopes = [s for s in selected_scopes if s not in scopes]

                if len(unrecognized_scopes):

                    self.log_trace(trace_message=f'ERROR: unrecognized tags: {unrecognized_scopes}')

                else:

                    toggle = False

        # Extract unique descriptions

        descriptions = list(set([cm[2] for cm in self.code_manifest if any(s in selected_scopes for s in cm[1].split('; '))]))

        # Filter found code based on description

        if not self.code_filtering:

            selected_descriptions = descriptions

        else:

            self.log_trace(trace_message='Filtering code by description')

            selected_descriptions = []

            while not len(selected_descriptions):
                
                selected_descriptions = self.rag_filter.forward(
                    model=self.text_embedding_model,
                    query=notes,
                    tags=descriptions,
                    score_threshold=self.rag_score_threshold
                )

                selected_descriptions = list(set([v[0] for val in selected_descriptions.values() for v in val]))

            self.log_trace(trace_message='  \u2714 Task completed')

        # Writing script for gathering useful code

        self.log_trace(trace_message='Writing script for gathering useful code')

        shell_script = self.gathering_code_generation(
            code_dir=self.code_base_path,
            manifest=self.code_manifest,
            scopes=selected_scopes,
            descriptions=selected_descriptions
        )

        with open('gather_useful_code.sh', 'w') as shell_script_out:

            shell_script_out.write(shell_script)

        self.log_trace(trace_message='  \u2714 Task completed')

        return shell_script

    ### ------------------------------------ ###

    def summarize_and_gather(self, notes: str) -> tuple[list[str], str]:

        # Summarize notes

        project_descriptors = self.summarize_notes(notes)

        # Gather useful code

        shell_script = self.gather_code('\n'.join(project_descriptors))

        return project_descriptors, shell_script

    ### ------------------------------------ ###

    @staticmethod
    def gathering_code_generation(
        code_dir: str,
        manifest: list[str],
        scopes: list[str],
        descriptions: list[str]
    ) -> str:

        # Init shell script

        shell_script = [
            '#!/bin/bash',
            'work_dir=$(pwd)'
        ]

        # Define directories and subdirectories to be created and transfer instructions

        output_dir = 'useful_code'

        shell_script.append(f'rm -r {output_dir}')

        mkdir_commands, cp_commands = [], []

        for (code, scope_of_work, description, location) in manifest:

            # Only gather code if scope in scopes and description in descriptions

            matching_scopes = [sow for sow in scope_of_work.split('; ') if sow in scopes]

            if len(matching_scopes) and description in descriptions:

                # Define mkdir command

                output_subdir = '_'.join(matching_scopes)

                full_output_path = f'{output_dir}/{output_subdir}'.replace(' ', '_')

                new_mkdir_command = f'mkdir -p {full_output_path}'

                if new_mkdir_command not in mkdir_commands:

                    mkdir_commands.append(new_mkdir_command)

                # Define copy command

                if location.lower() == 'local':

                    new_cp_command = f'cp {code_dir}/{code} {full_output_path}/'

                elif location.startswith('git@github.com'):

                    new_cp_command = f'cd {full_output_path}\n' + f'git clone {location}' + '\ncd ${work_dir}'

                else:

                    continue

                if new_cp_command not in cp_commands:

                    cp_commands.append(new_cp_command)

        # Assemble shell script
        
        shell_script += (
            [""] +
            mkdir_commands +
            [""] +
            cp_commands +
            [""]
        )

        shell_script = "\n".join(shell_script)

        return shell_script
