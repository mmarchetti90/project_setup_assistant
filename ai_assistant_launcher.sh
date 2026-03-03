#!/bin/bash

set -e
set -o pipefail

### Main paths and variables

venv_path="/Users/earthsea/Documents/hf_transformers_env/.env/bin/activate"
ai_assistant_path="/Users/earthsea/Documents/hf_transformers_env/ai_assistant_v2"
module_root="src"
config_path="${ai_assistant_path}/config/config.json"
notes_path="${1:-missing}"
code_retrieval_stringency="${1:-false}"

if [ "${code_retrieval_stringency}" = "true" ]
then

    code_retrieval_stringency="--code_filtering"

else

    code_retrieval_stringency=""

fi

export PYTHONPATH="${ai_assistant_path}:$PYTHONPATH"

### Activate virtual environment

source ${venv_path}

### Options selection

PS3='Choose what to do: '

select run_type in "Summarize notes" "Gather resources" "Summarize notes and gather resources" "Exit"; do

    case $run_type in

        "Summarize notes")
            
            python -m ${module_root} \
            --config ${config_path} \
            --notes ${notes_path} \
            --mode summarize_notes
            
            break
            ;;

        "Gather resources")
            
            python -m ${module_root} \
            --config ${config_path} \
            --notes ${notes_path} \
            --mode gather_code \
            ${code_retrieval_stringency}

            /bin/bash gather_useful_code.sh

            break
            ;;
        
        "Summarize notes and gather resources")
            
            python -m ${module_root} \
            --config ${config_path} \
            --notes ${notes_path} \
            --mode summarize_notes_and_gather_code

            /bin/bash gather_useful_code.sh

            break
            ;;

        "Exit")
            
            break
            ;;

        *)
            
            echo "Invalid option."
            ;;
    
    esac

done