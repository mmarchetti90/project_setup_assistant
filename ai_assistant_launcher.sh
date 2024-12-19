#!/bin/bash

set -e
set -o pipefail

### Main paths

venv_path="/Users/earthsea/Documents/hf_transformers_env/.env/bin/activate"
ai_assistant_path="/Users/earthsea/Documents/hf_transformers_env/ai_assistant"
config_path="${ai_assistant_path}/config.json"
notes_path="${1:-'missing'}"

### Activate virtual environment

source ${venv_path}

### Options selection

PS3='Choose what to do: '

select run_type in "Download models" "Summarize notes" "Gather resources" "Summarize notes and gather resources" "Exit"; do

    case $run_type in
        
        "Download models")
            
            python ${ai_assistant_path}/scripts/download_models.py \
            --config ${config_path}
            
            break
            ;;

        "Summarize notes")
            
            python ${ai_assistant_path}/scripts/summarize_notes.py \
            --config ${config_path} \
            --notes ${notes_path}
            
            break
            ;;

        "Gather resources")
            
            python ${ai_assistant_path}/scripts/gather_resources.py \
            --config ${config_path} \
            --notes ${notes_path}

            /bin/bash gather_useful_code.sh

            break
            ;;
        
        "Summarize notes and gather resources")
            
            python ${ai_assistant_path}/scripts/summarize_notes.py \
            --config ${config_path} \
            --notes ${notes_path}

            python ${ai_assistant_path}/scripts/gather_resources.py \
            --config ${config_path} \
            --notes ReadME.txt

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