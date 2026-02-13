#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

import numpy as np

from collections.abc import Callable
from sentence_transformers import SentenceTransformer
from sys import argv
from torch import no_grad as torch_no_grad
#from torch import float16 as torch_float16
#from torch.cuda import is_available as cuda_is_available
#from transformers import BitsAndBytesConfig

### CLASSES AND FUNCTIONS ------------------ ###

def init_text_embedding_model(model_checkpoint: str, device_map: str='auto', similarity_fn_name: str='cosine') -> Callable:

    """
    Embedding for RAG
    """

    # Init model

    """
    quantization_config = BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_compute_dtype=torch_float16
    ) if cuda_is_available else None

    embedder = SentenceTransformer(
        model_checkpoint,
        similarity_fn_name=similarity_fn_name,
        device=device_map,
        model_kwargs={
            "low_cpu_mem_usage" : True,
            "quantization_config" : quantization_config
        }
    )
    """

    embedder = SentenceTransformer(
        model_checkpoint,
        similarity_fn_name=similarity_fn_name,
        device=device_map,
        model_kwargs={
            "low_cpu_mem_usage" : True
        }
    )

    return embedder

### ---------------------------------------- ###

def parse_obo(obo_file: str) -> dict[str, str]:
    
    # Read obo file and create hierarchy
    
    name_to_id = {}
    for term in open(obo_file).read().split('\n\n'):
        
        if term.startswith('[Term]'):
            
            # Parse term info
            
            term = term.split('\n')
            
            term_ids = [t.replace('alt_id: ', '').replace('id: ', '') for t in term if t.startswith('id: ') or t.startswith('alt_id: ')]
            
            term_name = [t.replace('name: ', '') for t in term if t.startswith('name: ')]
            term_name = term_name[0] if len(term_name) else ''
            
            # Update name_to_id
            
            for t in term_ids:
                
                name_to_id[term_name] = t
    
    return name_to_id

### ---------------------------------------- ###

@torch_no_grad()
def get_candidate_hpo(model: Callable, pheno_description: list[str], hpo_descriptions: list[str], max_hits: int=3) -> dict[str, list[str]]:

    # Embed pheno_description and hpo_descriptions

    pheno_embedding = model.encode(pheno_description)

    hpo_embedding = model.encode(hpo_descriptions)

    # Compute similarity of pheno_description to hpo_descriptions
    
    similarity = model.similarity(pheno_embedding, hpo_embedding)
    
    # Find top candidates
    
    relevant_hpo = {}
    
    for pheno, sim in zip(pheno_description, similarity):
        
        pheno_hpo_candidates_idx = np.argsort(sim.tolist())[::-1][:max_hits]
        
        pheno_hpo_candidates = [hpo_descriptions[hpo_idx] for hpo_idx in pheno_hpo_candidates_idx]
    
        relevant_hpo[pheno] = pheno_hpo_candidates

    return relevant_hpo

### ------------------MAIN------------------ ###

### Model vars

MODEL_CHECKPOINT = 'Qwen/Qwen3-Embedding-0.6B'

DEVICE_MAP = 'cpu'

SIMILARITY_FN_NAME = 'cosine'

### Parse CLI args

hpo_obo_path = argv[argv.index('--hpo_obo') + 1]
=
phenotype_path = argv[argv.index('--phenotype') + 1]

if '--max_candidates' in argv:

    max_candidates = int(argv[argv.index('--max_candidates') + 1])

else:

    max_candidates = 3

### Parse HPO obo

print('Parsing HPO obo file')

hpo_name_to_id = parse_obo(hpo_obo_path)

### Load phenotype

print('Parsing phenotype descriptions')

phenotype_data = open(phenotype_path, 'r').read().split('\n')

### First rough filtering based on word similarity (to reduce hpo terms to test since there's ~20K of them)

prepositions = ['of', 'in', 'to', 'for', 'with', 'on', 'at', 'from', 'by', 'about', 'as', 'into', 'like', 'through', 'after', 'over', 'between', 'out', 'against', 'during', 'without', 'before', 'under', 'around', 'and', 'among']

pheno_keywords = set(' '.join(phenotype_data).split(' '))

pheno_keywords = [k.lower() for k in pheno_keywords if k not in prepositions]

hpo_to_keep = []

for hpo_name in hpo_name_to_id.keys():
    
    hpo_keywords = [k.lower() for k in hpo_name.split(' ') if k not in prepositions]
    
    if any(k in pheno_keywords for k in hpo_keywords):
        
        hpo_to_keep.append(hpo_name)

hpo_name_to_id = {key : val for key,val in hpo_name_to_id.items() if key in hpo_to_keep}

### Find closest matching terms

print('Finding closest matching HPO terms')

text_embedding_model = init_text_embedding_model(MODEL_CHECKPOINT, DEVICE_MAP, SIMILARITY_FN_NAME)

all_candidate_hpo_terms = get_candidate_hpo(text_embedding_model, phenotype_data, list(hpo_name_to_id.keys()), max_candidates)

# Structure report

print('Reporting findings')

results_text = ['\t'.join(['phenotype', 'closest_hpo_ids', 'closest_hpo_description'])]

for pheno in phenotype_data:

    candidate_hpo_terms = all_candidate_hpo_terms[pheno]
    
    candidate_hpo_ids = [hpo_name_to_id[hpo] for hpo in candidate_hpo_terms]

    results_text.append('\t'.join([pheno, '; '.join(candidate_hpo_ids), '; '.join(candidate_hpo_terms)]))

results_text = '\n'.join(results_text)

with open('candidate_hpo_terms.tsv', 'w') as results_out:
    
    results_out.write(results_text)
