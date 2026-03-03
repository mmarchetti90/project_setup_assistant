#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

from collections.abc import Callable
from torch import no_grad as torch_no_grad

### CLASSES AND FUNCTIONS ------------------ ###

class rag_filter:

    """
    Class for selecting tags based on a query
    """
    
    # Variables for guiding LMM choice

    tool_type = "llm"
    
    name = 'rag_filter'
    
    description = """
    This is a tool that filters tags based on a query
    """
    
    inputs = {
        'query': {
            'type': 'string',
            'description': 'Context to guide filtering'
        },
        'tags': {
            'type': 'string',
            'description': 'List of tags to be filtered'
        }
    }
    
    output_type = 'list of strings'
    
    ### ------------------------------------ ###
    
    def __init__(self):
        
        self.is_initialized = False # For compatibility with smolagents
    
    ### ------------------------------------ ###
    
    @torch_no_grad()
    def forward(self, model: Callable, query: list[str], tags: list[str], score_threshold: float=0.75) -> dict[str]:

        # Embed query and tags

        query_embedding = model.encode(query)

        tags_embedding = model.encode(tags)

        # Compute similarity of query and tags
        
        similarity = model.similarity(query_embedding, tags_embedding)
        
        # Find top tags
        
        relevant_tags = {}
        
        for q, sim in zip(query, similarity):

            q_tags_candidates = [(t, s) for t,s in zip(tags, sim.tolist()) if s >= score_threshold]
        
            relevant_tags[q] = q_tags_candidates

        return relevant_tags
        
### ---------------------------------------- ###

if __name__ == '__main__':
    
    query = [
        "Middle Earth is home to many races, including dwarves, humans, and elves.",
        "Luke Skywalker is a central character in Star Wars.",
        "Captain America is a known Marvel superhero.",
        "Sauron waged war in Middle Earth.",
        "Eru Ilúvatar is the core deity in Tolkien's works.",
        "Bilbo Baggins had a great adventure in Middle Earth.",
        "Tattoine is a remote planet in the Star Wars universe."
    ]
    
    tags = ["Tolkien's Middle Earth.", 'Star Wars']

    text_filter = rag_filter()

    from src.models.embedding_model import init_text_embedding_model
    
    model = init_text_embedding_model(model_checkpoint='Qwen/Qwen3-Embedding-0.6B', device_map='cpu')

    output = text_filter.forward(model=model, query=query, tags=tags, score_threshold=0.5)

    output_text = []

    for q,ts in output.items():

        output_text.append('-' * 40)
        output_text.append(f'# {q}')
        output_text.append('\n'.join([f'* {t} (score={s:.3f})' for t,s in ts]) if len(ts) else 'None')
        output_text.append('-' * 40)

    output_text = '\n'.join(output_text)

    with open('test_tag_filter.txt', 'w') as out:
    
        out.write(output_text)
