#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

from collections.abc import Callable
from torch import no_grad as torch_no_grad

### CLASSES AND FUNCTIONS ------------------ ###

class notes_summarizer:

    """
    Class for summarizing text
    """
    
    # Variables for guiding LMM choice

    tool_type = "llm"
    
    name = 'data_summarizer'
    
    description = """
    This is a tool that summarizes text
    """
    
    inputs = {
        'query': {
            'type': 'string',
            'description': 'Input data to summarize'
        },
        'context': {
            'type': 'string',
            'description': 'Context to guide the text generation'
        }
    }
    
    output_type = 'string'
    
    ### ------------------------------------ ###
    
    def __init__(self):
        
        self.is_initialized = False # For compatibility with smolagents

        self.output_block_tokens = ['<summary>', '</summary>']

        self.prompt = f"""
<|im_start|>system
You are an expert at extracting key points from a provided text.
[CONTEXT]
<im_end>
<|im_start|>user
Extract key information from following data: [QUERY]
Report the key information as bullet points wrapped within {self.output_block_tokens[0]} {self.output_block_tokens[1]} XML tags.
/no_think
<im_end>
<|im_start|>assistant

"""
    
    ### ------------------------------------ ###
    
    @torch_no_grad()
    def forward(self, model: Callable, query: str, context: str='', max_new_tokens: int=512) -> list[str]:

        # Update prompt

        updated_prompt = self.prompt.replace('[QUERY]', query)

        if context != '':

            updated_prompt = updated_prompt.replace('[CONTEXT]', f"Here's some context to guide the summary generation: {context}")

        else:

            updated_prompt = updated_prompt.replace('[CONTEXT]', '')

        # Run LLM

        llm_output = model(
            updated_prompt,
            return_full_text=False,
            max_new_tokens=max_new_tokens
        )[0]['generated_text']
        
        # Parse output

        parsed_llm_output = self.parse_llm_output(llm_output, self.output_block_tokens)

        return parsed_llm_output

    ### ------------------------------------ ###

    @staticmethod
    def parse_llm_output(text: str, delimiters: list[str]) -> str:

        parsed_text = text.strip()

        # Check start token
        
        if '<|im_start|>assistant' in parsed_text:

            parsed_text = parsed_text.split('<|im_start|>assistant')[-1].strip()
        
        if delimiters[0] in parsed_text:
                
            parsed_text = parsed_text.split(delimiters[0])[-1].strip()
                
            if delimiters[1] in parsed_text:
                
                parsed_text = parsed_text.split(delimiters[1])[0].strip()

        # Parse bullets

        info = [line.strip() for line in parsed_text.split('\n') if len(line)]

        return info
        
### ---------------------------------------- ###

if __name__ == '__main__':
    
    query = "The epidermal growth factor receptor (EGFR; ErbB-1; HER1 in humans) is a transmembrane protein that is a receptor for members of the epidermal growth factor family (EGF family) of extracellular protein ligands. The epidermal growth factor receptor is a member of the ErbB family of receptors, a subfamily of four closely related receptor tyrosine kinases: EGFR (ErbB-1), HER2/neu (ErbB-2), Her 3 (ErbB-3) and Her 4 (ErbB-4). In many cancer types, mutations affecting EGFR expression or activity could result in cancer. Epidermal growth factor and its receptor was discovered by Stanley Cohen of Vanderbilt University. Cohen shared the 1986 Nobel Prize in Medicine with Rita Levi-Montalcini for their discovery of growth factors. Deficient signaling of the EGFR and other receptor tyrosine kinases in humans is associated with diseases such as Alzheimer's, while over-expression  is associated with the development of a wide variety of tumors. Interruption of EGFR signalling, either by blocking EGFR binding sites on the extracellular domain of the receptor or by inhibiting intracellular tyrosine kinase activity, can prevent the growth of EGFR-expressing tumours and improve the patient's condition. == Function ==Epidermal growth factor receptor (EGFR) is a transmembrane protein that is activated by binding of its specific ligands, including epidermal growth factor and transforming growth factor alpha (TGF-\u03b1). ErbB2 has no known direct activating ligand, and may be in an activated state constitutively or become active upon heterodimerization with other family members such as EGFR.   Upon activation by its growth factor ligands, EGFR undergoes a transition from an inactive monomeric form to an active homodimer. \u2013 although there is some evidence that preformed inactive dimers may also exist before ligand binding."
    
    text_generator = notes_summarizer()

    from src.models.generation_model import init_text_generation_model
    
    model = init_text_generation_model(model_checkpoint='Qwen/Qwen3-0.6B', device_map='cpu')

    output = text_generator.forward(model=model, query=query)

    with open('test_summarizer.txt', 'w') as out:
    
        out.write('\n'.join(output))
