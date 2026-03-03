#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

from sentence_transformers import SentenceTransformer
#from torch import float16 as torch_float16
#from torch.cuda import is_available as cuda_is_available
#from transformers import BitsAndBytesConfig

### CLASSES AND FUNCTIONS ------------------ ###

def init_text_embedding_model(model_checkpoint: str, device_map: str='auto', similarity_fn_name: str='cosine'):

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
