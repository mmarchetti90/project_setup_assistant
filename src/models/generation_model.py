#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

#from torch import float16 as torch_float16
#from torch.cuda import is_available as cuda_is_available
from transformers import (
    AutoModelForCausalLM,
    AutoTokenizer,
    BitsAndBytesConfig,
    TextGenerationPipeline
)

### CLASSES AND FUNCTIONS ------------------ ###

def init_text_generation_model(model_checkpoint: str, device_map: str='auto'):

    """
    Initialize the LLM
    """

    # Tokenizer
    
    tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)

    # Init model

    """
    quantization_config = BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_compute_dtype=torch_float16
    ) if cuda_is_available else None

    model = AutoModelForCausalLM.from_pretrained(
        model_checkpoint,
        use_safetensors=True,
        low_cpu_mem_usage=True,
        quantization_config=quantization_config,
        device_map=device_map
    )
    """

    model = AutoModelForCausalLM.from_pretrained(
        model_checkpoint,
        use_safetensors=True,
        low_cpu_mem_usage=True,
        device_map=device_map
    )

    # Init text generator

    text_generator = TextGenerationPipeline(
        model,
        tokenizer,
        framework="pt",
        task="text-generation"
    )

    return text_generator
