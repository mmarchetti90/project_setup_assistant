#!/usr/bin/env python3

### IMPORTS -------------------------------- ###

from collections.abc import Callable
from torch import no_grad as torch_no_grad
from transformers import (
    AutoTokenizer,
    BitsAndBytesConfig,
    BartForSequenceClassification,
    pipeline
)

### CLASSES AND FUNCTIONS ------------------ ###

def init_zero_shot_model(model_checkpoint: str, device_map: str='auto'):

	# Tokenizer

	tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)

	# Init model

	model = BartForSequenceClassification.from_pretrained(
		model_checkpoint,
        use_safetensors=True,
        low_cpu_mem_usage=True,
        device_map=device_map
    )

	# Init summarizer

	zeroshot = pipeline(
		model=model,
		tokenizer=tokenizer, 
		framework="pt",
		task='zero-shot-classification'
	)

	return zeroshot
