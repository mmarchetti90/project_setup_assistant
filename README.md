# Project Setup Assistant

AI assistant for setting up an analysis directory for a bioinformatics project

---

## OVERVIEW

* Download a local copy of needed models from Hugging Face

* Init project README from notes or email threads using "facebook/bart-large-cnn" for summarization

* Zero-shot classification for finding most relevant tags using "facebook/bart-large-mnli"

* Gather useful code/pipelines for analysis based on a tagged list of user-generated functions and code snippets

---

## SETUP

* Clone repository

* Update code manifest by running **init_code_manifest.py** in **code_base** directory

* Update paths in **ai_assistant_launcher.sh** and **config.json**

* Run **ai_assistant_launcher.sh** and download models (option 1)

---

## DEPENDENCIES

* Python 3.9.16+ with
    * torch 2.5.1+
    * transformers 4.47.0+
