# Project Setup Assistant

AI agent for setting up an analysis directory for a bioinformatics project

---

## OVERVIEW

* Summarizes notes or email threads and outputs key information as bulletpoints describing the project

* Matches project descriptors to tags describing software in the code base

* Optionally further filters the selected software by matching the verbose software descrition to the project descriptors

* Gather useful code/pipelines for analysis

**N.B.** The project description can also be directly provided by the user as a bullet point of main tasks (one bullet per line).

---

## SETUP [UPDATE]

* Clone repository

* Update code manifest in **code_base** directory, if necessary

* Update paths and variables in **ai_assistant_launcher.sh** and **config.json**

* Run **ai_assistant_launcher.sh**

---

## DEPENDENCIES [UPDATE]

* Python 3.9.16+ with
    * sentence-transformers 5.1.2
    * torch 2.5.1+
    * transformers 4.56.2+
