# Variant extraction from papers

## Install

- Create a virtual environment with Python-3.9 (only once):

  `python3 -m venv venv`

- Activate the virtual environment each time:

  `source venv/bin/activate`

- Install dependencies (only once):

  `pip install -r requirements.txt`


## Train model (only once)

- Download nltk and related tools, prepare data for training and train model:

  `bash train_model.sh`


## Configure credentials

Edit:

   `utils/all_config.cfg`

With wbtools access credentials

## Tests

- To test regexes:

 `python test_regexes.py`

- To fetch some papers locally:

 `python fetch_papers.py `

- To test the variant extraction on some papers

  `python hybrid_extraction.py`


## Execution

  `python main.py`
