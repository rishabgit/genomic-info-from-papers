#!/usr/bin/env bash

#BSUB -J train_var_ext
#BSUB -q production
#BSUB -n 1
#BSUB -M 64000
#BSUB -R "rusage[mem=64GB]"
#BSUB -o /hps/scratch/flicek/wormbase/var_extraction/Output_%J.out
#BSUB -e /hps/scratch/flicek/wormbase/var_extraction/Error_%J.err


VAREXTPATH=/hps/software/users/wormbase/variant-extraction/genomic-info-from-papers
VAREXTLOGS=/hps/scratch/flicek/wormbase/var_extraction


# Install dependencies and download nltk data
# -------------------------------------------------------------
module purge
module load python-3.9.10-gcc-9.3.0-i56je3q
cd "${VAREXTPATH}"
wget -N https://s3.amazonaws.com/net.tagtog.public/resources/corpora/tagtog_IDP4%2B_anndoc.zip
unzip -qq tagtog_IDP4+_anndoc.zip -d data/nala
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
mkdir data/nltk
export NLTK_DATA="${VAREXTPATH}/data/nltk"
python nltk-download.py


# Training:
# -------------------------------------------------------------
python train_ner.py --model_name_or_path dmis-lab/biobert-base-cased-v1.1 --train_file data/nala/train_dev.json --validation_file data/nala/devel.json --text_column_name tokens --label_column_name tags --pad_to_max_length --max_length 192 --per_device_train_batch_size 8 --learning_rate 2e-5 --num_train_epochs 10 --output_dir models/nala --seed 1
