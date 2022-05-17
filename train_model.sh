#!/bin/bash


#1. Fetch IDP4 dataset from web
wget -N https://s3.amazonaws.com/net.tagtog.public/resources/corpora/tagtog_IDP4%2B_anndoc.zip
unzip -qq tagtog_IDP4+_anndoc.zip -d data/nala
python nltk-download.py
python prepareData.py
python train_ner.py --model_name_or_path dmis-lab/biobert-base-cased-v1.1 --train_file data/nala/train_dev.json --validation_file data/nala/devel.json --text_column_name tokens --label_column_name tags --pad_to_max_length --max_length 192 --per_device_train_batch_size 8 --learning_rate 2e-5 --num_train_epochs 10 --output_dir models/nala --seed 1
